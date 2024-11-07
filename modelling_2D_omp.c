#include <stdio.h>      // Input and output
#include <stdlib.h>     // Pointers
#include <string.h>
#include <math.h>		// Math functions
#include "omp.h"		// Parallelism OpenMP
#include "modelling_utils.h"
#define  PI 3.14159265358979323846264338328
#define dtrec 0.004
#define sec2ms 1.0e6
#define cfl 0.25
#define NDERIV 3
#define NTHREADS 4
#define NB 800
#define CACHE_SZ 112

void swap(float** a,float** b){
    float* tmp;
    tmp=(*a);
    *a=*b;
    *b=tmp;
}



int main(){    
	// OpenMP settings
	int nthreads;
	int max_threads;
	max_threads=omp_get_max_threads();
	printf("Maximum number of threads:%d\n",max_threads);
	omp_set_num_threads(NTHREADS);
    printf("OpenMP version:%d\n",_OPENMP);

	// wavefield 
	float *p1 ;
    float *p2 ;
	// Parameters on original velocity model 
    int nvparam,nvparam_read,nx,nz;    
	float dx,dz,x0,z0;
	float* vel;

	// Parameters on extended velocity model
	int nxx,nzz,nbx,nbz,nbzcache;
	struct border* borda=malloc(sizeof(struct border));
	float* velextend;
	int i0_model,i0_extended_model;
	
	// Aquisition parameters
	int nparam,nparam_read,it0;
	int isrc,isx,isz,idxs,nshots;
	float xsmax,xsmin,dxs,zs,xg_offset_min,xg_offset_max,dxg,zg ;
	float freq,ttotal;
	float* ricker;
	
	// Stability parameters
	int ns,nt,ndt;
	float maxv,minv;

	float dxmax,dt;
	float deriv[NDERIV];
	deriv[0]=(float)-5/2;
	deriv[1]=(float) 4/3;
	deriv[2]=(float)-1/12;
	float derivsum=0;


	// Directories for reading or writing files  
    char* param_dir;
	char* velbin;
    char* velhdr;
    char* velextend_path;
	//char* moviedir;
	param_dir="./param.txt" ;
    velbin = "./vel.bin";
	velhdr = "./vel.h";
	velextend_path="./velextend.bin";
	//moviedir="./wavemovie.bin";
	
	// File pointers to open files
	FILE* velfile;
    FILE* velfile_hdr;
    FILE* velwrite;
    FILE* param;
  //  FILE* wavemovie;
    
	//Velocity model reading: hdr e bin
	nparam=10;
	nvparam=6;
    velfile_hdr=fopen(velhdr,"r");    
    velfile=fopen(velbin,"rb");
	param=fopen(param_dir,"r");
	
	nparam_read=fscanf(param,"%f %f %f %f %f %f %f %f %f %f",&xsmax,&xsmin,&dxs,&zs,&xg_offset_min,&xg_offset_max,&dxg,&zg,&freq,&ttotal);
	if(nparam_read != nparam){
		printf("Too few acquisition parameters!!");
		exit(EXIT_FAILURE);
	}
	
	nvparam_read=fscanf(velfile_hdr,"%i %i %f %f %f %f",&nz,&nx,&z0,&x0,&dz,&dx);
	if(nvparam_read != nvparam){
		printf("Too few acquisition parameters!!");
		exit(EXIT_FAILURE);
	}


    //Print of external parameters

	printf("smax=%f\n",xsmax);
	printf("smin=%f\n",xsmin);
	printf("dxs=%f\n",dxs);
	printf("zs=%f\n",zs);
	printf("hmin=%f\n",xg_offset_min);
	printf("hmax=%f\n",xg_offset_max);
	printf("dh=%f\n",dxg);
	printf("zg=%f\n",zg);
	printf("freq=%f\n",freq);
	printf("time=%f\n",ttotal);


	vel=malloc(sizeof(float)*nx*nz);
	if(!fread(vel,sizeof(float),nx*nz,velfile)){
		printf("Error reading velocity model!!\n");
		exit(EXIT_FAILURE);
	}
    fclose(velfile);
    fclose(velfile_hdr);

	// Extending velocity model *******************************************
	nxx=2*NB + nx;
	nbx=nxx - 2*(NDERIV-1);
	nbz=((nz+2*(NB -(NDERIV-1)))/CACHE_SZ + 1)*CACHE_SZ;
	borda->nb=NB;
	borda->izb=(int) (nbz-nz)/2 - 1 + (NDERIV-1);
	borda->ize=borda->izb + 1 + nz;
	nbzcache=nbz/CACHE_SZ;
    nzz=nbz + 2*(NDERIV-1);
    borda->ixb=NB -1;
    borda->nxx=nxx;
    borda->nzz=nzz;
	
   


	//  Maximum and minimum of velocity model ******************************
	maxv=vel[0];
	minv=vel[0];
	for(int ii=0;ii<nx*nz;ii++){
		if(maxv < vel[ii]){
			maxv=vel[ii];
		}
		if(minv > vel[ii]){
			minv=vel[ii];
		}
	}

	// Setting parameters dt,dx for stability of FD
	ns= (int) (ttotal/dtrec) ;
	ns+=1;
	dxmax=minv/(6.0*freq);
	for(int ii=0;ii<NDERIV;ii++){
		derivsum+=abs(deriv[ii]);
	}
	dt=cfl*sqrt(2.0/derivsum)*dx/maxv;
	ndt=(int)(dtrec/dt);
	ndt+=1;

	// Chosing dt for modelling to be a a divisor of dtrec
	if(ndt > 1){
		dt=(dtrec/ndt);
	}else{
		dt=dtrec;
	}



	if(dx > dxmax){
		printf("Frequency is too high");
		return 0;
	}
	

	

	printf("Dt for stability: %.8f\n",dt);
	printf("Dx for stability:%f\n",dxmax);
	printf("ndt is:%i \n",ndt);
	printf("Maximum value of velocity: %f \n",maxv);
	printf("Minimum value of velocity: %f \n",minv);
	printf("Number of border points:%d\n",borda->nb);
	

	velextend=malloc(sizeof(float)*nxx*nzz);
	velextension(velextend,vel,nx,nz,borda);
    velwrite=fopen(velextend_path,"wb");
    fwrite(velextend,sizeof(float),nxx*nzz,velwrite);    


	// Scaling  extended velocity for computational resource economy
	for(int ii=0;ii<nxx;ii++){
		for(int jj=0;jj<nzz;jj++){
			velextend[jj + ii*nzz]= powf(velextend[jj + ii*nzz],2.0)*powf(dt/dx,2.0);
		}
	}

	// Wavefield temporal evolution common-shot situation
	p1=calloc(nxx*nzz,sizeof(float));
	p2=calloc(nxx*nzz,sizeof(float));
	i0_model=nzz*(borda->ixb + 1) + (borda->izb+1);
	i0_extended_model=nzz*(NDERIV-1) + (NDERIV-1);
	

	// Index of source in the model
	isx=(int) (xsmin/dx);   // starting at 0 sample
	isz=(int) (zs/dz);	   // starting at 0 sample
	isrc=i0_model + nzz*(isx) + isz;
	idxs=(int)(dxs/dx);
	idxs=idxs*nzz;
	nshots=(int) ((xsmax - xsmin)/dxs);
	nshots+=1;
	
	// Border parameters
	int ix;
	float gamma_x[nxx];
	float gamma_z[nzz];
	float gamma,beta;
	beta=PI*freq*dt;
	for(int iz = 0; iz < nzz; ++iz) {
		gamma_z[iz]=0.0;
	}
	for (ix = 0; ix < nxx; ++ix) {
		gamma_x[ix]=0.0;
	}
	for (ix = 0; ix < borda->nb; ++ix) {
		gamma=beta*pow((((float) (ix+1))/((float) borda->nb)),2.0);
		gamma_x[borda->ixb - ix]=gamma;
		gamma_x[(borda->ixb + nx + 1) + ix]=gamma;
	}	

	for(int iz=0;iz<(borda->izb+1);iz++){
		gamma=beta*pow((((float) (iz+1))/((float) borda->izb + 1)),2.0);
		gamma_z[borda->izb - iz]=gamma;
	}

	for (int iz = 0; iz < borda->nb; ++iz) {
		gamma=beta*pow((((float) (iz+1))/((float) borda->nb)),2.0);
		gamma_z[borda->ize + iz]=gamma;
	}	

	nt=(int)(ttotal/dt);
	nt+=1;
	printf("nt:%d\n",nt);
	ricker=malloc(sizeof(float)*nt);
	// Half duration of ricker number of indexes it0	
	it0 = floor(sqrt(6.0/PI)/(freq*dt)); 
	for(int it=0;it<nt;it++){
		ricker[it]=fricker((it-it0)*dt,freq);
	}
	printf("Requested number of threads on OpenMP:%d\n",NTHREADS);
	printf("nzz:%d\n",nzz);
	printf("nxx:%d\n",nxx);
    int it=0,ishot=0,srcindex;
    float *temp=NULL,tempcount;
    FILE* wavemovie=fopen("./wavemovie.bin","wb");
    omp_lock_t lock;
	#pragma omp parallel  default(none) shared(p1,p2,wavemovie,isrc,srcindex,ishot,nbz,idxs,deriv,ricker,dx,velextend,nbzcache,nbx,i0_extended_model,i0_model,gamma_x,gamma_z,nx,nxx,nzz,nz,nt,ndt,nshots,nthreads)
		{
            int iconv,ii,ix,it=0;
			int mm;
			int	imodel;
			int iz,niz;
			float gama,mgamma,invpgamma;
			float laplacian;
			float tempfield[CACHE_SZ];
            float *ip1=p1,*ip2=p2;
			nthreads=omp_get_num_threads();
            #pragma omp single
			{   
				printf("Number of threads used:%d\n",nthreads);
            }
            while(ishot<nshots){
                #pragma omp single
                {
                    it=0;
                    srcindex=isrc + idxs*ishot;
                }
                #pragma omp for        
                for(int iswap=0;iswap<nxx*nzz;iswap++){
                    p1[iswap]=0.0;
                    p2[iswap]=0.0;
                }
            
                for(it=0;it<nt;it++){
                    #pragma omp for schedule(static)
                    for(mm=0;mm<nbzcache*nbx;mm++){
                            ix=(int) mm/nbzcache;
                            niz=mm%nbzcache;
                            iz=niz*CACHE_SZ; 
                        for(ii=0;ii<CACHE_SZ;ii++){
                            imodel=i0_extended_model + ix*nzz +iz + ii;
                            gama=gamma_x[ix + NDERIV-1] + gamma_z[iz + NDERIV-1 +ii];
                            mgamma=-(1-gama)/(1+gama);
                            invpgamma=0.5*(1-mgamma);
                            laplacian=2*deriv[0]*p2[imodel];
                            for(iconv=1;iconv< NDERIV;iconv++){
                                laplacian+=deriv[iconv]*p2[imodel - iconv];
                                laplacian+=deriv[iconv]*p2[imodel + iconv];
                                laplacian+=deriv[iconv]*p2[imodel - nzz*iconv];      		
                                laplacian+=deriv[iconv]*p2[imodel + nzz*iconv];
                            }
                            tempfield[ii]=invpgamma*velextend[imodel]*laplacian;
                            tempfield[ii]+=invpgamma*2.0*p2[imodel];
                            tempfield[ii]+=mgamma*p1[imodel];
                        }
                        for(ii=0;ii<CACHE_SZ;ii++){
                            p1[imodel-(CACHE_SZ-1) + ii]=tempfield[ii];
                        }
                    }
                    #pragma omp single 
                    {
                        p1[srcindex]+= -pow(dx,2.0)*ricker[it]*velextend[srcindex];
                        swap(&p1,&p2);
                    }
                   #pragma omp single nowait
                    {
                        if(it%ndt==0){
                            for(ii=0;ii<nx;ii++){
                            fwrite(p2+i0_model + ii*nzz,sizeof(float),nz,wavemovie);
                            }
                        }
                    }
                }
                #pragma omp single
                {
                    ishot++;
                    printf("ishot:%d\n",ishot);
                    printf("nshots:%d\n",nshots);
                }
          }
    }
	        return 0;
}
