#include <math.h>
#include "modelling_utils.h"
#define  PI 3.14159265358979323846264338328
void velextension(float* vex,float* v,int nx,int nz,struct border* borda){
	int i0;
	i0=(borda->ixb+1)*borda->nzz + borda->izb+1;
	
	for(int ii=0;ii< nx;ii++){
            for(int jj=0;jj<nz;jj++){
                vex[(int) (i0 + ii*borda->nzz + jj)]=v[ii*nz + jj];
            }
    }
    // model at the left side
        for(int ii=0;ii<(borda->ixb + 1);ii++){
            for(int jj=0;jj<nz;jj++){
					vex[borda->izb +  1 + ii*borda->nzz + jj]=v[jj];
             }
        }
        
    // model at the right side
        for(int ii=0;ii<(borda->ixb+1);ii++){
          for(int jj=0;jj<nz;jj++){
            vex[(borda->nzz*( borda->ixb+1 + nx) +  borda->izb + 1) 
				+ ii*borda->nzz + jj]=v[jj];
          }
        }
   
	// model at the top side
		for (int ii=0;ii<borda->nxx;ii++){
			for(int jj=0;jj<(borda->izb+1);jj++){
				vex[jj + ii*borda->nzz] = vex[(borda->izb+1) + ii*borda->nzz];
			}
		}	
	// model at the bottom side
		for (int ii=0;ii< borda->nxx;ii++){
			for(int jj=borda->ize;jj < borda->nzz;jj++){
				vex[jj + ii*borda->nzz] = vex[(borda->ize-1) +  ii*borda->nzz];
			}
		}	
}

float fricker(float t,float freq){
		float beta=pow(PI*freq*t,2);
		float ricker;	
		ricker=(1.0 - 2*beta)*exp(-beta);
		return ricker;
}
