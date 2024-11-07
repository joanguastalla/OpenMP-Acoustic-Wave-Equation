typedef struct{
	int ixb;
	int izb;
	int ize;
	int nb;
	int nx;
	int nz;
	int nxx;
	int nzz;
	float x0;
	float z0;
	float dx;
	float dz;
}vmodel;
struct border{
    int ixb;
    int izb;
    int ize;
    int nb;
    int nxx;
    int nzz;
};
typedef struct{
	float xsmax;
	float xsmin;
	float dxs;
	float zs;
	float hmin;
	float hmax;
	float dh;
	float zg;
	float freq;
	float ttotal;
}aquisition;


void velextension(float* vex,float* v,int nx,int nz,struct border* borda);
float fricker(float t,float freq);
