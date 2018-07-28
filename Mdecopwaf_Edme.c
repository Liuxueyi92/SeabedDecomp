/*
 * Decomposition wavefiled of lateral varied velocity media
 * (without zero padding operation during DFT)
 * Author:liuxueyi    DATE:2018.4.15 
 * Tongji University
*/

#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <rsf.h>
#include <math.h>
#define PI  3.141592654

int     nx,nt;
float   dx,dt;


float absolutevalue(float in);
void  puthead2d(sf_file output,int n1,float d1,int n2,float d2,float o1,float o2);
void  ts2fs(float **vxdata,float **vzdata,float **phdata,sf_complex **vxdatafs,sf_complex **vzdatafs,sf_complex **phdatafs);
void  fs2fk(sf_complex **vxdatafs,sf_complex **vzdatafs,sf_complex **phdatafs,sf_complex **uppurefk,sf_complex **usppurefk);
void  fk2ts(sf_complex **uppurefk,sf_complex **uspurefk,float **uppurets,float **uspurets);

int main(int argc,char *argv[])
{
	sf_init(argc,argv);
	sf_file vxfile,vzfile,phfile;
	vxfile=sf_input("in");
	vzfile=sf_input("vzrecord");
	phfile=sf_input("precord");

	sf_file uppure,uspure;
	uppure=sf_output("out");
	uspure=sf_output("uspure");

	/*get the axis*/
	sf_axis As,At;
	At=sf_iaxa(vxfile,1);
	As=sf_iaxa(vxfile,2);

	/*get axis parameters from input files*/
	nx=sf_n(As); dx=sf_d(As);
	nt=sf_n(At); dt=sf_d(At);

	float **vxdata,**vzdata,**phdata;
	vxdata=sf_floatalloc2(nt,nx);
	vzdata=sf_floatalloc2(nt,nx);
	phdata=sf_floatalloc2(nt,nx);
    sf_floatread(vxdata[0],nx*nt,vxfile);
    sf_floatread(vzdata[0],nx*nt,vzfile);
    sf_floatread(phdata[0],nx*nt,phfile);

	sf_complex **vxdatafs,**vzdatafs,**phdatafs;
	vxdatafs=sf_complexalloc2(nt,nx);
	vzdatafs=sf_complexalloc2(nt,nx);
	phdatafs=sf_complexalloc2(nt,nx);
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt;j++){
			vxdatafs[i][j]=sf_cmplx(.0,.0);
			vzdatafs[i][j]=sf_cmplx(.0,.0);
			phdatafs[i][j]=sf_cmplx(.0,.0);
		}
	ts2fs(vxdata,vzdata,phdata,vxdatafs,vzdatafs,phdatafs);
/*
	for(int j=0;j<nt;j++)
		sf_warning("j=%d,real=%10f,image=%10f",j,creal(vxdatafs[100][j]),cimag(vxdatafs[100][j]));
	FILE *fp;
	if((fp=fopen("vx.txt","w"))!=NULL){
		for(int j=0;j<nt;j++)
			fprintf(fp,"%f\n",vxdata[100][j]);
		fclose(fp);
	}
*/

	sf_complex **uppurefk=sf_complexalloc2(nt,nx);
	sf_complex **uspurefk=sf_complexalloc2(nt,nx);
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt;j++){
			uppurefk[i][j]=sf_cmplx(.0,.0);
			uspurefk[i][j]=sf_cmplx(.0,.0);
		}
	fs2fk(vxdatafs,vzdatafs,phdatafs,uppurefk,uspurefk);

	float **uppurets,**uspurets;
	uppurets=sf_floatalloc2(nt,nx);
	uspurets=sf_floatalloc2(nt,nx);
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt;j++){
			uppurets[i][j]=0.0;
			uspurets[i][j]=0.0;
		}
	fk2ts(uppurefk,uspurefk,uppurets,uspurets);

	/*write the axis of output files*/
	puthead2d(uppure,nt,dt,nx,dx,0,0);
	puthead2d(uspure,nt,dt,nx,dx,0,0);
	sf_floatwrite(uppurets[0],nt*nx,uppure);
	sf_floatwrite(uspurets[0],nt*nx,uspure);

	sf_close();
	exit(0);
}
void puthead2d(sf_file output,int n1,float d1,int n2,float d2,float o1,float o2)
{
	sf_putint(output,"n1",n1);
	sf_putint(output,"n2",n2);
	sf_putfloat(output,"d1",d1);
	sf_putfloat(output,"d2",d2);
	sf_putfloat(output,"o1",o1);
	sf_putfloat(output,"o2",o2);
}
void ts2fs(float **vxdata,float **vzdata,float **phdata,sf_complex **vxdatafs,sf_complex **vzdatafs,sf_complex **phdatafs)
{
	fftwf_plan vxplan;
	fftwf_plan vzplan;
	fftwf_plan phplan;

	fftwf_complex *vxin,*vxout;
	fftwf_complex *vzin,*vzout;
	fftwf_complex *phin,*phout;

	vxin  = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt);
	vxout = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt);
	vzin  = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt);
	vzout = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt);
	phin  = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt);
	phout = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt);

	vxplan=fftwf_plan_dft_1d(nt, vxin, vxout, FFTW_FORWARD, FFTW_MEASURE);
	vzplan=fftwf_plan_dft_1d(nt, vzin, vzout, FFTW_FORWARD, FFTW_MEASURE);
	phplan=fftwf_plan_dft_1d(nt, phin, phout, FFTW_FORWARD, FFTW_MEASURE);

	for(int index=0;index<nx;index++){
		for(int j=0;j<nt;j++){
			vxin[j][0]=vxdata[index][j];
			vxin[j][1]=0.0;
			vzin[j][0]=vzdata[index][j];
			vzin[j][1]=0.0;
			phin[j][0]=phdata[index][j]/1.5;
			phin[j][1]=0.0;
		}

		fftwf_execute(vxplan);
		fftwf_execute(vzplan);
		fftwf_execute(phplan);

		for(int j=0;j<nt;j++){
			vxdatafs[index][j]=sf_cmplx(vxout[j][0],vxout[j][1]);
			vzdatafs[index][j]=sf_cmplx(vzout[j][0],vzout[j][1]);
			phdatafs[index][j]=sf_cmplx(phout[j][0],phout[j][1]);
		}
	}

    fftwf_destroy_plan(vxplan);
    fftwf_destroy_plan(vzplan);
    fftwf_destroy_plan(phplan);

	fftwf_free(vxin);fftwf_free(vxout);
	fftwf_free(vzin);fftwf_free(vzout);
	fftwf_free(phin);fftwf_free(phout);
}
void fs2fk(sf_complex **vxdatafs,sf_complex **vzdatafs,sf_complex **phdatafs,sf_complex **uppurefk,sf_complex **uspurefk)
{
	/*get velocity and density parameters*/
	float *vp=sf_floatalloc(nx);
	float *vs=sf_floatalloc(nx);
	float *dens=sf_floatalloc(nx);
	for(int i=0;i<nx;i++){
		dens[i]=2.0;
	}

	FILE *fp;
	if((fp=fopen("vpvel.dat","r"))!=NULL){
		for(int i=0;i<nx;i++)
			fread(&vp[i],sizeof(float),1,fp);
		fclose(fp);
	}
	if((fp=fopen("vsvel.dat","r"))!=NULL){
		for(int i=0;i<nx;i++)
			fread(&vs[i],sizeof(float),1,fp);
		fclose(fp);
	}
	//remove("vpvel.dat"); remove("vsvel.dat");

	float kx,w;
	float pp,qp,qs;
	float tmp1,tmp2;

	sf_complex fftc;
	fftc=sf_cmplx(.0,.0);
	for(int j=nt/2;j<nt;j++){
		for(int i=0;i<nx;i++){
			if(i<nx/2)  kx=2*PI*i/nx/dx; else  kx=2*PI*(i-nx)/nx/dx;
			w=2*PI*(j-nt)/nt/dt;
			pp=-kx/w;
			if(absolutevalue(pp)<0.6){
				for(int index=0;index<nx;index++){
					tmp1=powf(1.0/vp[index],2.0)-powf(pp,2.0);
					tmp2=powf(1.0/vs[index],2.0)-powf(pp,2.0);
					if(tmp1>0) qp=sqrt(tmp1); else qp=.0;
					if(tmp2>0) qs=sqrt(tmp2); else qs=.0;
					fftc=sf_cmplx(cos(2*PI*index*i/nx),-sin(2*PI*index*i/nx));
					uppurefk[i][j]+=(vzdatafs[index][j]-qp/(1.0-2.0*pp*pp*vs[index]*vs[index])*
						            (1.5/dens[index]*phdatafs[index][j]+2.0*pp*vs[index]*vs[index]*vxdatafs[index][j]))*fftc;
					uspurefk[i][j]+=(vxdatafs[index][j]-pp/(1.0-2.0*pp*pp*vs[index]*vs[index])*
									(1.5/dens[index]*phdatafs[index][j]-2.0*qs*vs[index]*vs[index]*vzdatafs[index][j]))*fftc;
				}
			}
		}
	}
	for(int j=0;j<nt/2;j++){
		for(int i=0;i<nx;i++){
			if(i<nx/2)  kx=2*PI*i/nx/dx; else  kx=2*PI*(i-nx)/nx/dx;
			w=2*PI*j/nt/dt;
			pp=-kx/w;
			if(absolutevalue(pp)<0.6){
				for(int index=0;index<nx;index++){
					tmp1=powf(1.0/vp[index],2.0)-powf(pp,2.0);
					tmp2=powf(1.0/vs[index],2.0)-powf(pp,2.0);
					if(tmp1>0) qp=sqrt(tmp1); else qp=.0;    //qp is positive or nagetive?
					if(tmp2>0) qs=sqrt(tmp2); else qs=.0;    //qs is positive or nagetive?
					fftc=sf_cmplx(cos(2*PI*index*i/nx),-sin(2*PI*index*i/nx));
					uppurefk[i][j]+=(vzdatafs[index][j]-qp/(1.0-2.0*pp*pp*vs[index]*vs[index])*
						            (1.5/dens[index]*phdatafs[index][j]+2.0*pp*vs[index]*vs[index]*vxdatafs[index][j]))*fftc;
					uspurefk[i][j]+=(vxdatafs[index][j]-pp/(1.0-2.0*pp*pp*vs[index]*vs[index])*
									(1.5/dens[index]*phdatafs[index][j]-2.0*qs*vs[index]*vs[index]*vzdatafs[index][j]))*fftc;
				}
			}
		}
	}
	
/*
	float vp=2.2,vs=0.7;
	float kx,w;
	float pp,qp,qs;
	for(int i=0;i<nx;i++)
		for(int j=nt/2;j<nt;j++){
			if(i<nx/2)  kx=2*PI*i/nx/dx; else  kx=2*PI*(i-nx)/nx/dx;
			w=2*PI*(j-nt)/nt/dt;
			pp=-kx/w;
			tmps1=powf(1.0/vp,2.0)-powf(pp,2.0);
			tmps2=powf(1.0/vs,2.0)-powf(pp,2.0);
			if(tmps1>0) qp=sqrt(tmps1); else qp=.0;
			if(tmps2>0) qs=sqrt(tmps2); else qs=.0;
			if(abs(pp)<0.6){
				uppurefk[i][j]=vzdatafk[i][j]-qp/(1.0-2.0*pp*pp*vs*vs)*(1.5/2*phdatafk[i][j]+2.0*pp*vs*vs*vxdatafk[i][j]);
				uspurefk[i][j]=vxdatafk[i][j]-pp/(1.0-2.0*pp*pp*vs*vs)*(1.5/2*phdatafk[i][j]-2.0*qs*vs*vs*vzdatafk[i][j]);
			}
		}
*/	
}
void fk2ts(sf_complex **uppurefk,sf_complex **uspurefk,float **uppurets,float **uspurets)
{
	fftwf_plan upplan;
	fftwf_plan usplan;
	fftwf_complex *upin=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt*nx);
	fftwf_complex *usin=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt*nx);
	fftwf_complex *upout=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt*nx);
	fftwf_complex *usout=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt*nx);

    upplan=fftwf_plan_dft_2d(nx, nt, upin, upout, FFTW_BACKWARD, FFTW_MEASURE);
    usplan=fftwf_plan_dft_2d(nx, nt, usin, usout, FFTW_BACKWARD, FFTW_MEASURE);

	int index=0;
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt;j++){
			upin[index][0]=creal(uppurefk[i][j]);
			upin[index][1]=cimag(uppurefk[i][j]);
			usin[index][0]=creal(uspurefk[i][j]);
			usin[index][1]=cimag(uspurefk[i][j]);
			index++;
		}
	fftwf_execute(upplan);
	fftwf_execute(usplan);

	index=0;
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt;j++){
			uppurets[i][j]=upout[index][0]/nx/nt;
			uspurets[i][j]=usout[index][0]/nx/nt;
			index++;
		}

	fftwf_destroy_plan(upplan);
	fftwf_destroy_plan(usplan);

	fftwf_free(upin);	fftwf_free(usin);
	fftwf_free(upout);	fftwf_free(usout);
}
float absolutevalue(float in)
{
	if(in<0)
		return -in;
	else
		return in;
}
