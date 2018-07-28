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
#include "free.h"
#define PI  3.141592654

int     nx,nt;
float   dx,dt;
FILE    *fp;

float absolutevalue(float in);
void  puthead2d(sf_file output,int n1,float d1,int n2,float d2,float o1,float o2);
void  ts2fk(float **vxdata,float **vzdata,float **phdata,sf_complex **vxdatafs,sf_complex **vzdatafs,sf_complex **phdatafs);
void  fk2fk(sf_complex **vxdatafs,sf_complex **vzdatafs,sf_complex **phdatafs,sf_complex **uppurefk,sf_complex **usppurefk);
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
	
	//nt=1800;

	float **vxdatar,**vzdatar,**phdatar;
	vxdatar=sf_floatalloc2(nt,nx);
	vzdatar=sf_floatalloc2(nt,nx);
	phdatar=sf_floatalloc2(nt,nx);
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt-3;j++){
			vxdatar[i][j]=vxdata[i][j];
			vzdatar[i][j]=vzdata[i][j]; //caution:calibration in terms of time
			phdatar[i][j]=phdata[i][j];
		}


	sf_complex **vxdatafk,**vzdatafk,**phdatafk;
	vxdatafk=sf_complexalloc2(nt,nx);
	vzdatafk=sf_complexalloc2(nt,nx);
	phdatafk=sf_complexalloc2(nt,nx);
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt;j++){
			vxdatafk[i][j]=sf_cmplx(.0,.0);
			vzdatafk[i][j]=sf_cmplx(.0,.0);
			phdatafk[i][j]=sf_cmplx(.0,.0);
		}
	ts2fk(vxdatar,vzdatar,phdatar,vxdatafk,vzdatafk,phdatafk);

	sf_complex **uppurefk=sf_complexalloc2(nt,nx);
	sf_complex **uspurefk=sf_complexalloc2(nt,nx);
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt;j++){
			uppurefk[i][j]=sf_cmplx(.0,.0);
			uspurefk[i][j]=sf_cmplx(.0,.0);
		}
	fk2fk(vxdatafk,vzdatafk,phdatafk,uppurefk,uspurefk);

	float **uppurets,**uspurets;
	uppurets=sf_floatalloc2(nt,nx);
	uspurets=sf_floatalloc2(nt,nx);
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt;j++){
			uppurets[i][j]=0.0;
			uspurets[i][j]=0.0;
		}
	fk2ts(uppurefk,uspurefk,uppurets,uspurets);

	//write files of one channel
	if((fp=fopen("pup.txt","w"))!=NULL){
		for(int i=0;i<nx;i++)
			for(int j=0;j<nt;j++)
				fprintf(fp,"%f ",uppurets[i][j]);
		fclose(fp);
	}
	if((fp=fopen("pdown.txt","w"))!=NULL){
		for(int i=0;i<nx;i++)
			for(int j=0;j<nt;j++)
				fprintf(fp,"%f ",uspurets[i][j]);
		fclose(fp);
	}
	if((fp=fopen("precord.txt","w"))!=NULL){
		for(int i=0;i<nx;i++)
			for(int j=0;j<nt;j++)
				fprintf(fp,"%f ",phdata[i][j]);
		fclose(fp);
	}
	if((fp=fopen("vzrecord.txt","w"))!=NULL){
		for(int i=0;i<nx;i++)
			for(int j=0;j<nt;j++)
				fprintf(fp,"%f ",vzdata[i][j]);
		fclose(fp);
	}


	/*write the axis of output files*/
	puthead2d(uppure,nt,dt,nx,dx,0,0);
	puthead2d(uspure,nt,dt,nx,dx,0,0);
	sf_floatwrite(uppurets[0],nt*nx,uppure);
	sf_floatwrite(uspurets[0],nt*nx,uspure);

	//==============free_storage_space============//
	free_float_2d(vxdata); free_float_2d(vzdata); free_float_2d(phdata);
	free_float_2d(vxdatar); free_float_2d(vzdatar); free_float_2d(phdatar);
	free_complex_2d(vxdatafk); free_complex_2d(vzdatafk); free_complex_2d(phdatafk);
	free_float_2d(uppurets); free_float_2d(uspurets);
	free_complex_2d(uppurefk); free_complex_2d(uspurefk);

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
void ts2fk(float **vxdata,float **vzdata,float **phdata,sf_complex **vxdatafs,sf_complex **vzdatafs,sf_complex **phdatafs)
{
	fftwf_plan vxplan;
	fftwf_plan vzplan;
	fftwf_plan phplan;

	fftwf_complex *vxin,*vxout;
	fftwf_complex *vzin,*vzout;
	fftwf_complex *phin,*phout;

	vxin  = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt*nx);
	vxout = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt*nx);
	vzin  = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt*nx);
	vzout = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt*nx);
	phin  = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt*nx);
	phout = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt*nx);

    vxplan=fftwf_plan_dft_2d(nx, nt, vxin, vxout, FFTW_FORWARD, FFTW_MEASURE);
    vzplan=fftwf_plan_dft_2d(nx, nt, vzin, vzout, FFTW_FORWARD, FFTW_MEASURE);
    phplan=fftwf_plan_dft_2d(nx, nt, phin, phout, FFTW_FORWARD, FFTW_MEASURE);

	int index=0;
	for(int i=0;i<nx;i++){
		for(int j=0;j<nt;j++){
			vxin[index][0]=vxdata[i][j];
			vxin[index][1]=0.0;
			vzin[index][0]=vzdata[i][j];
			vzin[index][1]=0.0;
			phin[index][0]=phdata[i][j];
			phin[index][1]=0.0;

			index++;
		}
	}
	fftwf_execute(vxplan);
	fftwf_execute(vzplan);
	fftwf_execute(phplan);

	index=0;
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt;j++){
			vxdatafs[i][j]=sf_cmplx(vxout[index][0],vxout[index][1]);
			vzdatafs[i][j]=sf_cmplx(vzout[index][0],vzout[index][1]);
			phdatafs[i][j]=sf_cmplx(phout[index][0],phout[index][1]);

			index++;
		}

    fftwf_destroy_plan(vxplan);
    fftwf_destroy_plan(vzplan);
    fftwf_destroy_plan(phplan);

	fftwf_free(vxin);fftwf_free(vxout);
	fftwf_free(vzin);fftwf_free(vzout);
	fftwf_free(phin);fftwf_free(phout);
}
void fk2fk(sf_complex **vxdatafs,sf_complex **vzdatafs,sf_complex **phdatafs,sf_complex **uppurefk,sf_complex **uspurefk)
{
	float vp=1.500,vs=0.00,dens=1.000;
	float kx,w;
	float pp,qp;
	float tmp1;

	for(int j=nt/2;j<nt;j++){
		for(int i=0;i<nx;i++){
			if(i<nx/2)  kx=2*PI*i/nx/dx; else  kx=2*PI*(i-nx)/nx/dx;
			w=2*PI*(j-nt)/nt/dt;
			pp=-kx/w;
			if(absolutevalue(pp)<0.66){
				tmp1=powf(1.0/1.500,2.0)-powf(pp,2.0);
				if(tmp1>0) qp=sqrt(tmp1); else qp=.0;
				//if (qp>0.0)
				{
					uppurefk[i][j]=(phdatafs[i][j]-vzdatafs[i][j]/qp);
					uspurefk[i][j]=(phdatafs[i][j]+vzdatafs[i][j]/qp);
				}
			}
		}
	}
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
