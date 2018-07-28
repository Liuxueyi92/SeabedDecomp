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
#include "zero.h"
#include "utils.h"
#define PI  3.141592654

int     nx,nt;
float   dx,dt;
FILE    *fp;

float absolutevalue(float in);
void  puthead2d(sf_file output,int n1,float d1,int n2,float d2,float o1,float o2);
void  ts2fk(float **vxdata,float **vzdata,float **phdata,sf_complex **vxdatafs,sf_complex **vzdatafs,sf_complex **phdatafs);

void  elastic_decomposition(sf_complex **vxdatafs,sf_complex **vzdatafs,sf_complex **phdatafs,sf_complex **uptxzfk,sf_complex **dotxzfk,
	  sf_complex **uptzzfk,	sf_complex **dotzzfk,sf_complex **uppfk,sf_complex **dopfk,sf_complex **upsfk,sf_complex **dosfk);
void  acoustic_decomposition(sf_complex **vxdatafk,sf_complex **vzdatafk,sf_complex **phdatafk,
	  float *pha_cpz,float *amp_cpz,sf_complex **dotxzfk,sf_complex **uptxzfk);
void  acoustic_decomposition_compare(sf_complex **vxdatafk,sf_complex **vzdatafk,sf_complex **phdatafk,sf_complex **dotzzfk,sf_complex **uptzzfk);
void  fk2ts(sf_complex **uppurefk,sf_complex **uspurefk,float **uppurets,float **uspurets);
void  phase_cpz(sf_complex **vzdatafk,sf_complex **phdatafk,float *pha_cpz);
void  amplitude_cpz(sf_complex **vzdatafk,sf_complex **phdatafk,float *pha_cpz,float *amp_cpz);


int main(int argc,char *argv[])
{
	sf_init(argc,argv);
	sf_file vxfile,vzfile,phfile;
	vxfile=sf_input("in");
	vzfile=sf_input("vzrecord0");
	phfile=sf_input("precord0");
	sf_file vxrecord1,vzrecord1,phrecord1;
	vxrecord1=sf_input("vxrecord1");
	vzrecord1=sf_input("vzrecord1");
	phrecord1=sf_input("precord1");
	sf_file vxrecord2,vzrecord2,phrecord2;
	vxrecord2=sf_input("vxrecord2");
	vzrecord2=sf_input("vzrecord2");
	phrecord2=sf_input("precord2");

	sf_file uptxzout,dotxzout,uptzzout,dotzzout;
	uptxzout=sf_output("uptxzout");
	dotxzout=sf_output("dotxzout");
	uptzzout=sf_output("uptzzout");
	dotzzout=sf_output("dotzzout");
	sf_file uppout,dopout,upsout,dosout;
	uppout=sf_output("out");
	dopout=sf_output("dopout");
	upsout=sf_output("upsout");
	dosout=sf_output("dosout");

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

/*   
    //======need to rewrite=======//
	for(int i=0;i<nx;i++)
		for(int j=absolutevalue(i-390)*4+1600;j<nt;j++){
			vxdata[i][j]=0.0;
			vzdata[i][j]=0.0;
			phdata[i][j]=0.0;
		}
*/

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
	ts2fk(vxdata,vzdata,phdata,vxdatafk,vzdatafk,phdatafk);


	float **vxdata1,**vzdata1,**phdata1;
	vxdata1=sf_floatalloc2(nt,nx);
	vzdata1=sf_floatalloc2(nt,nx);
	phdata1=sf_floatalloc2(nt,nx);
	zero_float_2d(vxdata1,nt,nx);
	zero_float_2d(vzdata1,nt,nx);
	zero_float_2d(phdata1,nt,nx);
    sf_floatread(vxdata1[0],nx*nt,vxrecord1);
    sf_floatread(vzdata1[0],nx*nt,vzrecord1);
    sf_floatread(phdata1[0],nx*nt,phrecord1);
	sf_complex **vxdata1fk,**vzdata1fk,**phdata1fk;
	vxdata1fk=sf_complexalloc2(nt,nx);
	vzdata1fk=sf_complexalloc2(nt,nx);
	phdata1fk=sf_complexalloc2(nt,nx);
	zero_complex_2d(vxdata1fk,nt,nx);
	zero_complex_2d(vzdata1fk,nt,nx);
	zero_complex_2d(phdata1fk,nt,nx);
	ts2fk(vxdata1,vzdata1,phdata1,vxdata1fk,vzdata1fk,phdata1fk);

	//============calculate P/Vz calibration=============//
	float *pha_cpz=sf_floatalloc(nt);
	float *amp_cpz=sf_floatalloc(nt);
	zero_float_1d(pha_cpz,nt);
	zero_float_1d(amp_cpz,nt);
	//PZ_calibration
	phase_cpz(vzdatafk,phdatafk,pha_cpz);
	amplitude_cpz(vzdata1fk,phdata1fk,pha_cpz,amp_cpz);

	//======time-domain P/Vz calibration======//
	float *cpz_real=sf_floatalloc(nt);
	float *cpz_imag=sf_floatalloc(nt);
	zero_float_1d(cpz_real,nt);
	zero_float_1d(cpz_imag,nt);
	for(int i=0;i<nt;i++){
		cpz_real[i]=amp_cpz[i]*cos(pha_cpz[i]);
		cpz_imag[i]=amp_cpz[i]*sin(pha_cpz[i]);
	}
	if((fp=fopen("cpz_real.txt","w"))!=NULL){
		for(int i=0;i<nt;i++)
			fprintf(fp,"%f\n",cpz_real[i]*10000);
		fclose(fp);
	}
	if((fp=fopen("cpz_imag.txt","w"))!=NULL){
		for(int i=0;i<nt;i++)
			fprintf(fp,"%f\n",cpz_imag[i]*10000);
		fclose(fp);
	}


	//free variables' storage space
	free_complex_2d(vxdatafk);		free_complex_2d(vzdatafk);
	free_complex_2d(phdatafk);		free_complex_2d(vxdata1fk);
	free_complex_2d(vzdata1fk);		free_complex_2d(phdata1fk);
	free_float_2d(vzdata);			free_float_2d(vxdata);		    free_float_2d(phdata);
	free_float_2d(vzdata1);			free_float_2d(vxdata1);			free_float_2d(phdata1);
	//===================================================//

	float **vxdata2,**vzdata2,**phdata2;
	vxdata2=sf_floatalloc2(nt,nx);
	vzdata2=sf_floatalloc2(nt,nx);
	phdata2=sf_floatalloc2(nt,nx);
	zero_float_2d(vxdata2,nt,nx);
	zero_float_2d(vzdata2,nt,nx);
	zero_float_2d(phdata2,nt,nx);
    sf_floatread(vxdata2[0],nx*nt,vxrecord2);
    sf_floatread(vzdata2[0],nx*nt,vzrecord2);
    sf_floatread(phdata2[0],nx*nt,phrecord2);
	sf_complex **vxdata2fk,**vzdata2fk,**phdata2fk;
	vxdata2fk=sf_complexalloc2(nt,nx);
	vzdata2fk=sf_complexalloc2(nt,nx);
	phdata2fk=sf_complexalloc2(nt,nx);
	zero_complex_2d(vxdata2fk,nt,nx);
	zero_complex_2d(vzdata2fk,nt,nx);
	zero_complex_2d(phdata2fk,nt,nx);
	ts2fk(vxdata2,vzdata2,phdata2,vxdata2fk,vzdata2fk,phdata2fk);

	sf_complex **uptxzfk=sf_complexalloc2(nt,nx);
	sf_complex **dotxzfk=sf_complexalloc2(nt,nx);
	sf_complex **uptzzfk=sf_complexalloc2(nt,nx);
	sf_complex **dotzzfk=sf_complexalloc2(nt,nx);
	sf_complex **uppfk=sf_complexalloc2(nt,nx);
	sf_complex **dopfk=sf_complexalloc2(nt,nx);
	sf_complex **upsfk=sf_complexalloc2(nt,nx);
	sf_complex **dosfk=sf_complexalloc2(nt,nx);
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt;j++){
			uptxzfk[i][j]=sf_cmplx(.0,.0);
			dotxzfk[i][j]=sf_cmplx(.0,.0);
			uptzzfk[i][j]=sf_cmplx(.0,.0);
			dotzzfk[i][j]=sf_cmplx(.0,.0);
			uppfk[i][j]=sf_cmplx(.0,.0);
			dopfk[i][j]=sf_cmplx(.0,.0);
			upsfk[i][j]=sf_cmplx(.0,.0);
			dosfk[i][j]=sf_cmplx(.0,.0);
		}


	//function of acoustic decomposition just above the seabed
	acoustic_decomposition(vxdata2fk,vzdata2fk,phdata2fk,pha_cpz,amp_cpz,dotxzfk,uptxzfk);
	acoustic_decomposition_compare(vxdata2fk,vzdata2fk,phdata2fk,dotzzfk,uptzzfk);
	//function of elastic decomposition just below the seabed
	//elastic_decompositon(vxdata2fk,vzdata2fk,phdata2fk,pha_cpz,amp_cpz,uptxzfk,dotxzfk,uptzzfk,dotzzfk,uppfk,dopfk,upsfk,dosfk);

	float **uptxzts,**dotxzts,**uptzzts,**dotzzts;
	uptxzts=sf_floatalloc2(nt,nx);
	dotxzts=sf_floatalloc2(nt,nx);
	uptzzts=sf_floatalloc2(nt,nx);
	dotzzts=sf_floatalloc2(nt,nx);
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt;j++){
			uptxzts[i][j]=0.0;
			dotxzts[i][j]=0.0;
			uptzzts[i][j]=0.0;
			dotzzts[i][j]=0.0;
		}
	fk2ts(uptxzfk,dotxzfk,uptxzts,dotxzts);
	fk2ts(uptzzfk,dotzzfk,uptzzts,dotzzts);

	float **uppts,**dopts,**upsts,**dosts;
	uppts=sf_floatalloc2(nt,nx);
	dopts=sf_floatalloc2(nt,nx);
	upsts=sf_floatalloc2(nt,nx);
	dosts=sf_floatalloc2(nt,nx);
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt;j++){
			uppts[i][j]=0.0;
			dopts[i][j]=0.0;
			upsts[i][j]=0.0;
			dosts[i][j]=0.0;
		}
	fk2ts(uppfk,dopfk,uppts,dopts);
	fk2ts(upsfk,dosfk,upsts,dosts);

	/*write the axis of output files*/
	puthead2d(uptxzout,nt,dt,nx,dx,0,0);
	puthead2d(dotxzout,nt,dt,nx,dx,0,0);
	puthead2d(uptzzout,nt,dt,nx,dx,0,0);
	puthead2d(dotzzout,nt,dt,nx,dx,0,0);
	sf_floatwrite(uptxzts[0],nt*nx,uptxzout);
	sf_floatwrite(dotxzts[0],nt*nx,dotxzout);
	sf_floatwrite(uptzzts[0],nt*nx,uptzzout);
	sf_floatwrite(dotzzts[0],nt*nx,dotzzout);
	puthead2d(uppout,nt,dt,nx,dx,0,0);
	puthead2d(dopout,nt,dt,nx,dx,0,0);
	puthead2d(upsout,nt,dt,nx,dx,0,0);
	puthead2d(dosout,nt,dt,nx,dx,0,0);
	sf_floatwrite(uppts[0],nt*nx,uppout);
	sf_floatwrite(dopts[0],nt*nx,dopout);
	sf_floatwrite(upsts[0],nt*nx,upsout);
	sf_floatwrite(dosts[0],nt*nx,dosout);

	//free storage space
	free_complex_2d(uptxzfk);	free_complex_2d(dotxzfk);
	free_complex_2d(uptzzfk);	free_complex_2d(dotzzfk);
	free_complex_2d(uppfk);		free_complex_2d(dopfk);
	free_complex_2d(upsfk);		free_complex_2d(dosfk);
	free_float_2d(uptxzts); 	free_float_2d(dotxzts); 
	free_float_2d(uptzzts); 	free_float_2d(dotzzts); 
	free_float_2d(uppts);	 	free_float_2d(dopts); 
	free_float_2d(upsts);		free_float_2d(dosts); 

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
void elastic_decomposition(sf_complex **vxdatafs,sf_complex **vzdatafs,sf_complex **phdatafs,sf_complex **uptxzfk,
	 sf_complex **dotxzfk,sf_complex **uptzzfk,sf_complex **dotzzfk,sf_complex **uppfk,sf_complex **dopfk,sf_complex **upsfk,sf_complex **dosfk)
{
	float vp=2.100,vs=1.400,dens=2.000;
	float kx,w;
	float pp,qp,qs;
	float tmp1,tmp2,tmp3;


	for(int j=nt/2;j<nt;j++){
		for(int i=0;i<nx;i++){
			if(i<nx/2)  kx=2*PI*i/nx/dx; else  kx=2*PI*(i-nx)/nx/dx;
			w=2*PI*(j-nt)/nt/dt;
			pp=-kx/w;
			if(absolutevalue(pp)<0.47){
				tmp1=powf(1.0/vp,2.0)-powf(pp,2.0);
				if(tmp1>0) qp=sqrt(tmp1); else qp=.0;
				tmp2=powf(1.0/vs,2.0)-powf(pp,2.0);
				if(tmp2>0) qs=sqrt(tmp2); else qs=.0;
				tmp3=4.0*vs*vs*pp*pp*qp*qs+powf(vs*(1.0/(vs*vs)-pp*pp),2.0);
				{
					uptxzfk[i][j]=0.5*(-phdatafs[i][j]*pp*vs*vs*(2*qs*qp-1.0/(vs*vs)+2*pp*pp)/qs-vxdatafs[i][j]*dens*vs*vs*tmp3/qs);
					dotxzfk[i][j]=0.5*(+phdatafs[i][j]*pp*vs*vs*(2*qs*qp-1.0/(vs*vs)+2*pp*pp)/qs+vxdatafs[i][j]*dens*vs*vs*tmp3/qs);
					uptzzfk[i][j]=0.5*(+phdatafs[i][j]-vzdatafs[i][j]*dens*vs*vs*tmp3/qp);
					dotzzfk[i][j]=0.5*(+phdatafs[i][j]+vzdatafs[i][j]*dens*vs*vs*tmp3/qp);
				}
				if(tmp3!=0.0)
				{
					uppfk[i][j]=(-2*pp*qs*uptxzfk[i][j]+(1.0/(vs*vs)-2*pp*pp)*uptzzfk[i][j])/tmp3;
					dopfk[i][j]=(+2*pp*qs*dotxzfk[i][j]+(1.0/(vs*vs)-2*pp*pp)*dotzzfk[i][j])/tmp3;
					upsfk[i][j]=((2*pp*pp-1.0/(vs*vs))*uptxzfk[i][j]-2*pp*qp*uptzzfk[i][j])/tmp3;
					dosfk[i][j]=((2*pp*pp-1.0/(vs*vs))*dotxzfk[i][j]+2*pp*qp*dotzzfk[i][j])/tmp3;
				}
			}
		}
	}
}
void acoustic_decomposition(sf_complex **vxdatafk,sf_complex **vzdatafk,sf_complex **phdatafk,
	 float *pha_cpz,float *amp_cpz,sf_complex **dotxzfk,sf_complex **uptxzfk)
{
	float vpw=1.5,densw=1.0;
	float kx,w;
	float pp,qp;
	float tmp1;
	sf_complex tmp2;

	for(int j=nt/2;j<nt;j++){
		for(int i=0;i<nx;i++){
			if(i<nx/2)  kx=2*PI*i/nx/dx; else  kx=2*PI*(i-nx)/nx/dx;
			w=2*PI*(j-nt)/nt/dt;
			pp=-kx/w;
			tmp2=sf_cmplx(cos(pha_cpz[j]),sin(pha_cpz[j]));
			if(absolutevalue(pp)<0.45){
				tmp1=powf(1.0/vpw,2.0)-powf(pp,2.0);
				if(tmp1>0) qp=sqrt(tmp1); else qp=.0;
				{
					uptxzfk[i][j]=(phdatafk[i][j]-amp_cpz[j]*tmp2*vzdatafk[i][j]*densw/qp);
					dotxzfk[i][j]=(phdatafk[i][j]+amp_cpz[j]*tmp2*vzdatafk[i][j]*densw/qp);
				}
			}
		}
	}
}
void acoustic_decomposition_compare(sf_complex **vxdatafk,sf_complex **vzdatafk,sf_complex **phdatafk,sf_complex **dotzzfk,sf_complex **uptzzfk)
{
	float vpw=1.5,densw=1.0;
	float kx,w;
	float pp,qp;
	float tmp1;

	for(int j=nt/2;j<nt;j++){
		for(int i=0;i<nx;i++){
			if(i<nx/2)  kx=2*PI*i/nx/dx; else  kx=2*PI*(i-nx)/nx/dx;
			w=2*PI*(j-nt)/nt/dt;
			pp=-kx/w;
			if(absolutevalue(pp)<0.45){
				tmp1=powf(1.0/vpw,2.0)-powf(pp,2.0);
				if(tmp1>0) qp=sqrt(tmp1); else qp=.0;
				{
					uptzzfk[i][j]=(phdatafk[i][j]-vzdatafk[i][j]*densw/qp);
					dotzzfk[i][j]=(phdatafk[i][j]+vzdatafk[i][j]*densw/qp);
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
void  phase_cpz(sf_complex **vzdatafk,sf_complex **phdatafk,float *pha_cpz)
{
	utils_print_title("calculate the phase spectrum");
	float **vzphase,**phphase;
	vzphase=sf_floatalloc2(nt,nx);
	phphase=sf_floatalloc2(nt,nx);
	zero_float_2d(vzphase,nt,nx);
	zero_float_2d(phphase,nt,nx);
	
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt;j++){
			vzphase[i][j]=atan(cimag(vzdatafk[i][j])/creal(vzdatafk[i][j]))/PI*180;
			phphase[i][j]=atan(cimag(phdatafk[i][j])/creal(phdatafk[i][j]))/PI*180;
		}
	for(int i=0;i<nt;i++){
		pha_cpz[i]=10.0;
	}
	float Epz_phase=0.0;
	double ee=0.0000000000000001;
	float *gradient_phase=sf_floatalloc(nt);
	zero_float_1d(gradient_phase,nt);
	float step=0.001; //step of one iteration (units is degree)

	for(int i=0;i<nx;i++)
		for(int j=nt/2;j<nt;j++){
			Epz_phase+=powf(phphase[i][j]-vzphase[i][j]-pha_cpz[j],2.0);
		}
	sf_warning("Epz_phase_initial=%f",Epz_phase);

	int index=0;
    int numiteration=1000;
	float kx,w;
	float pp;
	while(Epz_phase>=ee&&index<numiteration){
		Epz_phase=0.0;
		zero_float_1d(gradient_phase,nt);
		for(int iw=nt/2;iw<nt;iw++)
			for(int ik=0;ik<nx;ik++){
				if(ik<nx/2)  kx=2*PI*ik/nx/dx; else  kx=2*PI*(ik-nx)/nx/dx;
				w=2*PI*(iw-nt)/nt/dt;
				pp=-kx/w;
				if(absolutevalue(pp)<0.55){
					gradient_phase[iw]+=-2.0*(phphase[ik][iw]-vzphase[ik][iw]-pha_cpz[iw]);
				}
		}

//		if(index==0)
//		for(int iw=nt/2;iw<nt;iw++)
//			sf_warning("gradient=%f",gradient_phase[iw]);

		for(int iw=nt/2;iw<nt;iw++){
			pha_cpz[iw]=pha_cpz[iw]-step*gradient_phase[iw];
		}
		for(int i=0;i<nx;i++)
			for(int j=nt/2;j<nt;j++){
				if(i<nx/2)  kx=2*PI*i/nx/dx; else  kx=2*PI*(i-nx)/nx/dx;
				w=2*PI*(j-nt)/nt/dt;
				pp=-kx/w;
				if(absolutevalue(pp)<0.55){
					Epz_phase+=powf(phphase[i][j]-vzphase[i][j]-pha_cpz[j],2.0);
				}
			}
		index++;
		utils_loadbar(index,numiteration,1000,50);
	}
//	for(int iw=nt/2;iw<nt;iw++){
//		if(iw%4==0)
//			sf_warning("pha_cpz=%f",pha_cpz[iw]);
//	}

	sf_warning("phase_iteration_time=%d",index);
	sf_warning("Epz_phase_end=%f",Epz_phase);
}
void  amplitude_cpz(sf_complex **vzdatafk,sf_complex **phdatafk,float *pha_cpz,float *amp_cpz)
{
	for(int i=0;i<nt;i++){
		amp_cpz[i]=0.3;
	}
	utils_print_title("calculate the amplitude spectrum");
	for(int i=0;i<nt;i++)
		pha_cpz[i]=pha_cpz[i]/180*PI;
	sf_complex *phase=sf_complexalloc(nt);
	for(int i=0;i<nt;i++){
		phase[i]=sf_cmplx(cos(pha_cpz[i]),sin(pha_cpz[i]));
	}

	double Epz_amplitude=0.0;
	float vp=1.5,dens=1.0;
	float kx,w;
	float qp,pp,tmp1;
	sf_complex tmp2,tmp3,tmp4,tmp5;
	for(int j=nt/2;j<nt;j++)
		for(int i=0;i<nx;i++){
			if(i<nx/2)  kx=2*PI*i/nx/dx; else  kx=2*PI*(i-nx)/nx/dx;
			w=2*PI*(j-nt)/nt/dt;
			pp=-kx/w;
			tmp1=powf(1.0/vp,2.0)-powf(pp,2.0);
			if(tmp1>0) qp=sqrt(tmp1); else qp=.0;
			if(absolutevalue(pp)<0.55){
				tmp2=(-qp*phdatafk[i][j]+amp_cpz[j]*phase[j]*dens*vzdatafk[i][j]);
				tmp3=sf_cmplx(creal(tmp2),-cimag(tmp2));
				tmp4=(qp*phdatafk[i][j]+amp_cpz[j]*phase[j]*dens*vzdatafk[i][j]);
				tmp5=sf_cmplx(creal(tmp4),-cimag(tmp4));
				Epz_amplitude+=powf(tmp2*tmp3-tmp4*tmp5,2.0);
			}
		}
	sf_warning("Epz_amplitude_initial=%f",Epz_amplitude);

	float ee=0.0000005;
	int index=0,numiteration=500;
	float *gradient=sf_floatalloc(nt);
	float step=0.0000000000006;
	sf_complex tmp6,tmp7,tmp8,tmp9,tmp10;
	while(Epz_amplitude>ee&&index<numiteration){
		Epz_amplitude=0.0;
		zero_float_1d(gradient,nt);
		for(int j=nt/2;j<nt;j++)
			for(int i=0;i<nx;i++){
				if(i<nx/2)  kx=2*PI*i/nx/dx; else  kx=2*PI*(i-nx)/nx/dx;
				w=2*PI*(j-nt)/nt/dt;
				pp=-kx/w;
				tmp1=powf(1.0/vp,2.0)-powf(pp,2.0);
				if(tmp1>0) qp=sqrt(tmp1); else qp=.0;
				if(absolutevalue(pp)<0.55){
					tmp2=(-qp*phdatafk[i][j]+amp_cpz[j]*phase[j]*dens*vzdatafk[i][j]);
					tmp3=sf_cmplx(creal(tmp2),-cimag(tmp2));
					tmp4=(qp*phdatafk[i][j]+amp_cpz[j]*phase[j]*dens*vzdatafk[i][j]);
					tmp5=sf_cmplx(creal(tmp4),-cimag(tmp4));
					tmp10=phase[j]*dens*vzdatafk[i][j];
					tmp6=tmp10*tmp3;
					tmp7=tmp2*sf_cmplx(creal(tmp10),-cimag(tmp10));
					tmp8=tmp10*tmp5;
					tmp9=tmp4*sf_cmplx(creal(tmp10),-cimag(tmp10));
					gradient[j]+=creal(2*(tmp2*tmp3-tmp4*tmp5)*(tmp6+tmp7-tmp8-tmp9));
				}
				if(i==300&&index==0&&j%10==0)
					sf_warning("amp_gradient=%f",gradient[j]);
			}
		for(int j=nt/2;j<nt;j++){
			amp_cpz[j]=amp_cpz[j]-step*gradient[j];
		}

		for(int j=nt/2;j<nt;j++)
			for(int i=0;i<nx;i++){
				if(i<nx/2)  kx=2*PI*i/nx/dx; else  kx=2*PI*(i-nx)/nx/dx;
				w=2*PI*(j-nt)/nt/dt;
				pp=-kx/w;
				tmp1=powf(1.0/vp,2.0)-powf(pp,2.0);
				if(tmp1>0) qp=sqrt(tmp1); else qp=.0;
				if(absolutevalue(pp)<0.55){
					tmp2=(-qp*phdatafk[i][j]+amp_cpz[j]*phase[j]*dens*vzdatafk[i][j]);
					tmp3=sf_cmplx(creal(tmp2),-cimag(tmp2));
					tmp4=(qp*phdatafk[i][j]+amp_cpz[j]*phase[j]*dens*vzdatafk[i][j]);
					tmp5=sf_cmplx(creal(tmp4),-cimag(tmp4));
					Epz_amplitude+=powf(tmp2*tmp3-tmp4*tmp5,2.0);
				}
			}

		index++;
		utils_loadbar(index,numiteration,1000,50);
	}

	//output for checking the result
//	for(int i=nt-400;i<nt;i++)
//		sf_warning("amp_cpz=%f",amp_cpz[i]);
	sf_warning("amp_iteration_time=%d",index);
	sf_warning("Epz_amplitude_end=%f",Epz_amplitude);
}
