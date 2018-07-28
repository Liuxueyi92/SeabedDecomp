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
FILE    *fp;

float absolutevalue(float in);
void  puthead2d(sf_file output,int n1,float d1,int n2,float d2,float o1,float o2);
void  ts2fk(float **vxdata,float **vzdata,float **phdata,sf_complex **vxdatafs,sf_complex **vzdatafs,sf_complex **phdatafs);
void  fk2fk(sf_complex **vxdatafs,sf_complex **vzdatafs,sf_complex **phdatafs,sf_complex **uptxzfk,sf_complex **dotxzfk,sf_complex **uptzzfk,
		    sf_complex **dotzzfk,sf_complex **uppfk,sf_complex **dopfk,sf_complex **upsfk,sf_complex **dosfk);
void  fk2ts(sf_complex **uppurefk,sf_complex **uspurefk,float **uppurets,float **uspurets);

int main(int argc,char *argv[])
{
	sf_init(argc,argv);
	sf_file vxfile,vzfile,phfile;
	vxfile=sf_input("in");
	vzfile=sf_input("vzrecord");
	phfile=sf_input("precord");

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
	
	//nt=1800;

	float **vxdatar,**vzdatar,**phdatar;
	vxdatar=sf_floatalloc2(nt,nx);
	vzdatar=sf_floatalloc2(nt,nx);
	phdatar=sf_floatalloc2(nt,nx);
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt;j++){
			vxdatar[i][j]=vxdata[i][j];
			vzdatar[i][j]=vzdata[i][j];
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

/*
	//output---fkdomain data!!
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt;j++){
			vxdata[i][j]=creal(vxdatafk[i][j]);
			vzdata[i][j]=cimag(vxdatafk[i][j]);
		}
	if((fp=fopen("vxreal.bin","wb"))!=NULL){
		fwrite(vxdata[0],sizeof(float),nx*nt,fp);
		fclose(fp);
	}
	if((fp=fopen("vximag.bin","wb"))!=NULL){
		fwrite(vzdata[0],sizeof(float),nx*nt,fp);
		fclose(fp);
	}
*/

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
	fk2fk(vxdatafk,vzdatafk,phdatafk,uptxzfk,dotxzfk,uptzzfk,dotzzfk,uppfk,dopfk,upsfk,dosfk);

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
void fk2fk(sf_complex **vxdatafs,sf_complex **vzdatafs,sf_complex **phdatafs,sf_complex **uptxzfk,sf_complex **dotxzfk,sf_complex **uptzzfk,
	       sf_complex **dotzzfk,sf_complex **uppfk,sf_complex **dopfk,sf_complex **upsfk,sf_complex **dosfk)
{
/*
	//get velocity and density parameters
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
*/
	float vpw=1.5,densw=1.0;
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
