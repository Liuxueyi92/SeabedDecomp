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
void  pzcomb(float **uppurets,float **uspurets,float **vzdata,float **phdata)
{
	float densw=1.0,vpw=1.5;
	float dens=2.0,vp=2.1,vs=1.4;
	float R=(dens*vp-densw*vpw)-(dens*vp+densw*vpw);
	sf_warning("R=%f",R);
	for(int i=0;i<nx;i++) 
		for(int j=0;j<nt;j++){
			uppurets[i][j]=(1+R)/(1-R)*vzdata[i][j]+phdata[i][j]/(densw*vpw);
			uspurets[i][j]=phdata[i][j]/(densw*vpw)-(1+R)/(1-R)*vzdata[i][j];
		}
}

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
	

	float **uppurets,**uspurets;
	uppurets=sf_floatalloc2(nt,nx);
	uspurets=sf_floatalloc2(nt,nx);
	for(int i=0;i<nx;i++)
		for(int j=0;j<nt;j++){
			uppurets[i][j]=0.0;
			uspurets[i][j]=0.0;
		}
	//fk2ts(uppurefk,uspurefk,uppurets,uspurets);
	pzcomb(uppurets,uspurets,vzdata,phdata);

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
