/*
 *Description: define functions -- initialize variables
 *Author: xueyiliu, xueyilie@tongji.edu.cn
 *Data: 31th, May, 2018
*/

#include <rsf.h>

void zero_int_1d(int *p,int nx)
/*< zero 1d int variables  >*/
{
	for(int i=0;i<nx;i++)
	 	p[i]=0;
}

void zero_float_1d(float *p,int nx)
/*< zero 1d float variables  >*/
{
	for(int i=0;i<nx;i++)
	 	p[i]=0.0;
}

void zero_double_1d(double *p,int nx)
/*< zero 1d double variables  >*/
{
	for(int i=0;i<nx;i++)
	 	p[i]=0.0;
}

void zero_complex_1d(sf_complex *p,int nx)
/*< zero 1d double variables  >*/
{
	for(int i=0;i<nx;i++)
	 	p[i]=sf_cmplx(.0,.0);
}

void zero_int_2d(float **p,int nz,int nx)
/*< zero 2d int variables  >*/
{
	for(int i=0;i<nx;i++)
		for(int j=0;j<nz;j++){
			p[i][j]=0;
		}
}

void zero_float_2d(float **p,int nz,int nx)
/*< zero 2d float variables  >*/
{
	for(int i=0;i<nx;i++)
		for(int j=0;j<nz;j++){
			p[i][j]=0.0;
		}
}

void zero_complex_2d(sf_complex **p,int nz,int nx)
/*< zero 2d float variables  >*/
{
	for(int i=0;i<nx;i++)
		for(int j=0;j<nz;j++){
			p[i][j]=sf_cmplx(.0,.0);
		}
}

void zero_double_2d(double **p,int nz,int nx)
/*< zero 2d double variables  >*/
{
	for(int i=0;i<nx;i++)
		for(int j=0;j<nz;j++){
			p[i][j]=0.0;
		}
}

void zero_int_3d(int ***p,int nz,int ny,int nx)
/*< zero 2d int variables  >*/
{
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			for(int k=0;k<nz;k++){
				p[i][j][k]=0;
			}
}

void zero_float_3d(float ***p,int nz,int ny,int nx)
/*< zero 2d float variables  >*/
{
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			for(int k=0;k<nz;k++){
				p[i][j][k]=0.0;
			}
}

void zero_double_3d(double ***p,int nz,int ny,int nx)
/*< zero 2d double variables  >*/
{
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			for(int k=0;k<nz;k++){
				p[i][j][k]=0.0;
			}
}


