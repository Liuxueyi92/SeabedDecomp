/* Utilities library for modeling and calculations.
*  Author: Wiktor W. Weibull, wiktor.weibull@ntnu.no, NTNU, 2013
*/

#include <math.h>
#include <rsf.h>
#include <sys/types.h>
#include <sys/stat.h>

float utils_find3dmax(float *array, const int nx, const int ny, const int nz) 
/*<Finding the maximum number in a 3D array>*/
{
	int ix;
	int iy;
	int iz;
	float max = array[0];
	for(iy=0; iy<ny; iy++) {
		for(ix=0; ix<nx; ix++) {
			for(iz=0; iz<nz; iz++) {
				if(array[iy*nx*nz + ix*nz + iz] >= max){
					max = array[iy*nx*nz + ix*nz + iz];
				}
			}
		}
	}
	return max;
}

void utils_copy3darray(float *to,  float *from, const int nx, const int ny, const int nz) 
/*<Array x float multiplication in  3D array>*/
{
	int ix;
	int iy;
	int iz;
	for(iy=0; iy<ny; iy++) {
		for(ix=0; ix<nx; ix++) {
			for(iz=0; iz<nz; iz++) {
				to[iy*nx*nz + ix*nz + iz]=from[iy*nx*nz + ix*nz + iz];
			}
		}
	}
}

void utils_scale3darray(float *array, const float scale, const int nx, const int ny, const int nz) 
/*<Array x float multiplication in  3D array>*/
{
	int ix;
	int iy;
	int iz;
	for(iy=0; iy<ny; iy++) {
		for(ix=0; ix<nx; ix++) {
			for(iz=0; iz<nz; iz++) {
				array[iy*nx*nz + ix*nz + iz]=array[iy*nx*nz + ix*nz + iz]*scale;
			}
		}
	}
}

void utils_array_array_mult3d(float *array1, float *array2, const int nx, const int ny, const int nz) 
/*<Array x array multiplication in 3D>*/
{
	int ix;
	int iy;
	int iz;
	for(iy=0; iy<ny; iy++) {
		for(ix=0; ix<nx; ix++) {
			for(iz=0; iz<nz; iz++) {
				array1[iy*nx*nz + ix*nz + iz] *= array2[iy*nx*nz + ix*nz + iz];
			}
		}
	}
}

void utils_array_array_add3d(float *array1, float *array2, const int nx, const int ny, const int nz) 
/*<Array x array addition in 3D>*/
{
	int ix;
	int iy;
	int iz;
	for(iy=0; iy<ny; iy++) {
		for(ix=0; ix<nx; ix++) {
			for(iz=0; iz<nz; iz++) {
				array1[iy*nx*nz + ix*nz + iz] += array2[iy*nx*nz + ix*nz + iz];
			}
		}
	}
}

void utils_array_array_sub3d(float *array1, float *array2, const int nx, const int ny, const int nz) 
/*<Array x array subtruction in 3D>*/
{
	int ix;
	int iy;
	int iz;
	for(iy=0; iy<ny; iy++) {
		for(ix=0; ix<nx; ix++) {
			for(iz=0; iz<nz; iz++) {
				array1[iy*nx*nz + ix*nz + iz] -= array2[iy*nx*nz + ix*nz + iz];
			}
		}
	}
}

void utils_padmodel1d(float *model, float *padded, const int n1, const int pad1)
/*<Pad input 1d model>*/
{
	int ix;

	for(ix=pad1;ix<n1-pad1; ix++){
		padded[ix]=model[ix-pad1];	
	}

	for(ix=0; ix<pad1; ix++){
		padded[ix]=padded[pad1];
		padded[n1-pad1+ix]=padded[n1-pad1-1];
	}

}

void utils_padmodel2d(float **model, float **padded, const int n1, const int n2, const int pad1, const int pad2)
/*<Pad input 2d model>*/
{
	int ix,iy;

	for(ix=pad1;ix<n1-pad1; ix++){
		for(iy=pad2;iy<n2-pad2; iy++){
			padded[ix][iy]=model[ix-pad1][iy-pad2];	
		}
	}

	for(ix=0; ix<pad1; ix++){
		for(iy=0; iy<n2; iy++){
			padded[ix][iy]=padded[pad1][iy];
			padded[n1-pad1+ix][iy]=padded[n1-pad1-1][iy];
		}
	}
	for(ix=0; ix<n1; ix++){
		for(iy=0; iy<pad2; iy++){
			padded[ix][iy]=padded[ix][pad2];
			padded[ix][n2-pad2+iy]=padded[ix][n2-pad2-1];
		}
	}

}

void utils_padmodel3d(float ***model, float ***padded, const int n1, const int n2, const int n3, const int pad1, const int pad2,const int pad3) 
/*<Pads input 3d model>*/
{
	int ix,iy,iz;

	for(ix=pad1;ix<n1-pad1; ix++){
		for(iy=pad2;iy<n2-pad2; iy++){
			for(iz=pad3;iz<n3-pad3; iz++){
				padded[ix][iy][iz]=model[ix-pad1][iy-pad2][iz-pad3];	
			}
		}
	}

	for(ix=0; ix<pad1; ix++){
		for(iy=0; iy<n2; iy++){
			for(iz=0; iz<n3; iz++){
				padded[ix][iy][iz]=padded[pad1][iy][iz];
				padded[n1-pad1+ix][iy][iz]=padded[n1-pad1-1][iy][iz];
			}
		}
	}
	for(ix=0; ix<n1; ix++){
		for(iy=0; iy<pad2; iy++){
			for(iz=0; iz<n3; iz++){
				padded[ix][iy][iz]=padded[ix][pad2][iz];
				padded[ix][n2-pad2+iy][iz]=padded[ix][n2-pad2-1][iz];
			}
		}
	}
	for(ix=0; ix<n1; ix++){
		for(iy=0; iy<n2; iy++){
			for(iz=0; iz<pad3; iz++){
				padded[ix][iy][iz]=padded[ix][iy][pad3];
				padded[ix][iy][n3-pad3+iz]=padded[ix][iy][n3-pad3-1];
			}
		}
	}

}

void utils_add_central2(float **padded, float **center, const int nx, const int nz, const int padx, const int padz)
/*<add value to the central of a padded model>*/
{
	int ix,iy;

	for(ix=padx;ix<nx-padx; ix++){
		for(iy=padz;iy<nz-padz; iy++){
			padded[ix][iy] += center[ix-padx][iy-padz];	
		}
	}
}

void utils_sub_central2(float **padded, float **center, const int nx, const int nz, const int padx, const int padz)
/*<sub value to the central of a padded model>*/
{
	int ix,iy;

	for(ix=padx;ix<nx-padx; ix++){
		for(iy=padz;iy<nz-padz; iy++){
			padded[ix][iy] -= center[ix-padx][iy-padz];	
		}
	}
}

void utils_loadbar(int x, int n, int r, int w)
/*<Prints progress bar on terminal>*/
{
	if ( r > n) r = n;
	// Only update r times.
	if ( x % (n/r) != 0 ) return;

	// Calculuate the ratio of complete-to-incomplete.
	float ratio = x/(float)n;
	int   c     = ratio * w;

	// Show the percentage complete.
	fprintf(stderr,"%3.2f%% [", (ratio*100) );

	// Show the load bar.
	for (x=0; x<c; x++)
		fprintf(stderr,"=");

	for (x=c; x<w; x++)
		fprintf(stderr," ");

	fprintf(stderr,"]\r"); // Move to the first column
	fflush(stderr);
}

void utils_print_title(const char *title)
/*<Printing out program title to terminal.>*/
{
	int i;
	int len;
	len=strlen(title);
	len=len+16;
	fprintf(stderr,"\n");
	for (i=0; i<len+2; i++) fprintf(stderr, "*");
	fprintf(stderr,"\n");
	fprintf(stderr,"***");
	for (i=0; i<len-4; i++) fprintf(stderr, " ");
	fprintf(stderr,"***\n");
	fprintf(stderr,"***      %s      ***\n",title);
	fprintf(stderr,"***");
	for (i=0; i<len-4; i++) fprintf(stderr, " ");
	fprintf(stderr,"***\n");
	for (i=0; i<len+2; i++) fprintf(stderr, "*");
	fprintf(stderr,"\n");
}

float norm(const float *v, const int n, const int type)
/*<Calculates the norm of vector v. The norm type is specified as the last input. NB: type=99999 is the infinity/max norm.>*/
{	
	double norm = 0.0;
	int i;

	// Calculation norm
	switch(type) {
		case 99999:
			norm = fabs(v[0]);
			for(i=0; i<n; i++) {
				if(norm < fabs(v[i]))
					norm = fabs(v[i]);
			}
			break;
		case 1:
			for(i=0; i<n; i++) {
				norm += fabs(v[i]);
			}
			break;
		case 2:
		default:
			for(i=0; i<n; i++) {
				norm += v[i]*v[i];
			}
			norm = sqrtf(norm);
			break;
	}

	return (float) norm;
}

float max(const float *v, const int n)
/*<Finding the maximum of a vector v with size n.>*/
{
	float val = v[0];
	int i;

	for(i=0; i<n; i++) {
		if(v[i]>val) val=v[i];
	}
	return val;
}

float min(const float *v, const int n)
/*<Finding the minimum of a vector v with size n.>*/
{
	float val = v[0];
	int i;

	for(i=0; i<n; i++) {
		if(v[i]<val) val=v[i];
	}
	return val;
}

void utils_cutarray(float *to, float *from, const int nx, const int ny, const int nz, const int snap)
/*<cut specific array >*/
{
	int i;
	for (i=0; i<nx*ny*nz; i++)	
	{
		to[i] = from[i+snap*nx*ny*nz];
	}
}

void utils_insertarray(float *to, float *from, const int nx, const int ny, const int nz, const int snap)
/*<cut specific array >*/
{
	int i;
	for (i=0; i<nx*ny*nz; i++)	
	{
		to[i+snap*nx*ny*nz]=from[i];
	}
}
