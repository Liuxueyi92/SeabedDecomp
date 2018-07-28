/*      Copyright(c)2013  Tongji Univ
*           All rights reserved
* free.c
*
* Free the Points
* Version:1.0
* Author:Chenlong Wang 
* 2013.7.31
*/

#include <rsf.h>

void free1(void *p)
/*<free a 1-d array>*/
{
	free(p);
}

void free2(void **p)
/*<free a 2-d array>*/
{
	free(p[0]);
	free(p);
}

void free3(void ***p)
/*<free a 3-d array>*/
{
    free(p[0][0]);
    free(p[0]);
    free(p);
}

void free4(void ****p)
/*<free a 4-d array>*/
{
    free(p[0][0][0]);
    free(p[0][0]);
    free(p[0]);
    free(p);
}

void free_int_1d(int *p)
/*<free a 1-d array of ints>*/
{
	free1(p);
}

void free_int_2d(int **p)
/*<free a 2-d array of ints>*/
{
    free2((void**)p);
}

void free_int_3d(int ***p)
/*<free a 3-d array of ints>*/
{
    free3((void***)p);
}

void free_int_4d(int ****p)
/*<free a 4-d array of ints>*/
{
    free4((void****)p);
}

void free_float_1d(float *p)
/*<free a 1-d array of floats>*/
{
    free1(p);
}

void free_float_2d(float **p)
/*<free a 2-d array of floats>*/
{
    free2((void**)p);
}

void free_float_3d(float ***p)
/*<free a 3-d array of floats>*/
{
    free3((void***)p);
}

void free_float_4d(float ****p)
/*<free a 4-d array of floats>*/
{
    free4((void****)p);
}

void free_double_1d(double *p)
/*<free a 1-d array of doubles>*/
{
    free1(p);
}

void free_double_2d(double **p)
/*<free a 2-d array of doubles>*/
{
    free2((void**)p);
}

void free_double_3d(double ***p)
/*<free a 3-d array of doubles>*/
{
    free3((void***)p);
}

void free_double_4d(double ****p)
/*<free a 4-d array of doubles>*/
{
    free4((void****)p);
}

void free_complex_1d(sf_complex *p)
/*<free a 1-d array of complexs>*/
{
    free1(p);
}

void free_complex_2d(sf_complex **p)
/*<free a 2-d array of complexs>*/
{
    free2((void**)p);
}
