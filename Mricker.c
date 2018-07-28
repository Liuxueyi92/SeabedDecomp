#include <rsf.h>
#include <stdio.h>
#include <stdlib.h>
#define PI  3.1415926
 
int main ( int argc, char *argv[] )
{
  int tn;
  float F,dt;
  int i=0;  

  sf_file  source;

  sf_init(argc,argv);

  if(!sf_getfloat("F",&F))    F=20.0;
  if(!sf_getint("tn",&tn))    tn=1000;
  if(!sf_getfloat("dt",&dt))  dt=0.005;  
  
  source=sf_output("out");     //OUTPUT TARGETS--Variable
  sf_putint(source, "n1" , tn);  //Caution:input database coordinate axis number!
  sf_putfloat(source, "d1" , dt);  


  float *w;
  w=sf_floatalloc(tn);  
  

  //****************calculate wavelet*******************//
  int del=0;
  float Nk=PI*PI*F*F*dt*dt;
  int tmp=500;
  for(i=0;i<tmp;i++)
  {      
    w[i]=(1.0-2.0*Nk*(i-tmp)*(i-tmp))*exp(-Nk*(i-tmp)*(i-tmp));
    if(w[i]>0.000005||w[i]<-0.000005)
        {
          del=tmp-i;
          break;
        }
  }
  for(i=0;i<tn;i++)
    w[i]=(1.0-2.0*Nk*(i-del)*(i-del))*exp(-Nk*(i-del)*(i-del));
  //***completed.

  //write ricker wavelet as @rsf file*******//
  sf_floatwrite(w,tn,source); 

 
  /*FILE *dataFile;  
  dataFile = fopen("ricker.txt","w+");
  fprintf(dataFile,"%d %d\n",400,1);
  for(i=0; i<600; i++)
     fprintf(dataFile,"%f\n",w[i]);
  fclose(dataFile);  */

  sf_close();
  exit(0);
}
