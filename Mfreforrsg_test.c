#include <rsf.h>
#include <stdlib.h>
#include <stdio.h>
#include "free.h"
#include "utils.h"

#ifndef  PI
#define  PI  3.141592654
#endif

#ifndef  npower            //using in PML Boundary
#define  npower  3
#endif

#ifndef  N
#define  N 5
#endif
 
#ifndef c51
#define  c51   1.2112427
#endif
#ifndef c52
#define  c52   -8.972168e-2
#endif
#ifndef c53
#define  c53   1.3842773e-2
#endif
#ifndef c54
#define  c54   -1.7656599e-3
#endif
#ifndef c55
#define  c55   1.1867947e-4
#endif


#ifndef Oxv
#define Oxv(u,i,j)  \
    (c51*(u[i+0][j-1]-u[i-1][j+0])+ \
     c52*(u[i+1][j-2]-u[i-2][j+1])+ \
     c53*(u[i+2][j-3]-u[i-3][j+2])+ \
     c54*(u[i+3][j-4]-u[i-4][j+3])+ \
     c55*(u[i+4][j-5]-u[i-5][j+4]))
#endif

#ifndef Ozv
#define Ozv(u,i,j)  \
    (c51*(u[i+0][j+0]-u[i-1][j-1])+ \
     c52*(u[i+1][j+1]-u[i-2][j-2])+ \
     c53*(u[i+2][j+2]-u[i-3][j-3])+ \
     c54*(u[i+3][j+3]-u[i-4][j-4])+ \
     c55*(u[i+4][j+4]-u[i-5][j-5]))
#endif

#ifndef Oxt
#define Oxt(u,i,j)  \
   (c51*(u[i+1][j-0]-u[i-0][j+1])+ \
    c52*(u[i+2][j-1]-u[i-1][j+2])+ \
    c53*(u[i+3][j-2]-u[i-2][j+3])+ \
    c54*(u[i+4][j-3]-u[i-3][j+4])+ \
    c55*(u[i+5][j-4]-u[i-4][j+5]))
#endif

#ifndef Ozt
#define Ozt(u,i,j)  \
   (c51*(u[i+1][j+1]-u[i-0][j-0])+ \
    c52*(u[i+2][j+2]-u[i-1][j-1])+ \
    c53*(u[i+3][j+3]-u[i-2][j-2])+ \
    c54*(u[i+4][j+4]-u[i-3][j-3])+ \
    c55*(u[i+5][j+5]-u[i-4][j-4]))
#endif

float calmin(float a,float b)
{
  if(a>b)
    return b;
  else 
    return a;
}

void puthead(sf_file output,int n1,float d1,float o1,int n2,float d2,float o2,int n3,float d3,float o3)
{
  sf_putint(output, "n1" , n1);
  sf_putint(output, "n2" , n2);
  sf_putint(output, "n3" , n3);
  sf_putfloat(output, "d1" , d1);
  sf_putfloat(output, "d2" , d2);    
  sf_putfloat(output, "d3" , d3);
  sf_putfloat(output, "o1" , o1);
  sf_putfloat(output, "o2" , o2);
  sf_putfloat(output, "o3" , o3);
}


float pml_coefficent(float **vp,float dx,float dz,float dt,int xn,int zn,int pml)
{
  float *dxc=(float *)malloc(sizeof(float)*pml);
  float R=1e-1,delt=1.0e+10;
  float vpmax=0.0;  

  dx=calmin(dx,dz);

  //find the maximum velocity//
  for(int i=0;i<xn;i++)
    for(int j=0;j<zn;j++)
    {
      if(vpmax<vp[i][j])
        vpmax=vp[i][j];
    }

  float temp=0.0;
  for(int iflag=0;iflag<100;iflag++)
  {
    for(int i=0;i<pml;i++)
    {
      dxc[i]=0.5*(npower+1)*vpmax*logf(1.0/R)*powf(1.0*(pml-i)/pml,npower)/(pml*dx);
      temp=(1.0-dt*dxc[i])/(1.0+dt*dxc[i]);
      delt=calmin(temp,delt);
    }
    if(delt>0.6)
    {
      R=R/5.0;
      delt=1.0e+10;
      continue;
    }
    else
      break;
  }
  return R;
}

void vxpropagation(int xn,int zn,float dt,float dx,float dz,int t,int pml,float **vp,
                   float **p,float **vx,float **vx1,float **vx2,float **txx,float **txz,float R)
{    
    int i,j;
    
    int disx=0,disz=0;
    float dxac=0.0,dzac=0.0;     
    float attmpx=0.0,attmpz=0.0;
    attmpx=0.5*(npower+1)*logf(1.0/R)*powf(1.0/pml,npower)/(pml*dx);  
    attmpz=0.5*(npower+1)*logf(1.0/R)*powf(1.0/pml,npower)/(pml*dz); 

    for(i=N+pml-1;i<xn-pml-N;i++)
      for(j=N-1;j<zn-pml-N;j++)    
        vx[i][j]=vx[i][j]+dt/dx/p[i][j]/2*(Oxt(txx,i,j)+Ozt(txx,i,j))+dt/dz/p[i][j]/2*(Ozt(txz,i,j)-Oxt(txz,i,j));

    for(i=N-1;i<pml+N-1;i++) //left
      for(j=N-1;j<zn-pml-N;j++)
      { 
        disx=N+pml-1-i;
        dxac=attmpx*vp[i][j]*powf(1.0*disx,npower);
        dzac=0.0; 
        vx1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*vx1[i][j]+dt/dx/p[i][j]/2*(Oxt(txx,i,j)+Ozt(txx,i,j)));
        vx2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*vx2[i][j]+dt/dz/p[i][j]/2*(Ozt(txz,i,j)-Oxt(txz,i,j)));
        vx[i][j]=vx1[i][j]+vx2[i][j];
      }
    for(i=xn-pml-N;i<xn-N;i++)//right
      for(j=N-1;j<zn-pml-N;j++)
      {
        disx=i-xn+pml+N+1;
        dxac=attmpx*vp[i][j]*powf(1.0*disx,npower);
        dzac=0.0;
        vx1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*vx1[i][j]+dt/dx/p[i][j]/2*(Oxt(txx,i,j)+Ozt(txx,i,j)));
        vx2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*vx2[i][j]+dt/dz/p[i][j]/2*(Ozt(txz,i,j)-Oxt(txz,i,j)));
        vx[i][j]=vx1[i][j]+vx2[i][j];
      }
    for(i=N+pml-1;i<xn-pml-N;i++)//down
      for(j=zn-pml-N;j<zn-N;j++)
      {
        disz=j-zn+pml+N+1;
        dxac=0.0;
        dzac=attmpz*vp[i][j]*powf(1.0*disz,npower);
        vx1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*vx1[i][j]+dt/dx/p[i][j]/2*(Oxt(txx,i,j)+Ozt(txx,i,j)));
        vx2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*vx2[i][j]+dt/dz/p[i][j]/2*(Ozt(txz,i,j)-Oxt(txz,i,j)));
        vx[i][j]=vx1[i][j]+vx2[i][j];
      }
    for(i=N-1;i<pml+N-1;i++)//left-down
      for(j=zn-pml-N;j<zn-N;j++)
      {
        disx=N+pml-1-i;
        disz=j-zn+pml+N+1;
        dxac=attmpx*vp[i][j]*powf(1.0*disx,npower);
        dzac=attmpz*vp[i][j]*powf(1.0*disz,npower);
        vx1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*vx1[i][j]+dt/dx/p[i][j]/2*(Oxt(txx,i,j)+Ozt(txx,i,j)));
        vx2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*vx2[i][j]+dt/dz/p[i][j]/2*(Ozt(txz,i,j)-Oxt(txz,i,j)));
        vx[i][j]=vx1[i][j]+vx2[i][j];
      }
    for(i=xn-pml-N;i<xn-N;i++)//right-down
      for(j=zn-pml-N;j<zn-N;j++)
      {
        disx=i-xn+pml+N+1;
        disz=j-zn+pml+N+1;
        dxac=attmpx*vp[i][j]*powf(1.0*disx,npower);
        dzac=attmpz*vp[i][j]*powf(1.0*disz,npower);
        vx1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*vx1[i][j]+dt/dx/p[i][j]/2*(Oxt(txx,i,j)+Ozt(txx,i,j)));
        vx2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*vx2[i][j]+dt/dz/p[i][j]/2*(Ozt(txz,i,j)-Oxt(txz,i,j)));
        vx[i][j]=vx1[i][j]+vx2[i][j];
      }   
}    
void vzpropagation(int xn,int zn,float dt,float dx,float dz,int t,int pml,float **vp,
                   float **p,float **vz,float **vz1,float **vz2,float **txz,float **tzz,float R)
{      
    int i,j;
    int disx=0,disz=0;
    float dxac=0.0,dzac=0.0;
    float attmpx=0.0,attmpz=0.0;
    attmpx=0.5*(npower+1)*logf(1.0/R)*powf(1.0/pml,npower)/(pml*dx);  
    attmpz=0.5*(npower+1)*logf(1.0/R)*powf(1.0/pml,npower)/(pml*dz); 

    for(i=N+pml-1;i<xn-pml-N;i++)
      for(j=N-1;j<zn-pml-N;j++) 
        vz[i][j]=vz[i][j]+dt/dx/p[i][j]/2*(Oxt(txz,i,j)+Ozt(txz,i,j))+dt/dz/p[i][j]/2*(Ozt(tzz,i,j)-Oxt(tzz,i,j)); 

    for(i=N-1;i<pml+N-1;i++) //left
      for(j=N-1;j<zn-pml-N;j++)
      { 
        disx=N+pml-1-i;
        dxac=attmpx*vp[i][j]*powf(1.0*disx,npower);
        dzac=0.0; 
        vz1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*vz1[i][j]+dt/dx/p[i][j]/2*(Oxt(txz,i,j)+Ozt(txz,i,j)));
        vz2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*vz2[i][j]+dt/dz/p[i][j]/2*(Ozt(tzz,i,j)-Oxt(tzz,i,j)));
        vz[i][j]=vz1[i][j]+vz2[i][j];  
      }
    for(i=xn-pml-N;i<xn-N;i++)//right
      for(j=N-1;j<zn-pml-N;j++)
      {
        disx=i-xn+pml+N+1;
        dxac=attmpx*vp[i][j]*powf(1.0*disx,npower);
        dzac=0.0;
        vz1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*vz1[i][j]+dt/dx/p[i][j]/2*(Oxt(txz,i,j)+Ozt(txz,i,j)));
        vz2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*vz2[i][j]+dt/dz/p[i][j]/2*(Ozt(tzz,i,j)-Oxt(tzz,i,j)));
        vz[i][j]=vz1[i][j]+vz2[i][j];   
      }
    for(i=N+pml-1;i<xn-pml-N;i++)//down
      for(j=zn-pml-N;j<zn-N;j++)
      {
        disz=j-zn+pml+N+1;
        dxac=0.0;
        dzac=attmpz*vp[i][j]*powf(1.0*disz,npower);
        vz1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*vz1[i][j]+dt/dx/p[i][j]/2*(Oxt(txz,i,j)+Ozt(txz,i,j)));
        vz2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*vz2[i][j]+dt/dz/p[i][j]/2*(Ozt(tzz,i,j)-Oxt(tzz,i,j)));
        vz[i][j]=vz1[i][j]+vz2[i][j]; 
      }
    for(i=N-1;i<pml+N-1;i++)//left-down
      for(j=zn-pml-N;j<zn-N;j++)
      {
        disx=N+pml-1-i;
        disz=j-zn+pml+N+1;
        dxac=attmpx*vp[i][j]*powf(1.0*disx,npower);
        dzac=attmpz*vp[i][j]*powf(1.0*disz,npower);
        vz1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*vz1[i][j]+dt/dx/p[i][j]/2*(Oxt(txz,i,j)+Ozt(txz,i,j)));
        vz2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*vz2[i][j]+dt/dz/p[i][j]/2*(Ozt(tzz,i,j)-Oxt(tzz,i,j)));
        vz[i][j]=vz1[i][j]+vz2[i][j]; ; 
      }
    for(i=xn-pml-N;i<xn-N;i++)//right-down
      for(j=zn-pml-N;j<zn-N;j++)
      {
        disx=i-xn+pml+N+1;
        disz=j-zn+pml+N+1;
        dxac=attmpx*vp[i][j]*powf(1.0*disx,npower);
        dzac=attmpz*vp[i][j]*powf(1.0*disz,npower);
        vz1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*vz1[i][j]+dt/dx/p[i][j]/2*(Oxt(txz,i,j)+Ozt(txz,i,j)));
        vz2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*vz2[i][j]+dt/dz/p[i][j]/2*(Ozt(tzz,i,j)-Oxt(tzz,i,j)));
        vz[i][j]=vz1[i][j]+vz2[i][j];   
      }
}
void txxzzpropagation(int xn,int zn,float dt,float dx,float dz,int t,int pml,float **vp,float **vs,
                      float **p,float *walet,float **txx,float **txx1,float **txx2,float **tzz,
                      float **tzz1,float **tzz2,float **vx,float **vz,int slx,int slz,float R)
{
    int i,j;
    int disx=0,disz=0;
    float dxac=0.0,dzac=0.0;
    float attmpx=0.0,attmpz=0.0;
    attmpx=0.5*(npower+1)*logf(1.0/R)*powf(1.0/pml,npower)/(pml*dx);  
    attmpz=0.5*(npower+1)*logf(1.0/R)*powf(1.0/pml,npower)/(pml*dz); 

    for(i=N+pml;i<xn-pml-N+1;i++)
      for(j=N;j<zn-pml-N+1;j++){
        txx[i][j]+=dt/dx*p[i][j]*vp[i][j]*vp[i][j]/2*(Oxv(vx,i,j)+Ozv(vx,i,j))+
                   dt/dz*p[i][j]*(vp[i][j]*vp[i][j]-2*vs[i][j]*vs[i][j])/2*(Ozv(vz,i,j)-Oxv(vz,i,j));
        tzz[i][j]+=dt/dx*p[i][j]*(vp[i][j]*vp[i][j]-2*vs[i][j]*vs[i][j])/2*(Oxv(vx,i,j)+Ozv(vx,i,j))+
                   dt/dz*p[i][j]*vp[i][j]*vp[i][j]/2*(Ozv(vz,i,j)-Oxv(vz,i,j));  
        {//add source
           txx[i][j]+=walet[t]*exp(-pow(0.6,2.0)*(pow(i-slx,2.0)+pow(j-slz,2.0)));
           tzz[i][j]+=walet[t]*exp(-pow(0.6,2.0)*(pow(i-slx,2.0)+pow(j-slz,2.0)));
        }
    }


    for(i=N;i<pml+N;i++) //left
      for(j=N;j<zn-pml-N+1;j++)
      { 
        disx=N+pml-i;
        dxac=attmpx*vp[i][j]*powf(1.0*disx,npower);
        dzac=0.0; 
        txx1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*txx1[i][j]+
		   dt/dx*p[i][j]*vp[i][j]*vp[i][j]/2*(Oxv(vx,i,j)+Ozv(vx,i,j)));
        txx2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*txx2[i][j]+
		   dt/dz*p[i][j]*(vp[i][j]*vp[i][j]-2*vs[i][j]*vs[i][j])/2*(Ozv(vz,i,j)-Oxv(vz,i,j)));
        txx[i][j]=txx1[i][j]+txx2[i][j];

        tzz1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*tzz1[i][j]+
		   dt/dx*p[i][j]*(vp[i][j]*vp[i][j]-2*vs[i][j]*vs[i][j])/2*(Oxv(vx,i,j)+Ozv(vx,i,j)));
        tzz2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*tzz2[i][j]+
		   dt/dz*p[i][j]*vp[i][j]*vp[i][j]/2*(Ozv(vz,i,j)-Oxv(vz,i,j)));
        tzz[i][j]=tzz1[i][j]+tzz2[i][j];
      }
    for(i=xn-pml-N+1;i<xn-N+1;i++)//right
      for(j=N;j<zn-pml-N+1;j++)
      {
        disx=i-xn+pml+N;
        dxac=attmpx*vp[i][j]*powf(1.0*disx,npower);
        dzac=0.0;
        txx1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*txx1[i][j]+
		   dt/dx*p[i][j]*vp[i][j]*vp[i][j]/2*(Oxv(vx,i,j)+Ozv(vx,i,j)));
        txx2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*txx2[i][j]+
		   dt/dz*p[i][j]*(vp[i][j]*vp[i][j]-2*vs[i][j]*vs[i][j])/2*(Ozv(vz,i,j)-Oxv(vz,i,j)));
        txx[i][j]=txx1[i][j]+txx2[i][j];

        tzz1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*tzz1[i][j]+
		   dt/dx*p[i][j]*(vp[i][j]*vp[i][j]-2*vs[i][j]*vs[i][j])/2*(Oxv(vx,i,j)+Ozv(vx,i,j)));
        tzz2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*tzz2[i][j]+
		   dt/dz*p[i][j]*vp[i][j]*vp[i][j]/2*(Ozv(vz,i,j)-Oxv(vz,i,j)));
        tzz[i][j]=tzz1[i][j]+tzz2[i][j]; 
      }
    for(i=N+pml;i<xn-pml-N+1;i++)//down
      for(j=zn-pml-N+1;j<zn-N+1;j++)
      {
        disz=j-zn+pml+N;
        dxac=0.0;
        dzac=attmpz*vp[i][j]*powf(1.0*disz,npower);
        txx1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*txx1[i][j]+
		   dt/dx*p[i][j]*vp[i][j]*vp[i][j]/2*(Oxv(vx,i,j)+Ozv(vx,i,j)));
        txx2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*txx2[i][j]+
		   dt/dz*p[i][j]*(vp[i][j]*vp[i][j]-2*vs[i][j]*vs[i][j])/2*(Ozv(vz,i,j)-Oxv(vz,i,j)));
        txx[i][j]=txx1[i][j]+txx2[i][j];

        tzz1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*tzz1[i][j]+
		   dt/dx*p[i][j]*(vp[i][j]*vp[i][j]-2*vs[i][j]*vs[i][j])/2*(Oxv(vx,i,j)+Ozv(vx,i,j)));
        tzz2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*tzz2[i][j]+
		   dt/dz*p[i][j]*vp[i][j]*vp[i][j]/2*(Ozv(vz,i,j)-Oxv(vz,i,j)));
        tzz[i][j]=tzz1[i][j]+tzz2[i][j]; 
      }
    for(i=N;i<pml+N;i++)//left-down
      for(j=zn-pml-N+1;j<zn-N+1;j++)
      {
        disx=N+pml-i;
        disz=j-zn+pml+N;
        dxac=attmpx*vp[i][j]*powf(1.0*disx,npower);
        dzac=attmpz*vp[i][j]*powf(1.0*disz,npower);
        txx1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*txx1[i][j]+
		   dt/dx*p[i][j]*vp[i][j]*vp[i][j]/2*(Oxv(vx,i,j)+Ozv(vx,i,j)));
        txx2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*txx2[i][j]+
		   dt/dz*p[i][j]*(vp[i][j]*vp[i][j]-2*vs[i][j]*vs[i][j])/2*(Ozv(vz,i,j)-Oxv(vz,i,j)));
        txx[i][j]=txx1[i][j]+txx2[i][j];

        tzz1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*tzz1[i][j]+
		   dt/dx*p[i][j]*(vp[i][j]*vp[i][j]-2*vs[i][j]*vs[i][j])/2*(Oxv(vx,i,j)+Ozv(vx,i,j)));
        tzz2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*tzz2[i][j]+
		   dt/dz*p[i][j]*vp[i][j]*vp[i][j]/2*(Ozv(vz,i,j)-Oxv(vz,i,j)));
        tzz[i][j]=tzz1[i][j]+tzz2[i][j];
      }
    for(i=xn-pml-N+1;i<xn-N+1;i++)//right-down
      for(j=zn-pml-N+1;j<zn-N+1;j++)
      {
        disx=i-xn+pml+N;
        disz=j-zn+pml+N;
        dxac=attmpx*vp[i][j]*powf(1.0*disx,npower);
        dzac=attmpz*vp[i][j]*powf(1.0*disz,npower);
        txx1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*txx1[i][j]+
		   dt/dx*p[i][j]*vp[i][j]*vp[i][j]/2*(Oxv(vx,i,j)+Ozv(vx,i,j)));
        txx2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*txx2[i][j]+
		   dt/dz*p[i][j]*(vp[i][j]*vp[i][j]-2*vs[i][j]*vs[i][j])/2*(Ozv(vz,i,j)-Oxv(vz,i,j)));
        txx[i][j]=txx1[i][j]+txx2[i][j];

        tzz1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*tzz1[i][j]+
		   dt/dx*p[i][j]*(vp[i][j]*vp[i][j]-2*vs[i][j]*vs[i][j])/2*(Oxv(vx,i,j)+Ozv(vx,i,j)));
        tzz2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*tzz2[i][j]+
		   dt/dz*p[i][j]*vp[i][j]*vp[i][j]/2*(Ozv(vz,i,j)-Oxv(vz,i,j)));
        tzz[i][j]=tzz1[i][j]+tzz2[i][j]; 
      } 

   for(i=N;i<xn-N+1;i++) //freesurface
   {
     tzz[i][N-1]=-tzz[i][N+0];
     tzz[i][N-2]=-tzz[i][N+1];
     tzz[i][N-3]=-tzz[i][N+2];
     tzz[i][N-4]=-tzz[i][N+3];
     tzz[i][N-5]=-tzz[i][N+4];
   }
}
void txzpropagation(int xn,int zn,float dt,float dx,float dz,int t,int pml,float **vp,float **vs,
                    float **p,float **txz,float **txz1,float **txz2,float **vx,float **vz,float R)
{
    int i,j;
    int disx=0,disz=0;
    float dxac=0.0,dzac=0.0;
    float attmpx=0.0,attmpz=0.0;
    attmpx=0.5*(npower+1)*logf(1.0/R)*powf(1.0/pml,npower)/(pml*dx);  
    attmpz=0.5*(npower+1)*logf(1.0/R)*powf(1.0/pml,npower)/(pml*dz); 

    for(i=N+pml;i<xn-pml-N+1;i++)
      for(j=N;j<zn-pml-N+1;j++){
        txz[i][j]+=dt/dx*p[i][j]*vs[i][j]*vs[i][j]/2*(Oxv(vz,i,j)+Ozv(vz,i,j))+
                   dt/dz*p[i][j]*vs[i][j]*vs[i][j]/2*(Ozv(vx,i,j)-Oxv(vx,i,j)); 
    }

    for(i=N;i<pml+N;i++) //left
      for(j=N;j<zn-pml-N+1;j++)
      { 
        disx=N+pml-i;
        dxac=attmpx*vp[i][j]*powf(1.0*disx,npower);
        dzac=0.0; 
        txz1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*txz1[i][j]+dt/dx*p[i][j]*vs[i][j]*vs[i][j]/2*(Oxv(vz,i,j)+Ozv(vz,i,j)));
        txz2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*txz2[i][j]+dt/dz*p[i][j]*vs[i][j]*vs[i][j]/2*(Ozv(vx,i,j)-Oxv(vx,i,j)));
        txz[i][j]=txz1[i][j]+txz2[i][j]; 
      }
    for(i=xn-pml-N+1;i<xn-N+1;i++)//right
      for(j=N;j<zn-pml-N+1;j++)
      {
        disx=i-xn+pml+N;
        dxac=attmpx*vp[i][j]*powf(1.0*disx,npower);
        dzac=0.0;
        txz1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*txz1[i][j]+dt/dx*p[i][j]*vs[i][j]*vs[i][j]/2*(Oxv(vz,i,j)+Ozv(vz,i,j)));
        txz2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*txz2[i][j]+dt/dz*p[i][j]*vs[i][j]*vs[i][j]/2*(Ozv(vx,i,j)-Oxv(vx,i,j)));
        txz[i][j]=txz1[i][j]+txz2[i][j];  
      }
    for(i=N+pml;i<xn-pml-N+1;i++)//down
      for(j=zn-pml-N+1;j<zn-N+1;j++)
      {
        disz=j-zn+pml+N;
        dxac=0.0;
        dzac=attmpz*vp[i][j]*powf(1.0*disz,npower);
        txz1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*txz1[i][j]+dt/dx*p[i][j]*vs[i][j]*vs[i][j]/2*(Oxv(vz,i,j)+Ozv(vz,i,j)));
        txz2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*txz2[i][j]+dt/dz*p[i][j]*vs[i][j]*vs[i][j]/2*(Ozv(vx,i,j)-Oxv(vx,i,j)));
        txz[i][j]=txz1[i][j]+txz2[i][j];  
      }
    for(i=N;i<pml+N;i++)//left-down
      for(j=zn-pml-N+1;j<zn-N+1;j++)
      {
        disx=N+pml-i;
        disz=j-zn+pml+N;
        dxac=attmpx*vp[i][j]*powf(1.0*disx,npower);
        dzac=attmpz*vp[i][j]*powf(1.0*disz,npower);
        txz1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*txz1[i][j]+dt/dx*p[i][j]*vs[i][j]*vs[i][j]/2*(Oxv(vz,i,j)+Ozv(vz,i,j)));
        txz2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*txz2[i][j]+dt/dz*p[i][j]*vs[i][j]*vs[i][j]/2*(Ozv(vx,i,j)-Oxv(vx,i,j)));
        txz[i][j]=txz1[i][j]+txz2[i][j]; 
      }
    for(i=xn-pml-N+1;i<xn-N+1;i++)//right-down
      for(j=zn-pml-N+1;j<zn-N+1;j++)
      {
        disx=i-xn+pml+N;
        disz=j-zn+pml+N;
        dxac=attmpx*vp[i][j]*powf(1.0*disx,npower);
        dzac=attmpz*vp[i][j]*powf(1.0*disz,npower);
        txz1[i][j]=1.0/(1.0+0.5*dt*dxac)*((1.0-0.5*dt*dxac)*txz1[i][j]+dt/dx*p[i][j]*vs[i][j]*vs[i][j]/2*(Oxv(vz,i,j)+Ozv(vz,i,j)));
        txz2[i][j]=1.0/(1.0+0.5*dt*dzac)*((1.0-0.5*dt*dzac)*txz2[i][j]+dt/dz*p[i][j]*vs[i][j]*vs[i][j]/2*(Ozv(vx,i,j)-Oxv(vx,i,j)));
        txz[i][j]=txz1[i][j]+txz2[i][j];  
      } 

   for(i=N;i<xn-N+1;i++) //freesurface
   {
     txz[i][N-1]=-txz[i][N+0];
     txz[i][N-2]=-txz[i][N+1];
     txz[i][N-3]=-txz[i][N+2];
     txz[i][N-4]=-txz[i][N+3];
     txz[i][N-5]=-txz[i][N+4];
   }
}

int main ( int argc, char *argv[] )
{
  utils_print_title("Forward modeling using Rotated Staggered Grid");
	
  int xn,zn,tn,pml; //ts=write time start,ti=write time interval,tn=time number//
  float dx,dz,dt;
  int slx,slz,shotnum;
  int re_dep;
  int i=0,j=0,t=0;  
  FILE *fp;

  sf_init(argc,argv);  
      
  if(!sf_getint("xn",&xn))   xn=100;        
  if(!sf_getint("zn",&zn))   zn=200;
  if(!sf_getint("tn",&tn))   tn=2000;
  if(!sf_getfloat("dx",&dx)) dx=0.005;
  if(!sf_getfloat("dz",&dz)) dz=0.005;
  if(!sf_getfloat("dt",&dt)) dt=0.0005;
  if(!sf_getint("pml",&pml)) pml=30;
  if(!sf_getint("slx",&slx)) slx=(int)xn/2; //source xline position
  if(!sf_getint("slz",&slz)) slz=(int)zn/2; //source zline position
  if(!sf_getint("shotnum",&shotnum)) shotnum=1;
  if(!sf_getint("re_dep",&re_dep)) re_dep=104;

  sf_file ricker,vpvel,vsvel,pdens,vxrecord,vzrecord,precord;
 
  vxrecord=sf_output("out");  //seimic record file variables
  vzrecord=sf_output("vzrecord");
  precord=sf_output("precord");

  vpvel=sf_input("vpvel");
  vsvel=sf_input("vsvel");
  pdens=sf_input("pdens");

  float **vxdata,**vzdata,**pdata;   //seimic record CPU variables
  vxdata=sf_floatalloc2(tn,xn-2*(pml+N));
  vzdata=sf_floatalloc2(tn,xn-2*(pml+N));
  pdata=sf_floatalloc2(tn,xn-2*(pml+N));
  for(i=0;i<xn-2*(pml+N);i++)
    for(j=0;j<tn;j++)
    { vxdata[i][j]=0.0; vzdata[i][j]=0.0; pdata[i][j]=0.0; } 

  //*****allocating vp vs p model and initialization**********//
  float **vp,**vs,**p;
  vp=sf_floatalloc2(zn,xn);
  vs=sf_floatalloc2(zn,xn);
  p =sf_floatalloc2(zn,xn);
  for(i=0;i<xn;i++)
    for(j=0;j<zn;j++)
    { vp[i][j]=0.0;vs[i][j]=0.0;p[i][j]=0.0; } 

  //**write velocity model from command lines or files
  //**caurion: the order reading velocity model
  sf_floatread(vp[0],xn*zn,vpvel);
  sf_floatread(vs[0],xn*zn,vsvel);
  sf_floatread(p[0],xn*zn,pdens); 
/*
  for(i=250;i<xn;i++)
	  for(j=104;j<254;j++){
		  vp[i][j]=2.1;
		  vs[i][j]=0.6;
	  }
*/
/*
  for(i=0;i<xn;i++)
	  for(j=104;j<250;j++){
		  vp[i][j]=2.1;
		  vs[i][j]=0.6;
	  }
*/
/*
  int smoleg=60;
  float tmp1,tmp2;
  for(i=smoleg;i<xn-smoleg;i++)
	  for(j=104;j<254;j++){
		  tmp1=0;
		  tmp2=0;
		  for(int index=-smoleg;index<=smoleg;index++){
			tmp1+=vp[i+index][j]/(2*smoleg+1);
			tmp2+=vs[i+index][j]/(2*smoleg+1);
		  }
		  vp[i][j]=tmp1;
		  vs[i][j]=tmp2;
	  }
*/

/*
  if((fp=fopen("vp_new.dat","wb+"))!=NULL){
    for(i=pml+N;i<xn-pml-N;i++)
		for(j=N;j<zn-pml-N;j++)
		    fwrite(&vp[i][j],sizeof(float),1,fp);
	fclose(fp);
  }
*/

  /*output vp and vs velocity nearing the seabed*/ 
/*  if((fp=fopen("vpvel.dat","wb+"))!=NULL){
    for(i=pml+N;i<xn-pml-N;i++)
		fwrite(&vp[i][104],sizeof(float),1,fp);
	fclose(fp);
  }
  if((fp=fopen("vsvel.dat","wb+"))!=NULL){
    for(i=pml+N;i<xn-pml-N;i++){
		fwrite(&vs[i][104],sizeof(float),1,fp);
	}
	fclose(fp);
  }
*/

  float R=0.0;
  //calculate the pml attenuation coeffient !!R!!
  R=pml_coefficent(vp,dx,dz,dt,xn,zn,pml);

  //*****allocating vx vz txx tzz txz and initialing**********//
  float **vx,**vz,**txx,**tzz,**txz;
  vx =sf_floatalloc2(zn,xn);
  vz =sf_floatalloc2(zn,xn);
  txx=sf_floatalloc2(zn,xn);
  tzz=sf_floatalloc2(zn,xn);
  txz=sf_floatalloc2(zn,xn);
  for(i=0;i<xn;i++)
    for(j=0;j<zn;j++)
    { vx[i][j]=0.0;vz[i][j]=0.0;txx[i][j]=0.0;tzz[i][j]=0.0;txz[i][j]=0.0; } 

  float **vx1,**vz1,**txx1,**tzz1,**txz1;
  vx1 =sf_floatalloc2(zn,xn);
  vz1 =sf_floatalloc2(zn,xn);
  txx1=sf_floatalloc2(zn,xn);
  tzz1=sf_floatalloc2(zn,xn);
  txz1=sf_floatalloc2(zn,xn);
  for(i=0;i<xn;i++)
    for(j=0;j<zn;j++)
    { vx1[i][j]=0.0;vz1[i][j]=0.0;txx1[i][j]=0.0;tzz1[i][j]=0.0;txz1[i][j]=0.0; } 
  float **vx2,**vz2,**txx2,**tzz2,**txz2;
  vx2 =sf_floatalloc2(zn,xn);
  vz2 =sf_floatalloc2(zn,xn);
  txx2=sf_floatalloc2(zn,xn);
  tzz2=sf_floatalloc2(zn,xn);
  txz2=sf_floatalloc2(zn,xn);
  for(i=0;i<xn;i++)
    for(j=0;j<zn;j++)
    { vx2[i][j]=0.0;vz2[i][j]=0.0;txx2[i][j]=0.0;tzz2[i][j]=0.0;txz2[i][j]=0.0; } 

  //**********read source file********//
  //if (!sf_histint(Fw,"n1",&nt))   sf_error("No n1= in wav");//read coordinate parameters from input files.rsf
  //if (!sf_histfloat(Fw,"d1",&dt)) sf_error("No d1= in wav");
  ricker=sf_input("in"); 
  float *walet;
  walet=sf_floatalloc(tn);  sf_floatread(walet , tn, ricker);
  //if (!sf_histint(source, "n1",&tn))   sf_error("No n1= in wav");
  
  //***write coordinate axes of output files (snapshot)********// 
  puthead(vxrecord,tn,dt,0.0,xn-2*(pml+N),dx,0.0,shotnum,1.0,0.0);
  puthead(vzrecord,tn,dt,0.0,xn-2*(pml+N),dx,0.0,shotnum,1.0,0.0);
  puthead(precord,tn,dt,0.0,xn-2*(pml+N),dx,0.0,shotnum,1.0,0.0);

  //for(int shot=0;shot<shotnum;shot++) //multi-shot begin
  {
     //slx=200+8*shot;
     //initialization variables**************//
     for(i=0;i<xn;i++)
       for(j=0;j<zn;j++)
       { 
          vx[i][j]=0.0;vz[i][j]=0.0;txx[i][j]=0.0;tzz[i][j]=0.0;txz[i][j]=0.0; 
          vx1[i][j]=0.0;vz1[i][j]=0.0;txx1[i][j]=0.0;tzz1[i][j]=0.0;txz1[i][j]=0.0; 
          vx2[i][j]=0.0;vz2[i][j]=0.0;txx2[i][j]=0.0;tzz2[i][j]=0.0;txz2[i][j]=0.0;
       }
 
     for(t=0;t<tn;t++)
     {
       vxpropagation(xn,zn,dt,dx,dz,t,pml,vp,p,vx,vx1,vx2,txx,txz,R);
       vzpropagation(xn,zn,dt,dx,dz,t,pml,vp,p,vz,vz1,vz2,txz,tzz,R);
       txxzzpropagation(xn,zn,dt,dx,dz,t,pml,vp,vs,p,walet,txx,txx1,txx2,tzz,tzz1,tzz2,vx,vz,slx,slz,R);
       txzpropagation(xn,zn,dt,dx,dz,t,pml,vp,vs,p,txz,txz1,txz2,vx,vz,R);
  
       for(i=0;i<xn-2*(pml+N);i++) //write vx/vz record
       { 
         vxdata[i][t]=vx[i+pml+N][re_dep];
         vzdata[i][t]=vz[i+pml+N][re_dep];
//         pdata[i][t] -=(tzz[i+pml+N][re_dep-1]+tzz[i+pml+N+1][re_dep-1]+txx[i+pml+N][re_dep-1]+txx[i+pml+N+1][re_dep-1])/4/1.5;
         pdata[i][t] -=(tzz[i+pml+N][re_dep]+tzz[i+pml+N+1][re_dep]+txx[i+pml+N][re_dep]+txx[i+pml+N+1][re_dep])/4;
       }
     
	   utils_loadbar(t,tn,1000,50);
       //if(t%10==0)
         //sf_floatwrite(vx[0] ,zn*xn, vxrecord);     
     }

     sf_floatwrite(vxdata[0] ,tn*(xn-2*(pml+N)), vxrecord);     
     sf_floatwrite(vzdata[0] ,tn*(xn-2*(pml+N)), vzrecord); 
     sf_floatwrite(pdata[0]  ,tn*(xn-2*(pml+N)), precord ); 

     for(i=0;i<(xn-2*(pml+N));i++)
       for(j=0;j<tn;j++)
       { vxdata[i][j]=0.0; vzdata[i][j]=0.0; pdata[i][j]=0.0;}
  } // multi-shot end

  free_float_2d(vp); free_float_2d(vs); free_float_2d(p);
  free_float_2d(vx); free_float_2d(vz); free_float_2d(txx); free_float_2d(tzz); free_float_2d(txz);
  free_float_2d(vx1);free_float_2d(vz1);free_float_2d(txx1);free_float_2d(tzz1);free_float_2d(txz1);
  free_float_2d(vx2);free_float_2d(vz2);free_float_2d(txx2);free_float_2d(tzz2);free_float_2d(txz2);

  sf_close();
  exit(0);
}
