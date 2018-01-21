#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>
#include <string.h>
#include<time.h>
#include "phase_correction_ricker_decon.cpp"
#include "LSMOrCon_rec_2D.cpp"
#include "GPU_velocity_real.cpp"
#include "DisToTimeAndTimeToDis1D.cpp"
#include "segy.h"
#include "Resample.cpp"
#include "SGYWrite.cpp"
#include <cuda_runtime.h>
#define BLOCKSIZE 16
#define pi 3.1415926535898

__global__ void
	Equal(float *Dlf,float *Drt,float *Dup,float *Ddw,float *DFW0,float *DFW1,int NZ,int NX,int N2,int nfdmax)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;
	int NX1,NZ1,NX2,NZ2,mod_NX,mod_NZ;
	mod_NX=NX-2*N2;
	mod_NZ=NZ-2*N2;
	NX1=mod_NX+N2;
	NZ1=mod_NZ+N2;
	NX2=mod_NX*nfdmax;
	NZ2=mod_NZ*nfdmax;
	if(x<mod_NZ&&y<nfdmax) 
	{
		Dlf[x*nfdmax+y]=DFW0[(N2+x)*NX+N2-y-1];
		Drt[x*nfdmax+y]=DFW0[(N2+x)*NX+y+NX1];
		Dlf[x*nfdmax+y+NZ2]=DFW1[(N2+x)*NX+N2-y-1];
		Drt[x*nfdmax+y+NZ2]=DFW1[(N2+x)*NX+y+NX1];
	}
	if(x<nfdmax&&y<mod_NX)
	{
		Dup[x*mod_NX+y]=DFW0[(N2-x-1)*NX+y+N2];
		Ddw[x*mod_NX+y]=DFW0[(x+NZ1)*NX+y+N2];
		Dup[x*mod_NX+y+NX2]=DFW1[(N2-x-1)*NX+y+N2];
		Ddw[x*mod_NX+y+NX2]=DFW1[(x+NZ1)*NX+y+N2];
	}
	return;
}
__global__ void
	Add(float *Dv,float *DFW2,float *DFW0,float *DFW1,float *Dc,int *DIndex,float tao2,float h2,float vmin,float dv,int nfdmax,int NZ,int NX,int NXZ,int r_u,int r_x,float wavelet,float hzx2_1)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;
	int l,x1,x2,y1,y2,M_end,M_top,M_M;
	float w1;
	if(x<NZ&&y<NX)
	{
		M_end=DIndex[(int)((Dv[x*NX+y]-vmin)/dv+1.5)];
		M_top=DIndex[(int)((Dv[x*NX+y]-vmin)/dv+0.5)];					
		M_M=M_end-M_top;
		w1=(1.0+hzx2_1)*Dc[M_top]*DFW1[x*NX+y];
		for(l=1;l<M_M;l++)
		{
			x1=x-l;
			x2=x+l;
			y1=y-l;
			y2=y+l;
			if(x1<0) x1=-x1;
			if(x2>=NZ) x2=2*NZ-2-x2;
			if(y1<0) y1=-y1;
			if(y2>=NX) y2=2*NX-2-y2;

			w1+=Dc[l+M_top]*((DFW1[x1*NX+y]+DFW1[x2*NX+y])*hzx2_1+DFW1[x*NX+y1]+DFW1[x*NX+y2]);
		}

		DFW2[x*NX+y]=2*DFW1[x*NX+y]-DFW0[x*NX+y]+Dv[x*NX+y]*Dv[x*NX+y]*tao2*h2*w1;
		if(x==r_u&&y==r_x)
		{
			DFW2[x*NX+y]+=wavelet;
		}
	}
	return;
}

__global__ void
	Add_Con(float *Dv,float *DFW2,float *DFW0,float *DFW1,float *Dc,float tao2,float h2,int nfdmax,int NZ,int NX,int NXZ,int r_u,int r_x,float wavelet,float hzx2_1)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;
	int l,x1,x2,y1,y2;
	float w1;
	if(x<NZ&&y<NX)
	{
		w1=(1.0+hzx2_1)*Dc[0]*DFW1[x*NX+y];
		for(l=1;l<=nfdmax;l++)
		{
			x1=x-l;
			x2=x+l;
			y1=y-l;
			y2=y+l;
			if(x1<0) x1=-x1;
			if(x2>=NZ) x2=2*NZ-2-x2;
			if(y1<0) y1=-y1;
			if(y2>=NX) y2=2*NX-2-y2;

			w1+=Dc[l]*((DFW1[x1*NX+y]+DFW1[x2*NX+y])*hzx2_1+DFW1[x*NX+y1]+DFW1[x*NX+y2]);
		}

		DFW2[x*NX+y]=2.0*DFW1[x*NX+y]-DFW0[x*NX+y]+Dv[x*NX+y]*Dv[x*NX+y]*tao2*h2*w1;

		if(x==r_u&&y==r_x)
		{
			DFW2[x*NX+y]+=wavelet;
		}
	}
	return;
}

__global__ void
	Hybrid1(float *DFW1,float *DFW0,float *DFW2,float *DFWb,float *Dlf,float *Drt,float *Dup,float *Ddw,float *Dv,float *Dr_1,float *Dw,int NZ,int NX,int N2,int nfdmax,int k,float taoh,float taoh2)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;
	
	if(x>0&&x<=N2&&y>=(N2-x+2)&&y<(NX-N2-2+x))
	{
		DFWb[(N2-x)*NX+y]=1/(taoh*Dv[(N2-x)*NX+y]+1)*(taoh*Dv[(N2-x)*NX+y]*(DFW2[(N2-x+1)*NX+y]-DFW0[(N2-x+1)*NX+y]+DFW0[(N2-x)*NX+y])-(-2*DFW1[(N2-x)*NX+y]+DFW0[(N2-x)*NX+y]+DFW2[(N2-x+1)*NX+y]-2*DFW1[(N2-x+1)*NX+y]+DFW0[(N2-x+1)*NX+y])+taoh2*Dv[(N2-x)*NX+y]*Dv[(N2-x)*NX+y]*(DFW2[(N2-x+1)*NX+(y+1)]-2*DFW2[(N2-x+1)*NX+y]+DFW2[(N2-x+1)*NX+(y-1)]+DFW0[(N2-x)*NX+(y+1)]-2*DFW0[(N2-x)*NX+y]+DFW0[(N2-x)*NX+(y-1)]));//上边界
	}
	if(y>0&&y<=N2&&x>=(N2-y+2)&&x<(NZ-N2-2+y))
	{
		DFWb[x*NX+(N2-y)]=1/(taoh*Dv[x*NX+(N2-y)]+1)*(taoh*Dv[x*NX+(N2-y)]*(DFW2[x*NX+(N2-y+1)]-DFW0[x*NX+(N2-y+1)]+DFW0[x*NX+(N2-y)])-(-2*DFW1[x*NX+(N2-y)]+DFW0[x*NX+(N2-y)]+DFW2[x*NX+(N2-y+1)]-2*DFW1[x*NX+(N2-y+1)]+DFW0[x*NX+(N2-y+1)])+taoh2*Dv[(N2-y)*NX+x]*Dv[(N2-y)*NX+x]*(DFW2[(x+1)*NX+(N2-y+1)]-2*DFW2[x*NX+(N2-y+1)]+DFW2[(x-1)*NX+(N2-y+1)]+DFW0[(x+1)*NX+(N2-y)]-2*DFW0[x*NX+(N2-y)]+DFW0[(x-1)*NX+(N2-y)]));//左边界
	}
	if(x>0&&x<=N2&&y>=(N2-x+2)&&y<(NX-N2-2+x))
	{
		DFWb[(NZ-N2+x-1)*NX+y]=1/(taoh*Dv[(NZ-N2+x-1)*NX+y]+1)*(taoh*Dv[(NZ-N2+x-1)*NX+y]*(DFW2[(NZ-N2+x-2)*NX+y]-DFW0[(NZ-N2+x-2)*NX+y]+DFW0[(NZ-N2+x-1)*NX+y])-(-2*DFW1[(NZ-N2+x-1)*NX+y]+DFW0[(NZ-N2+x-1)*NX+y]+DFW2[(NZ-N2+x-2)*NX+y]-2*DFW1[(NZ-N2+x-2)*NX+y]+DFW0[(NZ-N2+x-2)*NX+y])+taoh2*Dv[(N2-x)*NX+y]*Dv[(N2-x)*NX+y]*(DFW2[(NZ-N2+x-2)*NX+(y+1)]-2*DFW2[(NZ-N2+x-2)*NX+y]+DFW2[(NZ-N2+x-2)*NX+(y-1)]+DFW0[(NZ-N2+x-1)*NX+(y+1)]-2*DFW0[(NZ-N2+x-1)*NX+y]+DFW0[(NZ-N2+x-1)*NX+(y-1)]));//下边界
	}
	if(y>0&&y<=N2&&x>=(N2-y+2)&&x<(NZ-N2-2+y))
	{
		DFWb[x*NX+NX-N2+y-1]=1/(taoh*Dv[x*NX+NX-N2+y-1]+1)*(taoh*Dv[x*NX+NX-N2+y-1]*(DFW2[x*NX+NX-N2+y-2]-DFW0[x*NX+NX-N2+y-2]+DFW0[x*NX+NX-N2+y-1])-(-2*DFW1[x*NX+NX-N2+y-1]+DFW0[x*NX+NX-N2+y-1]+DFW2[x*NX+NX-N2+y-2]-2*DFW1[x*NX+NX-N2+y-2]+DFW0[x*NX+NX-N2+y-2])+taoh2*Dv[(N2-y)*NX+x]*Dv[(N2-y)*NX+x]*(DFW2[(x+1)*NX+NX-N2+y-2]-2*DFW2[x*NX+NX-N2+y-2]+DFW2[(x-1)*NX+NX-N2+y-2]+DFW0[(x+1)*NX+NX-N2+y-1]-2*DFW0[x*NX+NX-N2+y-1]+DFW0[(x-1)*NX+NX-N2+y-1]));//右边界
	}
	if(x>0&&x<=N2&&y<NX)
	{
		DFWb[(N2-x+1)*NX+(N2-x)]=1/(2*Dr_1[(N2-x+1)*NX+(N2-x)]+1)*(DFW1[(N2-x+1)*NX+(N2-x)]+Dr_1[(N2-x+1)*NX+(N2-x)]*(DFW2[(N2-x+1)*NX+(N2-x+1)]+DFW2[(N2-x+2)*NX+(N2-x)]));
		DFWb[(N2-x)*NX+(N2-x+1)]=1/(2*Dr_1[(N2-x)*NX+(N2-x+1)]+1)*(DFW1[(N2-x)*NX+(N2-x+1)]+Dr_1[(N2-x)*NX+(N2-x+1)]*(DFW2[(N2-x)*NX+(N2-x+2)]+DFW2[(N2-x+1)*NX+(N2-x+1)]));// 左上角
		DFWb[(N2-x)*NX+(N2-x)]=1/(2*Dr_1[(N2-x)*NX+(N2-x)]+1)*(DFW1[(N2-x)*NX+(N2-x)]+Dr_1[(N2-x)*NX+(N2-x)]*(DFW2[(N2-x)*NX+(N2-x+1)]+DFW2[(N2-x+1)*NX+(N2-x)]));

		DFWb[(NZ-N2+x-2)*NX+NX-N2+x-1]=1/(2*Dr_1[(NZ-N2+x-2)*NX+NX-N2+x-1]+1)*(DFW1[(NZ-N2+x-2)*NX+NX-N2+x-1]+Dr_1[(NZ-N2+x-2)*NX+NX-N2+x-1]*(DFW2[(NZ-N2+x-2)*NX+NX-N2+x-2]+DFW2[(NZ-N2+x-3)*NX+NX-N2+x-1]));
		DFWb[(NZ-N2+x-1)*NX+NX-N2+x-2]=1/(2*Dr_1[(NZ-N2+x-1)*NX+NX-N2+x-2]+1)*(DFW1[(NZ-N2+x-1)*NX+NX-N2+x-2]+Dr_1[(NZ-N2+x-1)*NX+NX-N2+x-2]*(DFW2[(NZ-N2+x-1)*NX+NX-N2+x-3]+DFW2[(NZ-N2+x-2)*NX+NX-N2+x-2]));//右下角
		DFWb[(NZ-N2+x-1)*NX+NX-N2+x-1]=1/(2*Dr_1[(NZ-N2+x-1)*NX+NX-N2+x-1]+1)*(DFW1[(NZ-N2+x-1)*NX+NX-N2+x-1]+Dr_1[(NZ-N2+x-1)*NX+NX-N2+x-1]*(DFW2[(NZ-N2+x-1)*NX+NX-N2+x-2]+DFW2[(NZ-N2+x-2)*NX+NX-N2+x-1]));

		DFWb[(NZ-N2+x-1)*NX+(N2-x+1)]=1/(2*Dr_1[(NZ-N2+x-1)*NX+(N2-x+1)]+1)*(DFW1[(NZ-N2+x-1)*NX+(N2-x+1)]+Dr_1[(NZ-N2+x-1)*NX+(N2-x+1)]*(DFW2[(NZ-N2+x-2)*NX+(N2-x+1)]+DFW2[(NZ-N2+x-1)*NX+(N2-x+2)]));
		DFWb[(NZ-N2+x-2)*NX+(N2-x)]=1/(2*Dr_1[(NZ-N2+x-2)*NX+(N2-x)]+1)*(DFW1[(NZ-N2+x-2)*NX+(N2-x)]+Dr_1[(NZ-N2+x-2)*NX+(N2-x)]*(DFW2[(NZ-N2+x-3)*NX+(N2-x)]+DFW2[(NZ-N2+x-2)*NX+(N2-x+1)]));//左下角
		DFWb[(NZ-N2+x-1)*NX+(N2-x)]=1/(2*Dr_1[(NZ-N2+x-1)*NX+(N2-x)]+1)*(DFW1[(NZ-N2+x-1)*NX+(N2-x)]+Dr_1[(NZ-N2+x-1)*NX+(N2-x)]*(DFW2[(NZ-N2+x-2)*NX+(N2-x)]+DFW2[(NZ-N2+x-1)*NX+(N2-x+1)]));

		DFWb[(N2-x+1)*NX+NX-N2+x-1]=1/(2*Dr_1[(N2-x+1)*NX+NX-N2+x-1]+1)*(DFW1[(N2-x+1)*NX+NX-N2+x-1]+Dr_1[(N2-x+1)*NX+NX-N2+x-1]*(DFW2[(N2-x+1)*NX+NX-N2+x-2]+DFW2[(N2-x+2)*NX+NX-N2+x-1]));
		DFWb[(N2-x)*NX+NX-N2+x-2]=1/(2*Dr_1[(N2-x)*NX+NX-N2+x-2]+1)*(DFW1[(N2-x)*NX+NX-N2+x-2]+Dr_1[(N2-x)*NX+NX-N2+x-2]*(DFW2[(N2-x)*NX+NX-N2+x-3]+DFW2[(N2-x+1)*NX+NX-N2+x-2]));//右上角
		DFWb[(N2-x)*NX+NX-N2+x-1]=1/(2*Dr_1[(N2-x)*NX+NX-N2+x-1]+1)*(DFW1[(N2-x)*NX+NX-N2+x-1]+Dr_1[(N2-x)*NX+NX-N2+x-1]*(DFW2[(N2-x)*NX+NX-N2+x-2]+DFW2[(N2-x+1)*NX+NX-N2+x-1]));
	}

	return;
}

__global__ void
	Hybrid2(float *DFW1,float *DFW0,float *DFW2,float *DFWb,float *Dlf,float *Drt,float *Dup,float *Ddw,float *Dv,float *Dr_1,float *Dw,int NZ,int NX,int N2,int nfdmax,int k,float taoh,float taoh2)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;

	if (x>0&&x<=N2&&y>=(N2-x+1)&&y<(NX-N2+x))
	{
		DFW2[(N2-x)*NX+y]=(1-Dw[x])*DFW2[(N2-x)*NX+y]+Dw[x]*DFWb[(N2-x)*NX+y];//上边界
	}
	if(y>0&&y<=N2&&x>=(N2-y)&&x<(NZ-N2+y-1))
	{
		DFW2[x*NX+(N2-y)]=(1-Dw[y])*DFW2[x*NX+(N2-y)]+Dw[y]*DFWb[x*NX+(N2-y)];//左边界
	}
	if(x>0&&x<=N2&&y>=(N2-x)&&y<(NX-N2+x-1))
	{
		DFW2[(NZ-N2+x-1)*NX+y]=(1-Dw[x])*DFW2[(NZ-N2+x-1)*NX+y]+Dw[x]*DFWb[(NZ-N2+x-1)*NX+y];//下边界
	}
	if(y>0&&y<=N2&&x>=(N2-y+1)&&x<(NZ-N2+y))
	{
		DFW2[x*NX+NX-N2+y-1]=(1-Dw[y])*DFW2[x*NX+NX-N2+y-1]+Dw[y]*DFWb[x*NX+NX-N2+y-1];//右边界
	}
	return;
}
__global__ void
	Hybrid3(float *DFW1,float *DFW0,float *DFW2,float *DFWb,float *Dlf,float *Drt,float *Dup,float *Ddw,float *Dv,float *Dr_1,float *Dw,int NZ,int NX,int N2,int nfdmax,int k,float taoh,float taoh2)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;
	int NX1,NZ1,NX2,NZ2,mod_NX,mod_NZ;
	mod_NX=NX-2*N2;
	mod_NZ=NZ-2*N2;
	NX1=mod_NX+N2;
	NZ1=mod_NZ+N2;
	NX2=mod_NX*nfdmax;
	NZ2=mod_NZ*nfdmax;

	if(x<mod_NZ&&y<nfdmax) 
	{
		Dlf[x*nfdmax+y+k*NZ2]=DFW2[(N2+x)*NX+N2-y-1];
		Drt[x*nfdmax+y+k*NZ2]=DFW2[(N2+x)*NX+y+NX1];
	}
	if(x<nfdmax&&y<mod_NX)
	{
		Dup[x*mod_NX+y+k*NX2]=DFW2[(N2-x-1)*NX+y+N2];
		Ddw[x*mod_NX+y+k*NX2]=DFW2[(x+NZ1)*NX+y+N2];
	}
	return;
}

__global__ void
	Deliver(float *DFW1,float *DFW0,float *DFW2,int NZ,int NX)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;
	if(x<NZ&&y<NX)
	{
		DFW0[x*NX+y]=DFW1[x*NX+y];
		DFW1[x*NX+y]=DFW2[x*NX+y];
	}
	return;
}
__global__ void
	BKEqual(float *Dlf,float *Drt,float *Dup,float *Ddw,float *DFW1,int NZ,int NX,int N2,int nfdmax,int k)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;
	int NX1,NZ1,NX2,NZ2,mod_NX,mod_NZ;
	mod_NX=NX-2*N2;
	mod_NZ=NZ-2*N2;
	NX1=mod_NX+N2;
	NZ1=mod_NZ+N2;
	NX2=mod_NX*nfdmax;
	NZ2=mod_NZ*nfdmax;
	if(x<mod_NZ&&y<nfdmax) 
	{
		DFW1[(N2+x)*NX+N2-y-1]=Dlf[x*nfdmax+y+NZ2+k*NZ2];
		DFW1[(N2+x)*NX+y+NX1]=Drt[x*nfdmax+y+NZ2+k*NZ2];
	}
	if(x<nfdmax&&y<mod_NX)
	{
		DFW1[(N2-x-1)*NX+y+N2]=Dup[x*mod_NX+y+NX2+k*NX2];
		DFW1[(x+NZ1)*NX+y+N2]=Ddw[x*mod_NX+y+NX2+k*NX2];
	}
	return;
}
__global__ void
	BKAdd_EFF(float *Dv,float *DFW2,float *DFW0,float *DFW1,float *Dc,int *DIndex,float tao2,float h2,float vmin,float dv,int N2,int NZ,int NX,int NXZ,int r_u,int r_x,float wavelet,float hzx2_1)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;
	int l,x1,x2,y1,y2,M_end,M_top,M_M,NX1,NZ1;
	float w1;
	NX1=NX-N2;
	NZ1=NZ-N2;
	if(x>=N2&&x<NZ1&&y>=N2&&y<NX1)
	{
		M_end=DIndex[(int)((Dv[x*NX+y]-vmin)/dv+1.5)];
		M_top=DIndex[(int)((Dv[x*NX+y]-vmin)/dv+0.5)];					
		M_M=M_end-M_top;
		w1=(1.0+hzx2_1)*Dc[M_top]*DFW1[x*NX+y];
		for(l=1;l<M_M;l++)
		{
			x1=x-l;
			x2=x+l;
			y1=y-l;
			y2=y+l;
			if(x1<0) x1=-x1;
			if(x2>=NZ) x2=2*NZ-2-x2;
			if(y1<0) y1=-y1;
			if(y2>=NX) y2=2*NX-2-y2;

			w1+=Dc[l+M_top]*((DFW1[x1*NX+y]+DFW1[x2*NX+y])*hzx2_1+DFW1[x*NX+y1]+DFW1[x*NX+y2]);
		}

		DFW2[x*NX+y]=2.0*DFW1[x*NX+y]-DFW0[x*NX+y]+Dv[x*NX+y]*Dv[x*NX+y]*tao2*h2*w1;

		if(x==r_u&&y==r_x)
		{
			DFW2[x*NX+y]+=wavelet;
		}
	}
	return;
}


__global__ void
	BKAdd_EFF_Con(float *Dv,float *DFW2,float *DFW0,float *DFW1,float *Dc,float tao2,float h2,int N2,int NZ,int NX,int NXZ,int r_u,int r_x,int nfdmax,float wavelet,float hzx2_1)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;
	int l,x1,x2,y1,y2,NX1,NZ1;
	float w1;
	NX1=NX-N2;
	NZ1=NZ-N2;
	if(x>=N2&&x<NZ1&&y>=N2&&y<NX1)
	{
		w1=(1.0+hzx2_1)*Dc[0]*DFW1[x*NX+y];
		for(l=1;l<=nfdmax;l++)
		{
			x1=x-l;
			x2=x+l;
			y1=y-l;
			y2=y+l;
			if(x1<0) x1=-x1;
			if(x2>=NZ) x2=2*NZ-2-x2;
			if(y1<0) y1=-y1;
			if(y2>=NX) y2=2*NX-2-y2;

			w1+=Dc[l]*((DFW1[x1*NX+y]+DFW1[x2*NX+y])*hzx2_1+DFW1[x*NX+y1]+DFW1[x*NX+y2]);
		}

		DFW2[x*NX+y]=2.0*DFW1[x*NX+y]-DFW0[x*NX+y]+Dv[x*NX+y]*Dv[x*NX+y]*tao2*h2*w1;

		if(x==r_u&&y==r_x)
		{
			DFW2[x*NX+y]+=wavelet;
		}
	}
	return;
}


__global__ void
	Deliver_EFF(float *DFW1,float *DFW0,float *DFW2,int N2,int NZ,int NX)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;
	int NX1,NZ1;
	NX1=NX-N2;
	NZ1=NZ-N2;
	if(x>=N2&&x<NZ1&&y>=N2&&y<NX1)
	{
		DFW0[x*NX+y]=DFW1[x*NX+y];
		DFW1[x*NX+y]=DFW2[x*NX+y];
	}
	return;
}

__global__ void
	BKAdd(float *Dv,float *DBW2,float *DBW1,float *DBW0,float *Dc,int *DIndex,float *Dseis,float tao2,float h2,float vmin,float dv,int nfdmax, int NZ,int NX,int NXZ,int s_l,int s_r,int s_z,int ds,float hzx2_1)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;
	int l,judge,s,x1,x2,y1,y2,M_end,M_top,M_M;
	float w1;

	if(x<NZ&&y<NX)
	{		
		judge=(int)((y-s_l)%ds);
		s=(int)((y-s_l)/ds);
		if(x==s_z&&y>=s_l&&y<=s_r&&judge==0&&Dseis[s]!=0)
		{	
			DBW2[x*NX+y]=Dseis[s];
		}
		else
		{	
			M_end=DIndex[(int)((Dv[x*NX+y]-vmin)/dv+1.5)];
			M_top=DIndex[(int)((Dv[x*NX+y]-vmin)/dv+0.5)];					
			M_M=M_end-M_top;
			w1=(1.0+hzx2_1)*Dc[M_top]*DBW1[x*NX+y];
			for(l=1;l<M_M;l++)
			{
				x1=x-l;
				x2=x+l;
				y1=y-l;
				y2=y+l;
				if(x1<0) x1=-x1;
				if(x2>=NZ) x2=2*NZ-2-x2;
				if(y1<0) y1=-y1;
				if(y2>=NX) y2=2*NX-2-y2;

				w1=w1+Dc[l+M_top]*((DBW1[x1*NX+y]+DBW1[x2*NX+y])*hzx2_1+DBW1[x*NX+y1]+DBW1[x*NX+y2]);
			}

			DBW2[x*NX+y]=2*DBW1[x*NX+y]-DBW0[x*NX+y]+Dv[x*NX+y]*Dv[x*NX+y]*tao2*h2*w1;
		}

	}
	return;
}
__global__ void
	BKAdd_Con(float *Dv,float *DBW2,float *DBW1,float *DBW0,float *Dc,float *Dseis,float tao2,float h2,int nfdmax,int NZ,int NX,int NXZ,int s_l,int s_r,int s_z,int ds,float hzx2_1)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;
	int l,judge,s,x1,x2,y1,y2;
	float w1;

	if(x<NZ&&y<NX)
	{		
		judge=(int)((y-s_l)%ds);
		s=(int)((y-s_l)/ds);
		if(x==s_z&&y>=s_l&&y<=s_r&&judge==0&&Dseis[s]!=0)
		{			
			DBW2[x*NX+y]=Dseis[s];
		}
		else
		{	
			w1=(1.0+hzx2_1)*Dc[0]*DBW1[x*NX+y];
			for(l=1;l<=nfdmax;l++)
			{
				x1=x-l;
				x2=x+l;
				y1=y-l;
				y2=y+l;
				if(x1<0) x1=-x1;
				if(x2>=NZ) x2=2*NZ-2-x2;
				if(y1<0) y1=-y1;
				if(y2>=NX) y2=2*NX-2-y2;

				w1=w1+Dc[l]*((DBW1[x1*NX+y]+DBW1[x2*NX+y])*hzx2_1+DBW1[x*NX+y1]+DBW1[x*NX+y2]);
			}
			DBW2[x*NX+y]=2*DBW1[x*NX+y]-DBW0[x*NX+y]+Dv[x*NX+y]*Dv[x*NX+y]*tao2*h2*w1;
		}

	}
	return;
}


__global__ void
	BKHybrid1(float *DFW1,float *DFW0,float *DFW2,float *DFWb,float *Dv,float *Dr_1,float *Dw,int NZ,int NX,int N2,float taoh,float taoh2)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;

	if(x>0&&x<=N2&&y>=(N2-x+2)&&y<(NX-N2-2+x))
	{
		DFWb[(N2-x)*NX+y]=1/(taoh*Dv[(N2-x)*NX+y]+1)*(taoh*Dv[(N2-x)*NX+y]*(DFW2[(N2-x+1)*NX+y]-DFW0[(N2-x+1)*NX+y]+DFW0[(N2-x)*NX+y])-(-2*DFW1[(N2-x)*NX+y]+DFW0[(N2-x)*NX+y]+DFW2[(N2-x+1)*NX+y]-2*DFW1[(N2-x+1)*NX+y]+DFW0[(N2-x+1)*NX+y])+taoh2*Dv[(N2-x)*NX+y]*Dv[(N2-x)*NX+y]*(DFW2[(N2-x+1)*NX+(y+1)]-2*DFW2[(N2-x+1)*NX+y]+DFW2[(N2-x+1)*NX+(y-1)]+DFW0[(N2-x)*NX+(y+1)]-2*DFW0[(N2-x)*NX+y]+DFW0[(N2-x)*NX+(y-1)]));//上边界
	}
	if(y>0&&y<=N2&&x>=(N2-y+2)&&x<(NZ-N2-2+y))
	{
		DFWb[x*NX+(N2-y)]=1/(taoh*Dv[x*NX+(N2-y)]+1)*(taoh*Dv[x*NX+(N2-y)]*(DFW2[x*NX+(N2-y+1)]-DFW0[x*NX+(N2-y+1)]+DFW0[x*NX+(N2-y)])-(-2*DFW1[x*NX+(N2-y)]+DFW0[x*NX+(N2-y)]+DFW2[x*NX+(N2-y+1)]-2*DFW1[x*NX+(N2-y+1)]+DFW0[x*NX+(N2-y+1)])+taoh2*Dv[(N2-y)*NX+x]*Dv[(N2-y)*NX+x]*(DFW2[(x+1)*NX+(N2-y+1)]-2*DFW2[x*NX+(N2-y+1)]+DFW2[(x-1)*NX+(N2-y+1)]+DFW0[(x+1)*NX+(N2-y)]-2*DFW0[x*NX+(N2-y)]+DFW0[(x-1)*NX+(N2-y)]));//左边界
	}
	if(x>0&&x<=N2&&y>=(N2-x+2)&&y<(NX-N2-2+x))
	{
		DFWb[(NZ-N2+x-1)*NX+y]=1/(taoh*Dv[(NZ-N2+x-1)*NX+y]+1)*(taoh*Dv[(NZ-N2+x-1)*NX+y]*(DFW2[(NZ-N2+x-2)*NX+y]-DFW0[(NZ-N2+x-2)*NX+y]+DFW0[(NZ-N2+x-1)*NX+y])-(-2*DFW1[(NZ-N2+x-1)*NX+y]+DFW0[(NZ-N2+x-1)*NX+y]+DFW2[(NZ-N2+x-2)*NX+y]-2*DFW1[(NZ-N2+x-2)*NX+y]+DFW0[(NZ-N2+x-2)*NX+y])+taoh2*Dv[(N2-x)*NX+y]*Dv[(N2-x)*NX+y]*(DFW2[(NZ-N2+x-2)*NX+(y+1)]-2*DFW2[(NZ-N2+x-2)*NX+y]+DFW2[(NZ-N2+x-2)*NX+(y-1)]+DFW0[(NZ-N2+x-1)*NX+(y+1)]-2*DFW0[(NZ-N2+x-1)*NX+y]+DFW0[(NZ-N2+x-1)*NX+(y-1)]));//下边界
	}
	if(y>0&&y<=N2&&x>=(N2-y+2)&&x<(NZ-N2-2+y))
	{
		DFWb[x*NX+NX-N2+y-1]=1/(taoh*Dv[x*NX+NX-N2+y-1]+1)*(taoh*Dv[x*NX+NX-N2+y-1]*(DFW2[x*NX+NX-N2+y-2]-DFW0[x*NX+NX-N2+y-2]+DFW0[x*NX+NX-N2+y-1])-(-2*DFW1[x*NX+NX-N2+y-1]+DFW0[x*NX+NX-N2+y-1]+DFW2[x*NX+NX-N2+y-2]-2*DFW1[x*NX+NX-N2+y-2]+DFW0[x*NX+NX-N2+y-2])+taoh2*Dv[(N2-y)*NX+x]*Dv[(N2-y)*NX+x]*(DFW2[(x+1)*NX+NX-N2+y-2]-2*DFW2[x*NX+NX-N2+y-2]+DFW2[(x-1)*NX+NX-N2+y-2]+DFW0[(x+1)*NX+NX-N2+y-1]-2*DFW0[x*NX+NX-N2+y-1]+DFW0[(x-1)*NX+NX-N2+y-1]));//右边界
	}
	if(x>0&&x<=N2&&y<NX)
	{
		DFWb[(N2-x+1)*NX+(N2-x)]=1/(2*Dr_1[(N2-x+1)*NX+(N2-x)]+1)*(DFW1[(N2-x+1)*NX+(N2-x)]+Dr_1[(N2-x+1)*NX+(N2-x)]*(DFW2[(N2-x+1)*NX+(N2-x+1)]+DFW2[(N2-x+2)*NX+(N2-x)]));
		DFWb[(N2-x)*NX+(N2-x+1)]=1/(2*Dr_1[(N2-x)*NX+(N2-x+1)]+1)*(DFW1[(N2-x)*NX+(N2-x+1)]+Dr_1[(N2-x)*NX+(N2-x+1)]*(DFW2[(N2-x)*NX+(N2-x+2)]+DFW2[(N2-x+1)*NX+(N2-x+1)]));// 左上角
		DFWb[(N2-x)*NX+(N2-x)]=1/(2*Dr_1[(N2-x)*NX+(N2-x)]+1)*(DFW1[(N2-x)*NX+(N2-x)]+Dr_1[(N2-x)*NX+(N2-x)]*(DFW2[(N2-x)*NX+(N2-x+1)]+DFW2[(N2-x+1)*NX+(N2-x)]));

		DFWb[(NZ-N2+x-2)*NX+NX-N2+x-1]=1/(2*Dr_1[(NZ-N2+x-2)*NX+NX-N2+x-1]+1)*(DFW1[(NZ-N2+x-2)*NX+NX-N2+x-1]+Dr_1[(NZ-N2+x-2)*NX+NX-N2+x-1]*(DFW2[(NZ-N2+x-2)*NX+NX-N2+x-2]+DFW2[(NZ-N2+x-3)*NX+NX-N2+x-1]));
		DFWb[(NZ-N2+x-1)*NX+NX-N2+x-2]=1/(2*Dr_1[(NZ-N2+x-1)*NX+NX-N2+x-2]+1)*(DFW1[(NZ-N2+x-1)*NX+NX-N2+x-2]+Dr_1[(NZ-N2+x-1)*NX+NX-N2+x-2]*(DFW2[(NZ-N2+x-1)*NX+NX-N2+x-3]+DFW2[(NZ-N2+x-2)*NX+NX-N2+x-2]));//右下角
		DFWb[(NZ-N2+x-1)*NX+NX-N2+x-1]=1/(2*Dr_1[(NZ-N2+x-1)*NX+NX-N2+x-1]+1)*(DFW1[(NZ-N2+x-1)*NX+NX-N2+x-1]+Dr_1[(NZ-N2+x-1)*NX+NX-N2+x-1]*(DFW2[(NZ-N2+x-1)*NX+NX-N2+x-2]+DFW2[(NZ-N2+x-2)*NX+NX-N2+x-1]));

		DFWb[(NZ-N2+x-1)*NX+(N2-x+1)]=1/(2*Dr_1[(NZ-N2+x-1)*NX+(N2-x+1)]+1)*(DFW1[(NZ-N2+x-1)*NX+(N2-x+1)]+Dr_1[(NZ-N2+x-1)*NX+(N2-x+1)]*(DFW2[(NZ-N2+x-2)*NX+(N2-x+1)]+DFW2[(NZ-N2+x-1)*NX+(N2-x+2)]));
		DFWb[(NZ-N2+x-2)*NX+(N2-x)]=1/(2*Dr_1[(NZ-N2+x-2)*NX+(N2-x)]+1)*(DFW1[(NZ-N2+x-2)*NX+(N2-x)]+Dr_1[(NZ-N2+x-2)*NX+(N2-x)]*(DFW2[(NZ-N2+x-3)*NX+(N2-x)]+DFW2[(NZ-N2+x-2)*NX+(N2-x+1)]));//左下角
		DFWb[(NZ-N2+x-1)*NX+(N2-x)]=1/(2*Dr_1[(NZ-N2+x-1)*NX+(N2-x)]+1)*(DFW1[(NZ-N2+x-1)*NX+(N2-x)]+Dr_1[(NZ-N2+x-1)*NX+(N2-x)]*(DFW2[(NZ-N2+x-2)*NX+(N2-x)]+DFW2[(NZ-N2+x-1)*NX+(N2-x+1)]));

		DFWb[(N2-x+1)*NX+NX-N2+x-1]=1/(2*Dr_1[(N2-x+1)*NX+NX-N2+x-1]+1)*(DFW1[(N2-x+1)*NX+NX-N2+x-1]+Dr_1[(N2-x+1)*NX+NX-N2+x-1]*(DFW2[(N2-x+1)*NX+NX-N2+x-2]+DFW2[(N2-x+2)*NX+NX-N2+x-1]));
		DFWb[(N2-x)*NX+NX-N2+x-2]=1/(2*Dr_1[(N2-x)*NX+NX-N2+x-2]+1)*(DFW1[(N2-x)*NX+NX-N2+x-2]+Dr_1[(N2-x)*NX+NX-N2+x-2]*(DFW2[(N2-x)*NX+NX-N2+x-3]+DFW2[(N2-x+1)*NX+NX-N2+x-2]));//右上角
		DFWb[(N2-x)*NX+NX-N2+x-1]=1/(2*Dr_1[(N2-x)*NX+NX-N2+x-1]+1)*(DFW1[(N2-x)*NX+NX-N2+x-1]+Dr_1[(N2-x)*NX+NX-N2+x-1]*(DFW2[(N2-x)*NX+NX-N2+x-2]+DFW2[(N2-x+1)*NX+NX-N2+x-1]));
	}

	return;
}
__global__ void
	BKHybrid2(float *DFW1,float *DFW0,float *DFW2,float *DFWb,float *Dv,float *Dr_1,float *Dw,int NZ,int NX,int N2,float taoh,float taoh2)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;

	if(x>0&&x<=N2&&y>=(N2-x+1)&&y<(NX-N2+x))
	{
		DFW2[(N2-x)*NX+y]=(1-Dw[x])*DFW2[(N2-x)*NX+y]+Dw[x]*DFWb[(N2-x)*NX+y];//上边界
	}
	if(y>0&&y<=N2&&x>=(N2-y)&&x<(NZ-N2+y-1))
	{
		DFW2[x*NX+(N2-y)]=(1-Dw[y])*DFW2[x*NX+(N2-y)]+Dw[y]*DFWb[x*NX+(N2-y)];//左边界
	}
	if(x>0&&x<=N2&&y>=(N2-x)&&y<(NX-N2+x-1))
	{
		DFW2[(NZ-N2+x-1)*NX+y]=(1-Dw[x])*DFW2[(NZ-N2+x-1)*NX+y]+Dw[x]*DFWb[(NZ-N2+x-1)*NX+y];//下边界
	}
	if(y>0&&y<=N2&&x>=(N2-y+1)&&x<(NZ-N2+y))
	{
		DFW2[x*NX+NX-N2+y-1]=(1-Dw[y])*DFW2[x*NX+NX-N2+y-1]+Dw[y]*DFWb[x*NX+NX-N2+y-1];//右边界
	}
	return;
}

__global__ void
	Rel_NonCompen(float *DBW2,float *DFW2,float *Drel1,float *Drel2,int NZ,int NX)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;

	if(x<NZ&&y<NX)
	{
		Drel1[x*NX+y]+=DBW2[x*NX+y]*DFW2[x*NX+y];
		Drel2[x*NX+y]+=DFW2[x*NX+y]*DFW2[x*NX+y];
	}
	return;
}

__global__ void
	Rel_Compen(float *DBW2,float *DFW2,float *DBWADD,float *DFWADD,float *Drel1,float *Drel2,int NZ,int NX)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;

	if(x<NZ&&y<NX)
	{
		DFWADD[x*NX+y]+=DFW2[x*NX+y];
		DBWADD[x*NX+y]+=DBW2[x*NX+y];
		Drel1[x*NX+y]+=DBWADD[x*NX+y]*DFWADD[x*NX+y];
		Drel2[x*NX+y]+=DFW2[x*NX+y]*DFW2[x*NX+y];
	}
	return;
}

float f(float t1,float f0);  //ricker source

void velocity(char *OutNameVp,float *v,float *v_2,float *r,float *r_1,float *r_2,int NZ,int NX,int NXZ,int N2,float tao,float h,int choice) ;

void WriteSGY(float *Data,int NX,int NT,int tao3,float *SX,float *SY,float RX,float RY,float *DSR,char *FILENAME);

int main()
{
	cudaSetDevice(0);
	clock_t start, finish;
	float duration;
	start = clock (); 

	char name[100],name1[100],name3[100],OutNameseis[100],OutNameVp[100],OutNameDPR[100],OutPara[100],Result[100];
	int i,j,i1,i2,j1,j2,m,k,N,NC,nvel,mod_NZ,mod_NX,NZ,NX,Nmax,NT2,NT1,NT,NXZ,s_l,s_r,s_z,r_x,r_u,n,nfdmax,nfdmin,N2,nrec,dr,ds;
	int iNorm,iCompen,Nsmooth,iLSTE,nthita,ifv,size,sizec,sizew,sizeseis,sizein,sizeLF,sizeUP,NX_BG,NX_ED,NZ_BG,NZ_ED;
	int *M,*Index,*fdcoeneed,*DIndex;
	float h,hz,tao1,tao,taoh,tao2,h2,taoh2,f0,fmax,eps,df,vmax,vmin,dv,vel,wavelet,MIGmax,whitecoe,stable,hzx,hzx2_1,angle,wthite_phase;
	float *INRE,*FW2,*FW1,*FW0,*FWb,*BW2,*BW1,*BW0,*v,*v_2,*r,*r_1,*r_2,*c,**seis1,**seis,*Cseis,*w,*MIG1,*MIG2,*MIG3;
	float *DFWADD,*DFW2,*DFW0,*DFW1,*DFWb,*DBWADD,*DBW2,*DBW0,*DBW1,*Dv,*Dr,*Dr_1,*Dw,*Dc,*Drel1,*Drel2,*Dseis,*Dlf,*Drt,*Dup,*Ddw;

	FILE *inp,*out,*out1,*DepRec,*Rdseismic,*Af_la;

	inp=fopen("2D_Real_RVSP_RTM.txt","r");

	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%d\n",&nfdmax);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%d\n",&nfdmin);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%d\n",&N2);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%f\n",&f0);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%f\n",&fmax);	
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%f\n",&df);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%d\n",&nthita);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%f\n",&eps);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%f\n",&dv);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%d\n",&iLSTE);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%d\n",&ifv);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%f\n",&whitecoe);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%f\n",&hz);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%f\n",&tao);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%d\n",&iNorm);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%d\n",&iCompen);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%d\n",&Nsmooth);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%f\n",&wthite_phase);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%f\n",&angle);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%d\n",&NX_BG);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%d\n",&NX_ED);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%d\n",&NZ_BG);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%d\n",&NZ_ED);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%s\n",OutNameseis);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%s\n",OutNameVp);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%s\n",OutNameDPR);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%s\n",OutPara);
	fscanf(inp,"%[^\n]\n",name);
	fscanf(inp,"%s\n",Result);
	fclose(inp);

	inp=fopen(OutPara,"r");
	fscanf(inp,"%f \n%f \n%d \n%d \n%d \n%d \n%d \n%d \n%d \n%d \n%d \n%d",&h,&tao1,&mod_NZ,&mod_NX,&NT1,&s_l,&s_z,&n,&ds,&r_x,&nrec,&dr);
	fclose(inp);
	

	s_l=s_l+N2-1;
	r_x=r_x+N2-1;
	s_z=s_z+N2-1;
	s_r=(n-1)*ds+s_l;
	NZ=mod_NZ+2*N2;
	NX=mod_NX+2*N2;
	NT=(int)((NT1-1)*tao1/tao+1.5);
	taoh=tao/h;
	tao2=pow(tao,2);
	h2=1/pow(h,2);
	taoh2=tao2*h2/2;
	NXZ=NZ*NX;
	sizeLF=sizeof (float) *(nfdmax*mod_NZ*NT);
	sizeUP=sizeof (float) *(nfdmax*mod_NX*NT);
	size=sizeof (float) *(NXZ);
	sizew=sizeof(float)*(N2+1);
	sizeseis=sizeof(float)*(n);
	NT2=int(2.0/(f0*tao))+1;	
	hzx=hz/h;
	hzx2_1=1/(hzx*hzx);
	Nmax=NX;
	if(Nmax<NZ) Nmax=NZ;
	printf("ifv=%d\n",ifv);
	printf("Nmax=%d\n",Nmax);
	printf("hzx=%f,hzx2_1=%f\n",hzx,hzx2_1);
	printf("The maximum length of operator\nnfdmax=%d\n",nfdmax);
	printf("The minimum length of operator\nnfdmin=%d\n",nfdmin);
	printf("The hyrid absorbing boundary width\nN2=%d\n",N2);
	printf("space interval\nh=%f\n",h);
	printf("time interval\ntao=%f\n",tao);
	printf("z grid dimension\nmod_NZ=%d\n",mod_NZ);
	printf("x grid dimension\nmod_NX=%d\n",mod_NX);
	printf("actual grid number in z\nNZ=%d\n",NZ);
	printf("actual grid number in x\nNX=%d\n",NX);
	printf("Number of time number\nNT=%d\n",NT);
	printf("Source X\ns_x=%d\n",s_l);
	printf("Source Z\ns_z=%d\n",s_z);	
	printf("The number of sources\nn=%d\n",n);
	printf("The interval of sources\nds=%d\n",ds);
	printf("Dominant Frequency\nf0=%f\n",f0);
	printf("Maximum Frequency\nfmax=%f\n",fmax);
	printf("Interval of Frequency\ndf=%f\n",df);
	printf("Azimuth of the plane wave divide into the number\nnthita=%d\n",nthita);
	printf("dispersion value\neps=%f\n",eps);
	printf("velocity interval\ndv=%f\n",dv);
	printf("Interval of receiver\ndr=%d\n",dr);
	printf("The number of receivers\nnrec=%d\n",nrec);
	printf("Receiver Z\nr_x=%d\n",r_x);
	printf("LSM-0,TEM-1\niLSTE=%d\n",iLSTE);
	printf("hz=%f\n",hz);
	printf("iNorm=%d\n",iNorm);
	printf("iCompen=%d\n",iCompen);
	printf("Nsmooth=%d\n",Nsmooth);
	printf("wthite_phase=%f\n",wthite_phase);
	printf("whitecoe=%f\n",whitecoe);
	printf("angle=%f\n",angle);
	
	FW2=(float *)malloc(size);
	FW1=(float *)malloc(size);
	FW0=(float *)malloc(size);
	FWb=(float *)malloc(size);
	BW2=(float *)malloc(size);
	BW1=(float *)malloc(size);
	BW0=(float *)malloc(size);
	v=(float *)malloc(size);
	v_2=(float *)malloc(size);
	r=(float *)malloc(size);
	r_1=(float *)malloc(size);
	r_2=(float *)malloc(size);
	MIG1=(float *)malloc(size);
	MIG2=(float *)malloc(size);
	MIG3=(float *)malloc(size);
	Cseis=(float *)malloc(sizeseis);
	seis1=(float **)malloc(sizeof(float *)*(n));
	seis=(float **)malloc(sizeof(float *)*(n));
	for(i=0;i<n;i++)
	{
		seis1[i]=(float *)malloc(sizeof(float )*(NT1));
		seis[i]=(float *)malloc(sizeof(float )*(NT));
	}
	w=(float *)malloc(sizeof(float )*(N2+1));
	for (i=0;i<=N2;i++)
	{
		w[i]=(1.0*i)/(1.0*N2);
	}

	INRE=(float *)malloc(sizeof(float)*(nrec));
	
	DepRec=fopen(OutNameDPR,"r");
	for(i=0;i<nrec;i++)
	{
		fscanf(DepRec,"%f ",&INRE[i]);
	}
	fclose(DepRec);

	velocity(OutNameVp,v,v_2,r,r_1,r_2,NZ,NX,NXZ,N2,tao,h,ifv);

	vmin=v[0*NX+0];
	vmax=v[0*NX+0];

	for(i=0;i<NZ;i++)
	{
		for (j=0;j<NX;j++)
		{
			if(v[i*NX+j]<vmin) vmin=v[i*NX+j];
			if(v[i*NX+j]>vmax) vmax=v[i*NX+j];
		}
	}
	vel=(int(vmin/dv))*dv;
	if (vel>vmin) vmin=vel-dv;
	else vmin=vel;
	vel=(int(vmax/dv))*dv;
	if (vel<vmax) vmax=vel+dv;
	else vmax=vel;
	nvel=int((vmax-vmin)/dv+1.5);

	M=(int *)malloc(sizeof(int)*(nvel));
	fdcoeneed=(int *)malloc(sizeof(int)*(nvel-1));
	Index=(int *)malloc(sizeof(int)*(nvel+1));

	for (i=0;i<nvel;i++)
	{
		fdcoeneed[i]=0;
	}
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			k=int((v[i*NX+j]-vmin)/dv+0.5);
			fdcoeneed[k]=1;
		}
	}

	printf("vmin=%f\n",vmin);
	printf("vmax=%f\n",vmax);
	printf("nvel=%d\n",nvel);
	
	if(iLSTE==0)
	{
		NC=funMandC(nthita,nfdmax,nfdmin,nvel,tao,h,df,eps,fmax,vmin,vmax,dv,fdcoeneed,M,Index,iLSTE,&c,hzx);
	}
	else
	{
		NC=nfdmax+1;
		c=(float *)malloc(sizeof(float)*(NC));
		order(2*nfdmax,c);
	}
	
	sizec=sizeof (float) *(NC);
	sizein=sizeof(int)*(nvel+1);
	
	cudaMalloc (&DFWADD, size);
	cudaMalloc (&DFW2, size);
	cudaMalloc (&DFW1, size);
	cudaMalloc (&DFW0, size);
	cudaMalloc (&DFWb, size);
	cudaMalloc (&DBWADD, size);
	cudaMalloc (&DBW2, size);
	cudaMalloc (&DBW1, size);
	cudaMalloc (&DBW0, size);
	cudaMalloc (&Dv, size);
	cudaMalloc (&Dr, size);
	cudaMalloc (&Dr_1, size);
	cudaMalloc (&Dw, sizew);
	cudaMalloc (&Dc, sizec);
	cudaMalloc (&Drel1, size);
	cudaMalloc (&Drel2, size);
	cudaMalloc (&Dseis, sizeseis);
	cudaMalloc (&DIndex, sizein);
	cudaMalloc (&Dlf, sizeLF);
	cudaMalloc (&Drt, sizeLF);
	cudaMalloc (&Dup, sizeUP);
	cudaMalloc (&Ddw, sizeUP);

	cudaMemcpy (Dv, v, size, cudaMemcpyHostToDevice);
	cudaMemcpy (Dr, r, size, cudaMemcpyHostToDevice);
	cudaMemcpy (Dr_1, r_1, size, cudaMemcpyHostToDevice);
	cudaMemcpy (Dw, w, sizew, cudaMemcpyHostToDevice);
	cudaMemcpy (Dc, c, sizec, cudaMemcpyHostToDevice);
	cudaMemcpy (DIndex,Index, sizein, cudaMemcpyHostToDevice);
	
	dim3 blockDim (BLOCKSIZE, BLOCKSIZE);
	dim3 gridDim ((NZ+blockDim.x-1)/blockDim.x,(NX+blockDim.y-1)/blockDim.y);
	
	for(m=0;m<nrec;m++)
	{
		printf("/********************the number of %d receiver***********************/\n",m+1);
		N=int(INRE[m]);
		r_u=abs(N/hz)+N2-1;
		printf("r_u=%d r_x=%d N=%d\n",r_u,r_x,N);
	
		for(i=0;i<NXZ;i++)
		{
			FW0[i]=0.0; 
			FW1[i]=0.0;
		}
		FW1[r_u*NX+r_x]=f(0.0,f0)/2.0;

		cudaMemcpy (DFW0, FW0, size, cudaMemcpyHostToDevice);
		cudaMemcpy (DFW1, FW1, size, cudaMemcpyHostToDevice);

		Equal <<< gridDim, blockDim >>> (Dlf,Drt,Dup,Ddw,DFW0,DFW1,NZ,NX,N2,nfdmax);
		
		for(k=2;k<NT;k++)
		{
			if(k<NT2) wavelet=f((k-1)*tao,f0);
			else wavelet=0.0;
			if(iLSTE==0) Add <<< gridDim, blockDim >>> (Dv,DFW2,DFW0,DFW1,Dc,DIndex,tao2,h2,vmin,dv,nfdmax,NZ,NX,NXZ,r_u,r_x,wavelet,hzx2_1);
			else
			Add_Con<<< gridDim, blockDim >>> (Dv,DFW2,DFW0,DFW1,Dc,tao2,h2,nfdmax,NZ,NX,NXZ,r_u,r_x,wavelet,hzx2_1);
			Hybrid1 <<< gridDim, blockDim >>> (DFW1,DFW0,DFW2,DFWb,Dlf,Drt,Dup,Ddw,Dv,Dr_1,Dw,NZ,NX,N2,nfdmax,k,taoh,taoh2);
			Hybrid2 <<< gridDim, blockDim >>> (DFW1,DFW0,DFW2,DFWb,Dlf,Drt,Dup,Ddw,Dv,Dr_1,Dw,NZ,NX,N2,nfdmax,k,taoh,taoh2);
			Hybrid3 <<< gridDim, blockDim >>> (DFW1,DFW0,DFW2,DFWb,Dlf,Drt,Dup,Ddw,Dv,Dr_1,Dw,NZ,NX,N2,nfdmax,k,taoh,taoh2);
			Deliver <<< gridDim, blockDim >>> (DFW1,DFW0,DFW2,NZ,NX);
		}
		cudaMemcpy (BW0, DFW1, size, cudaMemcpyDeviceToHost);
		cudaMemcpy (BW1, DFW0, size, cudaMemcpyDeviceToHost);
		cudaMemcpy (DBW1, BW1, size, cudaMemcpyHostToDevice);
		cudaMemcpy (DBW0, BW0, size, cudaMemcpyHostToDevice);
	
		sprintf(name1,"NEW_L10-1932-X_%d.dat",N);
		strcpy(name3,OutNameseis);
		Rdseismic=fopen(strcat(name3,name1),"rb");
	
		for(i=0;i<n;i++)
		{
			for(k=0;k<NT1;k++)
			{
				fread(&seis1[i][k],sizeof(float),1,Rdseismic);
			}
		}
		fclose(Rdseismic);
		if(NT!=NT1)
		{
			for(i=0;i<n;i++)
			{
				resample (NT1, tao1, seis1[i], NT,tao,seis[i]);
			}		
		}
		else
		{
			for(i=0;i<n;i++)
			{
				for(k=0;k<NT;k++)
				{
					seis[i][k]=seis1[i][k];
				}
			
			}
		}
		printf("NT2=%d,NT1=%d,NT=%d,tao1=%f,tao=%f\n",NT2,NT1,NT,tao1,tao);
		
		for(i=0;i<NXZ;i++)
		{
			FWb[i]=0.0;
			if(iCompen==1)
			{
				BW2[i]=BW0[i]+BW1[i];
				FW2[i]=FW0[i]+FW1[i];
				MIG1[i]=FW2[i]*BW2[i]+FW0[i]*BW0[i];
				MIG2[i]=BW1[i]*BW1[i]+BW0[i]*BW0[i];
			}
			else
			{
				BW2[i]=0.0;	
				FW2[i]=0.0;	
				MIG1[i]=FW1[i]*BW1[i]+FW0[i]*BW0[i];
				MIG2[i]=BW1[i]*BW1[i]+BW0[i]*BW0[i];
			}	
		}
	
		cudaMemcpy (DFW2, FW2, size, cudaMemcpyHostToDevice);
		cudaMemcpy (DFWb, FWb, size, cudaMemcpyHostToDevice);
		cudaMemcpy (DFW0, FWb, size, cudaMemcpyHostToDevice);
		cudaMemcpy (DFW1, FWb, size, cudaMemcpyHostToDevice);
		cudaMemcpy (DFWADD, FW2, size, cudaMemcpyHostToDevice);
		cudaMemcpy (DBWADD, BW2, size, cudaMemcpyHostToDevice);
		cudaMemcpy (Drel1,MIG1, size, cudaMemcpyHostToDevice);
		cudaMemcpy (Drel2,MIG2, size, cudaMemcpyHostToDevice);
	
		for(k=NT-3;k>=0;k--)
		{
			if(k<NT2) wavelet=f((k+1)*tao,f0);
			else wavelet=0.0;
			
			BKEqual<<< gridDim, blockDim >>>(Dlf,Drt,Dup,Ddw,DBW1,NZ,NX,N2,nfdmax,k);
			
			if(iLSTE==0)
			{
				BKAdd_EFF <<< gridDim, blockDim >>> (Dv,DBW2,DBW0,DBW1,Dc,DIndex,tao2,h2,vmin,dv,N2,NZ,NX,NXZ,r_u,r_x,wavelet,hzx2_1);
			}
			else
			{
				BKAdd_EFF_Con<<< gridDim, blockDim >>> (Dv,DBW2,DBW0,DBW1,Dc,tao2,h2,N2,NZ,NX,NXZ,r_u,r_x,nfdmax,wavelet,hzx2_1);
			}
			
			Deliver_EFF <<< gridDim, blockDim >>> (DBW1,DBW0,DBW2,N2,NZ,NX);
			
			for(j=0;j<n;j++)
			{
				Cseis[j]=seis[j][k+1];
			}
			cudaMemcpy (Dseis, Cseis, sizeseis, cudaMemcpyHostToDevice);
			
			if(iLSTE==0)
			{
				BKAdd <<< gridDim, blockDim >>> (Dv,DFW2,DFW1,DFW0,Dc,DIndex,Dseis,tao2,h2,vmin,dv,nfdmax,NZ,NX,NXZ,s_l,s_r,s_z,ds,hzx2_1);
			}
			else
			{
				BKAdd_Con <<< gridDim, blockDim >>> (Dv,DFW2,DFW1,DFW0,Dc,Dseis,tao2,h2,nfdmax,NZ,NX,NXZ,s_l,s_r,s_z,ds,hzx2_1);
			}
			BKHybrid1 <<< gridDim, blockDim >>> (DFW1,DFW0,DFW2,DFWb,Dv,Dr_1,Dw,NZ,NX,N2,taoh,taoh2);
			BKHybrid2 <<< gridDim, blockDim >>> (DFW1,DFW0,DFW2,DFWb,Dv,Dr_1,Dw,NZ,NX,N2,taoh,taoh2);
			Deliver <<< gridDim, blockDim >>> (DFW1,DFW0,DFW2,NZ,NX);
		
			if(iCompen==1)
			{ 
				Rel_Compen <<< gridDim, blockDim >>>(DFW2,DBW2,DFWADD,DBWADD,Drel1,Drel2,NZ,NX); 
			}
			else
			{ 
				Rel_NonCompen <<< gridDim, blockDim >>>(DFW2,DBW2,Drel1,Drel2,NZ,NX);
			}
		}
		cudaMemcpy (MIG1, Drel1, size, cudaMemcpyDeviceToHost);
		cudaMemcpy (MIG2, Drel2, size, cudaMemcpyDeviceToHost);
		
		for(j=N2;j<NX-N2;j++)
		{
			for(i=N2;i<NZ-N2;i++)
			{		
				i1=i-1;
				i2=i+1;
				j1=j-1;
				j2=j+1;
				if(i1<N2) i1=2*N2-i1;
				if(j1<N2) j1=2*N2-j1;
				if(i2>=NZ-N2) i2=2*(NZ-N2-1)-i2;
				if(j2>=NX-N2) j2=2*(NX-N2-1)-j2;		
				MIG3[i*NX+j]=-1.0*(MIG1[i*NX+j2]+MIG1[i*NX+j1]+MIG1[i2*NX+j]+MIG1[i1*NX+j]-4*MIG1[i*NX+j])*v[(i)*NX+(j)]*v[(i)*NX+(j)]/(vmax*vmax);						
			}
		}
		
		sprintf(name1,"RVSP_RTM_up_%d.dat",m+1);
		strcpy(name3,Result);
		Af_la=fopen(strcat(name3,name1),"wb");	
		for(j=N2;j<NX-N2;j++)
		{
			for(i=N2;i<NZ-N2;i++)
			{
				fwrite(&MIG3[i*NX+j],sizeof(float),1,Af_la);
			}
		}
		fclose(Af_la);
		
		MIGmax=0.0;
		for(j=N2;j<NX-N2;j++)
		{
			for(i=N2;i<NZ-N2;i++)
			{	
				if(fabs(MIG2[i*NX+j])>MIGmax) MIGmax=fabs(MIG2[i*NX+j]);
			}
		}
		stable=MIGmax*whitecoe;
		printf("%0.16f\n",stable);
		for(j=N2;j<NX-N2;j++)
		{
			for(i=N2;i<NZ-N2;i++)
			{
				MIG2[i*NX+j]=MIG2[i*NX+j]+stable;
			}
		}
		sprintf(name1,"RVSP_RTM_down_%d.dat",m+1);
		strcpy(name3,Result);
		Af_la=fopen(strcat(name3,name1),"wb");	
		for(j=N2;j<NX-N2;j++)
		{
			for(i=N2;i<NZ-N2;i++)
			{
				fwrite(&MIG2[i*NX+j],sizeof(float),1,Af_la);
			}
		}
		fclose(Af_la);
	}
	for(j=0;j<NX;j++)
	{
		for(i=0;i<NZ;i++)
		{							
			FW0[i*NX+j]=0;
			FW1[i*NX+j]=0;
		}
	}
	for(m=0;m<nrec;m++)
	{
		sprintf(name1,"RVSP_RTM_up_%d.dat",m+1);
		strcpy(name3,Result);
		out=fopen(strcat(name3,name1),"rb");
		for(j=0;j<mod_NX;j++)
		{
			for(i=0;i<mod_NZ;i++)
			{
				fread(&MIG1[i*NX+j],sizeof(float),1,out);
			}
		}
		fclose(out);
		
		sprintf(name1,"RVSP_RTM_down_%d.dat",m+1);
		strcpy(name3,Result);
		out=fopen(strcat(name3,name1),"rb");
		for(j=0;j<mod_NX;j++)
		{
			for(i=0;i<mod_NZ;i++)
			{
				fread(&MIG2[i*NX+j],sizeof(float),1,out);
			}
		}
		fclose(out);

		for(j=0;j<mod_NX;j++)
		{
			for(i=0;i<mod_NZ;i++)
			{							
				FW0[i*NX+j]+=MIG1[i*NX+j];				
			}
		}
		for(j=0;j<mod_NX;j++)
		{
			for(i=0;i<mod_NZ;i++)
			{							
				FW1[i*NX+j]+=MIG2[i*NX+j];			
			}
		}
	}

	for(j=0;j<mod_NX;j++)
	{
		for(i=0;i<mod_NZ;i++)
		{
			MIG1[i*NX+j]=FW0[i*NX+j]/nrec;
			MIG2[i*NX+j]=FW1[i*NX+j]/nrec;
		}
	}
	if(iNorm==1)
	{
		for(j=0;j<mod_NX;j++)
		{
			for(i=0;i<mod_NZ;i++)
			{
				MIG1[i*NX+j]=MIG1[i*NX+j]/MIG2[i*NX+j];
			}
		}
	}

	if(ifv==1)
	{
		strcpy(name3,Result);
		out1=fopen(strcat(name3,"RVSP_Migration_Real_new2.dat"),"wb");
		for(j=NX_ED;j>NX_BG;j--)
		{
			for(i=NZ_BG;i<NZ_ED;i++)
			{
				fwrite(&MIG1[i*NX+j],sizeof(float),1,out1);
			}
		}
		fclose(out1);
		
		strcpy(name3,Result);
		out1=fopen(strcat(name3,"vnew.dat"),"wb");
		for(j=NX_ED+N2;j>NX_BG+N2;j--)
		{
			for(i=NZ_BG+N2;i<NZ_ED+N2;i++)
			{
				fwrite(&v[i*NX+j],sizeof(float),1,out1);
			}
		}
		fclose(out1);
	}
	else
	{
		strcpy(name3,Result);
		out1=fopen(strcat(name3,"RVSP_Migration_Real_new2.dat"),"wb");
		for(j=NX_BG;j<NX_ED;j++)
		{
			for(i=NZ_BG;i<NZ_ED;i++)
			{
				fwrite(&MIG1[i*NX+j],sizeof(float),1,out1);
			}
		}
		fclose(out1);
		
		strcpy(name3,Result);
		out1=fopen(strcat(name3,"vnew.dat"),"wb");
		for(j=NX_BG+N2;j<NX_ED+N2;j++)
		{
			for(i=NZ_BG;i<NZ_ED+N2;i++)
			{
				fwrite(&v[i*NX+j],sizeof(float),1,out1);
			}
		}
		fclose(out1);
	}
		
	int SGY_N;
	float *SGYseis,*DSR,*SX,*SY,RX,RY;
	char namesgy[100];
	SGY_N=NX_ED-NX_BG;
	mod_NX=SGY_N;
	float *MIG,*vDT,*MIG1_T,*MIG_T;
	int nt;

	MIG=(float *)malloc(sizeof (float) *(mod_NX*mod_NZ));
	vDT=(float *)malloc(sizeof (float) *(mod_NX*mod_NZ));
	
	
	strcpy(name3,Result);
	out1=fopen(strcat(name3,"RVSP_Migration_Real_new2.dat"),"rb");
	for(i=0;i<mod_NX;i++)
	{
		for(j=0;j<mod_NZ;j++)
		{		
			fread(&MIG[i*mod_NZ+j],sizeof(float),1,out1);
		}
	}
	fclose(out1);
	
	strcpy(name3,Result);
	out1=fopen(strcat(name3,"vnew.dat"),"rb");
	for(i=0;i<mod_NX;i++)
	{
		for(j=0;j<mod_NZ;j++)
		{		
			fread(&vDT[i*mod_NZ+j],sizeof(float),1,out1);
		}
	}
	fclose(out1);
		
	
	sprintf(name1,"RVSP_Migration_Real_T.dat",N);
	strcpy(name3,Result);
	nt=D2T(strcat(name3,name1), vDT, MIG, mod_NX, mod_NZ, 0, (mod_NX-1), hz, hz, tao);
	printf("nt=%d\n",nt);
	MIG_T=(float *)malloc(sizeof (float) *(mod_NX*nt));
	MIG1_T=(float *)malloc(sizeof (float) *(mod_NX*nt));
	
	strcpy(name3,Result);
	out1=fopen(strcat(name3,"RVSP_Migration_Real_T.dat"),"rb");
	for(i=0;i<mod_NX;i++)
	{
		for(j=0;j<nt;j++)
		{		
			fread(&MIG1_T[i*nt+j],sizeof(float),1,out1);
		}
	}
	fclose(out1);

	phase_correction(MIG1_T,MIG_T,mod_NX, nt,angle);
	
	strcpy(name3,Result);
	out1=fopen(strcat(name3,"RVSP_Migration_Real_T_phase.dat"),"wb");
	for(i=0;i<mod_NX;i++)
	{
		for(j=0;j<nt;j++)
		{		
			fwrite(&MIG_T[i*nt+j],sizeof(float),1,out1);
		}
	}
	fclose(out1);
	
	sprintf(name1,"RVSP_Migration_Real_D.dat",N);
	strcpy(name3,Result);
	nt=T2D(strcat(name3,name1),vDT, MIG_T, mod_NX, nt, mod_NZ, 0, (mod_NX-1), hz, tao,hz);
	printf("nt=%d\n",nt);
	
	SGYseis=(float *)malloc(sizeof(float )*(SGY_N*mod_NZ));

	SX=(float *)malloc(sizeof(float )*(SGY_N));
	SY=(float *)malloc(sizeof(float )*(SGY_N));
	DSR=(float *)malloc(sizeof(float )*(SGY_N));	

	strcpy(name3,Result);
	out1=fopen(strcat(name3,"RVSP_Migration_Real_new2.dat"),"rb");
	for(i=0;i<SGY_N;i++)
	{
		for(k=0;k<mod_NZ;k++)
		{
			fread(&SGYseis[i*mod_NZ+k],sizeof(float),1,out1);
		}
	}
	fclose(out1);
	for(i=0;i<SGY_N;i++)
	{
		SX[i]=1.0*i*10;
		SY[i]=-1.0*i*10;
		DSR[i]=1000-i;
	}
	RX=1.0;
	RY=1.0;
	
	strcpy(namesgy,Result);
	strcat(namesgy,"RVSP_migration_Real.sgy");

	WriteSGY(SGYseis,SGY_N,mod_NZ,(int)(hz),SX,SY,RX,RY,DSR,namesgy);
	
	finish = clock ();
	duration = (double) (finish - start) / CLOCKS_PER_SEC;
	printf ("%f seconds\n", duration);

	free(INRE);
	free(w);

	for(i=0;i<n;i++)
	{
		free(seis[i]);
		free(seis1[i]);
	}
	free(seis);
	free(Cseis);
	free(MIG2);
	free(MIG1);
	free(c);
	free(r_2);
	free(r_1);
	free(r);
	free(v_2);
	free(v);
	free(BW0);
	free(BW1);
	free(BW2);
	free(FWb);
	free(FW0);
	free(FW1);
	free(FW2);

	cudaFree(DIndex);
	cudaFree(Dseis);
	cudaFree(Drel2);
	cudaFree(Drel1);
	cudaFree(Dc);
	cudaFree(Dw);
	cudaFree(Dr_1);
	cudaFree(Dr);
	cudaFree(Dv);
	cudaFree(DBW0);
	cudaFree(DBW1);
	cudaFree(DBW2);
	cudaFree(DFWb);
	cudaFree(DFW0);
	cudaFree(DFW1);
	cudaFree(DFW2);
}


/**************************子波**************************/
float f(float t1,float f0)
{
	float t00=1/f0,y;
	y=(1-2*pow(pi*f0*(t1-t00),2))*exp(-pow(pi*f0*(t1-t00),2)); 
	return(y);
}
