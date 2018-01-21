#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include "string.h"
#define pai 3.1415926535898
float Ricker(float t1,float f0) ;
void ricker(float *r1,int nt,float fp,float ti) 
{
	float t1,t2,t3;
	int i;
	for(i=0;i<nt;i++)
	{  t1=(i-nt/2)*ti;
	t2=pai*fp*t1;
	t3=t2*t2;
	r1[i]=(1-2*t3)*exp(-t3); 
	}
}
void xcor(float *x,int nx,float *y,int ny,float *z,int nz)
{
	int i;
	int j;
	int k;
	float s;
	for (i=1-nx;i<ny;i++)
	{
		s=0;
		for (j=0;j<nx;j++)
		{
			k=i+j;
			if(k>=0&&k<=ny-1)
			{
				s=s+x[j]*y[k];
			}
		}
		z[i+nx-1]=s;
	}
}
void kbfft(double *pr,double *pi,int n,int k,double *fr,double *fi,int l,int il);
void phase_correction(float *din,float *dout,int ntr, int nt,float angle);
void SpectrumDivide(float *din,float *din2,float *dout,int ntr, int nt,int nt2,float whitecoe);
/*int main()
{
	FILE *fpe,*fpe1;	//debug file
	int i;
	int ni;
	int itr;
	int ntr;
	int ncor;
	double dt;
	float *f;
	float *fmat;
	float *fmat2;

	float *ric1;
	float *ric2;
	float *fph_c;
	float *fcor;

	ni=256;
	ntr=4;
	ncor=2*ni-1;
	dt=0.001;

	f=(float *)malloc(sizeof(float)*(ni));
	fph_c=(float *)malloc(sizeof(float)*(ni));

	fmat=(float *)malloc(sizeof(float)*(ncor)*ntr);   
	fcor=(float *)malloc(sizeof(float)*(ncor));
	ric1=(float *)malloc(sizeof(float)*(ni)); 
	ric2=(float *)malloc(sizeof(float)*(ni));

	for(i=0;i<ni;i++)
	{
		f[i]=0.0;
		fph_c[i]=0.0;
		ric1[i]=0.0;
		ric2[i]=0.0;
	}

	// make ricker
	ricker(f,ni,30,dt);
	ricker(ric1,ni,30,dt);
	ricker(ric2,ni,30,dt);
	// make xcor mat
	for(i=0;i<ncor;i++)
		fcor[i]=0;
	xcor(ric1,ni,ric2,ni,fcor,ncor);
	for(itr=0;itr<ntr;itr++)
	{
		for(i=0;i<ncor;i++)
			fmat[itr*ncor+i]=0.0;
	}
	for(itr=0;itr<ntr;itr++)
	{
		for(i=0;i<ncor;i++)
			fmat[itr*ncor+i]=fcor[i];
	}

	// phase correction
	fpe=fopen("debug.txt","w");
	phase_correction(f,fph_c,1,ni,90);
//					input
//						output
						//trace number
							//nt
							  //angle
	for(i=0;i<ni;i++)
	{
		fprintf(fpe,"%f %f\n",f[i],fph_c[i]);
	}
	fclose(fpe);

	// ricker decon
	fpe1=fopen("decon.dat","wb");
	for(itr=0;itr<ntr;itr++)
	{
		for(i=0;i<ncor;i++)
		{
				if(i<ni)
					fwrite(&ric1[i],sizeof(float),1,fpe1);
				else
					fwrite(&ric1[0],sizeof(float),1,fpe1);
		}
	}
	for(itr=0;itr<ntr;itr++)
	{
		for(i=0;i<ncor;i++)
			fwrite(&fmat[itr*ncor+i],sizeof(float),1,fpe1);
	}
	SpectrumDivide(fmat,ric2,fmat,ntr,ncor,ni,1e-3);
	//			input matrix( multi traces)
						// remove ricker( single trace )
							//output matrix( )
								// trace number 
									 // trace length -- nt
									 	  // ricker nt
									 	  	// whitecoe
	for(itr=0;itr<ntr;itr++)
	{
		for(i=0;i<ncor;i++)
			fwrite(&fmat[itr*ncor+i],sizeof(float),1,fpe1);
	}
	fclose(fpe1);

	free(f);
	free(fmat);
	free(fcor);
	free(ric1);
	free(ric2);
	return(0);
}*/
//*din---Input,*dout---Output,*ntr---trace number,nt---time number,angle---0-360
void phase_correction(float *din,float *dout,int ntr, int nt,float angle)
{
	int ntfft;
	int ntffth;

	int itr;
	int i;
	int k;

	double theta;
	double cost;
	double sint;

	double *pr;		
	double *pi;
	double *fr;		
	double *fi;

	theta=angle*pai/180;
	ntfft=2;
	k=1;
	while(ntfft<nt)
	{
		ntfft*=2;
		k++;
	}

	cost=cos(theta);
	sint=sin(theta);
	ntffth=ntfft/2;

	pr=(double *)malloc(sizeof(double)*(ntfft));
	pi=(double *)malloc(sizeof(double)*(ntfft));
	fr=(double *)malloc(sizeof(double)*(ntfft));
	fi=(double *)malloc(sizeof(double)*(ntfft));

	for(itr=0;itr<ntr;itr++)
	{
		for (i=0;i<nt;i++)
		{
			pr[i]=din[itr*nt+i];
			pi[i]=0.0;
		}
		for (i=nt;i<ntfft;i++)
		{
			pr[i]=0.0;
			pi[i]=0.0;
		}
		kbfft(pr,pi,ntfft,k,fr,fi,0,1);
		for(i=1;i<ntffth;i++)
		{
			pr[i]=fr[i]*cost-fi[i]*sint;
			pi[i]=fi[i]*cost+fr[i]*sint;
		}
		for(i=ntfft-1;i>ntffth;i--)
		{
			pr[i]=pr[ntfft-i];
			pi[i]=-pi[ntfft-i];
		}
		kbfft(pr,pi,ntfft,k,fr,fi,1,0);
		for(i=0;i<nt;i++)
		{
			dout[itr*nt+i]=fr[i];
		}
	}

	free(pr);
	free(pi);
	free(fr);
	free(fi);
}

//*din---input matrix( multi traces), din2-----remove ricker( single trace ), *dout----output matrix
// ntr---trace number, nt--trace length,nt2---ricker nt,whitecoe----1.0e-3,1.0e-4
void SpectrumDivide(float *din,float *din2,float *dout,int ntr, int nt,int nt2,float whitecoe)
{
	int ntfft;
	int ntffth;

	int itr;
	int i;
	int k;
	int ntmax;

	double amax;

	double *pr;		
	double *pr2;		
	double *pi;
	double *fr;		
	double *fi;

	if(nt>=nt2)
		ntmax=nt;
	else		
		ntmax=nt2;

	ntfft=2;
	k=1;
	while(ntfft<ntmax)
	{
		ntfft*=2;
		k++;
	}

	ntffth=ntfft/2;

	pr=(double *)malloc(sizeof(double)*(ntfft));
	pi=(double *)malloc(sizeof(double)*(ntfft));

	pr2=(double *)malloc(sizeof(double)*(ntfft));

	fr=(double *)malloc(sizeof(double)*(ntfft));
	fi=(double *)malloc(sizeof(double)*(ntfft));

	for (i=0;i<nt2;i++)
	{
		pr2[i]=din2[i];
		pi[i]=0.0;
	}
	for (i=nt2;i<ntfft;i++)
	{
		pr2[i]=0.0;
		pi[i]=0.0;
	}
	kbfft(pr2,pi,ntfft,k,fr,fi,0,1);

	for (amax=0.0,i=0;i<ntffth;i++)
		if(pr2[i]>amax) amax=pr2[i];
	amax=whitecoe*amax;
	printf("amax=%f\n",amax);
	for(i=0;i<ntffth;i++)
		pr2[i]=1.0/(pr2[i]+amax);

	for(itr=0;itr<ntr;itr++)
	{
		for (i=0;i<nt;i++)
		{
			pr[i]=din[itr*nt+i];
			pi[i]=0.0;
		}
		for (i=nt;i<ntfft;i++)
		{
			pr[i]=0.0;
			pi[i]=0.0;
		}
		kbfft(pr,pi,ntfft,k,fr,fi,0,1);

		for(i=1;i<ntffth;i++)
		{
			pr[i]=fr[i]*pr2[i];
			pi[i]=fi[i]*pr2[i];
		}
		for(i=ntfft-1;i>ntffth;i--)
		{
			pr[i]=pr[ntfft-i];
			pi[i]=-pi[ntfft-i];
		}
		kbfft(pr,pi,ntfft,k,fr,fi,1,0);
		for(i=0;i<nt;i++)
		{
			dout[itr*nt+i]=fr[i]/ntfft*25;
		}
	}

	free(pr);
	free(pr2);
	free(pi);
	free(fr);
	free(fi);
}

void kbfft(double *pr,double *pi,int n,int k,double *fr,double *fi,int l,int il)
{ 
	int it,m,is,i,j,nv,l0;
	double p,q,s,vr,vi,poddr,poddi;

	//排序
	for (it=0; it<=n-1; it++) 
	{ m=it; is=0;
	for (i=0; i<=k-1; i++)
	{ j=m/2; is=2*is+(m-2*j); m=j;
	fr[it]=pr[is]; fi[it]=pi[is];
	}
	}

	//蝶形运算
	pr[0]=1.0; pi[0]=0.0;
	p=6.283185306/(1.0*n);
	pr[1]=cos(p); pi[1]=-sin(p);
	if (l!=0) pi[1]=-pi[1];
	for (i=2; i<=n-1; i++)
	{ p=pr[i-1]*pr[1]; q=pi[i-1]*pi[1];
	s=(pr[i-1]+pi[i-1])*(pr[1]+pi[1]);
	pr[i]=p-q; pi[i]=s-p-q;
	}
	for (it=0; it<=n-2; it=it+2)
	{ vr=fr[it]; vi=fi[it];
	fr[it]=vr+fr[it+1]; fi[it]=vi+fi[it+1];
	fr[it+1]=vr-fr[it+1]; fi[it+1]=vi-fi[it+1];
	}
	m=n/2; nv=2;
	for (l0=k-2; l0>=0; l0--)
	{ m=m/2; nv=2*nv;
	for (it=0; it<=(m-1)*nv; it=it+nv)
		for (j=0; j<=(nv/2)-1; j++)
		{ p=pr[m*j]*fr[it+j+nv/2];
	q=pi[m*j]*fi[it+j+nv/2];
	s=pr[m*j]+pi[m*j];
	s=s*(fr[it+j+nv/2]+fi[it+j+nv/2]);
	poddr=p-q; poddi=s-p-q;
	fr[it+j+nv/2]=fr[it+j]-poddr;
	fi[it+j+nv/2]=fi[it+j]-poddi;
	fr[it+j]=fr[it+j]+poddr;
	fi[it+j]=fi[it+j]+poddi;
	}
	}
	if (l!=0)
		for (i=0; i<=n-1; i++)
		{ fr[i]=fr[i]/(1.0*n);
	fi[i]=fi[i]/(1.0*n);
	}
	if (il!=0)
		for (i=0; i<=n-1; i++)
		{ pr[i]=sqrt(fr[i]*fr[i]+fi[i]*fi[i]);
	pr[i]=(pr[i]/(n/2));						//各次谐波幅值，其中pr[1]为基波幅值
	if (fabs(fr[i])<0.000001*fabs(fi[i]))		//fabs()是取绝对值函数,浮点型的0 在内存中并不是严格等于0,可以认为当一个浮点数离原点足够近时,也就是f>0.00001 && f<-0.00001,认为f是0
	{ if ((fi[i]*fr[i])>0) pi[i]=90.0;
	else pi[i]=-90.0;
	}
	else
		pi[i]=atan(fi[i]/fr[i]);
	}
	return;
}
float Ricker(float t1,float f0)           //ricker source
{
	float t00=1/f0,y;
	y=(1-2*(float)pow(pai*f0*(t1-t00),2))*(float)exp(-(float)pow(pai*f0*(t1-t00),2)); 
	return(y);
}



