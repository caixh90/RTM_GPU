#include<stdio.h>
#include<math.h>
#include<malloc.h>
#define pi 3.1415926535898

void callenfd2d_ls(double *hk,double *disperr,int nfre,int nfdmin,int nfdmax,double vel,
                   double tao,double h,double df,int nthita,double eps,int *lenfd,int iLSTE,double hzx);

void calfdlen_ls(double *hk, double *disperr, int nfre, int nfdmax,int nfdmin, double tao, double h,
                  double df,int nthita,double eps, int *lenvar, int nvel, double vmin, double dv,int iok,int iLSTE,double hzx);

int Gauss(double **A,int n,double *b,double *x);

void CAL2DFDCOE_LSM(double *c,double r,double bb,int M,double hzx);

double fgaus(int n,int js[],int ai,int aj,double bb,double ar,int al,double hzx);

double fgausf(int n,double x[],int ai,int aj,double ar,int al,double hzx);

void fgauss(int j,int n,double bb,double x[],double y[]);

int funMandC(int nthita,int nfdmax,int nfdmin,int nvel,double tao,double h,double df,double eps,
			  double fmax,double vmin,double vmax,double dv,int *fdcoeneed,int *M, int *Index,int iLSTE,float **cc,double hzx)
{
	 int i,l,iok,nfre,NC;
	 double b,rat,ra,r;
	 double *c,*hk,*disperr;
	 
	 nfre=(int)(fmax/df)+1;
	 iok=1;

	 c=(double *)malloc(sizeof(double)*(nfdmax+1));
	 hk=(double *)malloc(sizeof(double)*nfre);
	 disperr=(double *)malloc(sizeof(double)*nfre);
	 
	 calfdlen_ls(hk,disperr,nfre,nfdmax,nfdmin,tao,h,df,nthita,eps,M,nvel,vmin,dv,iok,iLSTE,hzx);
	 
	 if(iok==0) printf("ERROR in calculate M");
	 
	 for (i=0;i<nvel;i++)
	 {
		 if(fdcoeneed[i]==0) M[i]=-1;
	 }
	 Index[0]=0;
	 for (i=1;i<=nvel;i++)
	 {
		 Index[i]=Index[i-1]+M[i-1]+1;
	 }
	 NC=Index[nvel];
	 *cc=(float *)malloc(sizeof(float)*NC);
	 
	 rat=tao/h;
	 ra=2.0*pi*fmax*tao;
	 for (i=0;i<nvel;i++)
	 {	
		 if(fdcoeneed[i]==1)
		 {
			 r=(vmin+i)*rat;
			 b=ra/r;
			 CAL2DFDCOE_LSM(c,r,b,M[i],hzx);
			
			 for (l=0;l<=M[i];l++)
			 {
				 (*cc)[l+Index[i]]=(float)c[l];
			 }
		 }
	 }
	
	 free(c);
	 free(hk);
	 free(disperr);
	 return (NC);
}
int Gauss(double **A,int n,double *b,double *x)  //高斯列主元消去法求解线性方程：Ax=b，A为n阶方阵。
{
	int i,j,k,p;
	double m,max1;
    for(i=0;i<n;i++)
	{
		max1=A[i][0];
		for(j=0;j<n;j++)
		{
			if(fabs(A[i][j])>fabs(max1))
				max1=A[i][j];
		}
		if(fabs(max1)<1e-10)
		{
			printf("error:A is not full");
            return(0);
			break;
		}
		for(j=0;j<n;j++)
		{
			A[i][j]=A[i][j]/max1;
		}
		b[i]=b[i]/max1;
	}
	for(i=0;i<n-1;i++)
	{
		k=i;
		for(j=i;j<n;j++)
		{
			if(fabs(A[j][i])>fabs(A[k][i]))
				k=j;
		}
		if(k!=i)
		{
			m=b[i];
            b[i]=b[k];
			b[k]=m;
			for(j=i;j<n;j++)
			{
				m=A[i][j];
                A[i][j]=A[k][j];
				A[k][j]=m;
			}
		}
        for(p=i+1;p<n;p++)
		{
            m=A[p][i]/A[i][i];
			b[p]=b[p]-m*b[i];
			for(j=i;j<n;j++)
			{
                A[p][j]=A[p][j]-m*A[i][j];
			}
		}
	}
	x[n-1]=b[n-1]/A[n-1][n-1];
	for(i=n-2;i>=0;i--)
	{
        m=0.0;
		for(j=i+1;j<n;j++)
            m=m+A[i][j]*x[j];
		x[i]=(b[i]-m)/A[i][i];
	}
	return(1);
}

double fgaus(int n,int js[],int ai,int aj,double bb,double ar,int al,double hzx)
{ 
	int m,j,k,q,l,*is;
	double y[2],p,s,*x,*a,*b;
	static double t[5]={-0.9061798459,-0.5384693101,0.0,
	0.5384693101,0.9061798459};
	static double c[5]={0.2369268851,0.4786286705,0.5688888889,
	0.4786286705,0.2369268851};
	is=(int *)malloc(2*(n+1)*sizeof(int));
	x=(double *)malloc(n*sizeof(double));
	a=(double *)malloc(2*(n+1)*sizeof(double));
	b=(double *)malloc((n+1)*sizeof(double));
	m=1; l=1;
	a[n]=1.0; a[2*n+1]=1.0;
	while (l==1)
	{ 
		for (j=m;j<=n;j++)
		{
			fgauss(j-1,n,bb,x,y);
			a[j-1]=0.5*(y[1]-y[0])/js[j-1];
			b[j-1]=a[j-1]+y[0];
			x[j-1]=a[j-1]*t[0]+b[j-1];
			a[n+j]=0.0;
			is[j-1]=1; is[n+j]=1;
		}
		j=n; q=1;
		while (q==1)
		{ 
			k=is[j-1];
			if (j==n) p=fgausf(n,x,ai,aj,ar,al,hzx);
			else p=1.0;
			a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];
			is[j-1]=is[j-1]+1;
			if (is[j-1]>5)
			if (is[n+j]>=js[j-1])
			{
				j=j-1; q=1;
				if (j==0)
				{ 
					s=a[n+1]*a[0]; free(is); free(x);
					free(a); free(b); return(s);
				}
			}
			else
			{ 
				is[n+j]=is[n+j]+1;
				b[j-1]=b[j-1]+a[j-1]*2.0;
				is[j-1]=1; k=is[j-1];
				x[j-1]=a[j-1]*t[k-1]+b[j-1];
				if (j==n) q=1;
				else q=0;
			}
			else
			{
				k=is[j-1];
				x[j-1]=a[j-1]*t[k-1]+b[j-1];
				if (j==n) q=1;
				else q=0;
			}
		}
		m=j+1;
	}
	return(1.0);
}

void fgauss(int j,int n,double bb,double x[],double y[])
{ 
	n=n;
	switch (j)
	{
		case 0: { y[0]=0.0; y[1]=2*pi; break;}
		case 1: { y[0]=0.0; y[1]=bb; break;}
		default: { }
	}
	return;
}

double fgausf(int n,double x[],int ai,int aj,double ar,int al,double hzx)
{ 
	double f1,f2,f,c,s,p,car,hzx2;
	hzx2=hzx*hzx;
	c=cos(x[0]);
	s=sin(x[0]);
	p=pow(ar,-2);
	car=cos(ar*x[1]);
	f1=((1+1/hzx2)-cos(ai*x[1]*c)-cos(hzx*ai*x[1]*s)/hzx2)/(p*(1-car));	
	if(al==1)
	{
		f2=((1+1/hzx2)-cos(aj*x[1]*c)-cos(hzx*aj*x[1]*s)/hzx2)/(p*(1-car));
		f=f1*f2;
		return(f);
	}
	else
	    return(f1);
 }


void CAL2DFDCOE_LSM(double *c,double r,double bb,int M,double hzx)		//c-要求解的差分系数，r-网比，fmax-最大频率，tao-时间步长，N-速度的变化范围（以1为单位），M-算子长度
{	
	int i,j;
	int js[2]={4,4};
	double **A,*x,*b,c0;												//b代表求解矩阵时的Ax=b中的b, bb代表kh的乘积
	
	b=(double *)malloc(sizeof(double )*(M));
	x=(double *)malloc(sizeof(double )*(M));
	A=(double **)malloc(sizeof(double *)*(M));
	for(i=0;i<M;i++)
	{
		A[i]=(double *)malloc(sizeof(double)*(M));
	}


	for(i=0;i<M;i++)
	{
		for(j=i;j<M;j++)
		{
			A[i][j]=fgaus(2,js,i+1,j+1,bb,r,1,hzx);
		}
		b[i]=fgaus(2,js,i+1,j+1,bb,r,2,hzx);
	}
	for(i=0;i<M;i++)
	{
		for(j=0;j<i;j++)
		{
			A[i][j]=A[j][i];
		}
	}

	
	Gauss(A,M,b,x);
	
	for(i=1;i<=M;i++)
	{
		c[i]=x[i-1];
    }
	c0=0.0;
	for(i=1;i<=M;i++)
	{
		c0+=c[i];
	}
	c[0]=-2.0*c0;
	for(i=0;i<M;i++)
	{
		free(A[i]);
	}
	free(A);
	free(x);
	free(b);
}

/********************N1——M的最大值*******************/
/********************nthita——角度采样个数*******************/
/********************nfre——最大频率的整数倍*******************/

void callenfd2d_ls(double *hk,double *disperr,int nfre,int nfdmin,int nfdmax,double vel,
                   double tao,double h,double df,int nthita,double eps,int *lenfd,int iLSTE,double hzx)
{
	int i,j,k,n,l;
	double *c,r,b,mid,rat,ra,nth,hzx2;
	hzx2=hzx*hzx;

	b=2.0*pi*df*h/vel;
	for(i=0;i<nfre;i++)
	{
	  hk[i]=b*i;
	}

	c=(double *)malloc(sizeof(double)*(nfdmax+1));
	
	for(i=0;i<=nfdmax;i++)
	{
		c[i]=0;
	}

	r=vel*tao/h;
	b=2*pi*df*(nfre-1)*tao/r;
	ra=h/vel;

	for (j=nfdmin;j<=nfdmax;j++)
	{	
		CAL2DFDCOE_LSM(c,r,b,j,hzx);	
		
		for (k=1;k<nfre;k++)    //对频率的循环
		{	
			rat=2/(r*hk[k]);
			for (n=0;n<=nthita;n++)	//对方位角的循环
			{
				nth=n*pi/(4*nthita);
				mid=0;
				for(l=1;l<=j;l++)
				{
					mid=mid+c[l]*(pow(sin(l*hzx*hk[k]*sin(nth)/2),2)/hzx2+pow(sin(l*hk[k]*cos(nth)/2),2));
				}
				disperr[k]=rat*asin(sqrt(r*r*mid));
				disperr[k]=fabs(ra*(1.0/disperr[k]-1.0));
				
				if(disperr[k]>eps) break;							
			}
			if(n<nthita) break;
		}		
		*lenfd=j;		
		nfdmin=j;
		
		if(k==nfre) break;
	}
	if(j==nfdmax+1) 	printf("M=%d is not enough",j);
	free(c);
}
void calfdlen_ls(double *hk, double *disperr, int nfre, int nfdmax,int nfdmin, double tao, double h,
				 double df,int nthita,double eps, int *lenvar, int nvel, double vmin, double dv,int iok,int iLSTE,double hzx)
 {
//    double hk[nfre],disperr[nfre];
//    int lenvar[nvel];
      int nsearchtime,num,iflag,ibeg,iend,lenfd,lbeg,lend,iinc;
	  int i,i1,i2,i3,i4,l1,l2,l3,l4;
	  double rat,vel;
      nsearchtime=0;
      num=0;
      iflag=0;
      ibeg=1;
      iend=nvel;

      i1=ibeg;            
      i=i1;
      vel=vmin+(i-1)*dv;
      rat=vel*tao/h;
      num++;


      callenfd2d_ls(hk,disperr,nfre,nfdmin,nfdmax,vel,tao,h,df,nthita,eps,&lenfd,iLSTE,hzx);

      nsearchtime=nsearchtime+lenfd-2+1;
      if(lenfd==0)
	  {
          printf("error===> fremain or eps is too lagre \n");
          printf("      or  nfdmax is too lagre !\n");
          goto a999;
	  }
      lbeg=lenfd;
      l1=lenfd;
      lenvar[i1-1]=l1;
      printf("%d, %d\n",i1,lenvar[i1-1]);

      i2=iend;            
      i=iend;
      vel=vmin+(i-1)*dv;
      rat=vel*tao/h;
      num++;
      callenfd2d_ls(hk,disperr,nfre,nfdmin,nfdmax,vel,tao,h,df,nthita,eps,&lenfd,iLSTE,hzx);

      nsearchtime+=lenfd-2+1;
      if(lenfd==0) 
	  {
         printf("error===> fremain or eps is too lagre !\n");
         printf("      or  nfdmax is too lagre ! \n");
         goto a999;
      }
      lend=lenfd;
      l2=lenfd;
      lenvar[i2-1]=l2;
      if(l1==lend)
	  {
         for(i=ibeg+1;i<=iend;i++) lenvar[i-1]=l1;
   //      printf("%d,%d,%d,%d\n",ibeg+1,l1,iend,l1);
         goto a300;
      }
      iinc=(iend-ibeg)/(2*(l1-lend));
      if(iinc<1) iinc=1;

a100: 
      i1=i2-iinc;
      if(i1<ibeg) 
	  {
         i1=ibeg;
         l1=lbeg;
         if(i1==i2) goto a300;
	  }
      else
	  {
         i=i1;
         vel=vmin+(i-1)*dv;
         rat=vel*tao/h;
         num++;
         callenfd2d_ls(hk,disperr,nfre,l2,nfdmax,vel,tao,h,df,nthita,eps,&lenfd,iLSTE,hzx);
         nsearchtime+=lenfd-l2+1;
         if(lenfd==0) 
		 {
             printf("error===> fremain or eps is too lagre\n");
             printf("      or  nfdmax is too lagre !\n");
             goto a999;
		 }
         l1=lenfd;
      }


a200:
	  if(l1==l2) 
	  {
          for(i=i1;i<=i2-1;i++)
              lenvar[i-1]=l1;
    //      if(i1<i2) printf("%d,%d,%d,%d\n",i1,lenvar[i1-1],i2-1,lenvar[i2-2]);
          i2=i1;
          l2=l1;
	  }
      else if(i1==i2-1) 
	  {
          lenvar[i1-1]=l1;
    //     printf("%d,%d\n",i1,lenvar[i1-1]);
          i2=i1;
          l2=l1;
	  }
      else
	  {
          i4=(i1+i2)/2;
          if(iflag==0) 
		  {
             iflag=1;
             i3=i1;
             l3=l1;
		  }
          i=i4;
          vel=vmin+(i-1)*dv;
          rat=vel*tao/h;
          num++;
          callenfd2d_ls(hk,disperr,nfre,l2,nfdmax,vel,tao,h,df,nthita,eps,&lenfd,iLSTE,hzx);
          nsearchtime+=lenfd-l2+1;
          if(lenfd==0) 
		  {
              printf("error===> fremain or eps is too lagre \n");
              printf("      or  nfdmax is too lagre !\n");
              goto a999;
          }
          l4=lenfd;
    //      printf("%d,%d,%d,%d,%d,%d\n",i1,l1,i4,l4,i2,l2);
          if(l4==l2) 
		  {
		      for(i=i4;i<=i2-1;i++)
                  lenvar[i-1]=l2;
    //         if(i4<i2) printf("%d,%d,%d,%d,%d\n",i4,lenvar[i4-1],i2-1,lenvar[i2-2]);
              i2=i4;
              l2=l4;
		  }
          else
		  {
              i1=i4;
              l1=l4;
          }
          goto a200;
	  }

      if(iflag==1) 
	  {
          iflag=0;
          if(l3==l2)
		  {
              for(i=i3;i<=i2-1;i++)
                  lenvar[i-1]=l2;
    //          if(i3<i2) printf("%d,%d,%d,%d\n",i3,lenvar[i3-1],i2-1,lenvar[i2-2]);
              i2=i3;
              l2=l3;
              goto a100;
		  }
          else
		  {
              i1=i3;
              l1=l3;
              goto a200;
		  }
	  }
      else if(l2==lbeg) 
	  {
          for(i=ibeg;i<=i2-1;i++) lenvar[i-1]=l2;
   //       if(ibeg<i2) printf("%d,%d,%d,%d\n",ibeg,lenvar[ibeg-1],i2-1,lenvar[i2-2]);
	  }
      else
          goto a100;
a300:     
	  iok=1;
   //   printf("effectiveness =%d / %d\n ",num,nvel);
      printf("effectiveness =%f\n ",((float)num)/nvel);
      printf("total search times for fd length : %d\n",nsearchtime);
      return;
a999:
	  iok=0;
      return;
}

void order(int N,float *c)	// any order difference coefficient
{
	int i,j;
	float *x;
	x=(float *)malloc(sizeof(float)*(N/2+1));
	c[0]=0.0;
	for(i=0;i<=N/2;i++)
	{
		x[i]=1.0;
	}

	for(i=1;i<=N/2;i++)
	{
		for(j=1;j<=N/2;j++)
		{
			if(j!=i)
			{
				x[i]=x[i]*fabs(pow(j,2)/(pow(j,2)-pow(i,2)));
			}
		}

			c[i]=pow(-1,i+1)/pow(i,2)*x[i];
			c[0]=(c[0]-2*c[i]);
	}
	free(x);
}
