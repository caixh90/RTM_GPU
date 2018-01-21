#include<stdio.h>
#include<math.h>
#include<malloc.h>

/**************************速度文件**************************/
void velocity(char *OutNameVp,float *v,float *v_2,float *r,float *r_1,float *r_2,int NZ,int NX,int NXZ,int N2,float tao,float h,int choice) 
{
	int i,j;
	FILE *RDV;
	
	RDV=fopen(OutNameVp,"rb");
	for(j=N2;j<NX-N2;j++)
	{
		for (i=N2;i<NZ-N2;i++)
		{
			fread(&v[i*NX+j],sizeof(float),1,RDV);
		}
	}
	fclose(RDV);
	//上边界
	for(i=0;i<N2;i++)
	{
		for (j=N2;j<NX-N2;j++)
		{
			v[i*NX+j]=v[N2*NX+j];
		}
	}
	//左边界
	for(j=0;j<N2;j++)
	{
		for (i=N2;i<NZ-N2;i++)
		{
			v[i*NX+j]=v[i*NX+N2];
		}
	}
	//下边界
	for(i=NZ-N2;i<NZ;i++)
	{
		for (j=N2;j<NX-N2;j++)
		{
			v[i*NX+j]=v[(NZ-N2-1)*NX+j];
		}
	}
	//右边界
	for(j=NX-N2;j<NX;j++)
	{
		for (i=N2;i<NZ-N2;i++)
		{
			v[i*NX+j]=v[i*NX+(NX-N2-1)];	
		}
	}
	//左上角
	for(i=0;i<N2;i++)
	{
		for (j=0;j<N2;j++)
		{
			v[i*NX+j]=v[N2*NX+N2];
		}
	}
	//右上角
	for(i=0;i<N2;i++)
	{
		for (j=NX-N2;j<NX;j++)
		{
			v[i*NX+j]=v[N2*NX+(NX-N2-1)];
		}
	}
	//左下角
	for(i=NZ-N2;i<NZ;i++)
	{
		for (j=0;j<N2;j++)
		{
			v[i*NX+j]=v[(NZ-N2-1)*NX+N2];
		}
	}
	//右下角
	for(i=NZ-N2;i<NZ;i++)
	{
		for (j=NX-N2;j<NX;j++)
		{
			v[i*NX+j]=v[(NZ-N2-1)*NX+(NX-N2-1)];
		}
	}
	if(choice==1)
	{
		for(i=0;i<NZ;i++)
		{	
			for(j=0;j<NX;j++)
			{
				v_2[i*NX+j]=v[i*NX+NX-1-j];
			}
		}
		for(i=0;i<NZ;i++)
		{	
			for(j=0;j<NX;j++)
			{
				v[i*NX+j]=v_2[i*NX+j];
			}
		}
	}
	
	for(i=0;i<NXZ;i++)
	{
		v_2[i]=pow(v[i],2);		
	}
	for(i=0;i<NXZ;i++)
	{
		r[i]=v[i]*tao/h;	
	}	
	for(i=0;i<NXZ;i++)
	{
		r_2[i]=pow(r[i],2)/2;	
	}
	for(i=0;i<NXZ;i++)
	{
		r_1[i]=sqrt(r_2[i]);
	}
}
