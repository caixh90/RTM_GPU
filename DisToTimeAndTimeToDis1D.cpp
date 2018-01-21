
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<malloc.h>


void interpolate(float *output,float *input_x,float *input_y,int N_out,int N_in,float dt)
{
	int i,j;
	int k1=0,k2=0;
	output[0]=input_y[0];
	for(i=1; i<N_in; i++)
	{
		k2=(int) (input_x[i]/dt);
		for(j=k1+1; j<=k2; j++)
		{
			if(j==N_out)
			{
				//printf("j=%d N_out=%d!!!!!!!!!",j, N_out);
				continue;
				//getchar();
			}
			output[j]=input_y[i-1]+(j*dt-input_x[i-1])/(input_x[i]-input_x[i-1])*(input_y[i]-input_y[i-1]);
		}
		k1=k2;
	}
}

//dz is dt; dt is dz
int T2D(const char *File_Name, float *V, float *D, int Nx, int Nz, int Nz_V, int Nx_d_st, int Nx_d_end, float dx, float dz, float dt)
{
	float temp=0, temp_max=0;
	int i, j, j_temp;
	int Nx_d=Nx_d_end-Nx_d_st+1;
	int Nt;
	float *t0;
	float *T;
	FILE *out;
	t0=(float *)malloc(Nx_d*Nz*sizeof(float ));
	
	for(i=0;i<Nx_d;i++)
	{
		for (j=1; j<Nz; j++)
		{
			t0[i*Nz+j]=0;
		}
	}
	
	for (i=Nx_d_st; i<=Nx_d_end; i++)
	{
		temp=0;
		for (j=1; j<Nz; j++)
		{
			/*if(i==178&&j==2405)
			{
				getchar();
			}*/
			/*if(j%500==0)
			{
				printf("%f\n",temp);
				getchar();
			}*/
			j_temp=(int)temp/dt;
			if(j_temp>=Nz_V)
			{
				j_temp=Nz_V-1;
			}
			temp+=dz*V[i*Nz_V+j_temp]/2;
			/*if(abs(temp)>6250)
			{
				printf("!!!!!!!\n");
				getchar();
			}*/
			t0[(i-Nx_d_st)*Nz+j]=temp;
		}
		if(temp_max<temp)
		{
			temp_max=temp;
		}
	}
	Nt=(int) (temp_max/dt);

	T=(float *)malloc(Nx_d*Nt*sizeof(float ));
	
	for(i=0;i<Nx_d;i++)
	{
		for (j=1; j<Nt; j++)
		{
			T[i*Nt+j]=0;
		}
	}

	for (i=0; i<Nx_d; i++)	
	{		
		interpolate(&T[i*Nt], &t0[i*Nz], &D[i*Nz], Nt, Nz, dt);
	}

	out=fopen(File_Name,"wb");
	for(i=0;i<Nx_d;i++)
	{
		for(j=0;j<Nt;j++)
		{			
			fwrite(&T[i*Nt+j],sizeof(float),1,out);
		}
	}

	free(T);
	free(t0);

	fclose(out);
	return Nt;
}

int D2T(const char *File_Name, float *V, float *D, int Nx, int Nz, int Nx_d_st, int Nx_d_end, float dx, float dz, float dt)
{
	float temp=0, temp_max=0;
	int i, j;
	int Nx_d=Nx_d_end-Nx_d_st+1;
	int Nt;
	float *t0;
	float *T;
	FILE *out;
	t0=(float *)malloc(Nz*Nx_d*sizeof(float ));
	
	
	for(i=0;i<Nx_d;i++)
	{
		for (j=1; j<Nz; j++)
		{
			t0[i*Nz+j]=0;
		}
	}
	
	for (i=Nx_d_st; i<=Nx_d_end; i++)
	{
		temp=0;
		for (j=1; j<Nz; j++)
		{
			temp+=2*dz/V[i*Nz+j];
			t0[(i-Nx_d_st)*Nz+j]=temp;
		}
		if(temp_max<temp)
		{
			temp_max=temp;
		}
	}
	Nt=(int) (temp_max/dt);
	
	T=(float *)malloc(Nx_d*Nt*sizeof(float ));
	
	
	for(i=0;i<Nx_d;i++)
	{
		for (j=1; j<Nt; j++)
		{
			T[i*Nt+j]=0;
		}
	}

	for (i=0; i<Nx_d; i++)
	{		
		interpolate(&T[i*Nt], &t0[i*Nz], &D[i*Nz], Nt, Nz, dt);
	}

	out=fopen(File_Name,"wb");
	for(i=0;i<Nx_d;i++)
	{
		for(j=0;j<Nt;j++)
		{			
			fwrite(&T[i*Nt+j],sizeof(float),1,out);
		}
	}

	free(T);
	free(t0);

	fclose(out);
	return Nt;
}
/*
int main()
{
	int i,j;
	int NX,NZ; 
	NX=202;NZ=415;
	float dt,dz;
	dt=0.001;dz=10;
	int nt;
	float *record,*V,*record1;
	FILE *Rdrecord,*RdV;
	record=(float *)malloc(NX*NZ*sizeof(float ));
	V=(float *)malloc(NX*NZ*sizeof(float ));
	
	Rdrecord=fopen("RVSP_Migration_f045Hz.dat","rb");
	
	for(i=0;i<NX;i++)
	{	
		for(j=0;j<NZ;j++)
		{		
			fread(&record[i*NZ+j],sizeof(float),1,Rdrecord);
		}
	
	}
	fclose(Rdrecord);
	
	RdV=fopen("L10-1932-X-Vp-202.dat","rb");
	
	for(i=0;i<NX;i++)
	{
		for(j=0;j<NZ;j++)
		{		
			fread(&V[i*NZ+j],sizeof(float),1,RdV);
		}
	}
	fclose(RdV);
	
	nt=D2T("RVSP_Migration_f045Hz_T.dat", V, record, NX, NZ, 0, (NX-1), dz, dz, dt);
	printf("%d\n",nt);



	record1=(float *)malloc(nt*NX*sizeof(float ));
	
	Rdrecord=fopen("RVSP_Migration_f045Hz_T.dat","rb");
	
	for(i=0;i<NX;i++)
	{	
		for(j=0;j<nt;j++)
		{		
			fread(&record1[i*nt+j],sizeof(float),1,Rdrecord);
		}
	
	}
	fclose(Rdrecord);
	nt=T2D("RVSP_Migration_f045Hz_D.dat", V, record1, NX, nt, NZ, 0, (NX-1), dz, dt, dz);
	printf("%d\n",nt);
}*/
