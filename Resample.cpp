#include<stdio.h>
#include<math.h>
#include<malloc.h>

#define LTABLE 8
#define NTABLE 513
#define PI 3.1415926535898

void intt8r (int ntable, float table[][8],
	int nxin, float dxin, float fxin, float yin[], float yinl, float yinr,
	int nxout, float dxout, float yout[]);
void resample (int nxin, float dxin, float yin[], 
	         int nxout, float dxout, float yout[]);
/*void main()
{
    float fm;
    float att;
    int nxin;
    float dxin;
    float *yin;
    int nxout;
    float dxout;
    float *yout;
    int i;
	float t;
    FILE *fdin,*fdout;

    fm=50;
    att=20;
    nxin=101;
    dxin=1.0;
    dxout=0.33;

    nxout=int((nxin-1)*dxin/dxout+1.5);

    yin=(float *)malloc(sizeof(float)*nxin);
    yout=(float *)malloc(sizeof(float)*nxout);

    for(i=0;i<nxin;i++)
	{
	   t=i*dxin*0.001;
	   yin[i]=sin(2*PI*fm*t)*exp(-att*att*t*t);
	}

     resample (nxin, dxin, yin, nxout, dxout,yout);


	 fdin=fopen("data_in.txt","w");
	 for (i=0;i<nxin;i++)
	 {
	    t=i*dxin*0.001;
		fprintf(fdin,"%f %f \n",t, yin[i]);
	 }
	 fclose(fdin);

	 fdout=fopen("data_out.txt","w");
	 for (i=0;i<nxout;i++)
	 {
	    t=i*dxout*0.001;
		fprintf(fdout,"%f %f \n",t, yout[i]);
	 }
	 fclose(fdout);

	 free(yin);
	 free(yout);

}*/

void intt8r (int ntable, float table[][8],
	int nxin, float dxin, float fxin, float yin[], float yinl, float yinr,
	int nxout, float dxout, float yout[])
{
	int ioutb,nxinm8,ixout,ixoutn,kyin,ktable,itable;
	float xoutb,xoutf,xouts,xoutn,frac,fntablem1,yini,sum,
		*yin0,*table00,*pyin,*ptable;
	float xout;

	/* compute constants */
	ioutb = -3-8;
	xoutf = fxin;
	xouts = 1.0/dxin;
	xoutb = 8.0-xoutf*xouts;
	fntablem1 = (float)(ntable-1);
	nxinm8 = nxin-8;
	yin0 = &yin[0];
	table00 = &table[0][0];

	/* loop over output samples */
	for (ixout=0; ixout<nxout; ixout++) {

		/* determine pointers into table and yin */
		xout=ixout*dxout;
		xoutn = xoutb+xout*xouts;
		ixoutn = (int)xoutn;
		kyin = ioutb+ixoutn;
		pyin = yin0+kyin;
		frac = xoutn-(float)ixoutn;
		ktable = frac>=0.0?frac*fntablem1+0.5:(frac+1.0)*fntablem1-0.5;
		ptable = table00+ktable*8;
		
		/* if totally within input array, use fast method */
		if (kyin>=0 && kyin<=nxinm8) {
			yout[ixout] = 
				pyin[0]*ptable[0]+
				pyin[1]*ptable[1]+
				pyin[2]*ptable[2]+
				pyin[3]*ptable[3]+
				pyin[4]*ptable[4]+
				pyin[5]*ptable[5]+
				pyin[6]*ptable[6]+
				pyin[7]*ptable[7];
		
		/* else handle end effects with care */
		} else {
	
			/* sum over 8 tabulated coefficients */
			for (itable=0,sum=0.0; itable<8; itable++,kyin++) {
				if (kyin<0)
					yini = yinl;
				else if (kyin>=nxin)
					yini = yinr;
				else
					yini = yin[kyin];
				sum += yini*(*ptable++);
			}
			yout[ixout] = sum;
		}
	}
}
void stoepd (int n, double r[], double g[], double f[], double a[])
{
	int i,j;
	double v,e,c,w,bot;

	if (r[0] == 0.0) return;

	a[0] = 1.0;
	v = r[0];
	f[0] = g[0]/r[0];

	for (j=1; j<n; j++) {
		
		/* solve Ra=v as in Claerbout, FGDP, p. 57 */
		a[j] = 0.0;
		f[j] = 0.0;
		for (i=0,e=0.0; i<j; i++)
			e += a[i]*r[j-i];
		c = e/v;
		v -= c*e;
		for (i=0; i<=j/2; i++) {
			bot = a[j-i]-c*a[i];
			a[i] -= c*a[j-i];
			a[j-i] = bot;
		}

		/* use a and v above to get f[i], i = 0,1,2,...,j */
		for (i=0,w=0.0; i<j; i++)
			w += f[i]*r[j-i];
		c = (w-g[j])/v;
		for (i=0; i<=j; i++)
			f[i] -= c*a[j-i];
	}
}
double dsinc (double x)
{
	double pix;

	if (x==0.0) {
		return 1.0;
	} else {
		pix = PI*x;
		return sin(pix)/pix;
	}
}
void mksinc (float d, int lsinc, float sinc[])
{
	int j;
	double s[20],a[20],c[20],work[20],fmax;

	/* compute auto-correlation and cross-correlation arrays */
	fmax = 0.066+0.265*log((double)lsinc);
	fmax = (fmax<1.0)?fmax:1.0;
	for (j=0; j<lsinc; j++) {
		a[j] = dsinc(fmax*j);
		c[j] = dsinc(fmax*(lsinc/2-j-1+d));
	}

	/* solve symmetric Toeplitz system for the sinc approximation */
	stoepd(lsinc,a,c,s,work);
	for (j=0; j<lsinc; j++)
		sinc[j] = s[j];
}
void resample (int nxin, float dxin, float yin[], 
	         int nxout, float dxout, float yout[])
{
	static float table[NTABLE][LTABLE];
	static int tabled=0;
	int jtable;
	float frac;
	float fxin;
    float yinl;
	float yinr;

    fxin=0.0;
	yinl=0.0;
	yinr=0.0;

	/* tabulate sinc interpolation coefficients if not already tabulated */
	if (!tabled) {
		for (jtable=1; jtable<NTABLE-1; jtable++) {
			frac = (float)jtable/(float)(NTABLE-1);
			mksinc(frac,LTABLE,&table[jtable][0]);
		}
		for (jtable=0; jtable<LTABLE; jtable++) {
			table[0][jtable] = 0.0;
			table[NTABLE-1][jtable] = 0.0;
		}
		table[0][LTABLE/2-1] = 1.0;
		table[NTABLE-1][LTABLE/2] = 1.0;
		tabled = 1;
	}

	/* interpolate using tabulated coefficients */
	intt8r(NTABLE,table,nxin,dxin,fxin,yin,yinl,yinr,nxout,dxout,yout);
}
