

void WriteSGY(float *Data,int NX,int NT,int tao3,float *SX,float *SY,float RX,float RY,float *DSR,char *FILENAME)
{
	int i,nsegy,format,ns,Ntao,itrace[SF_NKEYS];
	float dt;
	//bool verbose=1;
	char ahead[SF_EBCBYTES], bhead[SF_BNYBYTES],trhead[SF_HDRBYTES],*trace;

	FILE *RDSGYORG,*WTSGY;

	Ntao=tao3;

	RDSGYORG=fopen("SGY_Model.sgy","rb");
	fread(ahead,sizeof(char),SF_EBCBYTES,RDSGYORG);   //读卷头，3200字节
	if (SF_BNYBYTES!=fread(bhead,sizeof(char),SF_BNYBYTES,RDSGYORG))
		printf("Error reading binary header");
	format = segyformat (bhead);					//判断数据字符格式
	set_segyns(bhead,  NT);
	set_segydt(bhead, Ntao);

	ns=segyns(bhead);		//每道的样点数
	dt=segydt(bhead);       //采样间隔

	printf("ns=%d\n",ns);
	printf("dt=%0.6f\n",dt);

	nsegy =  ((3 == format)? ns*2: ns*4);  

	printf("nsegy=%u\n",nsegy);

	trace=(char*)malloc(nsegy*sizeof(char));
//	ftrace=(float*)malloc(NT*sizeof(float));
	WTSGY=fopen(FILENAME,"wb");
	fwrite(ahead,sizeof(char),SF_EBCBYTES,WTSGY);
	fwrite(bhead,sizeof(char),SF_BNYBYTES,WTSGY);				//写卷头
	fread(trhead, sizeof(char), SF_HDRBYTES, RDSGYORG);			//读道头	
	for (i=0;i<NX;i++)
	{
		segy2head(trhead,itrace,SF_NKEYS);
		itrace[11]=DSR[i];
		itrace[21]=SX[i];
		itrace[22]=SY[i];
		itrace[23]=RX;
		itrace[24]=RY;

		head2segy(trhead,itrace,SF_NKEYS);

		fwrite(trhead, sizeof(char), SF_HDRBYTES,WTSGY);	//写道头
		trace2segy(trace, &Data[i*NT], NT,format);				//数据转换
		fwrite(trace, sizeof(char), nsegy, WTSGY);			//写数据
	}
	fclose(RDSGYORG);
	fclose(WTSGY);
}

