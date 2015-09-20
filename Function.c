#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "Function.h"


float* fmatrix_allocate_1d(int hsize)
 {
  float* matrix;

  matrix=(float*)malloc(sizeof(float)*hsize); 
  if (matrix==NULL) printf("problem of allocating memory");

  return matrix; 
 }


float** fmatrix_allocate_2d(int vsize,int hsize)
 {
  int i;
  float** matrix;
  float *imptr;

  matrix=(float**)malloc(sizeof(float*)*vsize);
  if (matrix==NULL) printf("problem of allocating memory");

  imptr=(float*)malloc(sizeof(float)*hsize*vsize);
  if (imptr==NULL) printf("problem of allocating memory");
 
  for(i=0;i<vsize;i++,imptr+=hsize) matrix[i]=imptr;
  return matrix;
 }


void free_fmatrix_1d(float* pmat)
{ 
  free(pmat); 
}


void free_fmatrix_2d(float** pmat)
{ 
  free(pmat[0]);
  free(pmat);
}


void LoadImagePgm(char* name,float** mat,int length,int width)
{
  int i,j,k;
  unsigned char var;
  char buff[NBCHAR];

  char stringTmp1[NBCHAR],stringTmp2[NBCHAR],stringTmp3[NBCHAR];
 
  int ta1,ta2,ta3;
  FILE *fic;

  /*-----load image name-----*/
  strcpy(buff,name);
  strcat(buff,".pgm");
  printf("---> Open %s",buff);

  /*----open file----*/
  fic=fopen(buff,"r");
  if (fic==NULL)
    { printf("\n- There is an error when open image %s  -\n",buff);
      exit(-1); }

  /*--get info--*/
  fgets(stringTmp1,100,fic);
  fgets(stringTmp2,100,fic);
  fscanf(fic,"%d %d",&ta1,&ta2);
  fscanf(fic,"%d\n",&ta3);

  /*--display info--*/
  printf("\n\n--Info--");
  printf("\n----------");
  printf("\n%s%s%d %d \n%d\n",stringTmp1,stringTmp2,ta1,ta2,ta3);
   
  /*--load image data--*/
  for(i=0;i<length;i++)
      for(j=0;j<width;j++)  
        { fread(&var,1,1,fic);
          mat[i][j]=var; }

 
  fclose(fic);
 }


void SaveImagePgm(char* name,float** mat,int length,int width)
{
  int i,j,k;
  char buff[NBCHAR];
  FILE* fic;
  time_t tm;

  /*--extension--*/
  strcpy(buff,name);
  strcat(buff,".pgm");

  /*--open file--*/
  fic=fopen(buff,"w");
    if (fic==NULL) 
        { printf(" Problem of saving image %s",buff); 
          exit(-1); }
  printf("\n Save %s in format pgm\n",name);

  /*--save info--*/
  fprintf(fic,"P5");
  if (ctime(&tm)==NULL) fprintf(fic,"\n#\n");
  else fprintf(fic,"\n# IMG Module, %s",ctime(&tm));
  fprintf(fic,"%d %d",width,length);
  fprintf(fic,"\n255\n");

  /*--save data--*/
  for(i=0;i<length;i++)
      for(j=0;j<width;j++) 
        fprintf(fic,"%c",(char)mat[i][j]);
   
 
   fclose(fic); 
 } 

/*-----------------*/
/* FOURIER 2D -----*/
/*-----------------*/
void fourn(float data[], unsigned long nn[], int ndim, int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot=1,idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
						tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}



void FFTDD(float** mtxR,float** mtxI,int lgth, int wdth)
{
 int i,j;
 int posx,posy;

 float* data;
 float* ImgFreqR;
 float* ImgFreqI;
 unsigned long* nn;

 /*allocating memory*/
 data=(float*)malloc(sizeof(float)*(2*wdth*lgth)+1);
 ImgFreqR=(float*)malloc(sizeof(float)*(wdth*lgth));
 ImgFreqI=(float*)malloc(sizeof(float)*(wdth*lgth));
 nn=(unsigned long*)malloc(sizeof(unsigned long)*(FFT2D+1)); 

 /*fill nn*/
 nn[1]=lgth; nn[2]=wdth;

 /*fill data*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
   { data[2*(i*lgth+j)+1]=mtxR[i][j];
     data[2*(i*lgth+j)+2]=mtxI[i][j]; }

 /*FFTDD*/
 fourn(data,nn,FFT2D,FFT);

 /*fill data*/
 for(i=0;i<(wdth*lgth);i++)
  { ImgFreqR[i]=data[(2*i)+1];
    ImgFreqI[i]=data[(2*i)+2];  }

 /*matrix conversion*/
 for(i=0;i<(wdth*lgth);i++)
  { posy=(int)(i/wdth);
    posx=(int)(i%wdth);

    mtxR[posy][posx]=ImgFreqR[i]/(wdth*lgth);  
    mtxI[posy][posx]=ImgFreqI[i]/(wdth*lgth); }

 /*free memory*/
 free(data);
 free(ImgFreqR);
 free(ImgFreqI);
 free(nn);
}


void IFFTDD(float** mtxR,float**  mtxI,int lgth,int wdth)
{
 int i,j;
 int posx,posy;

 float* data;
 float* ImgFreqR;
 float* ImgFreqI;
 unsigned long* nn;

 /*allocating memory*/
 data=(float*)malloc(sizeof(float)*(2*wdth*lgth)+1);
 ImgFreqR=(float*)malloc(sizeof(float)*(wdth*lgth));
 ImgFreqI=(float*)malloc(sizeof(float)*(wdth*lgth));
 nn=(unsigned long*)malloc(sizeof(unsigned long)*(FFT2D+1));

 /*fill nn*/
 nn[1]=lgth; nn[2]=wdth;

 /*fill data*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
   { data[2*(i*lgth+j)+1]=mtxR[i][j];
     data[2*(i*lgth+j)+2]=mtxI[i][j]; }

 /*FFTDD*/
 fourn(data,nn,FFT2D,IFFT);

 /*fill data*/
 for(i=0;i<(wdth*lgth);i++)
  { ImgFreqR[i]=data[(2*i)+1];
    ImgFreqI[i]=data[(2*i)+2]; }

 /*matrix conversion*/
 for(i=0;i<(wdth*lgth);i++)
  { posy=(int)(i/wdth);
    posx=(int)(i%wdth);

   mtxR[posy][posx]=ImgFreqR[i];  
   mtxI[posy][posx]=ImgFreqI[i]; }

 /*free memory*/
 free(data);
 free(ImgFreqR);
 free(ImgFreqI);
 free(nn);
}
/*----------end fourier 2d----------*/


/*-----------------*/
/* FOURIER 1D -----*/
/*-----------------*/
void four1(float data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}


void FFT1D(float* vctR,float* vctI,int wdth)
{
 int i,j;

 float* data;
 unsigned long nn;

 /*allocating memory*/
 data=(float*)malloc(sizeof(float)*(2*wdth)+1);

 /*fill nn*/
 nn=wdth;

 /*fill data*/
 for(i=0;i<wdth;i++) 
   { data[(2*i)+1]=vctR[i];
     data[(2*i)+2]=vctI[i]; }

 /*FFT1D*/
 four1(data,nn,FFT);

 /*get data*/
 for(i=0;i<wdth;i++)
  { vctR[i]=data[(2*i)+1];
    vctI[i]=data[(2*i)+2];  }

 /*caliberation*/
 for(i=0;i<wdth;i++)
  { vctR[i]/=(wdth);
    vctI[i]/=(wdth);  }

 /*free memory*/
 free(data);
}


void IFFT1D(float* vctR,float* vctI,int wdth)
{
 int i,j;

 float* data;
 unsigned long nn;

 /*allocating memory*/
 data=(float*)malloc(sizeof(float)*(2*wdth)+1);

 /*fill nn*/
 nn=wdth;

 /*fill data*/
 for(i=0;i<wdth;i++) 
   { data[(2*i)+1]=vctR[i];
     data[(2*i)+2]=vctI[i]; }

 /*FFT1D*/
 four1(data,nn,IFFT);

 /*get data*/
 for(i=0;i<wdth;i++)
  { vctR[i]=data[(2*i)+1];
    vctI[i]=data[(2*i)+2];  }

 /*free memory*/
 free(data);
}
/*-------end fourier 1d-------------*/


void ReMkeVct(float* Vct,int wdth)
{
 int j;
 int cj;
 float*  Vct_tmp;

 /*initialization*/
 cj=(int)(wdth/2);

 /*allocating memory*/
 Vct_tmp=fmatrix_allocate_1d(wdth);

 /*adjust*/
 for(j=0;j<cj;j++) Vct_tmp[cj+j]=Vct[j];

 for(j=cj;j<wdth;j++) Vct_tmp[j-cj]=Vct[j];

 /*transfer*/
 for(j=0;j<wdth;j++) Vct[j]=Vct_tmp[j];

 /*free memory*/
 free(Vct_tmp);
}


void ModVct(float* VctM,float* VctR,float* VctI,int wdth)
{
 int j;

 /*calcul module*/
 for(j=0;j<wdth;j++)
 VctM[j]=sqrt((VctR[j]*VctR[j])+(VctI[j]*VctI[j]));
}



void Mod(float** matM,float** matR,float** matI,int lgth,int wdth)
{
 int i,j;

 /*calcul module*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
 matM[i][j]=sqrt((matR[i][j]*matR[i][j])+(matI[i][j]*matI[i][j]));
}


void ReMkeImg(float** mat,int lgth,int wdth)
{
 int i,j;
 int ci,cj;
 float** mattmp;

 /*initialization*/
 ci=(int)(lgth/2);
 cj=(int)(wdth/2);

 /*allocating memory*/
 mattmp=fmatrix_allocate_2d(lgth,wdth);

 /*adjust*/
 for(i=0;i<ci;i++) for(j=0;j<cj;j++)
 mattmp[ci+i][cj+j]=mat[i][j];

 for(i=ci;i<lgth;i++) for(j=cj;j<wdth;j++)
 mattmp[i-ci][j-cj]=mat[i][j];

 for(i=0;i<ci;i++) for(j=cj;j<wdth;j++)
 mattmp[ci+i][j-cj]=mat[i][j];

 for(i=ci;i<lgth;i++) for(j=0;j<cj;j++)
 mattmp[i-ci][cj+j]=mat[i][j];

 /*transfer*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
  mat[i][j]=mattmp[i][j];

 /*free memoire*/
 free_fmatrix_2d(mattmp);
}


void InverseSense(float** mat,int lgth,int wdth)
{
 int i,j;
 int ci,cj;
 float** mattmp;

 /*allocating memory*/
 mattmp=fmatrix_allocate_2d(lgth,wdth);

 /*transfer*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    mattmp[lgth-i-1][wdth-j-1]=mat[i][j];

 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    mat[i][wdth-j-1]=mattmp[i][j];

 /*free memory*/
 free_fmatrix_2d(mattmp);
}


void Mult(float** mat,float coef,int lgth,int wdth)
{
 int i,j;

 /*multiplication*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
    { mat[i][j]*=coef;
      if (mat[i][j]>GREY_LEVEL) mat[i][j]=GREY_LEVEL; }
}


void Recal(float** mat,int lgth,int wdth)
{
 int i,j;
 float max,min;

 /*initialization*/
 max=0.0;
 min=100000000;

 /*look for min*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    if (mat[i][j]<min) min=mat[i][j];

 /*plus min*/
   for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    mat[i][j]-=min;

 /*look for max*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
    if (mat[i][j]>max) max=mat[i][j];

 /*matrix caliberation*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   mat[i][j]*=(GREY_LEVEL/max);      
}


void MultMatrix(float** matRout,float** matIout,float** mat1Rin,float** mat1Iin,
float** mat2Rin,float** mat2Iin,int lgth,int wdth)
{
 int i,j;

 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   { matRout[i][j]=mat1Rin[i][j]*mat2Rin[i][j]-mat1Iin[i][j]*mat2Iin[i][j];
     matIout[i][j]=mat1Rin[i][j]*mat2Iin[i][j]+mat2Rin[i][j]*mat1Iin[i][j]; }
}


void SquareMatrix(float** matRout,float** matIout,float** matRin,float** matIin,int lgth,int wdth)
{
 int i,j;

 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   { matRout[i][j]=CARRE(matRin[i][j])-CARRE(matIin[i][j]);
     matIout[i][j]=2*matRin[i][j]*matIin[i][j]; }
}


/*--------------*/
/* Noise -------*/
/*--------------*/
float randomize(void)
{ return ((float)rand()/RAND_MAX); }

/*****************************************************/
/* gaussien noise                                    */
/*****************************************************/
float gaussian_noise(float var,float mean)
{
 float noise,theta;

 /*generation of noise*/
 noise=sqrt(-2*var*log(1.0-((float)rand()/RAND_MAX)));
 theta=(float)rand()*1.9175345E-4-PI;
 noise=noise*cos(theta);
 noise+=mean;
 if (noise>GREY_LEVEL) noise=GREY_LEVEL;
 if (noise<0) noise=0;
 return noise;
}

/*****************************************************/
/* add gaussian noise to matrix                      */
/*****************************************************/
void add_gaussian_noise(float** mat,int lgth,int wdth,float var)
{
 int i,j;

 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
 if (var!=0.0) mat[i][j]=gaussian_noise(var,mat[i][j]);
}



/*---------------*/
/* HISTO --------*/
/*---------------*/
/*****************************************************/
/* calcul of histogram                               */
/*****************************************************/
void compute_histo(float** mat,int lgth,int wdth,float* hist)
{
  int i,j;

  for(i=0;i<=GREY_LEVEL;i++) hist[i]=0.0;

  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    if ((mat[i][j]>=0)&&(mat[i][j]<=GREY_LEVEL))
       hist[(int)(mat[i][j])]++;
  
  for(i=0;i<=GREY_LEVEL;i++)  hist[i]/=(wdth*lgth);
}

/*******************************************************/
/* save histogram of image                             */
/*******************************************************/   
void SaveHistoPgm(char* name,float* hist)
{
 int i,j;
 float min,max;
 float** Mat;
 float* histtmp;

 /*allocating memory*/
 Mat=fmatrix_allocate_2d(256,256);
 histtmp=fmatrix_allocate_1d(256);
 
 /*initialization*/
 for(i=0;i<256;i++) for(j=0;j<256;j++) Mat[i][j]=BLACK; 

 /*adjust*/
 for(max=0.0,min=1000000.0,i=0;i<GREY_LEVEL;i++)
   { if (hist[i]>max) max=hist[i];
     if (hist[i]<min) min=hist[i]; }

 for(i=0;i<GREY_LEVEL;i++) histtmp[i]=((hist[i]-min)/max)*255.0; 

 /*display*/
 for(i=0;i<GREY_LEVEL;i++) for(j=0;j<histtmp[i];j++) Mat[GREY_LEVEL-j][i]=WHITE;

 /*save*/
 SaveImagePgm(name,Mat,256,256);

 /*free memory*/
 free_fmatrix_2d(Mat); 
 free_fmatrix_1d(histtmp);  
}

/*---------------*/
/* MARKOV -------*/
/*---------------*/
float funcgauss(float pix,float moy,float var)
{
  return (float)(1/(sqrt(2*PI*var)))*exp(-1*(CARRE((float)pix-moy))/(2*var));
}

/*---------------*/
/* GRAPH --------*/
/*---------------*/
float funcgauss2D(int row,int col,int row_moy,int col_moy,float var)
{
 float dist2;
 float amp;

 dist2=CARRE((float)row-(float)row_moy)+CARRE((float)col-(float)col_moy);
 amp=(1.0/(2.0*PI*var))*exp(-1.0*(dist2/(2*var)));

 return amp;
}


void Convol_U2DBlur(float** mat,int lgth, int wdth,int size)
 {
  int i,j,k,l,m;
  float** mattmp;
  float tmp;
  float nb;

  /*allocating memory*/
  mattmp=fmatrix_allocate_2d(lgth,wdth);

  /*uniform 2D Blur*/
  for(i=0;i<lgth;i++)  for(j=0;j<wdth;j++) 
   {
     nb=0.0; tmp=0.0;
     for(l=-(int)(size/2);l<=(int)(size/2);l++)
     for(m=-(int)(size/2);m<=(int)(size/2);m++)
       {
	if (((i+l)>0)&&((i+l)<lgth)&&((j+m)>0)&&((j+m)<wdth))
          {  nb++; tmp+=mat[i+l][j+m]; }
       
        else continue;
       }    
     mattmp[i][j]=(tmp/nb);
    }

  /*copy matrix*/
  for(i=0;i<lgth;i++)  for(j=0;j<wdth;j++) 
   mat[i][j]=mattmp[i][j];

  /*free memory*/
  free_fmatrix_2d(mattmp);
 }


void RecalMoy(float** matin,float** matout,int lgth,int wdth)
{
 int i,j;
 float moyin,moyout;

 /*look for two means*/
 moyin=moyout=0.0;
 for(i=0;i<lgth;i++)  for(j=0;j<wdth;j++)
   { moyin+=matin[i][j];
     moyout+=matout[i][j];}

 /*adjust value*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    matin[i][j]*=(moyout/moyin);

  /*Verification*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   {
    if (matin[i][j]<0.0)   matin[i][j]=0.0;
    if (matin[i][j]>255.0) matin[i][j]=255.0;
   }
}

void rotation(float** MatriceImgG, float** Mat1, float rotat, int LENGTH, int WIDTH)
{

  int i,j,k,alpha,x,y;

  float a,b,c,d;
  float** temp1;
  float** temp2;
	
  temp1=fmatrix_allocate_2d(2*LENGTH,2*WIDTH);
  temp2=fmatrix_allocate_2d(2*LENGTH,2*WIDTH);
		
  for(i=0;i<2*LENGTH;i++){
	for(j=0;j<2*WIDTH;j++){
		temp1[i][j]=0.0;
		temp2[i][j]=0.0;
	}	
  }
	
  for(i=0;i<2*LENGTH;i++){
	for(j=0;j<2*WIDTH;j++){
		if((i>=LENGTH/2&&i<3*LENGTH/2)&&(j>=WIDTH/2&&j<3*WIDTH/2)){
			temp1[i][j]=MatriceImgG[i-LENGTH/2][j-WIDTH/2];	
		}
		else{
			temp1[i][j]=0.0;			
		}
	}
  }
	
  for(i=0;i<2*LENGTH;i++){
	for(j=0;j<2*WIDTH;j++){
		if((i>=LENGTH/2&&i<3*LENGTH/2)&&(j>=WIDTH/2&&j<3*WIDTH/2)){
			x=(int)((i-LENGTH)*cos(rotat)+(j-WIDTH)*sin(rotat)+LENGTH);
			y=(int)((-1)*(i-LENGTH)*sin(rotat)+(j-WIDTH)*cos(rotat)+WIDTH);
			temp2[x][y]=temp1[i][j];
		}
	}	
  }
	
  /*interpolation bilineaire*/
  for(i=0;i<2*LENGTH;i++){
	for(j=0;j<2*WIDTH;j++){
		if((i>=1&&i<2*LENGTH-2)&&(j>=1&&j<2*WIDTH-2)){
			if(temp2[i][j]==0){
				a=(temp2[i-1][j-1]+temp2[i-1][j+1])/2;
				b=(temp2[i+1][j-1]+temp2[i+1][j+1])/2;	
				c=(temp2[i][j-1]+temp2[i][j+1])/2;
				d=(temp2[i-1][j]+temp2[i+1][j])/2;
				temp2[i][j]=(a+b+2*c+2*d)/6;
			}
		}		
	}	
  }

  for(i=0;i<LENGTH;i++){
	for(j=0;j<WIDTH;j++){
		Mat1[i][j]=temp2[i+LENGTH/2][j+WIDTH/2];		
	}	
  }

	
}


void high_filter(float** MatRFFT, int LENGTH, int WIDTH)
{
  int i,j,k;

  float** tempM;
  float** tempR;
  float** tempI;
  float** MatIFFT;
  float** tempOutR;
  float** tempOutI;
  float** temp1;

  float temp2;
  float temp3;	

  tempM=fmatrix_allocate_2d(LENGTH,WIDTH);
  tempR=fmatrix_allocate_2d(LENGTH,WIDTH);
  tempI=fmatrix_allocate_2d(LENGTH,WIDTH);
  MatIFFT=fmatrix_allocate_2d(LENGTH,WIDTH);
  tempOutR=fmatrix_allocate_2d(LENGTH,WIDTH);
  tempOutI=fmatrix_allocate_2d(LENGTH,WIDTH);
  temp1=fmatrix_allocate_2d(LENGTH,WIDTH);
		
  for(i=0;i<LENGTH;i++){
	for(j=0;j<WIDTH;j++){
	        tempM[i][j]=0.0;
		tempR[i][j]=0.0;
		tempI[i][j]=0.0;
		MatIFFT[i][j]=0.0;
		tempOutR[i][j]=0.0;
		tempOutI[i][j]=0.0;
		temp1[i][j]=0.0;
	}
  }	

  FFTDD(MatRFFT,MatIFFT,LENGTH,WIDTH);
  H(tempR,LENGTH,WIDTH);
  FFTDD(tempR,tempI,LENGTH,WIDTH);
  MultMatrix(tempOutR,tempOutI,tempR,tempI,MatRFFT,MatIFFT,LENGTH,WIDTH);
	
  IFFTDD(tempOutR,tempOutI,LENGTH,WIDTH);
	
  for(i=0;i<LENGTH;i++){
	for(j=0;j<WIDTH;j++){
		MatRFFT[i][j]=tempOutR[i][j];
	}	
  }

  for(i=0;i<LENGTH/2;i++){
	for(j=0;j<WIDTH/2;j++){
                temp2=MatRFFT[i][j];
		MatRFFT[i][j]=MatRFFT[i+LENGTH/2][j+WIDTH/2];
		MatRFFT[i+LENGTH/2][j+WIDTH/2]=temp2;
		temp3=MatRFFT[LENGTH/2+i][j];
		MatRFFT[LENGTH/2+i][j]=MatRFFT[i][WIDTH/2+j];
		MatRFFT[i][WIDTH/2+j]=temp3;				
	}	
  }
}

void H(float** tempR,int LENGTH_RADON,int WIDTH_RADON)
{
   int i,j;
	
   for(i=0;i<LENGTH_RADON;i++){
	 for(j=0;j<WIDTH_RADON;j++){
		if(i>=LENGTH_RADON/2-5&&i<LENGTH_RADON/2+5&&j>=WIDTH_RADON/2-5&&j<WIDTH_RADON/2+5){
			tempR[i][j]=-1;
		}
		else{
			tempR[i][j]=0.0;			
		}
	}		
   }
	
   tempR[LENGTH_RADON/2][WIDTH_RADON/2]=150;
}

void ramp_filter(float** MatriceRadonRFFT,float** MatriceRadonIFFT, int LENGTH_RADON,int WIDTH_RADON)
{
   int i,j,k;

   float*  VctR;
   float*  VctI;

   VctR=fmatrix_allocate_1d(WIDTH_RADON);
   VctI=fmatrix_allocate_1d(WIDTH_RADON);

   for(i=0;i<WIDTH_RADON;i++) 
   {    VctR[i]=0.0;
        VctI[i]=0.0;}

   for(i=0;i<LENGTH_RADON;i++){
	for(j=0;j<WIDTH_RADON;j++){
		MatriceRadonRFFT[i][j] = MatriceRadonRFFT[i][j]*abs(j);	
		MatriceRadonIFFT[i][j] = MatriceRadonIFFT[i][j]*abs(j);
	}
   }
	
}

void reconstruction(float** MatriceRadonRFFT,float** MatRFFT,int LENGTH,int WIDTH)
{
   int i,j,x,angle;
   float rotat;
   float** Mat1;
   float** Mat2;

   Mat1=fmatrix_allocate_2d(LENGTH,WIDTH);
   Mat2=fmatrix_allocate_2d(LENGTH,WIDTH);

   for(i=0;i<LENGTH;i++) for(j=0;j<WIDTH;j++) 
   {  	
     	Mat1[i][j]=0.0;
     	Mat2[i][j]=0.0; 
   }

   for(angle=0;angle<180;angle++)
   {
	rotat=angle*PI/180;
	for(i=0;i<128;i++)
        {
		for(j=0;j<128;j++)
                {
			Mat1[i][j]=MatriceRadonRFFT[angle][i];
		}
	}
	rotation(Mat1,Mat2,PI-rotat,LENGTH,WIDTH);
	for(i=0;i<128;i++)
        {
		for(j=0;j<128;j++)
                {
			MatRFFT[i][j]=MatRFFT[i][j]+Mat2[i][j];
		}
	}
   }
	
   rotation(MatRFFT,MatRFFT,PI/2,LENGTH,WIDTH);
	
}

