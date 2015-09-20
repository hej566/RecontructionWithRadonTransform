#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "Function.h"

#define NAME_IMG_IN  "lenna"
#define NAME_IMG_OUT0 "ImgOut0"
#define NAME_IMG_OUT1 "ImgOut1" 
#define NAME_IMG_OUT2 "ImgOut2" 
#define NAME_IMG_OUT3 "ImgOut3"  
#define NAME_IMG_OUT4 "ImgOut4"
#define NAME_IMG_OUT5 "ImgOut5" 

#define NB_PROJECTIONS 180
#define LENGTH 128
#define WIDTH  128

#define LENGTH_RADON NB_PROJECTIONS
#define WIDTH_RADON  WIDTH


int main(int argc,char** argv)
{
  int i,j,k;
  float rotat,angle;
  float sum,a,b,c,d;
  
  float** MatriceImgG;
  float** MatriceRadon;
  float** MatriceRadonRFFT;
  float** MatriceRadonIFFT;
  float** MatriceRadonMFFT;
  float** MatRFFT;
  float** MatIFFT;
  float** MatMFFT;
  float** Mat1;
  float** Mat2;
  float*  VctR;
  float*  VctI;
  
  /*allocating memory*/
  MatriceImgG=fmatrix_allocate_2d(LENGTH,WIDTH);
  MatriceRadon=fmatrix_allocate_2d(LENGTH_RADON,WIDTH_RADON);
  MatriceRadonRFFT=fmatrix_allocate_2d(LENGTH_RADON,WIDTH_RADON);
  MatriceRadonIFFT=fmatrix_allocate_2d(LENGTH_RADON,WIDTH_RADON);
  MatriceRadonMFFT=fmatrix_allocate_2d(LENGTH_RADON,WIDTH_RADON);
  MatRFFT=fmatrix_allocate_2d(LENGTH,WIDTH);
  MatIFFT=fmatrix_allocate_2d(LENGTH,WIDTH);
  MatMFFT=fmatrix_allocate_2d(LENGTH,WIDTH);
  Mat1=fmatrix_allocate_2d(LENGTH,WIDTH);
  Mat2=fmatrix_allocate_2d(LENGTH,WIDTH);

  /*allocating memory*/
  VctR=fmatrix_allocate_1d(WIDTH);
  VctI=fmatrix_allocate_1d(WIDTH);
  
  /*initialization*/
  for(i=0;i<LENGTH;i++) for(j=0;j<WIDTH;j++) 
  {  MatriceImgG[i][j]=0.0;
     MatRFFT[i][j]=0.0;
     MatIFFT[i][j]=0.0;
     MatMFFT[i][j]=0.0;
     Mat1[i][j]=0.0;
     Mat2[i][j]=0.0; }

  for(i=0;i<LENGTH_RADON;i++) for(j=0;j<WIDTH_RADON;j++)
  {   MatriceRadon[i][j]=0.0;
      MatriceRadonRFFT[i][j]=0.0;
      MatriceRadonIFFT[i][j]=0.0;
      MatriceRadonMFFT[i][j]=0.0; }

  /*initialization*/
  for(i=0;i<WIDTH;i++) 
  {    VctR[i]=0.0;
       VctI[i]=0.0;}
    
  LoadImagePgm(NAME_IMG_IN,MatriceImgG,LENGTH,WIDTH);

  for(k=0;k<180;k++){
      rotat=k*PI/180;
      rotation(MatriceImgG,Mat1,rotat,LENGTH,WIDTH);
      if(rotat>24*PI/180&&rotat<=25*PI/180){
		Recal(Mat1,LENGTH,WIDTH);	
   		SaveImagePgm("ROTAT_25",Mat1,LENGTH,WIDTH);
   		system("display ROTAT_25.pgm&");
      }
      for(j=0;j<WIDTH;j++){
	  sum=0;
	  for(i=0;i<LENGTH;i++){
		sum=sum+Mat1[i][j];
	  }
	  MatriceRadon[k][j]=sum;
      }
   }
	
    
   for(i=0;i<LENGTH_RADON;i++){
	for(j=0;j<WIDTH_RADON;j++){
		VctR[j]=MatriceRadon[i][j];	
	}
	FFT1D(VctR,VctI,WIDTH);
	for(k=0;k<WIDTH_RADON;k++){
		MatriceRadonRFFT[i][k]=VctR[k];
		MatriceRadonIFFT[i][k]=VctI[k];
	}
   }
  
  
   ramp_filter(MatriceRadonRFFT,MatriceRadonIFFT, LENGTH_RADON,WIDTH_RADON);
  
   for(i=0;i<LENGTH_RADON;i++){
	for(j=0;j<WIDTH_RADON;j++){
		VctR[j]=MatriceRadonRFFT[i][j];	
		VctI[j]=MatriceRadonIFFT[i][j];
	}
	IFFT1D(VctR,VctI,WIDTH);
	for(k=0;k<WIDTH_RADON;k++){
		MatriceRadonRFFT[i][k]=VctR[k];
		MatriceRadonIFFT[i][k]=VctI[k];
	}
   }
  
   /*-------- FIN ---------------------------------------------*/
   /*----------------------------------------------------------*/
   
   reconstruction(MatriceRadonRFFT,MatRFFT,LENGTH,WIDTH);
   high_filter(MatRFFT,LENGTH,WIDTH);

   Recal(MatriceRadon,LENGTH_RADON,WIDTH_RADON);	
   SaveImagePgm("Transform_of_Radon",MatriceRadon,LENGTH_RADON,WIDTH_RADON);
   system("display Transform_of_Radon.pgm&");

   Mod(MatMFFT,MatRFFT,MatIFFT,LENGTH,WIDTH);
   Recal(MatMFFT,LENGTH,WIDTH);	
   SaveImagePgm("reconstruction",MatMFFT,LENGTH,WIDTH);
   system("display reconstruction.pgm&");  

   /*free matrix*/
   free_fmatrix_2d(MatriceImgG); 
   free_fmatrix_2d(MatriceRadon);
   free_fmatrix_2d(MatriceRadonRFFT);
   free_fmatrix_2d(MatriceRadonIFFT);
   free_fmatrix_2d(MatriceRadonMFFT);
   free_fmatrix_2d(MatRFFT);
   free_fmatrix_2d(MatIFFT);
   free_fmatrix_2d(MatMFFT);
   free_fmatrix_2d(Mat1); 
   free_fmatrix_2d(Mat2); 

   /*free vecteurs*/
   free(VctR);
   free(VctI);   

   printf("\n Ending ... \n\n\n");
   return 0; 	 
}
