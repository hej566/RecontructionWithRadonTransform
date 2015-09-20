#ifndef FONCTIONDEMO_H
#define FONCTIONDEMO_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define CARRE(X) ((X)*(X))

#define NBCHAR 200

#define FFT   1
#define IFFT -1
#define FFT2D 2

#define GREY_LEVEL 255
#define PI 3.141592654

#define WHITE 255
#define BLACK 0


float*  fmatrix_allocate_1d(int);
float** fmatrix_allocate_2d(int,int);
void    free_fmatrix_1d(float*);
void    free_fmatrix_2d(float**);
void    LoadImagePgm(char*,float**,int,int);
void    SaveImagePgm(char*,float**,int,int);

void    fourn(float*,unsigned long*,int,int);
void    FFTDD(float**,float**,int,int);
void    IFFTDD(float**,float**,int,int);
void    Mod(float**,float**,float**,int,int);
void    ReMkeImg(float**,int,int);

void    four1(float*,unsigned long,int);
void    FFT1D(float*,float*,int);
void    IFFT1D(float*,float*,int);
void    ModVct(float*,float*,float*,int);
void    ReMkeVct(float*,int);

void    Recal(float**,int,int);
void    InverseSense(float**,int,int);
void    Mult(float**,float,int,int);

float   randomize(void);
float   gaussian_noise(float,float);
void    add_gaussian_noise(float**,int,int,float);

float   ISNR(float**,float**,float**,int,int);
float   BSNR(float**,int,int,float);

void    MultMatrix(float**,float**,float**,float**,float**,float**,int,int);
void    SquareMatrix(float**,float**,float**,float**,int,int);

void    compute_histo(float**,int,int,float*);
void    SaveHistoPgm(char*,float*);

float   funcgauss(float,float,float);
float   funcgauss2D(int,int,int,int,float);

void    ConvolFreq(float**,float**,float**,int,int);
void    Convol_U2DBlur(float**,int,int,int);
void    RecalMoy(float**,float**,int,int);

void    rotation(float**, float**, float, int, int);
void    high_filter(float**, int, int);
void    H(float**,int,int);
void    ramp_filter(float**, float**, int, int);
void    reconstruction(float**,float**,int, int);

#endif
