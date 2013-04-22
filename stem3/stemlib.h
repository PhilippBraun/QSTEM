#ifndef STEMLIB_H
#define STEMLIB_H

// #define WIN



#include "stemtypes_fftw3.h"


/**********************************************
 * This function creates a incident STEM probe 
 * at position (dx,dy)
 * with parameters given in muls
 *********************************************/
// int probe(MULS *muls,double dx, double dy);
void probe(MULS *muls, WAVEFUNC *wave, double dx, double dy);
void probePlot(MULS *muls,WAVEFUNC *wave);

void initSTEMSlices(MULS *muls, int nlayer);
void interimWave(MULS *muls,WAVEFUNC *wave,int slice);
void interimCollect(MULS *muls, WAVEFUNC *wave, int slices,int finalSlice);
void detectorCollect(MULS *muls, WAVEFUNC *wave);

void make3DSlices(MULS *muls,int nlayer,char *fileName,atom *center);
void make3DSlicesFFT(MULS *muls,int nlayer,char *fileName,atom *center);
void createAtomBox(MULS *muls, int Znum, atomBox *aBox);
void transmit(void **wave,void **trans,int nx, int ny,int posx,int posy);
void propagate_slow(void** wave,int nx, int ny,MULS *muls);
fftwf_complex *getAtomPotential3D_3DFFT(int Znum, MULS *muls,double B);
fftwf_complex *getAtomPotential3D(int Znum, MULS *muls,double B,int *nzSub,int *Nr,int*Nz_lut);
fftwf_complex *getAtomPotentialOffset3D(int Znum, MULS *muls,double B,int *nzSub,int *Nr,int*Nz_lut,float q);
fftwf_complex *getAtomPotential2D(int Znum, MULS *muls,double B);

WAVEFUNC initWave(int nx, int ny);
void readStartWave(MULS *muls, WAVEFUNC *wave);
/******************************************************************
 * runMulsSTEM() - do the multislice propagation in STEM mode
 *
 * waver, wavei are expected to contain incident wave function 
 * the will be updated at return
 *****************************************************************/
int runMulsSTEM_old(MULS *muls,int lstart);
int runMulsSTEM(MULS *muls, WAVEFUNC *wave);
void writePix(char *outFile,fftw_complex **pict,MULS *muls,int iz);
void fft_normalize(void **array,int nx, int ny);
void showPotential(fftw_complex ***pot,int nz,int nx,int ny,
		   double dx,double dy,double dz);
void atomBoxLookUp(fftw_complex *vlu,MULS *muls,int Znum,double x,double y,
			   double z,double B);
void writeBeams(MULS *muls, WAVEFUNC *wave,int ilayer);

/***********************************************************************************
 * old image read/write functions, may soon be outdated
 */
void readRealImage_old(fftw_real **pix, int nx, int ny, real *t, char *fileName);
void readImage_old(fftw_complex **pix, int nx, int ny, real *t, char *fileName);
void writeRealImage_old(fftw_real **pix, int nx, int ny, real t,char *fileName);
void writeImage_old(fftw_complex **pix, int nx, int ny, real t,char *fileName);

#endif /* STEMLIB_H */