#ifndef IMAGELIB_H
#define IMAGELIB_H

// #include "fftw.h"
// #include "floatdef.h"
#include "stemtypes_fftw3.h"

#include "boost\shared_ptr.hpp"
#include "boost\shared_array.hpp"

typedef struct imageStructType {
  int headerSize;  // first byte of image will be size of image header (in bytes)
                   // This is the size without the data and comment pointers!!!
  int paramSize;   // number of additional parameters
  int commentSize; // length of comment string
  int nx,ny;
  int complexFlag;
  int dataSize;    // size of one data element in bytes (e.g. complex double: 16)
  int version;     // The version flag will later help to find out how to 
                   // distinguish between images produced by different versions of stem
  double t;        // thickness
  double dx,dy;    // size of one pixel
  double *params;  // array for additional parameters
  boost::shared_array<char>comment;   // comment of prev. specified length
} imageStruct;


void getImageHeader(boost::shared_ptr<imageStruct>header,FILE * fp);
boost::shared_ptr<imageStruct>makeNewHeader(int nx,int ny);
boost::shared_ptr<imageStruct>makeNewHeaderCompact(int cFlag,int nx,int ny,double t,double dx,double dy,
				  int paramSize, double *params,char *comment); 
void setHeaderComment(boost::shared_ptr<imageStruct>header, char *comment);
  
boost::shared_ptr<imageStruct>readImage(void ***pix,int nx,int ny,char *fileName);
void writeImage(void **pix, boost::shared_ptr<imageStruct>header, char *fileName);
void writeRealImage(void **pix, boost::shared_ptr<imageStruct>header, char *fileName, int dataSize);
/*
void writeRealImage(fftw_real **pix, int nx, int ny, float_t dx, 
		   float_t dy, float_t t,char *fileName);
*/

void readRealImage(fftw_real **pix, int nx, int ny,real *dx, 
		   real *dy, real *t, char *fileName);

// old image I/O functions:
// void readRealImage_old(fftw_real **pix, int nx, int ny,float_t *t, char *fileName);
// void readImage_old(fftw_complex **pix, int nx, int ny,float_t *t, char *fileName);
// void writeRealImage_old(fftw_real **pix, int nx, int ny, float_t t,char *fileName);
// void writeImage_old(fftw_complex **pix, int nx, int ny, float_t t,char *fileName);


/**************************************************************
 * Here is how to use the new image writing routines
 *
 * static boost::shared_ptr<imageStruct>header = boost::shared_ptr<imageStruct>();
 *
 * if (header == NULL) header = makeNewHeaderCompact(cFlag,Nx,Ny,t,dx,dy,0,NULL,comment);
 * writeImage(cimage,header,filename);
 * or : writeRealImage(rimage,header,filename,sizeof(float));
 *
 **************************************************************
 * Reading an image works like this:
 * 
 * boost::shared_ptr<imageStruct>header;
 * header = readImage((void ***)(&pix),nx,ny,fileName);
 *
 * [This function will read an image.  It reuses the same header 
 * struct over and over.  Therefore, values must be copied from 
 * the header members before calling this function again.
 *
 * The image pointer may also be NULL, in which case memory will be
 * allocated for it, and its size will be returned in the header struct
 * members nx, and ny.]
 **************************************************************/

#endif
