/*
  QSTEM - image simulation for TEM/STEM/CBED
  Copyright (C) 2000-2010  Christoph Koch
  Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CONFIG_READER_INTERFACE_H
#define CONFIG_READER_INTERFACE_H

#include "stemtypes_fftw3.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include <vector>
#include <glog/logging.h>
#include <boost/property_tree/ptree.hpp>
using namespace std;
using boost::property_tree::ptree;
namespace QSTEM
{

class IConfigReader;
typedef boost::shared_ptr<IConfigReader> ConfigReaderPtr;
typedef ConfigReaderPtr (*CreateReaderFn)(boost::filesystem::path &filename);

enum class ExperimentType {STEM = 3, CBED = 1, TEM = 4, NBED = 2};
enum class SliceThicknessCalculation {Auto = 1, Thickness = 2, NumberOfSlices = 3};
enum class StructureFactorType {WeickKohl = 1, Rez=2};
enum class SaveLevel{Everything = 1, Something = 2, Results=3};

class IPropertyTreeReader{
public:
	virtual void Read(ptree& t) = 0;
protected:
	virtual ~IPropertyTreeReader(){};
};

class StructureConfig : IPropertyTreeReader{
public:
	string structureFilename;
	int nCellX;
	int nCellY;
	int nCellZ;
	float_tt temperatureK;

	virtual void Read(ptree& t);
};
class BeamConfig : IPropertyTreeReader{
public:
	float_tt EnergykeV, SourceDiameterAngstrom, BeamCurrentpA, DwellTimeMsec;

	virtual void Read(ptree& t);
};
class ModelConfig : IPropertyTreeReader{
public:
	bool UseTDS, TiltBack,CenterSlices;
	int TDSRuns,nPixels,nSlices;
	QSTEM::SliceThicknessCalculation SliceThicknessCalculation;
	float_tt sliceThicknessAngstrom,xOffset,yOffset,zOffset,resolutionXAngstrom,resolutionYAngstrom,crystalTiltX,
				crystalTiltY,crystalTiltZ,beamTiltX,beamTiltY;

	virtual void Read(ptree& t);
};
class PotentialConfig : IPropertyTreeReader{
public:
	bool Use3D, UseFFT, BandlimitTransmissionFunction,SavePotential,SaveProjectedPotential,OneTimeIntegration, PlotVrr;
	QSTEM::StructureFactorType StructureFactorType;
	float_tt AtomRadiusAngstrom;
	int NSubSlabs;

	virtual void Read(ptree& t);
};
class WaveConfig : IPropertyTreeReader{
public:
	float_tt Cs,C5,Cc,dV_V,alpha,Defocus,Astigmatism,AstigmatismAngle,
	a_33,a_31,a_44,a_42,a_55,a_53,a_51,a_66,a_64,a_62,phi_33,phi_31,phi_44,phi_42,phi_55,phi_53,phi_51,phi_66,phi_64,phi_62,gaussScale,
	dI_I,dE_E,AISaperture;
	bool Smooth,Gaussian;

	virtual void Read(ptree& t);
};
class OutputConfig : IPropertyTreeReader{
public:
	int LogLevel,SaveSliceAfterIterations,PropagationProgressInterval,PotentialProgressInterval;
	string SaveFolder;
	QSTEM::SaveLevel SaveLevel;
	bool ShowProbe,PendelloesungPlot;

	virtual void Read(ptree& t);
};

class Config{
public:
	Config(ptree& t);
	int nThreads;
	QSTEM::ExperimentType ExperimentType;
	StructureConfig Structure;
	ModelConfig Model;
	PotentialConfig Potential;
	OutputConfig Output;
	WaveConfig Wave;
	BeamConfig Beam;
//	mode 					cbed
//
//	structure
//	{
//		structure_filename		/home/philipp/QSTEM/bin/Al.cfg
//		ncellx 					10
//		ncelly 					10
//		ncellz 					1
//		temperatureK			300.000000
//	}
//
//	beam
//	{
//		energy_keV				200.000000
//		sourceDiameterAngstrom	0
//		beamCurrentpA         	1
//		dwellTimeMsec			1.6021773e-4
//	}
//
//	model
//	{
//		tds							false
//		{
//			tdsRuns						1
//		}
//		sliceThicknessCalculation	1			; 1=auto , 2=thickness, 3=number of slices
//		sliceThicknessAngstrom		1.012375
//		slices						80
//		centerSlices				false
//		xOffset						0.000000   ;%  x-position offset in cartesian coords
//		yOffset						0.000000   ;%  y-position offset in cartesian coords
//		zOffset						0.976000   ;%  slize z-position offset in cartesian coords
//		resolutionXAngstrom			0.050000
//		resolutionYAngstrom			0.050000
//		crystalTiltX				0.000000	;% tilt in rad
//		crystalTiltY				0.000000
//		crystalTiltZ				0.000000
//		beamTiltX					0.000000 ;deg	% beam tilt deg
//		beamTiltY					0.000000 ;deg
//		tiltBack					false
//		nPixels						1024
//		potential
//		{
//			3D								true
//			FFT								true
//			atomRadiusAngstrom				5.0
//			plotVr_r         				false		;% will create a plot for V_proj(r)*r vs. r  for all slices and 1st atom
//			bandlimitTransmissionFunction	false	;% indicate whether to band limit transmission function or not
//			savePotential					false	;% whether we want to save the projected potential in files
//			saveProjectedPotential        	false		;% whether we want to save the total projected potential
//			oneTimeIntegration 				true  	;% calculate V_proj once and then copy (much faster)
//			structureFactors				1		; 1=WK 2= ...
//		}
//	}
//	wave
//	{
//		Cs							0.000000	; Spherical abberation in mm
//		C5							0.000000	; C_5 abberation in mm
//		Cc							1.000000	; Chromatic abberation in mm
//		dV/V						0.000003	;energy spread in eV (FWHM)
//		alpha						0.00000	 	;Illumination angle in mrad
//		defocus						1.000000
//		astigmatism					0.000000
//		astigmatismAngle: 			0.000000
//		smooth						true		;smoothen edge in rec. space
//		gaussian					false
//	}
//
//	output
//	{
//		loglevel					1;  1=info 			2=warning 	3=error
//		savelevel  					3;  1=everything	2=			3=s
//		saveSliceAfter				1
//		folder						CBED
//		showProbe					false		;% displays a graph of the crosssection of the inc. e-beam
//		propagationProgressInterval	10 			;% show progress every N_prop_prog beam positions
//		potentialProgressInterval	1000 		;% show progress every N_pot_prog atoms
//		pendelloesungPlot			false		;% flag indicates whether to store Pendeloesung plot
//	}
//
//	stem
//	{
//		scan_x_start				10.570000  % X position of top left corner of scan window
//		scan_x_stop					10.570000	% X position of bottom right corner of scan window
//		scan_x_pixels				1	% number of pixels in X-direction
//		scan_y_start				11.300000
//		scan_y_stop					11.300000
//		scan_y_pixels				1
//	}
// from old c code
	//virtual ExperimentType ExperimentType() = 0;
	//virtual PrintLevel PrintLevel() = 0;
	//virtual int SaveLevel() = 0;

//	  real czOffset;
//	  real xOffset;
//	  real yOffset;
//
//	  char cin2[1024];				/* stacking sequence */
//	  char fileBase[512];
//	  char **filein;			/* array of input potential files */
//
//	  char fileWaveIn[512];  // RAM: input .IMG file for entry wavefunction
//	  int lpartl, lstartl;	                /* flags indicating partial
//						   coherence */
//	  char atomPosFile[512];
//	                                        /* and start wavefunction */
//	  float_tt v0;				/* inc. beam energy */
//	  float_tt resolutionX;                  /* real space pixelsize for wave function and potential */
//	  float_tt resolutionY;                  /* real space pixelsize for wave function and potential */
//	  float_tt ctiltx,ctilty,ctiltz;	        /* crystal tilt in mrad */
//	  char cfgFile[512];                        /* file name for writing tilted atomic configuration */
//	  float_tt cubex,cubey,cubez;            /* dimension of crystal cube, if zero, then nx,ny,nz *
//						 * will be used */
//	  int adjustCubeSize;
//	  float_tt btiltx,btilty;   	        /* beam tilt in mrad*/
//	  int tiltBack;               /* tilt back the wave below the specimen */
//	  int *hbeam,*kbeam;		        /* arrays to hold recorded
//						   beam indicies */
//	  int lbeams;				/* flag indicating, whether
//						   to record beams */
//	  char filebeam[512];		 	/* file, that beams get recorded in */
//	  int nbout;				/* number of recorded beams */
//
//	  //int nslic0;				/* slice counter */
//	  int mulsRepeat1;                      /* # of times to repeat structure */
//	  int mulsRepeat2;                      /* for REFINE mode # of mulsRun repeats */
//	  int slices;                           /* number of different slices */
//	  int centerSlices;                     /* flag indicating how to cut the sample */
//	  float_tt **pendelloesung;              /* pendelloesung plot for REFINE mode */
//	  float_tt ax,by,c;	                /* lattice parameters */
//	  float_tt cAlpha,cBeta,cGamma;
//	  double **Mm;                          /* metric matrix Mm(ax,by,cz,alpha,beta,gamma) */
//	  int nCellX,nCellY,nCellZ;             /* number of unit cells in x-y-z dir*/
//	  int natom;				/* number of atoms in "atoms" */
//	  atom *atoms;				/* 3D atoms array */
//	  float_tt atomRadius;                   /* for atom potential boxes */
//	  float_tt potOffsetX,potOffsetY;        /* offset of potential array from zero */
//	  float_tt potSizeX,potSizeY;            /* real space dimensions of potential array in A */
//	  int potNx,potNy;                      /* size of projected potential in pixels */
//	  int nx,ny;				/* size of wave function arrays */
//	  int avgCount;
//	  //float_tt thickness;
//
//	  float_tt C5;
//	  float_tt dE_E;
//	  float_tt dV_V;
//	  float_tt dI_I;
//	  float_tt alpha;
//	  float_tt sourceRadius;
//	  float_tt Cc;
//	  float_tt df0;				/* defocus */
//	  float_tt astigMag;				/* astigmatism*/
//	  float_tt astigAngle;				/* angle of astigmatism */
//
//	  int ismoth;                          /* smoothen the probe wave function */
//	  int gaussFlag;
//	  float_tt gaussScale;
//	  int showProbe;
//	  int displayProgInterval;             /* show progress every .. beam positions */
//	  int displayPotCalcInterval;             /* show progress every .. beam positions */
//
//	  float_tt beamCurrent;  // pico Ampere
//	  float_tt dwellTime;    // msec
//	  float_tt electronScale;  // number of electrons
//
//
//	  int totalSliceCount;
//	  int outputInterval;    // output results every n slices
//
//	  float_tt aobj;				/* obj aperture */
//	  float_tt aAIS;                         /* condensor aperture in A (projected size, */
//	                                        /* for Koehler illumination)                */
//	  // float_tt areaAIS;                      /* fractional area illuminated by AIS (def=1) */
//	  float_tt Cs;			      	/* spher. aberration */
//	  /////////////////////////////////////////////
//	  // more aberrations:
//	  float_tt a33;
//	  float_tt a31;
//	  float_tt a44;
//	  float_tt a42;
//	  float_tt a55;
//	  float_tt a53;
//	  float_tt a51;
//	  float_tt a66;
//	  float_tt a64;
//	  float_tt a62;
//	  float_tt phi33;
//	  float_tt phi31;
//	  float_tt phi44;
//	  float_tt phi42;
//	  float_tt phi55;
//	  float_tt phi53;
//	  float_tt phi51;
//	  float_tt phi66;
//	  float_tt phi64;
//	  float_tt phi62;
//
//	  float_tt acmax,acmin;
//	  float_tt sigmaf;
//	  float_tt dfdelt;
//	  float_tt dfa2,dfa3;
//	  float_tt dfa2phi,dfa3phi;
//	  float_tt chi,phi;
//	  float_tt *sparam;
//
//	  int saveFlag;			/* flag indicating, whether to save the result */
//	  float_tt rmin,rmax;		/* min and max of real part */
//	  float_tt aimin,aimax;		/* min and max of imag part */
//	  float_tt *kx2,*ky2,k2max,*kx,*ky;
//
//	  int nlayer;
//	  float_tt *cz;
//	  float_tt sliceThickness;
//	  int onlyFresnel;
//	  int startSpherical;
//	  float_tt startDistance;
//	  float_tt maxAngle;
//	  float_tt gaussWidth;
//	  float_tt defInfinity;
//	  int accumulateIntensity;
//	  int deconvolute;
//	  int showPhaseplate;
//	  int normHolog;
//	  int gaussianProp;    /* convolute fresnel propagator with gaussian or not */
//	  int nonPeriodZ;      /* for slicecell (make non periodic in Z */
//	  int nonPeriod;       /* for slicecell (make non periodic in x,y */
//	  int bandlimittrans;  /* flag for bandwidth limiting transmission function */
//	  int fftpotential;    /* flag indicating that we should use FFT for V_proj calculation */
//	  int plotPotential;
//	  int storeSeries;
//	  int tds;
//	  int Einstein;        /* if set (default=set), the Einstein model will be used */
//	  char phononFile[512];    /* file name for detailed phonon modes */
//	  int atomKinds;
//	  int *Znums;
//	  double **rPotential;   /* array containing real space potential LUT for each atom kind present */
//	  double *sfkArray;
//	  double **sfTable;
//	  int sfNk;              /* number of k-points in sfTable and sfkArray */
//	  double *u2,*u2avg;     /* (current/averaged) rms displacement of atoms */
//	  float_tt tds_temp;
//	  int savePotential;
//	  int saveTotalPotential;
//	  int readPotential;
//	  float_tt scanXStart,scanXStop,scanYStart,scanYStop;
//	  int scanXN,scanYN;
//	  float_tt intIntensity;
//	  double imageGamma;
//	  char folder[1024];
//	  int avgRuns; // RAM: What is this?
//	  int potential3D;
//	  int scatFactor;
//	  int Scherzer;
//	  std::vector<double> chisq;
//	  int webUpdate;
//	  int cellDiv;
//	  int equalDivs;           // this flag indicates whether we can reuse already pre-calculated potential data
//
//	  /* Parameters for STEM-detectors */
//	  int detectorNum;
//	  /* we will alow as many detector
//				   definitions as the user wants */
//	  std::vector<std::vector<DetectorPtr> > detectors;
//	  //DETECTOR *detectors;
//	  int save_output_flag;
//
//	  double *dE_EArray;
//
//	  // Tomography parameters:
//	  double tomoTilt;  // current tilt in tomography series
//	  double tomoStart; // in rad
//	  double tomoStep;  // in rad
//	  int    tomoCount;  // number of diffraction patterns.
//	  double zoomFactor; // increases the size of the super-box in x,y, in order to
//	                     // make full use of atoms present, creates vacuum edge around sample.
//	std::string &mode;
//	unsigned &printLevel;
//	unsigned &saveLevel;
//	unsigned &displayPotCalcInterval;
//	unsigned &displayProgInterval;
//	boost::filesystem::path &fileOrFolderName;
//	unsigned &nCellX, unsigned &nCellY, unsigned &nCellZ;
//	unsigned &cellDiv;
//	float_tt &btiltx, float_tt &btilty, bool &tiltBack;
//	float_tt &tiltx, float_tt &tilty, float_tt &tiltz,
//	                                      float_tt &cubex, float_tt &cubey, float_tt &cubez,
//	                                      bool &adjustCubeSize;
//	                                      bool &doTDS, float_tt &tdsTemperature,
//	                                                                         boost::filesystem::path &phononFile, bool &useEinstein;
//	                                                                         float_tt &xOffset, float_tt &yOffset;
};

class IConfigReader
{
public:
  // Any ConfigReader class should define this as a private static member.  See read_qsc.hpp for example.
  // virtual ConfigReaderPtr Create()=0;


  virtual void ReadMode(std::string &mode)=0;
  virtual void ReadPrintLevel(unsigned &printLevel)=0;
  virtual void ReadSaveLevel(unsigned &saveLevel)=0;
  virtual void ReadPotentialOutputInterval(unsigned &displayPotCalcInterval)=0;
  virtual void ReadSTEMProgressInterval(unsigned &displayProgInterval)=0;
  virtual void ReadOutputName(boost::filesystem::path &fileOrFolderName)=0;
  virtual void ReadNCells(unsigned &nCellX, unsigned &nCellY, unsigned &nCellZ)=0;
  virtual void ReadNSubSlabs(unsigned &cellDiv)=0;
  virtual void ReadBeamTilt(float_tt &btiltx, float_tt &btilty, bool &tiltBack)=0;
  virtual void ReadCrystalCubeAndTilt(float_tt &tiltx, float_tt &tilty, float_tt &tiltz, 
                                      float_tt &cubex, float_tt &cubey, float_tt &cubez,
                                      bool &adjustCubeSize)=0;
  virtual void ReadTemperatureData(bool &doTDS, float_tt &tdsTemperature, 
                                   boost::filesystem::path &phononFile, bool &useEinstein)=0;
  virtual void ReadSliceOffset(float_tt &xOffset, float_tt &yOffset)=0;
  virtual void ReadProbeArraySize(unsigned &nx, unsigned &ny)=0;
  virtual void ReadResolution(float_tt &resolutionX, float_tt &resolutionY)=0;
  virtual void ReadVoltage(float_tt &voltage)=0;
  virtual void ReadSliceParameters(bool &centerSlices, float_tt &sliceThickness, 
                                   unsigned &nslices, unsigned &outputInterval,
                                   float_tt &zOffset)=0;
  virtual void ReadPeriodicParameters(bool &periodicXY, bool &periodicZ)=0;
  virtual void ReadBandLimitTrans(bool &limit)=0;
  virtual void ReadLoadPotential(bool &loadPotential, boost::filesystem::path &filename)=0;
  virtual void ReadPotentialOutputParameters(bool &savePotential, bool &saveProjectedPotential, 
                                             bool &plotPotential)=0;
  virtual void ReadPotentialCalculationParameters(bool &fftPotential, bool &potential3D)=0;
  virtual void ReadAtomRadius(float_tt &radius)=0;
  virtual void ReadStructureFactorType(std::string &type)=0;
  virtual void ReadPendelloesungParameters(std::vector<int> &hbeams, std::vector<int> &kbeams,
                                           bool &lbeams, unsigned &nbout)=0;
  virtual void ReadAverageParameters(unsigned &avgRuns, bool &storeSeries)=0;
  virtual void ReadScanParameters(float_tt &scanXStart, float_tt &scanXStop, unsigned &scanXN,
                                  float_tt &scanYStart, float_tt &scanYStop, unsigned &scanYN)=0;
  // CBED uses this because stop and npixels are irrelevant
  inline void ReadScanParameters(float_tt &scanXStart, float_tt &scanYStart)
  {
    float_tt dummy_f;
    unsigned dummy_u;
    ReadScanParameters(scanXStart, dummy_f, dummy_u, scanYStart, dummy_f, dummy_u);
  }
  virtual void ReadStructureFileName(boost::filesystem::path &structureFile)=0;

  virtual void ReadNumberOfDetectors(int &numDetectors)=0;
  virtual void ReadDetectorParameters(int det_idx, float_tt &rInside, float_tt &rOutside, std::string &name, 
                              float_tt &shiftX, float_tt &shiftY)=0;
  //void ReadDetectors(std::vector<std::vector<DetectorPtr> > &detectors, std::vector<float_tt> &thicknesses,
  //                 DetectorPtr &detector_to_copy)=0;
  virtual void ReadDoseParameters(float_tt &beamCurrent, float_tt &dwellTimeMs)=0;
  virtual void ReadProbeParameters(float_tt &dE_E, float_tt &dI_I, float_tt &dV_V, float_tt &alpha, float_tt &aAIS)=0;
  virtual void ReadSourceRadius(float_tt &sourceRadius)=0;
  virtual void ReadSmoothingParameters(bool &smooth, float_tt &gaussScale, bool &gaussFlag)=0;

  virtual void ReadTomoParameters(float_tt &tomoTilt, float_tt &tomoStart, float_tt &tomoStep, int &tomoCount,
                     float_tt &zoomFactor)=0;
  virtual void ReadAberrationAmplitudes(float_tt &Cs, float_tt &C5, float_tt &Cc,
                           float_tt &df0, std::string &Scherzer, float_tt &astig,
                           float_tt &a33, float_tt &a31,
                           float_tt &a44, float_tt &a42,
                           float_tt &a55, float_tt &a53, float_tt &a51,
                           float_tt &a66, float_tt &a64, float_tt &a62)=0;
  virtual void ReadAberrationAngles(float_tt &astig,
                       float_tt &phi33, float_tt &phi31,
                       float_tt &phi44, float_tt &phi42,
                       float_tt &phi55, float_tt &phi53, float_tt &phi51,
                       float_tt &phi66, float_tt &phi64, float_tt &phi62)=0;
  inline bool IsValid() {return m_isValid;}
protected:
  bool m_isValid;
};

}

#endif
