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

#ifndef READ_QSC_H
#define READ_QSC_H

#include "read_interface.hpp"
#include <boost/filesystem.hpp>

namespace QSTEM 
{

class CQscReader : public IConfigReader
{
public:
  CQscReader(boost::filesystem::path &filename);
  ~CQscReader();

  void ReadMode(std::string &mode);
  void ReadPrintLevel(unsigned &printLevel);
  void ReadSaveLevel(unsigned &saveLevel);
  void ReadPotentialOutputInterval(unsigned &displayPotCalcInterval);
  void ReadSTEMProgressInterval(unsigned &displayProgInterval);
  void ReadOutputName(boost::filesystem::path &fileOrFolderName);
  void ReadNCells(unsigned &nCellX, unsigned &nCellY, unsigned &nCellZ);
  void ReadNSubSlabs(unsigned &cellDiv);
  void ReadBeamTilt(float_tt &btiltx, float_tt &btilty, bool &tiltBack);
  void ReadCrystalCubeAndTilt(float_tt &tiltx, float_tt &tilty, float_tt &tiltz, 
                                      float_tt &cubex, float_tt &cubey, float_tt &cubez,
                                      bool &adjustCubeSize);
  void ReadTemperatureData(bool &doTDS, float_tt &tdsTemperature, 
                           boost::filesystem::path &phononFile, bool &useEinstein);
  void ReadSliceOffset(float_tt &xOffset, float_tt &yOffset);
  void ReadProbeArraySize(unsigned &nx, unsigned &ny);
  void ReadResolution(float_tt &resolutionX, float_tt &resolutionY);
  void ReadVoltage(float_tt &voltage);
  void ReadSliceParameters(bool &centerSlices, float_tt &sliceThickness, 
                                   unsigned &nslices, unsigned &outputInterval,
                                   float_tt &zOffset);
  void ReadPeriodicParameters(bool &periodicXY, bool &periodicZ);
  void ReadBandLimitTrans(bool &limit);
  void ReadLoadPotential(bool &loadPotential, boost::filesystem::path &filename);
  void ReadPotentialOutputParameters(bool &savePotential, bool &saveProjectedPotential, 
                                             bool &plotPotential);
  void ReadPotentialCalculationParameters(bool &fftPotential, bool &potential3D);
  void ReadAverageParameters(unsigned &avgRuns, bool &storeSeries);
  void ReadScanParameters(float_tt &scanXStart, float_tt &scanXStop, unsigned &scanXN,
                                  float_tt &scanYStart, float_tt &scanYStop, unsigned &scanYN);  
  void ReadAtomRadius(float_tt &radius);
  void ReadStructureFactorType(std::string &type);
  void ReadPendelloesungParameters(std::vector<int> &hbeams, std::vector<int> &kbeams,
                                   bool &lbeams, unsigned &nbout);
  void ReadStructureFileName(boost::filesystem::path &structure_file); // 
  
  void ReadNumberOfDetectors(int &numDetectors);
  void ReadDetectorParameters(int det_idx, float_tt &rInside, float_tt &rOutside, std::string &name, 
                              float_tt &shiftX, float_tt &shiftY);
  //void ReadDetectors(std::vector<std::vector<DetectorPtr> > &detectors, std::vector<float_tt> &thicknesses,
  //                           DetectorPtr &detector_to_copy)
  void ReadDoseParameters(float_tt &beamCurrent, float_tt &dwellTimeMs);
  void ReadProbeParameters(float_tt &dE_E, float_tt &dI_I, float_tt &dV_V, float_tt &alpha, float_tt &aAIS);
  void ReadSourceRadius(float_tt &sourceRadius);
  void ReadSmoothingParameters(bool &ismoth, float_tt &gaussScale, bool &gaussFlag);
  void ReadTomoParameters(float_tt &tomoTilt, float_tt &tomoStart, float_tt &tomoStep, int &tomoCount,
                     float_tt &zoomFactor);
  void ReadAberrationAmplitudes(float_tt &Cs, float_tt &C5, float_tt &Cc,
                           float_tt &df0, std::string &Scherzer, float_tt &astig,
                           float_tt &a33, float_tt &a31,
                           float_tt &a44, float_tt &a42,
                           float_tt &a55, float_tt &a53, float_tt &a51,
                           float_tt &a66, float_tt &a64, float_tt &a62);
  void ReadAberrationAngles(float_tt &astig,
                       float_tt &phi33, float_tt &phi31,
                       float_tt &phi44, float_tt &phi42,
                       float_tt &phi55, float_tt &phi53, float_tt &phi51,
                       float_tt &phi66, float_tt &phi64, float_tt &phi62);
protected:
  std::string m_buf;
  char answer[256];
  FILE *m_fp;

  bool IsBufferYes(std::string &buf);
private:
  friend class CConfigReaderFactory;
  // Create an instance of this class, wrapped in a shared ptr
  //     This should not be inherited - any subclass needs its own implementation.
  static ConfigReaderPtr Create(boost::filesystem::path &filename){return ConfigReaderPtr(new CQscReader(filename));}
};

class QSTEM_HELPER_DLL_EXPORT StructureConfig : IPropertyTreeReader{
public:
	string structureFilename;
	int nCellX;
	int nCellY;
	int nCellZ;
	float_tt temperatureK;

	virtual void Read(ptree& t);
};
class QSTEM_HELPER_DLL_EXPORT BeamConfig : IPropertyTreeReader{
public:
	float_tt EnergykeV, SourceDiameterAngstrom, BeamCurrentpA, DwellTimeMsec;

	virtual void Read(ptree& t);
};
class QSTEM_HELPER_DLL_EXPORT ModelConfig : IPropertyTreeReader{
public:
	bool UseTDS, TiltBack, CenterSlices;
	int TDSRuns, nx, ny, nSlices;
	QSTEM::SliceThicknessCalculation SliceThicknessCalculation;
	float_tt sliceThicknessAngstrom, xOffset, yOffset, zOffset, dx, dy, crystalTiltX,
		crystalTiltY, crystalTiltZ, beamTiltX, beamTiltY;

	virtual void Read(ptree& t);
};
class QSTEM_HELPER_DLL_EXPORT PotentialConfig : IPropertyTreeReader{
public:
	bool Use3D, UseFFT, BandlimitTransmissionFunction,  OneTimeIntegration, PlotVrr,periodicXY,periodicZ;
	QSTEM::StructureFactorType StructureFactorType;
	float_tt AtomRadiusAngstrom;
	int NSubSlabs;

	virtual void Read(ptree& t);
};
class QSTEM_HELPER_DLL_EXPORT WaveConfig : IPropertyTreeReader{
public:
	float_tt Cs, C5, Cc, dV_V, alpha, Defocus, Astigmatism, AstigmatismAngle,
		a_33, a_31, a_44, a_42, a_55, a_53, a_51, a_66, a_64, a_62, phi_33, phi_31, phi_44, phi_42, phi_55, phi_53, phi_51, phi_66, phi_64, phi_62, gaussScale,
		dI_I, dE_E, AISaperture;
	bool Smooth, Gaussian;

	virtual void Read(ptree& t);
};
class QSTEM_HELPER_DLL_EXPORT OutputConfig : IPropertyTreeReader{
public:
	int LogLevel, SaveSliceAfterIterations, PropagationProgressInterval, PotentialProgressInterval,SavePotential,
	SaveProjectedPotential;
	string SaveFolder;
	QSTEM::SaveLevel SaveLevel;
	bool ShowProbe, PendelloesungPlot;

	virtual void Read(ptree& t);
};

class QSTEM_HELPER_DLL_EXPORT Config{
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
typedef boost::shared_ptr<Config> ConfigPtr;
} // end namespace QSTEM

#endif
