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

#ifndef WAVE_BASE_H
#define WAVE_BASE_H

#include "stemtypes_fftw3.hpp"
#include "imagelib_fftw3.hpp"
#include "memory_fftw3.hpp"
#include "config_IO/config_reader_factory.hpp"

#include "wave_interface.hpp"
#include "omp.h"
#include "boost/format.hpp"
#include <boost/log/trivial.hpp>

#include <cvmlcpp/signal/Fourier>

void CreateWaveFunctionDataSets(unsigned x, unsigned y, std::vector<unsigned> positions, std::string output_ext);

static std::string waveFilePrefix="wave";
static std::string dpFilePrefix="diff";
static std::string probeFilePrefix="probe_wave";
static std::string waveIntensityFilePrefix="wave_int";

namespace QSTEM
{

class QSTEM_HELPER_DLL_EXPORT CBaseWave : public IWave
{
public:
	CBaseWave(const ConfigPtr c);
	CBaseWave( const CBaseWave& other );
	virtual ~CBaseWave();

	virtual WavePtr Clone()=0;

	void Resize(unsigned x, unsigned y);
	void CreateDataSets();
	virtual void FormProbe()=0;

	void DisplayParams();

	void ToRealSpace();
	void ToFourierSpace();
	bool IsRealSpace(){return m_realSpace;}

	//void CopyDPToAvgArray(float_tt *avgArray);
	//void AddDPToAvgArray(unsigned avgCount);

	void GetSizePixels(unsigned &x, unsigned &y) const ;
	unsigned GetTotalPixels() const {return m_nx*m_ny;}
	void GetResolution(float_tt &x, float_tt &y) const ;
	void GetPositionOffset(unsigned &x, unsigned &y) const ;
	float_tt GetK2(unsigned ix, unsigned iy) const ;
	inline float_tt GetKX2(unsigned ix) const {return m_kx2[ix];}
	inline float_tt GetKY2(unsigned iy) const {return m_ky2[iy];}
	inline float_tt GetK2Max() const {return m_k2max;}

	inline float_tt GetVoltage()  const {return m_v0;}
	inline float_tt GetWavelength()  const {return m_wavlen;}
	inline ComplexArray2DPtr GetWave() const {return m_wave;}

	inline float_tt GetPixelIntensity(unsigned i) const {return (m_wave.data())[i].real()*(m_wave.data())[i].real() + (m_wave.data())[i].imag()*(m_wave.data())[i].imag();}
	inline float_tt GetPixelIntensity(unsigned x, unsigned y) const  {return GetPixelIntensity(x+m_nx*y);}
	inline float_tt GetDiffPatPixel(unsigned i)  const {return m_diffpat[i];}
	inline float_tt GetDiffPatPixel(unsigned x, unsigned y) const  { return m_diffpat[x+m_nx*y];}
	inline void SetDiffPatPixel(unsigned i, float_tt value) {m_diffpat[i]=value;}
	inline void SetDiffPatPixel(unsigned x, unsigned y, float_tt value) {m_diffpat[x+m_nx*y]=value;}

	void ApplyTransferFunction(boost::shared_array<complex_tt> &wave);

	void WriteBeams(unsigned absoluteSlice);

	inline void Save(std::string filename){
		_WriteWave(filename);
	}

	inline void WriteProbe()
	{
		_WriteWave(probeFilePrefix);
	}

	inline void WriteWave(std::string comment="Wavefunction")
	{
		m_position.clear();
		std::map<std::string, double> params;
		_WriteWave(waveFilePrefix, comment, params);
	}
	inline void WriteWave(unsigned navg, std::string comment="Wavefunction",
			std::map<std::string, double>params = std::map<std::string, double>())
	{
		SetWavePosition(navg);
		_WriteWave(waveFilePrefix, comment, params);
	}
	inline void WriteWave(unsigned posX, unsigned posY, std::string comment="Wavefunction",
			std::map<std::string, double>params = std::map<std::string, double>())
	{
		SetWavePosition(posX, posY);
		_WriteWave(waveFilePrefix, comment, params);
	}
	inline void WriteWave(unsigned posX, unsigned posY, unsigned posZ, std::string comment="Wavefunction",
			std::map<std::string, double>params = std::map<std::string, double>())
	{
		SetWavePosition(posX, posY, posZ);
		_WriteWave(waveFilePrefix, comment, params);
	}
	inline void WriteWave(string filename, std::string comment,
			std::map<std::string, double>params = std::map<std::string, double>())
	{
		_WriteWave(filename, comment, params);
	}
	inline void WriteDiffPat(std::string comment="Diffraction Pattern",
			std::map<std::string, double>params = std::map<std::string, double>())
	{
		m_position.clear();
		_WriteDiffPat(dpFilePrefix, comment, params);
	}
	inline void WriteDiffPat(unsigned navg, std::string comment="Diffraction Pattern",
			std::map<std::string, double>params = std::map<std::string, double>())
	{
		SetWavePosition(navg);
		_WriteDiffPat(dpFilePrefix, comment, params);
	}
	inline void WriteDiffPat(unsigned posX, unsigned posY, std::string comment="Diffraction Pattern",
			std::map<std::string, double>params = std::map<std::string, double>())
	{
		SetWavePosition(posX, posY);
		_WriteDiffPat(dpFilePrefix, comment, params);
	}
	inline void WriteDiffPat(unsigned posX, unsigned posY, unsigned posZ, std::string comment="Diffraction Pattern",
			std::map<std::string, double>params = std::map<std::string, double>())
	{
		SetWavePosition(posX, posY, posZ);
		_WriteDiffPat(dpFilePrefix, comment, params);
	}

	// People can change the wavefunction - for example, that's what we have to do when we
	//    transmit the wave through the sample's potential.
	complex_tt *GetWavePointer(){return m_wave.data();}
	// People should not directly change the diffraction pattern, since we'll re-calculate it when
	//   the wavefunction changes.
	//   They can, however, access it.
	const float_tt *GetDPPointer(){return &m_diffpat[0];}

	float_tt GetIntegratedIntensity() const ;

	void ReadWave();
	void ReadWave(unsigned navg);
	void ReadWave(unsigned posX, unsigned posY);
	void ReadDiffPat();
	void ReadDiffPat(unsigned navg);
	void ReadDiffPat(unsigned posX, unsigned posY);

protected:
	ImageIOPtr m_imageIO;
	ConfigPtr _config;

	bool m_realSpace;  // If true, the m_wave is in real space.  Else, it's in Fourier space.

	std::string m_fileStart;
	std::string m_avgName;
	std::string m_fileout;
	unsigned m_detPosX, m_detPosY;
	unsigned m_nx, m_ny;		      /* size of wavefunc and diffpat arrays */
	RealVector m_diffpat;
	//float_tt **m_avgArray;
	//float_tt m_thickness;
	//float_tt m_intIntensity;
	//float_tt m_electronScale;
	//float_tt m_beamCurrent;
	//float_tt m_dwellTime;
	float_tt m_v0;
	std::vector<unsigned> m_position;
	std::map<std::string, double> m_params;

	// defocus mode: 1 = Scherzer, 2 = ???
	int m_Scherzer;
	int m_printLevel;
	float_tt m_dx, m_dy;  // physical pixel size of wavefunction array
	ComplexArray2D m_wave; /* complex wave function */
	RealVector m_kx2,m_ky2,m_kx,m_ky;
	float_tt m_k2max;

	cvmlcpp::DFT<float_tt, 2> _fft;
	cvmlcpp::DFT<float_tt, 2> _ifft;

#if FLOAT_PRECISION == 1
	fftwf_plan m_fftPlanWaveForw,m_fftPlanWaveInv;
#else
	fftw_plan m_fftPlanWaveForw,m_fftPlanWaveInv;
#endif

protected:
	void Initialize(std::string input_ext, std::string output_ext);
	void InitializeKVectors();

	void _WriteWave(std::string &prefix, std::string comment="Wavefunction",
			std::map<std::string, double>params = std::map<std::string, double>());
	void _WriteDiffPat(std::string &prefix, std::string comment="Diffraction Pattern",
			std::map<std::string, double>params = std::map<std::string, double>());
	//void _WriteAvgArray(std::string &prefix, std::string comment="Average Array",
	//                   std::map<std::string, double>params = std::map<std::string, double>());

	// For CBED ( &TEM? )
	void SetWavePosition(unsigned navg);
	// For STEM
	void SetWavePosition(unsigned posX, unsigned posY);
	// If you need to save things wrt 3 dimensions...
	void SetWavePosition(unsigned posX, unsigned posY, unsigned posZ);

	float_tt Wavelength(float_tt keV);
	float_tt m_wavlen;

	// m_transferFunction  // The transfer function - optionally applied (used by TEM mode)
};

}

#endif
