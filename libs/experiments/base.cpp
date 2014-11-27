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

#include "base.hpp"
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
using boost::format;

namespace QSTEM
{

CExperimentBase::CExperimentBase(ConfigPtr c) : IExperiment()
{
	_config = c;
	m_sample = StructurePtr(new CCrystal(c));

	float_tt max_x, min_x, max_y, min_y, max_z, min_z, zTotal;
	m_sample->GetCrystalBoundaries(min_x, max_x, min_y, max_y, min_z, max_z);
	zTotal = max_z - min_z;

	switch (c->Model.SliceThicknessCalculation) {
	case SliceThicknessCalculation::Auto:
		m_dz = c->Model.sliceThicknessAngstrom = zTotal/((int)zTotal);
		c->Model.nSlices = (int)zTotal;
		break;
	case SliceThicknessCalculation::NumberOfSlices:
		m_dz = c->Model.sliceThicknessAngstrom = zTotal/c->Model.nSlices;
		break;
	case SliceThicknessCalculation::Thickness:
		m_dz = c->Model.sliceThicknessAngstrom;
		c->Model.nSlices = (int)(zTotal / c->Model.sliceThicknessAngstrom);
		break;
	default:
		break;
	}

	int atomRadiusSlices = ceil(_config->Potential.AtomRadiusAngstrom / c->Model.sliceThicknessAngstrom);
	if(_config->Potential.Use3D) c->Model.nSlices += 2*atomRadiusSlices;


	atomRadiusSlices = ceil(_config->Potential.AtomRadiusAngstrom / _config->Model.dx);
	_config->Model.nx = ceil((max_x - min_x) / _config->Model.dx) ;
	if(_config->Potential.periodicXY == false)  _config->Model.nx += atomRadiusSlices;


	atomRadiusSlices = ceil(_config->Potential.AtomRadiusAngstrom / _config->Model.dy);
	_config->Model.ny = ceil((max_y - min_y) / _config->Model.dy) ;
	if(_config->Potential.periodicXY == false)  _config->Model.ny += atomRadiusSlices;

	m_wave = CWaveFactory::Get()->GetWave("Convergent", c);

	m_potential = CPotFactory::Get()->GetPotential(c);
	m_potential->SetStructure(m_sample);

	m_equalDivs = true;
	m_saveLevel = static_cast<unsigned>(c->Output.SaveLevel);

	DisplayParams();
	m_sample->DisplayParams();
	m_potential->DisplayParams();
	m_wave->DisplayParams();

	m_avgArray = RealVector();
}

void CExperimentBase::DisplayParams() {
	FILE *fpDir;
	char systStr[64];
	double k2max,temp;
	int i,j;
	static char Date[16],Time[16];
	time_t caltime;
	struct tm *mytime;
	const double pi=3.1415926535897;

	/*
  if (wave->printLevel < 1) {
    if ((fpDir = fopen(muls.folder.c_str(),"r"))) {
      fclose(fpDir);
      // printf(" (already exists)\n");
    }
    else {
      sprintf(systStr,"mkdir %s",muls.folder.c_str());
      system(systStr);
      // printf(" (created)\n");
    }	  
    return;
  }
	 */
	caltime = time( NULL );
	mytime = localtime( &caltime );
	strftime( Date, 12, "%Y:%m:%d", mytime );
	strftime( Time, 9, "%H:%M:%S", mytime );

	BOOST_LOG_TRIVIAL(info) << "**************************************************************************************************";
	BOOST_LOG_TRIVIAL(info) << format("* Running program STEM3 (version %s) in %d mode") % VERSION % static_cast<int>(_config->ExperimentType);
	BOOST_LOG_TRIVIAL(info) << format("* Date: %s, Time: %s") % Date%Time;
	BOOST_LOG_TRIVIAL(info) << format("* Output file/folder:          ./%s/ ") % _config->Output.SaveFolder.c_str();
	BOOST_LOG_TRIVIAL(info) << format("* Super cell divisions: %d (in z direction) %s") % _config->Potential.NSubSlabs%
			(m_equalDivs ? "equal" : "non-equal");
	BOOST_LOG_TRIVIAL(info) << format("* Output every:         %d slices") %
			_config->Output.PropagationProgressInterval;
	BOOST_LOG_TRIVIAL(info) << format("* TDS:                  %d runs)") % _config->Model.TDSRuns;
	BOOST_LOG_TRIVIAL(info) << "******************************************************************************************";
}

void CExperimentBase::DisplayProgress(int flag)
{
	// static double timer;
	static double timeAvg = 0;
	static double intensityAvg = 0;
	static time_t time0,time1;
	double curTime;
	int jz;

	if (flag < 0) {
		time(&time0);
		// timer = cputim();
		return;
	}
	time(&time1);
	curTime = difftime(time1,time0);
	/*   curTime = cputim()-timer;
       if (curTime < 0) {
       printf("timer: %g, curr. time: %g, diff: %g\n",timer,cputim(),curTime);
       }
	 */
	if (_config->Output.LogLevel > 0) {
		if (m_sample->GetTDS()) {
			timeAvg = ((m_avgCount)*timeAvg+curTime)/(m_avgCount+1);
			intensityAvg = ((m_avgCount)*intensityAvg+m_intIntensity)/(m_avgCount+1);
			BOOST_LOG_TRIVIAL(info) << format("********************** run %3d ************************") % (m_avgCount+1);

			std::map<unsigned, float_tt> displacements(m_sample->GetU2());
			std::map<unsigned, float_tt>::iterator disp=displacements.begin(), end=displacements.end();

			BOOST_LOG_TRIVIAL(info) << format("* <u>: %3d |") % (*disp++).first;
			while(disp!=end) BOOST_LOG_TRIVIAL(info) << format(" %8d |") % (*disp++).first;

			BOOST_LOG_TRIVIAL(info) << " intensity | time(sec) |    chi^2  |";
			// }
			/*
        printf("* %9g | %9g | %9g \n",muls.u2,muls.intIntensity,curTime);  
        }
        else {
			 */
			BOOST_LOG_TRIVIAL(info) << "*";

			//ComputeAverageU2();

			disp = displacements.begin();
			while (disp!=end) BOOST_LOG_TRIVIAL(info) << format(" %8f |") % (*disp++).second;
			BOOST_LOG_TRIVIAL(info) << format(" %9f | %9f | %9f |") % m_intIntensity%curTime%(m_avgCount > 0 ? m_chisq[m_avgCount-1] : 0);
			BOOST_LOG_TRIVIAL(info) << format("*");

			/*
        // TODO: averaging should be handled on this class, not on lower level crystal class.
      atom = atomTypes.begin();
      while (atom!=end) printf(" %8f |",(float)(m_crystal->GetU2avg((*atom++))));  
			 */
			BOOST_LOG_TRIVIAL(info) << format(" %9f | %9f ") %intensityAvg%timeAvg;
		}
		else {
			BOOST_LOG_TRIVIAL(info) << format("**************** finished after %.1f sec ******************") % curTime;
		}
	}  // end of printLevel check.

	time(&time0);
	//  timer = cputim();
}

/*
// TODO: this is very broken.  displaced here from Crystal because I want to handle averages outside of the 
//     lower level classes.
void CExperimentBase::ComputeAverageU2()
{

  (*z)->second /= u2Count[(*z)->first];
  if (runCount > 0) 
    m_u2avg[(*z)] = sqrt(((runCount-1)*(m_u2avg[(*z)]*m_u2avg[(*z)])+u2[(*z)])/runCount);
  else
    m_u2avg[(*z)] = sqrt(m_u2[(*z)]);
}
 */

////////////////////////////////////////////////////////////////
// save the current wave function at this intermediate thickness:
void CExperimentBase::InterimWave(int slice) {
	int t;
	char fileName[256];
	std::map<std::string, double> params;

	if ((slice < _config->Model.nSlices*_config->Potential.NSubSlabs-1) && ((slice+1) % _config->Output.PropagationProgressInterval != 0)) return;

	t = (int)((slice)/_config->Output.PropagationProgressInterval);

	// produce the following filename:
	// wave_avgCount_thicknessIndex.img or
	// wave_thicknessIndex.img if tds is turned off
	if (_config->Model.UseTDS) m_wave->WriteWave(m_avgCount, t, "Wave Function", params);
	else m_wave->WriteWave(t, "Wave Function", params);
}

void CExperimentBase::InitializePropagators(WavePtr wave){
	unsigned nx, ny;

	wave->GetSizePixels(nx, ny);
	m_propxr.resize(nx);
	m_propxi.resize(nx);
	m_propyr.resize(ny);
	m_propyi.resize(ny);

	float_tt scale = m_dz*PI;

#pragma omp parallel for
	for(int ixa=0; ixa<nx; ixa++) {
		float_tt t = scale * (wave->GetKX2(ixa)*wave->GetWavelength());
		m_propxr[ixa] = (float_tt)  cos(t);
		m_propxi[ixa] = (float_tt) -sin(t);
	}
#pragma omp parallel for
	for(int iya=0; iya<ny; iya++) {
		float_tt t = scale * (wave->GetKY2(iya)*wave->GetWavelength());
		m_propyr[iya] = (float_tt)  cos(t);
		m_propyi[iya] = (float_tt) -sin(t);
	}
}

int CExperimentBase::RunMultislice(WavePtr wave) 
{
	int printFlag = 0;
	int showEverySlice=1;
	int islice,i,ix,iy,mRepeat;
	float_tt cztot=0.0;
	float_tt wavlen,sum=0.0; //,zsum=0.0
	// static int *layer=NULL;
	float_tt x,y;
	int absolute_slice;

	char outStr[64];
	double fftScale;

	unsigned nx, ny;

	wave->GetSizePixels(nx, ny);

	printFlag = (_config->Output.LogLevel > 3);
	fftScale = 1.0/(nx*ny);

	wavlen = wave->GetWavelength();

	m_avgArray.resize(wave->GetTotalPixels());
	m_imageIO = ImageIOPtr(new CImageIO(nx,ny));

	/*  calculate the total specimen thickness and echo */
	cztot=0.0;

	if (printFlag){
		for( islice=0; islice<_config->Model.nSlices; islice++) {
			cztot += _config->Model.sliceThicknessAngstrom;
		}
		BOOST_LOG_TRIVIAL(info) << format("Specimen thickness: %g Angstroms\n") % cztot;
	}

	InitializePropagators(wave);

	for( islice=0; islice < _config->Model.nSlices; islice++ )
	{
		absolute_slice = (m_totalSliceCount+islice);
		Transmit(wave, islice);

		std::stringstream s;
		s << boost::format("after_transmit_%d") % islice;
		wave->Save(s.str());

		//remember: prop must be here to anti-alias propagate is a simple multiplication of wave with prop but it also takes care of the bandwidth limiting
		wave->ToFourierSpace();

		Propagate(wave, islice);
		CollectIntensity(absolute_slice);

		wave->ToRealSpace();
		fft_normalize(wave);


		// Call any additional saving/post-processing that should occur on a per-slice basis
		PostSliceProcess(absolute_slice);
		BOOST_LOG_TRIVIAL(info) << format("Slice %d of %d finished.")
				% (islice+1) % _config->Model.nSlices;
	} /* end for(islice...) */
	// collect intensity at the final slice
	//collectIntensity(muls, wave, m_totalSliceCount+m_slices*(1+mRepeat));
	BOOST_LOG_TRIVIAL(info) << "***************************************";
	if ((m_saveLevel > 1) || (_config->Potential.NSubSlabs > 1)) {
		wave->WriteWave();
	}
	return 0;
}  // end of runMulsSTEM

/******************************************************************
 * propagate_slow()
 * Propagates a wave
 *****************************************************************/
void CExperimentBase::Propagate(WavePtr wave, float_tt dz)
{
	int ixa, iya;
	float_tt wr, wi, tr, ti;
	float_tt scale,t;
	float_tt dzs=0;

	float_tt dx, dy;
	unsigned nx, ny, px;

	wave->GetResolution(dx, dy);
	wave->GetSizePixels(nx, ny);

	complex_tt *w=wave->GetWavePointer();

	px=nx*ny;
#pragma omp parallel for private(wr, wi, tr, ti)
	for (int i=0; i<px; i++)
	{
		try {
			ixa=i%nx;
			iya=i/nx;
			if( wave->GetKX2(ixa) < wave->GetK2Max() ) {
				if( (wave->GetKX2(ixa) + wave->GetKY2(iya)) < wave->GetK2Max() ) {
					wr = w[i].real();
					wi = w[i].imag();
					tr = wr*m_propyr[iya] - wi*m_propyi[iya];
					ti = wr*m_propyi[iya] + wi*m_propyr[iya];
					w[i] = complex_tt( tr*m_propxr[ixa] - ti*m_propxi[ixa], tr*m_propxi[ixa] + ti*m_propxr[ixa]);

				} else
					w[i] = 0.0F;
			} /* end for(iy..) */

			else w[i] = 0.0F;
		} catch (const std::exception &e) {
			std::cerr << e.what();
		}
	} /* end for(ix..) */
} /* end propagate */

/*------------------------ transmit() ------------------------*/
/*
transmit the wavefunction thru one layer 
(simply multiply wave by transmission function)

waver,i[ix][iy]  = real and imaginary parts of wavefunction
transr,i[ix][iy] = real and imag parts of transmission functions

nx, ny = size of array

on entrance waver,i and transr,i are in real space

only waver,i will be changed by this routine
 */
void CExperimentBase::Transmit(WavePtr wave, unsigned sliceIdx) {
	double wr, wi, tr, ti;

	complex_tt *w;
	unsigned nx, ny;

	w = wave->GetWavePointer();
	wave->GetSizePixels(nx, ny);

	/*  trans += posx; */
	for(unsigned ix=0; ix<nx; ix++)
		for(unsigned iy=0; iy<ny; iy++) {
			unsigned offset=ix+nx*iy;
			complex_tt t = m_potential->GetSlicePixel(sliceIdx, ix+m_iPosX, iy+m_iPosY);

			wr = w[offset].real();
			wi = w[offset].imag();
			tr = t.real();
			ti = t.imag();
			w[offset] = complex_tt(wr*tr - wi*ti,wr*ti + wi*tr);

			if(tr != 0 || ti != 0)
				BOOST_LOG_TRIVIAL(trace) << boost::format("w=(%2.3f,%2.3f) t=(%2.3f,%2.3f)") % wr % wi % tr % ti;
		} /* end for(iy.. ix .) */
} /* end transmit() */

void CExperimentBase::AddDPToAvgArray(const WavePtr &wave)
{
	unsigned px=wave->GetTotalPixels();
	// get the pointer to the first data element, and do 1D addressing (it's faster)
	float_tt chisq;

	const float_tt *dp = wave->GetDPPointer();

	for (unsigned i=0; i<px; i++)
	{
		float_tt t=m_avgArray[i]*m_avgCount+dp[i]/(m_avgCount+1);
		chisq+=(m_avgArray[i]-t)*(m_avgArray[i]-t);
		m_avgArray[i]=t;
	}
#pragma omp atomic
	m_chisq[m_avgCount]+=chisq/px;
}

void CExperimentBase::_WriteAvgArray(std::string &fileName, std::string &comment, 
		std::map<std::string, double> &params,
		std::vector<unsigned> &position)
{
	//params["dx"]=1.0/(m_nx*m_dx);
	//params["dy"]=1.0/(m_ny*m_dy);
	params["Thickness"]=m_thickness;
	m_imageIO->WriteImage(m_avgArray, fileName, params, comment, position);
}

void CExperimentBase::ReadAvgArray()
{
	std::vector<unsigned> position;
	m_imageIO->ReadImage(m_avgArray, avgFilePrefix, position);
}

void CExperimentBase::ReadAvgArray(unsigned navg)
{
	std::vector<unsigned> position(1);
	position[0]=navg;
	m_imageIO->ReadImage(m_avgArray, avgFilePrefix, position);
}

void CExperimentBase::ReadAvgArray(unsigned positionx, unsigned positiony)
{
	std::vector<unsigned>position(2);
	position[0]=positionx;
	position[1]=positiony;
	m_imageIO->ReadImage(m_avgArray, avgFilePrefix, position);
}

void CExperimentBase::fft_normalize(WavePtr wave) 
{
	complex_tt *w = wave->GetWavePointer();
	unsigned px = wave->GetTotalPixels();

	float_tt fftScale = 1.0/px;
	for (unsigned i=0; i<px; i++)
	{
		w[i] = complex_tt(w[i].real()*fftScale,w[i].imag() * fftScale);
	}
}

} // end namespace QSTEM
