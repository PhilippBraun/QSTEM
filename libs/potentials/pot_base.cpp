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

#include "pot_base.hpp"
#include <omp.h>
#include <algorithm>
#include <boost/format.hpp>
using boost::format;

const int BUF_LEN = 256;

namespace QSTEM {

CPotential::CPotential() :	IPotential(), m_dkz(0), m_dkx(0) {

}

CPotential::CPotential(const ConfigPtr c) :	CPotential() {
	_config = c;
	m_saveProjectedPotential = c->Output.SaveProjectedPotential;
	m_plotPotential= c->Potential.PlotVrr;
	m_centerSlices = c->Model.CenterSlices;
	m_sliceThickness = c->Model.sliceThicknessAngstrom;
	m_outputInterval=c->Output.PotentialProgressInterval;
	m_zOffset = c->Model.zOffset;
	m_offsetY=c->Model.yOffset;
	m_sliceThicknesses = std::vector<float_tt>();
	// TODO: Read if we should load the pot from a file
	//	configReader->ReadLoadPotential(m_readPotential, m_potFileBase);
	m_ddx = _config->Model.dx / (double) OVERSAMPLING;
	m_ddy = _config->Model.dy / (double) OVERSAMPLING;
	m_ddz = m_dz / (double) OVERSAMPLINGZ;

	/* For now we don't care, if the box has only small
	 * prime factors, because we will not fourier transform it
	 * especially not very often.
	 */
	m_boxNx = (int) (_config->Potential.AtomRadiusAngstrom / m_ddx + 2.0);
	m_boxNy = (int) (_config->Potential.AtomRadiusAngstrom / m_ddy + 2.0);

	m_c = m_sliceThickness * _config->Model.nSlices;
	m_dr = _config->Model.dx / OVERSAMPLING; // define step width in which radial V(r,z) is defined
	m_iRadX = (int) ceil(_config->Potential.AtomRadiusAngstrom / _config->Model.dx);
	m_iRadY = (int) ceil(_config->Potential.AtomRadiusAngstrom / _config->Model.dy);
	m_iRadZ = (int) ceil(_config->Potential.AtomRadiusAngstrom / m_sliceThickness);
	m_iRad2 = m_iRadX * m_iRadX + m_iRadY * m_iRadY;
	m_atomRadius2 = _config->Potential.AtomRadiusAngstrom * _config->Potential.AtomRadiusAngstrom;

	m_sliceThicknesses.resize(_config->Model.nSlices);
	m_slicePos.resize(_config->Model.nSlices);

	if (m_sliceThickness == 0)
		m_sliceThicknesses[0] = m_c / (float_tt) _config->Model.nSlices;
	else
		m_sliceThicknesses[0] = m_sliceThickness;

	m_slicePos[0] = m_zOffset;
	m_divCount = -1;

	m_imageIO = ImageIOPtr(new CImageIO(_config->Model.nx,_config->Model.ny,".img",".img"));
}

CPotential::~CPotential() {

}

void CPotential::SetStructure(StructurePtr structure) {
	m_crystal = structure;
	m_atoms = &(m_crystal->m_atoms);
}

void CPotential::DisplayParams() {
	BOOST_LOG_TRIVIAL(info)<< "*****************************************************";
	BOOST_LOG_TRIVIAL(info)<<"********* Potential Parameters **********************";
	BOOST_LOG_TRIVIAL(info)<< "*****************************************************";
	BOOST_LOG_TRIVIAL(info)<<format("* Print level:          %d") % _config->Output.LogLevel;
	BOOST_LOG_TRIVIAL(info)<<format("* Save level:           %d") % static_cast<int>(_config->Output.SaveLevel);
	if (_config->Output.SavePotential)
		BOOST_LOG_TRIVIAL(info)<<format("* Potential file name:  %s") % m_fileBase.c_str();
	BOOST_LOG_TRIVIAL(info)<<format("* Model Sampling:  %g x %g x %g A")
													% _config->Model.dx% _config->Model.dy% m_sliceThickness;

	BOOST_LOG_TRIVIAL(info)<<format("* Pot. array offset:    (%g,%g,%g)A")
													% _config->Model.xOffset%m_offsetY%	m_zOffset;
	BOOST_LOG_TRIVIAL(info)<<format("* Potential periodic:   (x,y): %s, z: %s")
													%((m_periodicXY) ? "yes" : "no") % ((m_periodicZ) ? "yes" : "no");
	if (k_fftMeasureFlag == FFTW_MEASURE)
		BOOST_LOG_TRIVIAL(info)<<format("* Potential array:      %d x %d (optimized)"), _config->Model.nx, _config->Model.ny;
	else
		BOOST_LOG_TRIVIAL(info)<<format("* Potential array:      %d x %d (estimated)")% _config->Model.nx% _config->Model.ny;
	BOOST_LOG_TRIVIAL(info)<<format("*                       %g x %gA")
											% (_config->Model.nx * _config->Model.dx)%	(_config->Model.ny * _config->Model.dy);
	BOOST_LOG_TRIVIAL(info)<<format("* Scattering factors:   %d")% m_scatFactor;

	BOOST_LOG_TRIVIAL(info)<<format("* Slices per division:  %d (%gA thick slices [%scentered])")
													% _config->Model.nSlices% m_sliceThickness% ((m_centerSlices) ? "" : "not ");
}

/****************************************************************************
 * function: atomBoxLookUp - looks up potential at position x, y, z, relative to atom center
 *
 * Znum = element
 * x,y,z = real space position (in A)
 * B = Debye-Waller factor, B=8 pi^2 <u^2>
 ***************************************************************************/
void CPotential::AtomBoxLookUp(complex_tt &val, int Znum, float_tt x,
		float_tt y, float_tt z, float_tt B) {
	int boxNx, boxNy, boxNz;
	float_tt dx, dy, dz, ddx, ddy, ddz;
	float_tt maxRadius2;
	char fileName[256], systStr[256];
	int tZ, tnx, tny, tnz, tzOversample;
	float_tt tdx, tdy, tdz, tv0, tB;
	FILE *fp;
	int numRead = 0, dummy;

	/* initialize all the atoms to non-used */
	if (!m_atomBoxes.count(Znum)) {
		m_atomBoxes[Znum] = atomBoxPtr(new atomBox());
		m_atomBoxes[Znum]->potential = NULL;
		m_atomBoxes[Znum]->rpotential = NULL;
		m_atomBoxes[Znum]->B = -1.0;

		BOOST_LOG_TRIVIAL(trace)<<format("Atombox has real space resolution of %g x %g x %gA (%d x %d x %d pixels)")
															%ddx% ddy% ddz% boxNx% boxNy% boxNz;
	}
	// printf("Debugging: %d %g %g: %g",Znum,m_atomBoxes[Znum]->B,B,fabs(m_atomBoxes[Znum]->B - B));

	/* Creating/Reading a atombox for every new kind of atom, but only as needed */
	if (fabs(m_atomBoxes[Znum]->B - B) > 1e-6) {
		//  printf("Debugging 1 (%d: %.7g-%.7g= %.7g), %d",
		//	   Znum,m_atomBoxes[Znum]->B,B,fabs(m_atomBoxes[Znum]->B - B),fabs(m_atomBoxes[Znum]->B - B) > 1e-6);
		m_atomBoxes[Znum]->B = B;
		/* Open the file with the projected potential for this particular element
		 */
		sprintf(fileName, "potential_%d_B%d.prj", Znum, (int) (100.0 * B));
		if ((fp = fopen(fileName, "r")) == NULL) {
			sprintf(systStr, "scatpot %s %d %g %d %d %d %g %g %g %d %g",
					fileName, Znum, B, boxNx, boxNy, boxNz, ddx, ddy, ddz,
					OVERSAMPLINGZ, _config->Beam.EnergykeV);
			BOOST_LOG_TRIVIAL(info)<<format("Could not find precalculated potential for Z=%d, will calculate now.")% Znum;
			BOOST_LOG_TRIVIAL(info)<<format("Calling: %s") %systStr;
			system(systStr);
			for (dummy = 0; dummy < 10000; dummy++)
				;
			if ((fp = fopen(fileName, "r")) == NULL) {
				BOOST_LOG_TRIVIAL(error)<<format("cannot calculate projected potential using scatpot - exit!");
				exit(0);
			}

		}
		fgets(systStr, 250, fp);
		sscanf(systStr, "%d %f %d %d %d %f %f %f %d %f", &tZ, &tB, &tnx, &tny,
				&tnz, &tdx, &tdy, &tdz, &tzOversample, &tv0);
		/* If the parameters in the file don't match the current ones,
		 * we need to create a new potential file
		 */
		if ((tZ != Znum) || (fabs(tB - B) > 1e-6) || (tnx != boxNx)
				|| (tny != boxNy) || (tnz != boxNz) || (fabs(tdx - ddx) > 1e-5)
				|| (fabs(tdy - ddy) > 1e-5) || (fabs(tdz - ddz) > 1e-5)
				|| (tzOversample != OVERSAMPLINGZ) || (tv0 != _config->Beam.EnergykeV)) {
			BOOST_LOG_TRIVIAL(info)<<format("Potential input file %s has the wrong parameters") %fileName;
			BOOST_LOG_TRIVIAL(info) << format(	"Parameters:");
			BOOST_LOG_TRIVIAL(info) << format("file:    Z=%d, B=%.3f A^2 (%d, %d, %d) (%.7f, %.7f %.7f) nsz=%d V=%g")
												%tZ% tB% tnx% tny% tnz% tdx% tdy% tdz% tzOversample% tv0;
			BOOST_LOG_TRIVIAL(info) << format("program: Z=%d, B=%.3f A^2 (%d, %d, %d) (%.7f, %.7f %.7f) nsz=%d V=%g")
												%Znum% B% boxNx% boxNy% boxNz% ddx% ddy% ddz%(OVERSAMPLINGZ)% _config->Beam.EnergykeV;
			BOOST_LOG_TRIVIAL(info) << format("will create new potential file, please wait ...");

			/* Close the old file, Create a new potential file now
			 */
			fclose(fp);
			sprintf(systStr, "scatpot %s %d %g %d %d %d %g %g %g %d %g",
					fileName, Znum, B, boxNx, boxNy, boxNz, ddx, ddy, ddz,
					OVERSAMPLINGZ, _config->Beam.EnergykeV);
			system(systStr);
			if ((fp = fopen(fileName, "r")) == NULL) {
				BOOST_LOG_TRIVIAL(error)<<format("cannot calculate projected potential using scatpot - exit!");
				exit(0);
			}
			fgets(systStr, 250, fp);
		}

		/* Finally we can read in the projected potential
		 */
		if (B == 0) {
			m_atomBoxes[Znum]->rpotential = float3D(boxNz, boxNx, boxNy,
					"atomBox");
			numRead = fread(m_atomBoxes[Znum]->rpotential[0][0],
					sizeof(float_tt), (size_t) (boxNx * boxNy * boxNz), fp);
		} else {
			m_atomBoxes[Znum]->potential = complex3D(boxNz, boxNx, boxNy,
					"atomBox");
			numRead = fread(m_atomBoxes[Znum]->potential[0][0],
					sizeof(complex_tt), (size_t) (boxNx * boxNy * boxNz), fp);
		}

		/* writeImage_old(m_atomBoxes[Znum]->potential[0],boxNx,boxNy, 0.0,"potential.img");
		 system("showimage potential.img");
		 */
		fclose(fp);

		if (numRead == boxNx * boxNy * boxNz) {
			BOOST_LOG_TRIVIAL(info)<<format("Sucessfully read in the projected potential");
		} else {
			BOOST_LOG_TRIVIAL(error)<<format("error while reading potential file %s: read %d of %d values")
													%fileName% numRead%( boxNx * boxNy * boxNz);
			exit(0);
		}
	}
}

/** Shuffle the structure with random offsets and recompute potential.
 Only for when you're computing potential, not loading it from files.
 */
void CPotential::Refresh() {
	//
}

void CPotential::ReadPotential(std::string &fileName, unsigned subSlabIdx) {
	/*************************************************************************
	 * read the potential that has been created externally!
	 */
	boost::filesystem::path path = fileName;
	unsigned slice_idx = 0;
	// Set the potential size based on the contents of the first slice
	DataReaderPtr reader = CDataReaderFactory::Get()->GetReader(
			path.extension().string());

	//	TODO: FIX THIS reader->ReadSize(path.stem().string(), slice_idx, _config->Model.nx, _config->Model.ny);
	ResizeSlices();
	for (unsigned i = (subSlabIdx + 1) * _config->Model.nSlices - 1;
			i >= (subSlabIdx) * _config->Model.nSlices; i--, slice_idx++) {
		ReadSlice(path.stem().string(),
				m_trans1[boost::indices[slice_idx][range(0, _config->Model.ny)][range(0,
						_config->Model.nx)]], i);
	}
	return;
}

void CPotential::ReadSlice(const std::string &fileName,ComplexArray2DView slice, unsigned idx) {
}

void CPotential::SliceSetup() {
	FILE *sliceFp = NULL;
	char buf[BUF_LEN];
	/**************************************************************
	 *        setup the slices with their start and end positions
	 *        then loop through all the atoms and add their potential to
	 *        the slice that their potential reaches into (up to RMAX)
	 *************************************************************/

	for (unsigned i = 1; i < _config->Model.nSlices; i++) {
		if (sliceFp == NULL)
			m_sliceThicknesses[i] = m_sliceThicknesses[0];
		/* need to all be the same for fast 3D-FFT method, otherwise OK to be different */
		else {
			fgets(buf, BUF_LEN, sliceFp);
			m_sliceThicknesses[i] = atof(buf);
		}
		// Slice position is last slice position + half width of last slice + half width of this slice
		m_slicePos[i] = m_slicePos[i - 1] + m_sliceThicknesses[i - 1] / 2.0
				+ m_sliceThicknesses[i] / 2.0;
	}
	//	m_trans1.resize(boost::extents[_config->Model.nSlices][_config->Model.ny][_config->Model.nx]);
	// If we are going to read the potential, then we need to size the slices according to the first read pot slice
	if (m_readPotential) {
	}
	// TODO If we are going to calculate the potential, then we need to size the slices according to the size of the
	//    structure and the corresponding resolution.
	else {

	}
	m_trans1.resize(boost::extents[_config->Model.nSlices][_config->Model.ny][_config->Model.nx]);
	//	m_trans.resize(_config->Model.nSlices * _config->Model.nx * _config->Model.ny);
}

complex_tt CPotential::GetSlicePixel(unsigned iz, unsigned ix, unsigned iy) {
	return m_trans1[iz][iy][ix];
}

void CPotential::CenterAtomZ(std::vector<atom>::iterator &atom, float_tt &z) {

	/*
	 * Since cellDiv is always >=1, and divCount starts at 0, the actual position
	 * of this atom with the super cell is given by:
	 */
	/* c*(float_tt)((*muls).cellDiv-divCount-1) will pick the right super-cell
	 * division in the big super-cell
	 * The z-offset 0.5*cz[0] will position atoms at z=0 into the middle of the first
	 * slice.
	 */
	z = atom->z;
	z -= m_c * (float_tt) ((int) m_cellDiv - (int) m_divCount - 1);
	z += m_zOffset;
	z -= (0.5 * m_sliceThickness * (1 - (int) m_centerSlices));
}

/*****************************************************
 * This function will create a 3D potential from whole
 * unit cell, slice it, and make transr/i and propr/i
 * Call this function with center = NULL, if you don't
 * want the array to be shifted.
 ****************************************************/
void CPotential::MakeSlices(int nlayer, StructurePtr crystal) {
	time_t time0, time1;
	/* return, if there is nothing to do */
	if (nlayer < 1)
		return;

	/* we need to keep track of which subdivision of the unit cell we are in
	 * If the cell is not subdivided, then m_cellDiv-1 = 0.
	 */
	//if ((m_divCount == 0) || (m_equalDivs))
	//m_divCount = m_cellDiv;
	//m_divCount--;
	/* we only want to reread and shake the atoms, if we have finished the
	 * current unit cell.
	 */
	//if (m_divCount == m_cellDiv-1) {
	//ReadAtoms();
	//} /* end of if divCount==cellDiv-1 ... */
	SliceSetup();

	std::fill( m_trans1.origin(), m_trans1.origin() + m_trans1.size(), complex_tt(0,0) );

#pragma omp parallel for
	for (std::vector<atom>::iterator atom = m_crystal->m_uniqueAtoms.begin();	atom < m_crystal->m_uniqueAtoms.end(); atom=atom+1) {
		ComputeAtomPotential(atom);
	}

	time(&time0);
	int atomsAdded = 0;


#pragma omp parallel for shared(atomsAdded)
	for (std::vector<atom>::iterator atom = m_crystal->m_atoms.begin(); atom < m_crystal->m_atoms.end(); atom = atom + 1) {


		// make sure we skip vacancies:
		while (atom->Znum == 0)	atom++;
		size_t iatom = atom - m_crystal->m_atoms.begin();

		if ((m_printLevel >= 4) &&
				(m_displayPotCalcInterval > 0)&&
				((iatom + 1) % (m_displayPotCalcInterval)) == 0) {
			BOOST_LOG_TRIVIAL(info)<<format("Adding potential for atom %d (Z=%u, pos=[%.1f, %.1f, %.1f])")%
					(iatom + 1)% atom->Znum% atom->x% atom->y% atom->z;
		}

		// atom coordinates in cartesian coords.  The x- and y-position will be offset by the starting point
		// of the actually needed array of projected potential

		float_tt atomX = atom->x - _config->Model.xOffset;
		float_tt atomY = atom->y - m_offsetY;
		float_tt atomZ = atom->z;

		CenterAtomZ(atom, atomZ);
		AddAtomToSlices(atom, atomX, atomY, atomZ);

#pragma omp critical
		atomsAdded++;
		if(atomsAdded % 50 == 0)
			atomZ = 1;

		if(atomsAdded % 100 == 0) BOOST_LOG_TRIVIAL(info)<<format("%2.1f percent of atoms added") % ((float)atomsAdded/m_atoms->size()*100);

	} // for iatom =0 ... 

	MakePhaseGratings();

	if(_config->Potential.BandlimitTransmissionFunction) BandlimitTransmissionFunction();

	time(&time1);
	BOOST_LOG_TRIVIAL(info) << format(	"%g sec used for real space potential calculation (%g sec per atom)")
							%	difftime(time1, time0)%(	difftime(time1, time0) / m_crystal->m_atoms.size());


	/*************************************************/
	/* Save the potential slices                                         */

	if (_config->Output.SavePotential) {
		for (unsigned iz = 0; iz < nlayer; iz++) {
			// find the maximum value of each layer:
			float_tt potVal = m_trans1[iz][0][0].real();

			if (_config->Output.LogLevel < 2) {
				float_tt ddx = potVal;
				float_tt ddy = potVal;
				for (unsigned ix = 0; ix <  _config->Model.nx; ix++)
					for (unsigned iy = 0; iy < _config->Model.ny ; iy++){
						potVal = m_trans1[iz][iy][ix].real();
						if (ddy < potVal)
							ddy = potVal;
						if (ddx > potVal)
							ddx = potVal;
						BOOST_LOG_TRIVIAL(trace)<<format("m_trans1[%d][%d][%d] = %g")%iz%iy%ix%potVal;
					}
				BOOST_LOG_TRIVIAL(info)<<format("Saving (complex) potential layer %d to file (r: %g..%g)")%iz% ddx% ddy;
			}

			WriteSlice(iz);

			//imageIO->WriteImage((void **)m_trans[iz], params, std::string(buf));
		} // loop through all slices
	} /* end of if savePotential ... */
	if (m_saveProjectedPotential) {
		WriteProjectedPotential();
	}
} // end of make3DSlices
void CPotential::BandlimitTransmissionFunction(){

}
void CPotential::MakePhaseGratings(){
	float_tt mm0 = 1.0F + _config->Beam.EnergykeV/511.0F;
	float_tt scale = mm0*_config->Beam.wavelength;

	BOOST_LOG_TRIVIAL(info)<<format("Making phase gratings for %d layers (scale=%g rad/VA, gamma=%g, sigma=%g) ... ")
			%m_trans1.shape()[0]%scale%mm0%_config->Beam.sigma;

#pragma omp parallel for
	for(int iz =0; iz < m_trans1.shape()[0];iz++){
		for(int iy =0; iy < m_trans1.shape()[1];iy++){
			for(int ix=0; ix < m_trans1.shape()[2];ix++){
				complex_tt t = m_trans1[iz][iy][ix];
				m_trans1[iz][iy][ix] = complex_tt(cos(scale*t.real()),sin(scale*t.imag()));
			}
		}
	}
}



void CPotential::WriteSlice(unsigned idx) {
	char buf[255];
	std::map<std::string, double> params;
	params["Thickness"] = m_sliceThickness;
	params["dx"] = _config->Model.dx;
	params["dy"] = _config->Model.dy;
	sprintf(buf, "Projected Potential (slice %d)", idx);
	std::string comment = buf;
	std::stringstream filename;
	filename << "pot_slice_" << idx;
	m_imageIO->WriteImage(m_trans1[boost::indices[idx][range(0,_config->Model.nx)][range(0,_config->Model.ny)]], filename.str().c_str(), params, comment);
}

void CPotential::WriteProjectedPotential() {
	std::map<std::string, double> params;
	char buf[255];
	RealVector tempPot(_config->Model.ny*_config->Model.nx);
	float_tt potVal = 0;

	for (unsigned idy = 0; idy < _config->Model.ny; idy++)
		for (unsigned idx = 0; idx < _config->Model.nx ; idx++){
			tempPot[idx+_config->Model.nx*idy] = 0;
			for (unsigned iz = 0; iz < _config->Model.nSlices; iz++)
				tempPot[idx+_config->Model.nx*idy] += m_trans1[iz][idy][idx].real();
			//			BOOST_LOG_TRIVIAL(trace)<<format("temppot = %g") % tempPot[idx+_config->Model.nx*idy];
		}

	float_tt ddx = tempPot[0], ddy = potVal;
	for (unsigned ix = 0; ix < _config->Model.ny * _config->Model.nx; potVal = tempPot[++ix]) {
		if (ddy < potVal)
			ddy = potVal;
		if (ddx > potVal)
			ddx = potVal;
	}
	BOOST_LOG_TRIVIAL(info)<<format("Saving total projected potential to file (r: %g..%g)")%ddx%ddy;
	params["Thickness"] = m_sliceThickness;
	params["dx"] = _config->Model.dx;
	params["dy"] = _config->Model.dy;
	sprintf(buf, "Projected Potential (sum of %d slices)", _config->Model.nSlices);
	std::string comment = buf;
	std::string fileName = "pot_projected";
	m_imageIO->WriteImage(tempPot, fileName, params, comment);
}

/*
 // TODO: this was taken from stemutils.  It seems to be used only in customslice, which then isn't used anywhere.
 float_tt CPotential::sfLUT(float_tt s,int atKind)
 {
 int i;
 double sf;
 static double *splinx=NULL;
 static double **spliny=NULL;
 static double **splinb=NULL;
 static double **splinc=NULL;
 static double **splind=NULL;
 static int sfSize = 0;
 static int atKinds = 0;
 static double maxK = 0;

 if(splinx == NULL) {
 // splinx = s;
 // spliny = sfC;
 sfSize = m_sfNk;
 splinx = m_sfkArray;
 spliny = m_sfTable;
 atKinds = m_atoms.size();
 splinb = double2D(atKinds,sfSize, "splinb" );
 splinc = double2D(atKinds,sfSize, "splinc" );
 splind = double2D(atKinds,sfSize, "splind" );
 maxK = splinx[sfSize-1];

 for (i=0;i<atKinds;i++)
 splinh(splinx,spliny[i],splinb[i],splinc[i],splind[i],sfSize);
 }


 // now that everything is set up find the
 //   scattering factor by interpolation in the table

 if (s > maxK) return 0.0;
 if (atKind < atKinds) {
 sf = seval(splinx,spliny[atKind],splinb[atKind],splinc[atKind],splind[atKind],sfSize,s);
 if (sf < 0) return 0.0;
 return(sf);
 }
 printf("sfLUT: invalid atom kind (%d) - exit!\n",atKind);
 exit(0);
 }  // end sfLUT()
 */

/*------------------ splinh() -----------------------------*/
/*
 fit a quasi-Hermite cubic spline

 [1] Spline fit as in H.Akima, J. ACM 17(1970)p.589-602
 'A New Method of Interpolation and Smooth
 Curve Fitting Based on Local Procedures'

 [2] H.Akima, Comm. ACM, 15(1972)p.914-918

 E. Kirkland 4-JUL-85
 changed zero test to be a small nonzero number 8-jul-85 ejk
 converted to C 24-jun-1995 ejk

 The inputs are:
 x[n] = array of x values in ascending order, each X(I) must
 be unique
 y[n] = array of y values corresponding to X(N)
 n = number of data points must be 2 or greater

 The outputs are (with z=x-x(i)):
 b[n] = array of spline coeficients for (x-x[i])
 c[n] = array of spline coeficients for (x-x[i])**2
 d[n] = array of spline coeficients for (x-x[i])**3
 ( x[i] <= x <= x[i+1] )
 To interpolate y(x) = yi + bi*z + c*z*z + d*z*z*z

 The coeficients b[i], c[i], d[i] refer to the x[i] to x[i+1]
 interval. NOTE that the last set of coefficients,
 b[n-1], c[n-1], d[n-1] are meaningless.
 */
void CPotential::splinh(float_tt x[], float_tt y[],std::vector<float_tt>& b,std::vector<float_tt>& c,
		std::vector<float_tt>& d, int n) {
#define SMALL 1.0e-25

	int i, nm1, nm4;
	float_tt m1, m2, m3, m4, m5, t1, t2, m54, m43, m32, m21, x43;

	if (n < 4)
		return;

	/* Do the first end point (special case),
	 and get starting values */

	m5 = (y[3] - y[2]) / (x[3] - x[2]); /* mx = slope at pt x */
	m4 = (y[2] - y[1]) / (x[2] - x[1]);
	m3 = (y[1] - y[0]) / (x[1] - x[0]);

	m2 = m3 + m3 - m4; /* eq. (9) of reference [1] */
	m1 = m2 + m2 - m3;

	m54 = fabs(m5 - m4);
	m43 = fabs(m4 - m3);
	m32 = fabs(m3 - m2);
	m21 = fabs(m2 - m1);

	if ((m43 + m21) > SMALL)
		t1 = (m43 * m2 + m21 * m3) / (m43 + m21);
	else
		t1 = 0.5 * (m2 + m3);

	/* Do everything up to the last end points */

	nm1 = n - 1;
	nm4 = n - 4;

	for (i = 0; i < nm1; i++) {

		if ((m54 + m32) > SMALL)
			t2 = (m54 * m3 + m32 * m4) / (m54 + m32);
		else
			t2 = 0.5 * (m3 + m4);

		x43 = x[i + 1] - x[i];
		b[i] = t1;
		c[i] = (3.0 * m3 - t1 - t1 - t2) / x43;
		d[i] = (t1 + t2 - m3 - m3) / (x43 * x43);

		m1 = m2;
		m2 = m3;
		m3 = m4;
		m4 = m5;
		if (i < nm4) {
			m5 = (y[i + 4] - y[i + 3]) / (x[i + 4] - x[i + 3]);
		} else {
			m5 = m4 + m4 - m3;
		}

		m21 = m32;
		m32 = m43;
		m43 = m54;
		m54 = fabs(m5 - m4);
		t1 = t2;
	}

	return;

} /* end splinh() */

/*----------------------- seval() ----------------------*/
/*
 Interpolate from cubic spline coefficients

 E. Kirkland 4-JUL-85
 modified to do a binary search for efficiency 13-Oct-1994 ejk
 converted to C 26-jun-1995 ejk
 fixed problem on end-of-range 16-July-1995 ejk

 The inputs are:
 x[n] = array of x values in ascending order, each x[i] must
 be unique
 y[n] = array of y values corresponding to x[n]
 b[n] = array of spline coeficients for (x-x[i])
 c[n] = array of spline coeficients for (x-x[i])**2
 d[n] = array of spline coeficients for (x-x[i])**3
 n = number of data points
 x0 = the x value to interpolate at
 (x[i] <= x <= x[i+1]) and all inputs remain unchanged

 The value returned is the interpolated y value.

 The coeficients b[i], c[i], d[i] refer to the x[i] to x[i+1]
 interval. NOTE that the last set of coefficients,
 b[n-1], c[n-1], d[n-1] are meaningless.
 */
float_tt CPotential::seval(float_tt *x, float_tt *y,std::vector<float_tt>& b,std::vector<float_tt>& c,
		std::vector<float_tt>& d, int n, float_tt x0) {
	int i, j, k;
	float_tt z, seval1;

	/* exit if x0 is outside the spline range */
	if (x0 <= x[0])
		i = 0;
	else if (x0 >= x[n - 2])
		i = n - 2;
	else {
		i = 0;
		j = n;
		do {
			k = (i + j) / 2;
			if (x0 < x[k])
				j = k;
			else if (x0 >= x[k])
				i = k;
		} while ((j - i) > 1);
	}

	z = x0 - x[i];
	seval1 = y[i] + (b[i] + (c[i] + d[i] * z) * z) * z;

	return (seval1);

} /* end seval() */

void CPotential::GetSizePixels(unsigned int &nx, unsigned int &ny) const {
	nx = _config->Model.nx;
	ny = _config->Model.ny;
}

void CPotential::ResizeSlices() {}

} // end namespace QSTEM
