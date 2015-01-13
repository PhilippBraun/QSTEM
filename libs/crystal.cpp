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


#include <stdlib.h>
#include "crystal.hpp"
#include "memory_fftw3.hpp"
#include "random.hpp"
#include "structure_factories.hpp"
#include <string>

#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/array.hpp>
using boost::format;

#define PI180 1.7453292519943e-2
#define THZ_AMU_HBAR 0.15745702964189 /* A°^2*THz*amu/(hbar) */
// 4.46677327584453 /* 1e10/sqrt(THz*amu/(pi*hbar)) */
#define THZ_HBAR_KB 1.90963567802059 /* THz*hbar/kB */


namespace QSTEM {

static const float_tt k_wobScale = 1.0 / (8 * M_PI * M_PI);
static const float_tt k_sq3 = 1.0 / sqrt(3.0); /* sq3 is an additional needed factor which stems from
 * int_-infty^infty exp(-x^2/2) x^2 dx = sqrt(pi)
 * introduced in order to match the wobble factor with <u^2>
 */

CCrystal::CCrystal() :
		m_minX(0),
		m_maxX(0),
		m_minY(0),
		m_maxY(0),
		m_minZ(0),
		m_maxZ(0),
		m_Einstein(true),
		m_phononFile(boost::filesystem::path()),
		m_cubex(0),
		m_cubey(0),
		m_cubez(0),
		m_ax(0)
		{}

CCrystal::CCrystal(const ConfigPtr c) : CCrystal()  {
	_config = c;

	m_wobble_temp_scale = sqrt(_config->Structure.temperatureK / 300.0);
	// Get the object that we'll use to read in the array of atoms
	m_structureReader = CStructureReaderFactory::Get()->GetReader(c->Structure.structureFilename);
	m_structureWriter = CStructureWriterFactory::Get()->GetWriter(c->Structure.structureFilename,m_ax, m_by, m_cz);
	m_Mm = float2D(3, 3, "");
	m_structureReader->ReadCellParams(m_Mm);
	m_structureReader->ReadAtoms(m_baseAtoms,m_uniqueAtoms);
	BOOST_LOG_TRIVIAL(info)<<format("Read %d atoms, tds: %d") % m_baseAtoms.size()% _config->Model.UseTDS;

	CalculateCellDimensions();
	MakeCrystal(true);
	CalculateCrystalBoundaries();
}

CCrystal::CCrystal(unsigned ncx, unsigned ncy, unsigned ncz, float_tt tx, float_tt ty, float_tt tz){}

CCrystal::~CCrystal() {}

void CCrystal::CalculateCrystalBoundaries() {
	int largestZindex = -1;
#pragma omp parallel for
	for (int i = 0; i < m_atoms.size(); i++) {
#pragma omp critical
		{
			if (m_atoms[i].x < m_minX)
				m_minX = m_atoms[i].x;
		}
#pragma omp critical
		{
			if (m_atoms[i].x > m_maxX)
				m_maxX = m_atoms[i].x;
		}
#pragma omp critical
		{
			if (m_atoms[i].y < m_minY)
				m_minY = m_atoms[i].y;
		}
#pragma omp critical
		{
			if (m_atoms[i].y > m_maxY)
				m_maxY = m_atoms[i].y;
		}
#pragma omp critical
		{
			if (m_atoms[i].z < m_minZ)
				m_minZ = m_atoms[i].z;
		}
#pragma omp critical
		{
			if (m_atoms[i].z > m_maxZ){
				m_maxZ = m_atoms[i].z;
				largestZindex = i;
			}
		}
	}
	BOOST_LOG_TRIVIAL(trace)<<format("largest Z at index %d : (%f,%f,%f)")
												%largestZindex%m_atoms[largestZindex].x%m_atoms[largestZindex].y%m_atoms[largestZindex].z;
}

void CCrystal::Init(unsigned run_number) {
	BOOST_LOG_TRIVIAL(trace)<<format("range of thermally displaced atoms (%d atoms): ") %m_atoms.size();
	BOOST_LOG_TRIVIAL(trace)<<format("X: %g .. %g")% m_minX% m_maxX;
	BOOST_LOG_TRIVIAL(trace)<<format("Y: %g .. %g")% m_minY% m_maxY;
	BOOST_LOG_TRIVIAL(trace)<<format("Z: %g .. %g")% m_minZ% m_maxZ;

	qsort(&m_atoms[0], m_atoms.size(), sizeof(atom),&CCrystal::AtomCompareZnum);
	WriteStructure(run_number);
}

void CCrystal::DisplayParams() {
	BOOST_LOG_TRIVIAL(info)<<format("* Input file:           %s") %  m_structureFile.c_str();

	if ((m_cubex == 0) || (m_cubey == 0) || (m_cubez == 0))
		BOOST_LOG_TRIVIAL(info)<<format("* Unit cell:            ax=%g by=%g cz=%g")% m_ax% m_by% m_cz;
	else {
		BOOST_LOG_TRIVIAL(info)<<format("* Size of Cube:         ax=%g by=%g cz=%g")% m_cubex% m_cubey%	m_cubez;
		BOOST_LOG_TRIVIAL(info)<<format("* Cube size adjusted:   %s")%( m_adjustCubeSize ? "yes" : "no");
	}
	BOOST_LOG_TRIVIAL(info)<<format("* Super cell:           %d x %d x %d unit cells")
			%_config->Structure.nCellX%		_config->Structure.nCellY% _config->Structure.nCellZ;
	BOOST_LOG_TRIVIAL(info)<<format("* Number of atoms:      %d (super cell)")% m_atoms.size();
	BOOST_LOG_TRIVIAL(info)<<format("* Crystal tilt:         x=%g deg, y=%g deg, z=%g deg")
			%(_config->Model.crystalTiltX* RAD2DEG)% (_config->Model.crystalTiltY * RAD2DEG)%( _config->Model.crystalTiltZ * RAD2DEG);
	BOOST_LOG_TRIVIAL(info)<<format("* Model dimensions:     ax=%gA, by=%gA, cz=%gA (after tilt)")
			%m_ax% m_by% m_cz;
	BOOST_LOG_TRIVIAL(info)<<format("* Temperature:          %gK") %_config->Structure.temperatureK;
	if (_config->Model.UseTDS)
		BOOST_LOG_TRIVIAL(info)<<format("* TDS:                  yes");
	else
		BOOST_LOG_TRIVIAL(info)<<format("* TDS:                  no");
}

void CCrystal::SetCellParameters(float_tt ax, float_tt by, float_tt cz) {
	m_ax = ax;
	m_by = by;
	m_cz = cz;
}

void CCrystal::SetNCells(unsigned nx, unsigned ny, unsigned nz) {
	_config->Structure.nCellX = nx;
	_config->Structure.nCellY = ny;
	_config->Structure.nCellZ = nz;
	// TODO: should this ever be false?
	bool handleVacancies = true;
	ReplicateUnitCell(handleVacancies);
}

void CCrystal::GetCellAngles(float_tt &alpha, float_tt &beta, float_tt &gamma) {
	alpha = m_cAlpha;
	beta = m_cBeta;
	gamma = m_cGamma;
}

void CCrystal::GetCellParameters(float_tt &ax, float_tt &by, float_tt &cz) {
	ax = m_ax;
	by = m_by;
	cz = m_cz;
}

void CCrystal::OffsetCenter(atom &center) {
	/*
	 define the center of our unit cell by moving the atom specified
	 by "center" at position (0.5,0.5,0.0)
	 */

	float_tt dx = m_ax / 2.0f - center.x;
	float_tt dy = m_by / 2.0f - center.y;
	float_tt dz = -center.z;
	for (size_t i = 0; i < m_atoms.size(); i++) {
		m_atoms[i].x += dx;
		if (m_atoms[i].x < 0.0f)
			m_atoms[i].x += m_ax;
		else if (m_atoms[i].x > m_ax)
			m_atoms[i].x -= m_ax;
		m_atoms[i].y += dy;
		if (m_atoms[i].y < 0.0f)
			m_atoms[i].y += m_by;
		else if (m_atoms[i].y > m_by)
			m_atoms[i].y -= m_by;
		m_atoms[i].z += dz;
		if (m_atoms[i].z < 0.0f)
			m_atoms[i].z += m_cz;
		else if (m_atoms[i].z > m_cz)
			m_atoms[i].z -= m_cz;
	}
}

// This function reads the atomic positions from fileName and also adds 
// Thermal displacements to their positions, if _config->Model.UseTDS is turned on.
void CCrystal::MakeCrystal(bool handleVacancies) {
	int ncoord = m_baseAtoms.size(), ncx, ncy, ncz, icx, icy, icz, jz;
	// float_t dw,occ,dx,dy,dz,r;
	int i2, j, ix, iy, iz = 0;
	// char s1[16],s2[16],s3[16];
	float_tt boxXmin = 0, boxXmax = 0, boxYmin = 0, boxYmax = 0, boxZmin = 0,
			boxZmax = 0;
	float_tt boxCenterX, boxCenterY, boxCenterZ, boxCenterXrot, boxCenterYrot,
	boxCenterZrot, bcX, bcY, bcZ;
	float_tt totOcc;
	float_tt choice, lastOcc;
	float_tt *u = NULL;
	float_tt **Mm = m_Mm;
	static int ncoord_old = 0;
	u = (float_tt *) malloc(3 * sizeof(float_tt));
	ncx = _config->Structure.nCellX;
	ncy = _config->Structure.nCellY;
	ncz = _config->Structure.nCellZ;

	BOOST_LOG_TRIVIAL(trace)<<format("Lattice parameters: ax=%g by=%g cz=%g (%d atoms)")
							% m_ax % m_by % m_cz % ncoord;
	if ((m_cubex == 0) || (m_cubey == 0) || (m_cubez == 0))
		BOOST_LOG_TRIVIAL(trace)<<format("Size of Super-lattice: ax=%g by=%g cz=%g (%d x %d x %d)")
		% (m_ax * ncx)% (m_by * ncy)%( m_cz * ncz)% ncx% ncy% ncz;
	else
		BOOST_LOG_TRIVIAL(trace)<<format("Size of Cube: ax=%g by=%g cz=%g") % m_cubex% m_cubey%	m_cubez;

	unsigned natom = m_baseAtoms.size() * ncx * ncy * ncz;
	m_atoms.resize(natom);

	if (handleVacancies) {
		qsort(&m_baseAtoms[0], ncoord, sizeof(atom),&CCrystal::AtomCompareZYX);
	}
	/***********************************************************
	 * Read actual Data
	 ***********************************************************/
	for (int i = m_baseAtoms.size() - 1; i >= 0; i--) {
		if (_config->Model.UseTDS) {
			m_u2[m_baseAtoms[i].Znum] = 0;
		}
	}
	if ((m_cubex > 0) && (m_cubey > 0) && (m_cubez > 0)) {
		/* at this point the atoms should have fractional coordinates */
		// printf("Entering tiltBoxed");
		TiltBoxed(ncoord, handleVacancies);
		// printf("ncoord: %d, natom: %d\n",ncoord,*natom);
	} else {  // work in NCell mode
		// atoms are in fractional coordinates so far, we need to convert them to
		// add the phonon displacement in this condition, because there we can
		// actually do the correct Eigenmode treatment.
		// but we will probably just do Einstein vibrations anyway:
		ReplicateUnitCell(handleVacancies);
		natom = ncoord * ncx * ncy * ncz;

		BOOST_LOG_TRIVIAL(trace) << "Atoms after replication of unit cells";

		for (int j = 0; j < natom; j++) {
			//			LOG(INFO) << format("atom %d: (%3.3f, %3.3f, %3.3f)\n") % j % m_atoms[j].x % m_atoms[j].y % m_atoms[j].z;
			// This converts also to cartesian coordinates
			float_tt x = Mm[0][0] * m_atoms[j].x + Mm[1][0] * m_atoms[j].y + Mm[2][0] * m_atoms[j].z;
			float_tt y = Mm[0][1] * m_atoms[j].x + Mm[1][1] * m_atoms[j].y + Mm[2][1] * m_atoms[j].z;
			float_tt z = Mm[0][2] * m_atoms[j].x + Mm[1][2] * m_atoms[j].y + Mm[2][2] * m_atoms[j].z;

			m_atoms[j].x = x;
			m_atoms[j].y = y;
			m_atoms[j].z = z;

			BOOST_LOG_TRIVIAL(trace) << format("atom %d: (%3.3f, %3.3f, %3.3f)") % j % m_atoms[j].x % m_atoms[j].y % m_atoms[j].z;
		}

		/***************************************************************
		 * Now let us tilt around the center of the full crystal
		 */

		bcX = ncx / 2.0;
		bcY = ncy / 2.0;
		bcZ = ncz / 2.0;
		u[0] = Mm[0][0] * bcX + Mm[1][0] * bcY + Mm[2][0] * bcZ;
		u[1] = Mm[0][1] * bcX + Mm[1][1] * bcY + Mm[2][1] * bcZ;
		u[2] = Mm[0][2] * bcX + Mm[1][2] * bcY + Mm[2][2] * bcZ;
		boxCenterX = u[0];
		boxCenterY = u[1];
		boxCenterZ = u[2];

		// rotateVect(u,u,m_ctiltx,m_ctilty,m_ctiltz);  // simply applies rotation matrix
		// boxCenterXrot = u[0]; boxCenterYrot = u[1];	boxCenterZrot = u[2];
		// Determine the size of the (rotated) super cell
		for (icx = 0; icx <= ncx; icx += ncx)
			for (icy = 0; icy <= ncy; icy += ncy)
				for (icz = 0; icz <= ncz; icz += ncz) {
					u[0] = Mm[0][0] * (icx - bcX) + Mm[1][0] * (icy - bcY)+ Mm[2][0] * (icz - bcZ);
					u[1] = Mm[0][1] * (icx - bcX) + Mm[1][1] * (icy - bcY)+ Mm[2][1] * (icz - bcZ);
					u[2] = Mm[0][2] * (icx - bcX) + Mm[1][2] * (icy - bcY)+ Mm[2][2] * (icz - bcZ);
					if ((_config->Model.crystalTiltX != 0) || (_config->Model.crystalTiltY != 0) || (_config->Model.crystalTiltZ != 0))
						RotateVect(u, u, _config->Model.crystalTiltX, _config->Model.crystalTiltY, _config->Model.crystalTiltZ); // simply applies rotation matrix
					// x = u[0]+boxCenterXrot; y = u[1]+boxCenterYrot; z = u[2]+boxCenterZrot;
					float_tt x = u[0] + boxCenterX;
					float_tt y = u[1] + boxCenterY;
					float_tt z = u[2] + boxCenterZ;
					if ((icx == 0) && (icy == 0) && (icz == 0)) {
						boxXmin = boxXmax = x;
						boxYmin = boxYmax = y;
						boxZmin = boxZmax = z;
					} else {
						boxXmin = boxXmin > x ? x : boxXmin;
						boxXmax = boxXmax < x ? x : boxXmax;
						boxYmin = boxYmin > y ? y : boxYmin;
						boxYmax = boxYmax < y ? y : boxYmax;
						boxZmin = boxZmin > z ? z : boxZmin;
						boxZmax = boxZmax < z ? z : boxZmax;
					}

					BOOST_LOG_TRIVIAL(info) << format("u=(%2.3f,%2.3f,%2.3f) x-(%2.3f,%2.3f) y-(%2.3f,%2.3f) z-(%2.3f,%2.3f)") % u[0] % u[1] % u[2] % boxXmin % boxXmax % boxYmin % boxYmax % boxZmin % boxZmax;

				}

		//		printf("(%f, %f, %f): %f .. %f, %f .. %f, %f .. %f\n",m_ax,m_by,m_cz,boxXmin,boxXmax,boxYmin,boxYmax,boxZmin,boxZmax);
		if ((_config->Model.crystalTiltX != 0) || (_config->Model.crystalTiltY != 0) || (_config->Model.crystalTiltZ != 0)) {
			for (int j =  0; j < (natom); j++) {
				u[0] = m_atoms[j].x - boxCenterX;
				u[1] = m_atoms[j].y - boxCenterY;
				u[2] = m_atoms[j].z - boxCenterZ;
				RotateVect(u, u, _config->Model.crystalTiltX, _config->Model.crystalTiltY, _config->Model.crystalTiltZ); // simply applies rotation matrix
				u[0] += boxCenterX;
				u[1] += boxCenterY;
				u[2] += boxCenterZ;
				m_atoms[j].x = u[0];
				m_atoms[j].y = u[1];
				m_atoms[j].z = u[2];
			}
		} /* if tilts != 0 ... */

		// rebase to some atom at (0,0,0)
		BOOST_LOG_TRIVIAL(trace) << "Atoms after tilting and going to cartesian";
		for (int j = 0; j < natom; j++) {
			//			BOOST_LOG_TRIVIAL(trace) << format("atom %d: (%3.3f, %3.3f, %3.3f)\n") % j % m_atoms[j].x % m_atoms[j].y % m_atoms[j].z;
			m_atoms[j].x -= boxXmin;
			m_atoms[j].y -= boxYmin;
			m_atoms[j].z -= boxZmin;
			BOOST_LOG_TRIVIAL(trace) << format("atom %d: (%3.3f, %3.3f, %3.3f)") % j % m_atoms[j].x % m_atoms[j].y % m_atoms[j].z;
		}

		// Offset the atoms in x- and y-directions:
		// Do this after the rotation!
		if ((_config->Model.xOffset != 0) || (_config->Model.yOffset != 0)) {
			for (int j = 0 ; j < natom; j++) {
				m_atoms[j].x += _config->Model.xOffset;
				m_atoms[j].y += _config->Model.yOffset;
			}
		}
	} // end of Ncell mode conversion to cartesian coords and tilting.
	// end of loop over atoms

	_sizeX = boxXmax - boxXmin;
	_sizeY = boxYmax - boxYmin;
	_sizeZ = boxZmax - boxZmin;

	free(u);
}

// Uses m_Mm to calculate ax, by, cz, and alpha, beta, gamma
void CCrystal::CalculateCellDimensions() {
	m_ax = sqrt(m_Mm[0][0] * m_Mm[0][0] + m_Mm[0][1] * m_Mm[0][1]+ m_Mm[0][2] * m_Mm[0][2]);
	m_by = sqrt(m_Mm[1][0] * m_Mm[1][0] + m_Mm[1][1] * m_Mm[1][1]+ m_Mm[1][2] * m_Mm[1][2]);
	m_cz = sqrt(m_Mm[2][0] * m_Mm[2][0] + m_Mm[2][1] * m_Mm[2][1]+ m_Mm[2][2] * m_Mm[2][2]);
	m_cGamma = atan2(m_Mm[1][1], m_Mm[1][0]);
	m_cBeta = acos(m_Mm[2][0] / m_cz);
	m_cAlpha = acos(m_Mm[2][1] * sin(m_cGamma) / m_cz + cos(m_cBeta) * cos(m_cGamma));
	m_cGamma /= (float) PI180;
	m_cBeta /= (float) PI180;
	m_cAlpha /= (float) PI180;
}

// TODO: old version return a pointer to new atom positions.  Can we do this in place?
//		If not, just set atoms to the new vector.
void CCrystal::TiltBoxed(int ncoord, bool handleVacancies) {
	int atomKinds = 0;
	int iatom, jVac, jequal, jChoice, i2, ix, iy, iz, atomCount = 0, atomSize;
	static float_tt *axCell, *byCell, *czCell = NULL;
	static float_tt **Mm = NULL, **Mminv = NULL, **MpRed = NULL, **MpRedInv =
			NULL;
	static float_tt **MbPrim = NULL, **MbPrimInv = NULL, **MmOrig = NULL,
			**MmOrigInv = NULL;
	//static float_tt **a = NULL,**aOrig = NULL,**b= NULL,**bfloor=NULL,**blat=NULL;
	std::vector<float_tt> a(3, 0), aOrig(3, 0), b(3, 0), bfloor(3, 0), blat(3, 0);
	//static float_tt *uf;
	static int oldAtomSize = 0;
	double x, y, z, dx, dy, dz;
	double totOcc, lastOcc, choice;
	atom newAtom;
	int Ncells;
	std::vector<float_tt> u(3, 0), uf(3, 0);
	//static float_tt *u;

	//unsigned jz;

	Ncells = _config->Structure.nCellX * _config->Structure.nCellY * _config->Structure.nCellZ;

	if (Mm == NULL) {
		MmOrig = float2D(3, 3, "MmOrig");
		MmOrigInv = float2D(3, 3, "MmOrigInv");
		MbPrim = float2D(3, 3, "MbPrim");// double version of primitive lattice basis
		MbPrimInv = float2D(3, 3, "MbPrim"); // double version of inverse primitive lattice basis
		MpRed = float2D(3, 3, "MpRed"); /* conversion lattice to obtain red. prim. coords
		 * from reduced cubic rect.
		 */
		MpRedInv = float2D(3, 3, "MpRedInv"); /* conversion lattice to obtain red. cub. coords
		 * from reduced primitive lattice coords
		 */
		Mm = float2D(3, 3, "Mm");
		Mminv = float2D(3, 3, "Mminv");
		axCell = Mm[0];
		byCell = Mm[1];
		czCell = Mm[2];
		//a			= float2D(1,3,"a");
		//aOrig		= float2D(1,3,"aOrig");
		//b			= float2D(1,3,"b");
		//bfloor		= float2D(1,3,"bfloor");
		//blat		= float2D(1,3,"blat");
		//uf			= (float_tt *)malloc(3*sizeof(float_tt));
		//u			= (float_tt *)malloc(3*sizeof(float_tt));
	}

	dx = 0;
	dy = 0;
	dz = 0;
	dx = m_offsetX;
	dy = m_offsetY;
	// We need to copy the transpose of m_Mm to Mm.
	// we therefore cannot use the following command:
	// memcpy(Mm[0],m_Mm[0],3*3*sizeof(double));
	for (ix = 0; ix < 3; ix++)
		for (iy = 0; iy < 3; iy++)
			Mm[ix][iy] = m_Mm[iy][ix];

	memcpy(MmOrig[0], Mm[0], 3 * 3 * sizeof(double));
	Inverse_3x3(MmOrigInv[0], MmOrig[0]);
	/* remember that the angles are in rad: */
	RotateMatrix(Mm[0], Mm[0], _config->Model.crystalTiltX, _config->Model.crystalTiltY, _config->Model.crystalTiltZ);
	Inverse_3x3(Mminv[0], Mm[0]);  // computes Mminv from Mm!
	/* find out how far we will have to go in unit of unit cell vectors.
	 * when creating the supercell by checking the number of unit cell vectors
	 * necessary to reach every corner of the supercell box.
	 */
	// matrixProduct(a,1,3,Mminv,3,3,b);
	MatrixProduct(Mminv, 3, 3, (float_tt **) &a[0], 3, 1, (float_tt **) &b[0]);
	// showMatrix(Mm,3,3,"M");
	// showMatrix(Mminv,3,3,"M");
	unsigned nxmax, nymax, nzmax;
	unsigned nxmin = nxmax = (int) floor(b[0] - dx);
	unsigned nymin = nymax = (int) floor(b[1] - dy);
	unsigned nzmin = nzmax = (int) floor(b[2] - dz);
	for (ix = 0; ix <= 1; ix++)
		for (iy = 0; iy <= 1; iy++)
			for (iz = 0; iz <= 1; iz++) {
				a[0] = ix * m_cubex - dx;
				a[1] = iy * m_cubey - dy;
				a[2] = iz * m_cubez - dz;

				// matrixProduct(a,1,3,Mminv,3,3,b);
				MatrixProduct(Mminv, 3, 3, (float_tt **) &a[0], 3, 1,
						(float_tt **) &b[0]);

				if (nxmin > (int) floor(b[0]))
					nxmin = (int) floor(b[0]);
				if (nxmax < (int) ceil(b[0]))
					nxmax = (int) ceil(b[0]);
				if (nymin > (int) floor(b[1]))
					nymin = (int) floor(b[1]);
				if (nymax < (int) ceil(b[1]))
					nymax = (int) ceil(b[1]);
				if (nzmin > (int) floor(b[2]))
					nzmin = (int) floor(b[2]);
				if (nzmax < (int) ceil(b[2]))
					nzmax = (int) ceil(b[2]);
			}

	atomSize =	(1 + (nxmax - nxmin) * (nymax - nymin) * (nzmax - nzmin) * ncoord);
	if (atomSize != oldAtomSize) {
		m_atoms.resize(atomSize);
		oldAtomSize = atomSize;
	}
	// showMatrix(Mm,3,3,"Mm");
	// showMatrix(Mminv,3,3,"Mminv");
	// printf("Range: (%d..%d, %d..%d, %d..%d)\n",
	// nxmin,nxmax,nymin,nymax,nzmin,nzmax);

	atomCount = 0;
	jVac = 0;
	for (iatom = 0; iatom < ncoord;) {
		// printf("%d: (%g %g %g) %d\n",iatom,unitAtoms[iatom].x,unitAtoms[iatom].y,
		//   unitAtoms[iatom].z,unitAtoms[iatom].Znum);
		memcpy(&newAtom, &m_baseAtoms[iatom], sizeof(atom));
		//for (jz=0;jz<m_Znums.size();jz++)	if (m_Znums[jz] == newAtom.Znum) break;
		// allocate more memory, if there is a new element
		/*
		 if (jz == atomKinds) {
		 atomKinds++;
		 if (atomKinds > m_atomKinds) {
		 m_Znums = (int *)realloc(m_Znums,atomKinds*sizeof(int));
		 m_atomKinds = atomKinds;
		 // printf("%d kinds (%d)\n",atomKinds,atoms[i].Znum);
		 }
		 m_Znums[jz] = newAtom.Znum;
		 }
		 */
		/////////////////////////////////////////////////////
		// look for atoms at equal position
		if ((handleVacancies) && (newAtom.Znum > 0)) {
			totOcc = newAtom.occ;
			for (jequal = iatom + 1; jequal < ncoord; jequal++) {
				// if there is anothe ratom that comes close to within 0.1*sqrt(3) A we will increase
				// the total occupany and the counter jequal.
				if ((fabs(newAtom.x - m_baseAtoms[jequal].x) < 1e-6)
						&& (fabs(newAtom.y - m_baseAtoms[jequal].y) < 1e-6)
						&& (fabs(newAtom.z - m_baseAtoms[jequal].z) < 1e-6)) {
					totOcc += m_baseAtoms[jequal].occ;
				} else
					break;
			} // jequal-loop
		} else {
			jequal = iatom + 1;
			totOcc = 1;
		}

		unsigned atomCount = 0;

		// printf("%d: %d\n",atomCount,jz);
		for (ix = nxmin; ix <= nxmax; ix++) {
			for (iy = nymin; iy <= nymax; iy++) {
				for (iz = nzmin; iz <= nzmax; iz++) {
					// atom position in cubic reduced coordinates:
					aOrig[0] = ix + newAtom.x;
					aOrig[1] = iy + newAtom.y;
					aOrig[2] = iz + newAtom.z;

					// Now is the time to remove atoms that are on the same position or could be vacancies:
					// if we encountered atoms in the same position, or the occupancy of the current atom is not 1, then
					// do something about it:
					// All we need to decide is whether to include the atom at all (if totOcc < 1
					// of which of the atoms at equal positions to include
					jChoice = iatom;  // This will be the atom we wil use.
					if ((totOcc < 1) || (jequal > iatom + 1)) { // found atoms at equal positions or an occupancy less than 1!
						// ran1 returns a uniform random deviate between 0.0 and 1.0 exclusive of the endpoint values.
						//
						// if the total occupancy is less than 1 -> make sure we keep this
						// if the total occupancy is greater than 1 (unphysical) -> rescale all partial occupancies!
						if (totOcc < 1.0)
							choice = ran1();
						else
							choice = totOcc * ran1();
						// printf("Choice: %g %g %d, %d %d\n",totOcc,choice,j,i,jequal);
						lastOcc = 0;
						for (i2 = iatom; i2 < jequal; i2++) {
							// atoms[atomCount].Znum = unitAtoms[i2].Znum;
							// if choice does not match the current atom:
							// choice will never be 0 or 1(*totOcc)
							if ((choice < lastOcc)
									|| (choice >= lastOcc + m_baseAtoms[i2].occ)) {
								// printf("Removing atom %d, Z=%d\n",jCell+i2,atoms[jCell+i2].Znum);
								// atoms[atomCount].Znum =  0;  // vacancy
								jVac++;
							} else {
								jChoice = i2;
							}
							lastOcc += m_baseAtoms[i2].occ;
						}
					}
					// here we select the index of our chosen atom
					//if (jChoice != iatom) {
					//std::vector<unsigned>::iterator item = std::find(m_Znums.begin(), m_Znums.end(), m_baseAtoms[jChoice].Znum);
					//if (item!=m_Znums.end())
					//jz=item-m_Znums.begin();
					//for (jz=0;jz<m_Znums.size();jz++)	if (m_Znums[jz] == m_baseAtoms[jChoice].Znum) break;
					//}

					// here we need to call phononDisplacement:
					// phononDisplacement(u,muls,iatom,ix,iy,iz,atomCount,atoms[i].dw,*natom,atoms[i].Znum);
					if (m_Einstein) {
						if (_config->Model.UseTDS) {
							// 20131218 MCS
							// final "false" is whether to print report in PhononDisplacement.  I think it's vestigial code...
							// This function creates the vector offset u
							PhononDisplacement(u, jChoice, ix, iy, iz,
									m_baseAtoms[jChoice], false);
							// here we add the vector u to the original position
							a[0] = aOrig[0] + u[0];
							a[1] = aOrig[1] + u[1];
							a[2] = aOrig[2] + u[2];
						} else {
							a[0] = aOrig[0];
							a[1] = aOrig[1];
							a[2] = aOrig[2];
						}
					} else {
						BOOST_LOG_TRIVIAL(fatal)<<format("Cannot handle phonon-distribution mode for boxed sample yet. Exiting...");
						exit(0);
					}
					// matrixProduct(aOrig,1,3,Mm,3,3,b);
					MatrixProduct(Mm, 3, 3, (float_tt **) &aOrig[0], 3, 1,
							(float_tt **) &b[0]);

					// if (atomCount < 2) {showMatrix(a,1,3,"a");showMatrix(b,1,3,"b");}
					// b now contains atom positions in cartesian coordinates */
					x = b[0] + dx;
					y = b[1] + dy;
					z = b[2] + dz;

					// include atoms that are within the box
					if ((x >= 0) && (x <= m_cubex) && (y >= 0) && (y <= m_cubey)
							&& (z >= 0) && (z <= m_cubez)) {
						// matrixProduct(a,1,3,Mm,3,3,b);
						MatrixProduct(Mm, 3, 3, (float_tt **) &a[0], 3, 1,
								(float_tt **) &b[0]);
						m_atoms[atomCount].x = b[0] + dx;
						m_atoms[atomCount].y = b[1] + dy;
						m_atoms[atomCount].z = b[2] + dz;
						m_atoms[atomCount].dw = m_baseAtoms[jChoice].dw;
						m_atoms[atomCount].occ = m_baseAtoms[jChoice].occ;
						m_atoms[atomCount].q = m_baseAtoms[jChoice].q;
						m_atoms[atomCount].Znum = m_baseAtoms[jChoice].Znum;


						atomCount++;
						/*
						 if (m_baseAtoms[jChoice].Znum > 22)
						 printf("Atomcount: %d, Z = %d\n",atomCount,m_baseAtoms[jChoice].Znum);
						 */
					}
				} /* iz ... */
			} /* iy ... */
		} /* ix ... */
		iatom = jequal;
	} /* iatom ... */
	BOOST_LOG_TRIVIAL(trace)<<format("Removed %d atoms because of multiple occupancy or occupancy < 1") % jVac;
	m_ax = m_cubex;
	m_by = m_cubey;
	m_cz = m_cubez;
	// call phononDisplacement again to update displacement data:
	PhononDisplacement(u, iatom, ix, iy, iz, newAtom, false);
}  // end of 'tiltBoxed(...)'

////////////////////////////////////////////////////////////////////////
// replicateUnitCell
// 
// Replicates the unit cell NcellX x NCellY x NCellZ times
// applies phonon displacement and removes vacancies and atoms appearing
// on same position:
// ncoord is the number of atom positions that has already been read.
// memory for the whole atom-array of size natom has already been allocated
// but the sites beyond natom are still empty.
void CCrystal::ReplicateUnitCell(int handleVacancies) {
	int i, j, i2, jChoice, ncx, ncy, ncz, icx, icy, icz;
	int jequal; // Number of atoms that share a position (mixed position if > 1)
	int jVac;  // Number of vacancies total in replicated supercell
	int jCell; // Offset (in number of coordinates) from origin cell
	int atomKinds = 0;
	double totalOccupancy;
	double choice, lastOcc;

	ncx = _config->Structure.nCellX;
	ncy = _config->Structure.nCellY;
	ncz = _config->Structure.nCellZ;
	std::vector<float_tt> u(3, 0);
	//u = (double *)malloc(3*sizeof(double));

	m_atoms.resize(ncx * ncy * ncz * m_baseAtoms.size());

	//////////////////////////////////////////////////////////////////////////////
	// Look for atoms which share the same position:
	jVac = 0;  // no atoms have been removed yet
	for (i = m_baseAtoms.size() - 1; i >= 0;) {

		////////////////
		if ((handleVacancies) && (m_atoms[i].Znum > 0)) {
			totalOccupancy = m_atoms[i].occ;
			for (jequal = i - 1; jequal >= 0; jequal--) {
				// if there is anothe ratom that comes close to within 0.1*sqrt(3) A we will increase
				// the total occupany and the counter jequal.
				if ((fabs(m_atoms[i].x - m_atoms[jequal].x) < 1e-6)
						&& (fabs(m_atoms[i].y - m_atoms[jequal].y) < 1e-6)
						&& (fabs(m_atoms[i].z - m_atoms[jequal].z) < 1e-6)) {
					totalOccupancy += m_atoms[jequal].occ;
				} else
					break;
			} // jequal-loop
		} else {
			jequal = i - 1;
			totalOccupancy = 1;
		}

		////////////////
		/* replicate unit cell ncx,y,z times: */
		/* We have to start with the last atoms first, because once we added the displacements
		 * to the original unit cell (icx=icy=icz=0), we cannot use those positions
		 * as unit cell coordinates for the other atoms anymore
		 */
		BOOST_LOG_TRIVIAL(trace) << "in RepliceteUnitCell";
		for (icx = ncx - 1; icx >= 0; icx--) {
			for (icy = ncy - 1; icy >= 0; icy--) {
				for (icz = ncz - 1; icz >= 0; icz--) {
					jCell = (icz + icy * ncz + icx * ncy * ncz)
																																			* m_baseAtoms.size();
					j = jCell + i;
					/* We will also add the phonon displacement to the atomic positions now: */
					m_atoms[j].dw = m_baseAtoms[i].dw;
					m_atoms[j].occ = m_baseAtoms[i].occ;
					m_atoms[j].q = m_baseAtoms[i].q;
					m_atoms[j].Znum = m_baseAtoms[i].Znum;

					// Now is the time to remove atoms that are on the same position or could be vacancies:
					// if we encountered atoms in the same position, or the occupancy of the current atom is not 1, then
					// do something about it:
					jChoice = i;
					if ((totalOccupancy < 1) || (jequal < i - 1)) { // found atoms at equal positions or an occupancy less than 1!
						// ran1
						//
						// if the total occupancy is less than 1 -> make sure we keep this
						// if the total occupancy is greater than 1 (unphysical) -> rescale all partial occupancies!
						if (totalOccupancy < 1.0)
							choice = ran1();
						else
							choice = totalOccupancy * ran1();
						// printf("Choice: %g %g %d, %d %d\n",totOcc,choice,j,i,jequal);
						lastOcc = 0;
						for (i2 = i; i2 > jequal; i2--) {
							m_atoms[jCell + i2].dw = m_baseAtoms[i2].dw;
							m_atoms[jCell + i2].occ = m_baseAtoms[i2].occ;
							m_atoms[jCell + i2].q = m_baseAtoms[i2].q;
							m_atoms[jCell + i2].Znum = m_baseAtoms[i2].Znum;

							// if choice does not match the current atom:
							// choice will never be 0 or 1(*totOcc)
							if ((choice < lastOcc)
									|| (choice >= lastOcc + m_baseAtoms[i2].occ)) {
								// printf("Removing atom %d, Z=%d\n",jCell+i2,atoms[jCell+i2].Znum);
								m_atoms[jCell + i2].Znum = 0;  // vacancy
								jVac++;
							} else {
								jChoice = i2;
							}
							lastOcc += m_baseAtoms[i2].occ;
						}
					}
					// this function does nothing, if _config->Model.UseTDS == 0
					PhononDisplacement(u, jChoice, icx, icy, icz, m_atoms[jChoice], false);

					for (i2 = i; i2 > jequal; i2--) {
						m_atoms[jCell + i2].x = m_baseAtoms[i2].x + icx + u[0];
						m_atoms[jCell + i2].y = m_baseAtoms[i2].y + icy + u[1];
						m_atoms[jCell + i2].z = m_baseAtoms[i2].z + icz + u[2];
						BOOST_LOG_TRIVIAL(trace) << format("atom %d: (%3.3f, %3.3f, %3.3f)\n")
																		% (jCell + i2) % m_atoms[jCell + i2].x % m_atoms[jCell + i2].y % m_atoms[jCell + i2].z;
					}
				}  // for (icz=ncz-1;icz>=0;icz--)
			} // for (icy=ncy-1;icy>=0;icy--)
		} // for (icx=ncx-1;icx>=0;icx--)
		i = jequal;
	} // for (i=ncoord-1;i>=0;)
	if ((jVac > 0) )
		BOOST_LOG_TRIVIAL(trace)<<format("Removed %d atoms because of occupancies < 1 or multiple atoms in the same place") %jVac;
}

void CCrystal::DisplaceAtoms() {
	std::vector<atom>::iterator at = m_atoms.begin(), end = m_atoms.end();
	std::vector<float_tt> u(3, 0);

	for (at; at != end; ++at) {
		EinsteinDisplacement(u, (*at));
		// Add obtained u to atomic position
		(*at).x += u[0];
		(*at).y += u[1];
		(*at).z += u[2];
	}
}

void CCrystal::EinsteinDisplacement(std::vector<float_tt>&u, atom &_atom) {
	/* convert the Debye-Waller factor to sqrt(<u^2>) */
	float_tt wobble = m_wobble_temp_scale * sqrt(_atom.dw * k_wobScale);
	u[0] = (wobble * k_sq3 * gasdev());
	u[1] = (wobble * k_sq3 * gasdev());
	u[2] = (wobble * k_sq3 * gasdev());
	///////////////////////////////////////////////////////////////////////
	// Book keeping:
	m_u2[_atom.Znum] += u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
	//ux += u[0]; uy += u[1]; uz += u[2];
	m_u2Count[_atom.Znum]++;

	/* Finally we must convert the displacement for this atom back into its fractional
	 * coordinates so that we can add it to the current position in vector a
	 */
	std::vector<float_tt> uf(3, 0);
	MatrixProduct((float_tt **) &u[0], 1, 3, m_MmInv, 3, 3,
			(float_tt **) &uf[0]);
	memcpy(&u[0], &uf[0], 3 * sizeof(float_tt));
}

/*******************************************************************************
 * int phononDisplacement:
 * This function will calculate the phonon displacement for a given atom i of the
 * unit cell, which has been replicated to the larger cell (icx,icy,icz)
 * The phonon displacement is either defined by the phonon data file, or,
 * the Einstein model, if the appropriate flags in the muls struct are set
 * The displacement will be given in fractional coordinates of a single unit cell.
 *
 * Input parameters:
 * Einstein-mode:
 * need only: dw, Znum, atomCount
 * atomCount: give statistics report, if 0, important only for non-Einstein mode
 * maxAtom: total number of atoms (will be called first, i.e. atomCount=maxAtoms-1:-1:0)
 *
 * Phonon-file mode:
 * ...
 *
 ********************************************************************************/
//  phononDisplacement(u,muls,jChoice,icx,icy,icz,j,atoms[jChoice].dw,*natom,jz);
//  j == atomCount
void CCrystal::PhononDisplacement(std::vector<float_tt> &u, int id, int icx,
		int icy, int icz, atom &atom, bool printReport) {
	int ix, iy, idd; // iz;
	static FILE *fpPhonon = NULL;
	static int Nk, Ns;  // number of k-vectors and atoms per primitive unit cell
	static float_tt *massPrim;   // masses for every atom in primitive basis
	static float_tt **omega;     // array of eigenvalues for every k-vector
	static complex_tt ***eigVecs;  // array of eigenvectors for every k-vector
	static float_tt **kVecs;     // array for Nk 3-dim k-vectors
	static float_tt **q1 = NULL, **q2 = NULL;
	int ik, lambda, icoord; // Ncells, nkomega;
	double kR, kRi, kRr, wobble;
	static float_tt *u2T, ux = 0, uy = 0, uz = 0; // u2Collect=0; // Ttotal=0;
	std::map<unsigned, float_tt> u2;
	std::map<unsigned, unsigned> u2Count;
	// static double uxCollect=0,uyCollect=0,uzCollect=0;
	static int *u2CountT, runCount = 1, u2Size = -1;
	static long iseed = 0;
	static float_tt **Mm = NULL, **MmInv = NULL;
	// static float_tt **MmOrig=NULL,**MmOrigInv=NULL;
	static float_tt *axCell, *byCell, *czCell;
	static float_tt wobScale = 0, sq3, scale = 0;
	std::vector<float_tt> uf(3), b(3);

	float_tt dw = atom.dw;

	if (_config->Model.UseTDS == 0)
		return;

	if (Mm == NULL) {
		// MmOrig = float2D(3,3,"MmOrig");
		// MmOrigInv = float2D(3,3,"MmOrigInv");
		Mm = float2D(3, 3, "Mm");
		MmInv = float2D(3, 3, "Mminv");
		//uf = float1D(3,"uf");
		//b = float1D(3,"uf");

		axCell = Mm[0];
		byCell = Mm[1];
		czCell = Mm[2];
		// memcpy(Mm[0],m_Mm[0],3*3*sizeof(double));
		// We need to copy the transpose of m_Mm to Mm.
		// we therefore cannot use the following command:
		// memcpy(Mm[0],m_Mm[0],3*3*sizeof(double));
		// or Mm = m_Mm;
		for (ix = 0; ix < 3; ix++)
			for (iy = 0; iy < 3; iy++)
				Mm[ix][iy] = m_Mm[iy][ix];

		// makeCellVectMuls(muls, axCell, byCell, czCell);
		// memcpy(MmOrig[0],Mm[0],3*3*sizeof(double));
		Inverse_3x3(MmInv[0], Mm[0]);
	}

	/***************************************************************************
	 * Thermal Diffuse Scattering according to accurate phonon-dispersion
	 *
	 * Information in the phonon file will be stored in binary form as follows:
	 * Nk (number of k-points: 32-bit integer)
	 * Ns (number of atomic species 32-bit integer)
	 * M_1 M_2 ... M_Ns  (32-bit floats)
	 * kx(1) ky(1) kz(1) (32-bit floats)
	 * w_1(1) q_11 q_21 ... q_(3*Ns)1    (32-bit floats)
	 * w_2(1) q_12 q_22 ... q_(3*Ns)2
	 * :
	 * w_(3*Ns)(1) q_1Ns q_2Ns ... q_(3*Ns)Ns
	 * kx(2) ky(2) kz(2)
	 * :
	 * kx(Nk) ky(Nk) kz(Nk)
	 * :
	 *
	 *
	 * Note: only k-vectors in half of the Brillouin zone must be given, since
	 * w(k) = w(-k)
	 * also: 2D arrays will be read slowly varying index = first index (i*m+j)
	 **************************************************************************/

	if (wobScale == 0) {
		wobScale = 1.0 / (8 * M_PI * M_PI);
		sq3 = 1.0 / sqrt(3.0); /* sq3 is an additional needed factor which stems from
		 * int_-infty^infty exp(-x^2/2) x^2 dx = sqrt(pi)
		 * introduced in order to match the wobble factor with <u^2>
		 */
		scale = (float) sqrt(_config->Structure.temperatureK / 300.0);
	}

	if (fpPhonon == NULL) {
		if ((fpPhonon = fopen(m_phononFile.string().c_str(), "r")) != NULL) {
			if (2 * sizeof(float) != sizeof(fftwf_complex)) {
				printf(
						"phononDisplacement: data type mismatch: fftw_complex != 2*float!\n");
				exit(0);
			}
			fread(&Nk, sizeof(int), 1, fpPhonon);
			fread(&Ns, sizeof(int), 1, fpPhonon);
			massPrim = (float_tt *) malloc(Ns * sizeof(float)); // masses for every atom in primitive basis
			fread(massPrim, sizeof(float), Ns, fpPhonon);
			kVecs = float2D(Nk, 3, "kVecs");
			omega = float2D(Nk, 3 * Ns, "omega"); /* array of eigenvalues for every k-vector
			 * omega is given in THz, but the 2pi-factor
			 * is still there, i.e. f=omega/2pi
			 */
			eigVecs = complex3D(Nk, 3 * Ns, 3 * Ns, "eigVecs"); // array of eigenvectors for every k-vector
			for (ix = 0; ix < Nk; ix++) {
				fread(kVecs[ix], sizeof(float), 3, fpPhonon);  // k-vector
				for (iy = 0; iy < 3 * Ns; iy++) {
					fread(omega[ix] + iy, sizeof(float), 1, fpPhonon);
					fread(eigVecs[ix][iy], 2 * sizeof(float), 3 * Ns, fpPhonon);
				}
			}

			/* convert omega into q scaling factors, since we need those, instead of true omega:    */
			/* The 1/sqrt(2) term is from the dimensionality ((q1,q2) -> d=2)of the random numbers */
			for (ix = 0; ix < Nk; ix++) {
				for (idd = 0; idd < Ns; idd++)
					for (iy = 0; iy < 3; iy++) {
						// quantize the energy distribution:
						// tanh and exp give different results will therefore use exp
						// nkomega = (int)(1.0/tanh(THZ_HBAR_2KB*omega[ix][iy+3*id]/_config->Structure.temperatureK));
						// wobble  =      (1.0/tanh(THZ_HBAR_2KB*omega[ix][iy+3*id]/_config->Structure.temperatureK)-0.5);
						// nkomega = (int)(1.0/(exp(THZ_HBAR_KB*omega[ix][iy+3*id]/_config->Structure.temperatureK)-1)+0.5);
						if (omega[ix][iy + 3 * idd] > 1e-4) {
							wobble =
									_config->Structure.temperatureK > 0 ?
											(1.0
													/ (exp(
															THZ_HBAR_KB
															* omega[ix][iy
															            + 3
															            * idd]
															            / _config->Structure.temperatureK)
															- 1)) :
															0;
							// if (ix == 0) printf("%g: %d %g\n",omega[ix][iy+3*id],nkomega,wobble);
							wobble = sqrt(
									(wobble + 0.5)
									/ (2 * M_PI * Nk * 2 * massPrim[idd]
									                                * omega[ix][iy + 3 * idd]
									                                            * THZ_AMU_HBAR));
						} else
							wobble = 0;
						/* Ttotal += 0.25*massPrim[id]*((wobble*wobble)/(2*Ns))*
						 omega[ix][iy+3*id]*omega[ix][iy+3*id]*AMU_THZ2_A2_KB;
						 */
						omega[ix][iy + 3 * idd] = wobble;
					}  // idd
				// if (ix == 0) printf("\n");
			}
			// printf("Temperature: %g K\n",Ttotal);
			// printf("%d %d %d\n",(int)(0.4*(double)Nk/11.0),(int)(0.6*(double)Nk),Nk);
			q1 = float2D(3 * Ns, Nk, "q1");
			q2 = float2D(3 * Ns, Nk, "q2");

		}
		fclose(fpPhonon);
	}  // end of if phononfile

	//
	// in the previous bracket: the phonon file is only read once.
	/////////////////////////////////////////////////////////////////////////////////////
	if (Nk > 800)
		printf("Will create phonon displacements for %d k-vectors - please wait ...\n",Nk);
	for (lambda = 0; lambda < 3 * Ns; lambda++)
		for (ik = 0; ik < Nk; ik++) {
			q1[lambda][ik] = (omega[ik][lambda] * gasdev());
			q2[lambda][ik] = (omega[ik][lambda] * gasdev());
		}
	// printf("Q: %g %g %g\n",q1[0][0],q1[5][8],q1[0][3]);
	// id seems to be the index of the correct atom, i.e. ranges from 0 .. Natom
	printf("created phonon displacements %d, %d, %d %d (eigVecs: %d %d %d)!\n",
			atom.Znum, Ns, Nk, id, Nk, 3 * Ns, 3 * Ns);
	/* loop over k and lambda:  */
	memset(&u[0], 0, 3 * sizeof(double));
	for (lambda = 0; lambda < 3 * Ns; lambda++)
		for (ik = 0; ik < Nk; ik++) {
			// if (kVecs[ik][2] == 0){
			kR = 2 * M_PI
					* (icx * kVecs[ik][0] + icy * kVecs[ik][1]
					                                        + icz * kVecs[ik][2]);
			//  kR = 2*M_PI*(blat[0][0]*kVecs[ik][0]+blat[0][1]*kVecs[ik][1]+blat[0][2]*kVecs[ik][2]);
			kRr = cos(kR);
			kRi = sin(kR);
			for (icoord = 0; icoord < 3; icoord++) {
				u[icoord] +=
						q1[lambda][ik]
						           * (eigVecs[ik][lambda][icoord + 3 * id].real() * kRr
						        		   - eigVecs[ik][lambda][icoord + 3 * id].imag()
						        		   * kRi)
						        		   - q2[lambda][ik]
						        		                * (eigVecs[ik][lambda][icoord + 3 * id].real()
						        		                		* kRi
						        		                		+ eigVecs[ik][lambda][icoord
						        		                		                      + 3 * id].imag() * kRr);
			}
		}
	// printf("u: %g %g %g\n",u[0],u[1],u[2]);
	/* Convert the cartesian displacements back to reduced coordinates
	 */
	///////////////////////////////////////////////////////////////////////
	// Book keeping:
	u2[atom.Znum] += u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
	ux += u[0];
	uy += u[1];
	uz += u[2];
	u[0] /= m_ax;
	u[1] /= m_by;
	u[2] /= m_cz;
	u2Count[atom.Znum]++;

	return;
}

int CCrystal::AtomCompareZnum(const void *atPtr1, const void *atPtr2) {
	int comp = 0;
	atom *atom1, *atom2;

	atom1 = (atom *) atPtr1;
	atom2 = (atom *) atPtr2;
	/* Use the fact that z is the first element in the atom struct */
	comp = (atom1->Znum == atom2->Znum) ?
			0 : ((atom1->Znum > atom2->Znum) ? -1 : 1);
	return comp;
}

int CCrystal::AtomCompareZYX(const void *atPtr1, const void *atPtr2) {
	int comp = 0;
	atom *atom1, *atom2;

	atom1 = (atom *) atPtr1;
	atom2 = (atom *) atPtr2;
	/* Use the fact that z is the first element in the atom struct */
	comp = (atom1->z == atom2->z) ? 0 : ((atom1->z > atom2->z) ? 1 : -1);
	if (comp == 0) {
		comp = (atom1->y == atom2->y) ? 0 : ((atom1->y > atom2->y) ? 1 : -1);
		if (comp == 0) {
			comp = (atom1->x == atom2->x) ?
					0 : ((atom1->x > atom2->x) ? 1 : -1);
		}
	}
	return comp;
}

void CCrystal::WriteStructure(unsigned run_number) {
	// to_string is C++0x - may not work on older compilers
	m_structureWriter->Write(m_atoms, std::to_string(run_number));
	/*
	 if (m_cfgFile != NULL) {
	 sprintf(buf,"%s/%s",m_folder,m_cfgFile);
	 // append the TDS run number
	 if (strcmp(buf+strlen(buf)-4,".cfg") == 0) *(buf+strlen(buf)-4) = '\0';
	 if (_config->Model.UseTDS) sprintf(buf+strlen(buf),"_%d.cfg",m_avgCount);
	 else sprintf(buf+strlen(buf),".cfg");

	 // printf("Will write CFG file <%s> (%d)\n",buf,_config->Model.UseTDS)
	 writeCFG(atoms,natom,buf,muls);
	 }
	 */
}

void CCrystal::GetCrystalBoundaries(float_tt &min_x, float_tt &max_x, float_tt &min_y, float_tt &max_y, float_tt &min_z, float_tt &max_z) {
	min_x = m_minX;
	max_x = m_maxX;
	min_y = m_minY;
	max_y = m_maxY;
	min_z = m_minZ;
	max_z = m_maxZ;
}

} // end namespace QSTEM

// *******************  Matrix manipulation ***********************
//    For now, use our own internal routines as has been done always.
//    For the future, consider using a linear algebra library instead - Eigen, BLAS/ATLAS/MKL/GOTO

#include "matrixlib.hpp"

namespace QSTEM {

void CCrystal::Inverse_3x3(float_tt *res, const float_tt *a) {
	// use function from matrixlib for now
	return inverse_3x3(res, a);
}

void CCrystal::RotateVect(float_tt *vectIn, float_tt *vectOut, float_tt phi_x,
		float_tt phi_y, float_tt phi_z) {
	return rotateVect(vectIn, vectOut, phi_x, phi_y, phi_z);
}

void CCrystal::MatrixProduct(float_tt **a, int Nxa, int Nya, float_tt **b,
		int Nxb, int Nyb, float_tt **c) {
	return matrixProduct(a, Nxa, Nya, b, Nxb, Nyb, c);
}

void CCrystal::MatrixProduct(float_tt **a, int Nxa, int Nya, float_tt **b,
		int Nxb, int Nyb, float_tt *c) {
	return matrixProduct(a, Nxa, Nya, b, Nxb, Nyb, &c);
}

void CCrystal::RotateMatrix(float_tt *matrixIn, float_tt *matrixOut,
		float_tt phi_x, float_tt phi_y, float_tt phi_z) {
	return rotateMatrix(matrixIn, matrixOut, phi_x, phi_y, phi_z);
}

// ******************  end matrix manipulation

}// end namespace QSTEM
