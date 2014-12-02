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

#include "../stemtypes_fftw3.hpp"
#include "../imagelib_fftw3.hpp"
//#include "config_reader_factory.hpp"
#include "../memory_fftw3.hpp"
#include "../crystal.hpp"
#include "scatfactsRez.hpp"
#include <map>

#include "pot_interface.hpp"


#ifndef POTENTIAL_BASE_H
#define POTENTIAL_BASE_H

#define OVERSAMPLING 2
#define OVERSAMPLINGZ 3*OVERSAMPLING

namespace QSTEM
{

class CPotential : public IPotential
{
public:
	CPotential();
	CPotential(unsigned nx, unsigned ny, unsigned nz, float_tt dx, float_tt dy, float_tt dz, float_tt atomRadius, float_tt v0);
	CPotential(const ConfigPtr c);

	void GetSizePixels(unsigned &nx, unsigned &ny) const;
	void WriteSlice(unsigned idx);
	void WriteProjectedPotential();
	inline ComplexArray2DView GetSlice(unsigned idx){return m_trans1[boost::indices[idx][range(0,_config->Model.ny)][range(0,_config->Model.nx)]];}
	complex_tt GetSlicePixel(unsigned iz, unsigned ix, unsigned iy);

	virtual ~CPotential();
	void AtomBoxLookUp(complex_tt &val, int Znum, float_tt x, float_tt y, float_tt z, float_tt B);
	virtual void CenterAtomZ(std::vector<atom>::iterator &atom, float_tt &z);
	virtual void DisplayParams();
	virtual void MakeSlices(int nlayer,StructurePtr crystal);
	virtual void Refresh();
	virtual void ReadPotential(std::string &fileName, unsigned subSlabIdx);
	virtual void SetStructure(StructurePtr structure);



	void SetSliceThickness(float_tt thickness_Angstroms);
	void SetSliceThickness(std::vector<float_tt> thickness_Angstroms);
	void SetNSlices(unsigned slices);
	atom GetAtom(unsigned idx){return m_crystal->m_atoms[idx];}

	// public members (should be moved to protected, and given getters/setters...
	bool m_potential3D;

protected:
	virtual void SliceSetup();
	void ResizeSlices();
	void MakePhaseGratings();
	void BandlimitTransmissionFunction();
	void ReadSlice(const std::string &fileName, ComplexArray2DView slice, unsigned idx);

	virtual void AddAtomToSlices(std::vector<atom>::iterator &atom, float_tt atomX, float_tt atomY, float_tt atomZ)=0;
	virtual void ComputeAtomPotential(std::vector<atom>::iterator &atom)=0;

	ImageIOPtr m_imageIO;
	StructurePtr m_crystal;
	ConfigPtr _config;
	ComplexArray3D m_trans1;

	float_tt m_dz;   // resolutions
	float_tt m_ddx, m_ddy, m_ddz;   // oversampled resolutions
	float_tt m_dkx,m_dky, m_dkz,m_kmax2;
	float_tt m_c; // the thickness of the current sub-slab (in A)
	float_tt m_dr, m_atomRadius2;
	float_tt m_offsetX, m_offsetY;
	int m_iRadX, m_iRadY, m_iRadZ, m_iRad2;

	int _nx,_ny;
	int m_boxNx, m_boxNy, m_boxNz;
	// ********* multi-slab parameters ********
	int m_cellDiv; // How many sub-slabs the model is divided into
	int m_divCount; // How many sub-slabs we've already processed

	int m_scatFactor;  // TODO: NOT USED YET The scattering factor type.  One of: 0 (Doyle-Turner); 1 (Wieck-Kohl); 2 (Custom)

	std::map<int, atomBoxPtr> m_atomBoxes;
	std::vector<atom> *m_atoms;
	std::vector<float_tt> m_sliceThicknesses;
	std::vector<float_tt> m_slicePos;

	std::string m_fileBase; // base filename for saving potential.  Will have slice index appended.
	boost::filesystem::path m_potFileBase; //base filename for loading potential from files.  Will have slice index appended.


	float_tt sfLUT(float_tt s,int atKind);

	void splinh( float_tt x[], float_tt y[],  std::vector<float_tt>& b,std::vector<float_tt>& c,std::vector<float_tt>& d, int n);
	float_tt seval( float_tt *x, float_tt *y,std::vector<float_tt>& b,std::vector<float_tt>& c,
			std::vector<float_tt>& d, int n, float_tt x0 );
};

}

#endif
