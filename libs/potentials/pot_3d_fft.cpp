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

#include "pot_3d_fft.hpp"

// for memset
#include <cstring>
#include <cvmlcpp/signal/Fourier>
#include <boost/format.hpp>
using boost::format;
namespace QSTEM {

C3DFFTPotential::C3DFFTPotential() :C3DPotential() {}

C3DFFTPotential::C3DFFTPotential(const ConfigPtr configReader) :	C3DPotential(configReader){
	for (unsigned i = 0; i < _config->Model.nSlices; i++) {
		if (m_sliceThicknesses[0] != m_sliceThicknesses[i]) {
			BOOST_LOG_TRIVIAL(warning) << format("Warning: slice thickness not constant, will give wrong results (iz=%d)!") % i;
		}
	}
}

void C3DFFTPotential::DisplayParams() {
	CPotential::DisplayParams();
	BOOST_LOG_TRIVIAL(info) << "* Potential calculation: 3D (FFT method)";
}

void C3DFFTPotential::AddAtomNonPeriodic(std::vector<atom>::iterator &atom,	float_tt atomBoxX, unsigned int iAtomX, float_tt atomBoxY, unsigned int iAtomY, float_tt atomZ) {

	unsigned iAtomZ = (int) floor(atomZ / _config->Model.sliceThicknessAngstrom + 0.5);
	unsigned nzSub, nRadius, Nz_lut;
	// FIXME array indexing
	ComplexArray2D atPotPtr = m_atPot[atom->Znum];

	nRadius = 2 * OVERSAMPLING * (int) ceil(_config->Potential.AtomRadiusAngstrom / _config->Model.dx) /2;
	nzSub = (int) floor(OVERSAMPLING * _config->Model.sliceThicknessAngstrom / _config->Model.dx);
	if (2.0 * (nzSub >> 1) == nzSub)	nzSub += 1;
	Nz_lut = (2 * (int) ceil(_config->Potential.AtomRadiusAngstrom / _config->Model.sliceThicknessAngstrom)) * nzSub / 2;

	int xstart = (int) iAtomX - (int) m_iRadX < 0 ? 0 : iAtomX - m_iRadX;
	int xend =(int) iAtomX + (int) m_iRadX >= _config->Model.nx ? _config->Model.nx - 1 : iAtomX + m_iRadX;
	int ystart = (int) iAtomY - (int) m_iRadY < 0 ? 0 : iAtomY - m_iRadY;
	int yend =	(int) iAtomY + (int) m_iRadY >= _config->Model.ny ? _config->Model.ny - 1 : iAtomY + m_iRadY;
	int zstart =  (int)iAtomZ -  (int)m_iRadZ < 0 ? -iAtomZ : -m_iRadZ;
	int zend =	 (int)iAtomZ +  (int)m_iRadZ >= _config->Model.nSlices ?	_config->Model.nSlices - iAtomZ - 1 : m_iRadZ;

	// if within the potential map range:
	if ((xstart < _config->Model.nx) && (xend >= 0) && (ystart < _config->Model.ny) && (yend >= 0)) {
		// define range of sampling from atomZ-/+atomRadius
		if (((int) iAtomZ + (int) zstart < (int) _config->Model.nSlices)	&& ((int) iAtomZ + (int) zend >= 0)) {
			// retrieve the pointer for this atom

#if USE_Q_POT_OFFSETS
			// retrieve the pointer to the array of charge-dependent potential offset
			// This function will return NULL; if the charge of this atom is zero:
			ComplexVector atPotOffsPtr;
			GetAtomPotentialOffset3D(atom->Znum,muls,m_tds ? 0 : atoms[iatom].dw,nzSub,nRadius,Nz_lut,atom->q, atPotOffsPtr);
#endif // USE_Q_POT_OFFSETS

			int iOffsetUpperLimit = nRadius * (Nz_lut - 1);
			int iOffsetLowerLimit = -nRadius * (Nz_lut - 1);
			int iOffsStep = nzSub * nRadius;

			// Slices around the slice that this atom is located in must be affected by this atom:
			// iaz must be relative to the first slice of the atom potential box.
			float_tt potVal = 0;
			int iay = 0;
			int iaz = 0;
			for (int iax = xstart; iax <= xend; iax++) {
				//////////////////////////////////////////////////////////////////////
				// Computation of Radius must be made faster by using pre-calculated ddx
				// and LUT for sqrt:
				// use of sqrt slows down from 120sec to 180 sec.
				float_tt x2 = iax * _config->Model.dx - atomBoxX;
				x2 *= x2;
				for (iay = ystart; iay <= yend; iay++) {
					// printf("iax=%d, iay=%d\n",iax,iay);
					float_tt y2 = iay * _config->Model.dy - atomBoxY;
					y2 *= y2;
					float_tt r = sqrt(x2 + y2);
					// r = (x2+y2);
					float_tt ddr = r / m_dr;
					int ir = (int) floor(ddr);
					// add in different slices, once r has been defined
					if (ir < nRadius - 1) {
						ddr = ddr - (double) ir;
						//						complex_tt *ptr = potPtr;

						float_tt dOffsZ = (iAtomZ + zstart
								- atomZ / _config->Model.sliceThicknessAngstrom) * nzSub;
#if Z_INTERPOLATION
						unsigned iOffsetZ = (unsigned)dOffsZ;
						ddz = fabs(dOffsZ - (float_tt)iOffsetZ);
#else
						unsigned iOffsetZ = (int) (dOffsZ + 0.5);
#endif
						iOffsetZ *= nRadius;

						for (iaz = zstart; iaz <= zend; iaz++) {
							potVal = 0;
							if (iOffsetZ < 0) {
								if (iOffsetZ > iOffsetLowerLimit) {
									// do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
									potVal = (1-ddz)*((1-ddr)*atPotPtr[ir-iOffsetZ+nRadius].real()+ddr*atPotPtr[ir+1-iOffsetZ+nRadius].real())+
											ddz *((1-ddr)*atPotPtr[ir-iOffsetZ ].real()+ddr*atPotPtr[ir+1-iOffsetZ ].real());
#if USE_Q_POT_OFFSETS
									// add the charge-dependent potential offset
									if (atPotOffsPtr != NULL) {
										potVal += atoms[iatom].q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir-iOffsetZ+nRadius][0]+ddr*atPotOffsPtr[ir+1-iOffsetZ+nRadius][0])+
												ddz *((1-ddr)*atPotOffsPtr[ir-iOffsetZ ][0]+ddr*atPotOffsPtr[ir+1-iOffsetZ ][0]));
									}
#endif // USE_Q_POT_OFFSETS
#else // Z_INTERPOLATION
									// FIXME indices
									potVal = (1 - ddr)* atPotPtr[ir - iOffsetZ + nRadius][1000].real()+ ddr* atPotPtr[ir + 1 - iOffsetZ + nRadius][1000].real();
#if USE_Q_POT_OFFSETS
									// add the charge-dependent potential offset
									if (atPotOffsPtr != NULL) {
										potVal += atoms[iatom].q*((1-ddr)*atPotOffsPtr[ir-iOffsetZ+nRadius].real()+ddr*atPotOffsPtr[ir+1-iOffsetZ+nRadius].real());
									}
#endif // USE_Q_POT_OFFSETS
#endif // Z_INTERPOLATION
								}
							} // if iOffZ < 0
							else {
								// select the pointer to the right layer in the lookup box
								if (iOffsetZ < iOffsetUpperLimit) {
									// do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
									potVal = (1-ddz)*((1-ddr)*atPotPtr[ir+iOffsetZ][0]+ddr*atPotPtr[ir+1+iOffsetZ][0])+
											ddz *((1-ddr)*atPotPtr[ir+iOffsetZ+nRadius][0]+ddr*atPotPtr[ir+1+iOffsetZ+nRadius][0]);
#if USE_Q_POT_OFFSETS
									// add the charge-dependent potential offset
									if (atPotOffsPtr != NULL) {
										potVal += atom->q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir+iOffsetZ][0]+ddr*atPotOffsPtr[ir+1+iOffsetZ][0])+
												ddz *((1-ddr)*atPotOffsPtr[ir+iOffsetZ+nRadius][0]+ddr*atPotOffsPtr[ir+1+iOffsetZ+nRadius][0]));
									}
#endif // USE_Q_POT_OFFSETS
#else // Z_INTERPOLATION
									// FIXME indices
									potVal =(1 - ddr)* atPotPtr[ir + iOffsetZ][1000].real()+ ddr* atPotPtr[ir + 1+ iOffsetZ][1000].real();
#if USE_Q_POT_OFFSETS
									// add the charge-dependent potential offset
									if (atPotOffsPtr != NULL) {
										potVal += atom->q*((1-ddr)*atPotOffsPtr[ir+iOffsetZ][0]+ddr*atPotOffsPtr[ir+1+iOffsetZ][0]);
									}
#endif // USE_Q_POT_OFFSETS
#endif // Z_INTERPOLATION

								}
							} // if iOffsZ >=0

#pragma omp critical
							m_trans1[iAtomZ + iaz][iay][iax] += potVal;

							iOffsetZ += iOffsStep;
						} // end of iaz-loop
					} // if ir < Nr
				} // iay=iay0 .. iay1
				//				if(potVal!=0)
				//					printf("_trans1[%-4d + %-4d][%-4d][%-4d] += %f\n",iAtomZ,iaz,iay,iax,potVal);
			} // iax=iax0 .. iax1
		} // iaz0+iAtomZ < m_slices
	} // if within bounds
}

void C3DFFTPotential::AddAtomToSlices(std::vector<atom>::iterator &atom,
		float_tt atomX, float_tt atomY, float_tt atomZ) {
	unsigned iAtomX = (int) floor(atomX / _config->Model.dx);
	unsigned iAtomY = (int) floor(atomY / _config->Model.dy);

	if (_config->Potential.periodicXY) {
		AddAtomPeriodic(atom, atomX, iAtomX, atomY, iAtomY, atomZ);
	} else {
		AddAtomNonPeriodic(atom, atomX, iAtomX, atomY, iAtomY, atomZ);
	}
}

void C3DFFTPotential::AddAtomPeriodic(std::vector<atom>::iterator &atom,float_tt atomBoxX, unsigned int iAtomX, float_tt atomBoxY,unsigned int iAtomY, float_tt atomZ) {
	int iAtomZ = (int) rint(atomZ / _config->Model.sliceThicknessAngstrom)+m_iRadZ;
	int iax0 = iAtomX - m_iRadX;
	int iax1 = iAtomX + m_iRadX;
	int iay0 = iAtomY - m_iRadY;
	int iay1 = iAtomY + m_iRadY;

//	int iaz0 = iAtomZ + atomRadiusZOffset - m_iRadZ < 0 ? -(iAtomZ+atomRadiusZOffset) : -m_iRadZ;
//	int iaz1 =	iAtomZ + m_iRadZ >= _config->Model.nSlices ? _config->Model.nSlices - (iAtomZ+atomRadiusZOffset) - 1 : m_iRadZ;

	// if (iatom < 2) printf("iatomZ: %d, %d cz=%g, %g: %d, %d\n",iAtomZ,iaz0,_config->Model.sliceThicknessAngstrom,atomZ,(int)(-2.5-atomZ),(int)(atomZ+2.5));
	// printf("%d: iatomZ: %d, %d cz=%g, %g\n",iatom,iAtomZ,iaz0,_config->Model.sliceThicknessAngstrom,atomZ);
//	bool isAtomZWithinSlices = (iAtomZ + iaz0 < _config->Model.nSlices) && (iAtomZ + iaz1 >= 0);
//	if (isAtomZWithinSlices) {
		ComplexArray2D pot = m_atPot[atom->Znum];

#if USE_Q_POT_OFFSETS
		// retrieve the pointer to the array of charge-dependent potential offset
		// This function will return NULL; if the charge of this atom is zero:
		ComplexVector atPotOffsPtr;
		BOOST_LOG_TRIVIAL(fatal) << "getAtomPotentialOffset3D not implemented yet";
		getAtomPotentialOffset3D(atoms[iatom].Znum,muls,m_tds ? 0 : atoms[iatom].dw,&_nzPerSlice,&_nx/2,&_nz/2,atoms[iatom].q, atPotOffsPtr);
#endif // USE_Q_POT_OFFSETS

		int iOffsLimHi = _nx/2 * (_nz/2 - 1);
		int iOffsLimLo = -_nx/2 * (_nz/2 - 1);
		int iOffsStep = _nzPerSlice * _nx/2;

		// Slices around the slice that this atom is located in must be affected by this atom:
		// iaz must be relative to the first slice of the atom potential box.
		for (int iax = iax0; iax < iax1; iax++) {
			for (int iay = iay0; iay < iay1; iay++) {
				float_tt x2 = iax * _config->Model.dx - atomBoxX;
				x2 *= x2;
				float_tt y2 = iay * _config->Model.dy - atomBoxY;
				y2 *= y2;
				float_tt ddr = sqrt(x2 + y2) / m_dr;
				int iCurrentRadius = (int) floor(ddr);
				// add in different slices, once r has been defined
				if (iCurrentRadius < _nx/2 - 1) {
					float fractionAboveCurrentRadius = ddr - (double) iCurrentRadius;
					// Include interpolation in z-direction as well (may do it in a very smart way here !):

					float_tt dOffsZ = (iAtomZ -m_iRadZ - atomZ / _config->Model.sliceThicknessAngstrom)* _nzPerSlice;
#if Z_INTERPOLATION
					int iOffsZ = (int)dOffsZ;
					float_tt ddz = fabs(dOffsZ - (double)iOffsZ);
#else // Z_INTERPOLATION
					int iOffsZ = (int) (dOffsZ + 0.5);
#endif // Z_INTERPOLATION
					iOffsZ *= _nx/2;

					for (int iaz = -m_iRadZ; iaz <= m_iRadZ; iaz++) {
						float_tt potVal = 0;
						// iOffsZ = (int)(fabs(iAtomZ+iaz-atomZ/_config->Model.sliceThicknessAngstrom)*nzSub+0.5);
						if (iOffsZ < 0) {
							if (iOffsZ > iOffsLimLo) {
								// TODO fix this to use ComplexArray2D, fix the confusing indices
								//do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
								potVal = (1-ddz)*((1-fractionAboveCurrentRadius)*pot[iCurrentRadius-iOffsZ+_nx/2][0]
								                     +fractionAboveCurrentRadius*pot[iCurrentRadius+1-iOffsZ+_nx/2][0])+
										ddz *((1-fractionAboveCurrentRadius)*pot[iCurrentRadius-iOffsZ][0]+
												fractionAboveCurrentRadius*pot[iCurrentRadius+1-iOffsZ][0]);
#if USE_Q_POT_OFFSETS
								// add the charge-dependent potential offset
								if (atPotOffsPtr != NULL) {
									potVal += atoms[iatom].q*((1-ddz)*((1-fractionAboveCurrentRadius)*atPotOffsPtr[iCurrentRadius-iOffsZ+_nx/2][0]+
											fractionAboveCurrentRadius*atPotOffsPtr[iCurrentRadius+1-iOffsZ+_nx/2][0])+
											ddz *((1-fractionAboveCurrentRadius)*atPotOffsPtr[iCurrentRadius-iOffsZ ][0]+
													fractionAboveCurrentRadius*atPotOffsPtr[iCurrentRadius+1-iOffsZ ][0]));
								}
#endif // USE_Q_POT_OFFSETS
#else // Z_INTERPOLATION
								//FIXME
								potVal =(1 - fractionAboveCurrentRadius)* pot[iCurrentRadius - iOffsZ + _nx/2][1000].real()+
										fractionAboveCurrentRadius* pot[iCurrentRadius + 1- iOffsZ + _nx/2][1000].real();
#if USE_Q_POT_OFFSETS
								// add the charge-dependent potential offset
								if (atPotOffsPtr != NULL) {
									potVal += atoms[iatom].q*((1-fractionAboveCurrentRadius)*atPotOffsPtr[iCurrentRadius-iOffsZ+_nx/2][0]+
											fractionAboveCurrentRadius*atPotOffsPtr[iCurrentRadius+1-iOffsZ+_nx/2][0]);
								}
#endif // USE_Q_POT_OFFSETS
#endif // Z_INTERPOLATION
							}
						} else {
							// select the pointer to the right layer in the lookup box
							// printf("%4d: iOffsZ: %d, iaz: %d (slice: %d, pos: %g [%d .. %d])\n",iatom,iOffsZ,iaz,iAtomZ+iaz,atomZ,iaz0,iaz1);
							if (iOffsZ < iOffsLimHi) {
								// do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
								potVal = (1-ddz)*((1-fractionAboveCurrentRadius)*pot[iCurrentRadius+iOffsZ][0]+
										fractionAboveCurrentRadius*pot[iCurrentRadius+1+iOffsZ][0])+
										ddz *((1-fractionAboveCurrentRadius)*pot[iCurrentRadius+iOffsZ+_nx][0]+
												fractionAboveCurrentRadius*pot[iCurrentRadius+1+iOffsZ+_nx][0]);
#if USE_Q_POT_OFFSETS
								// add the charge-dependent potential offset
								if (atPotOffsPtr != NULL) {
									potVal += atoms[iatom].q*((1-ddz)*((1-fractionAboveCurrentRadius)*atPotOffsPtr[iCurrentRadius+iOffsZ][0]+
											fractionAboveCurrentRadius*atPotOffsPtr[iCurrentRadius+1+iOffsZ][0])+
											ddz *((1-fractionAboveCurrentRadius)*atPotOffsPtr[iCurrentRadius+iOffsZ+_nx/2][0]+
													fractionAboveCurrentRadius*atPotOffsPtr[iCurrentRadius+1+iOffsZ+_nx/2][0]));
								}
#endif // USE_Q_POT_OFFSETS
#else // Z_INTERPOLATION
								//FIXME
								potVal =(1 - fractionAboveCurrentRadius) * pot[iCurrentRadius + iOffsZ][1000].real()+
										fractionAboveCurrentRadius* pot[iCurrentRadius + 1+ iOffsZ][1000].real();
#if USE_Q_POT_OFFSETS
								// add the charge-dependent potential offset
								if (atPotOffsPtr != NULL) {
									potVal += atoms[iatom].q*((1-fractionAboveCurrentRadius)*atPotOffsPtr[iCurrentRadius+iOffsZ][0]+
											fractionAboveCurrentRadius*atPotOffsPtr[iCurrentRadius+1+iOffsZ][0]);
								}
#endif // USE_Q_POT_OFFSETS
#endif // Z_INTERPOLATION
							}
						}
						int ix = (iAtomX + iax+_config->Model.nx) % _config->Model.nx;
						int iy = (iAtomY + iay+_config->Model.ny) % _config->Model.ny;
						m_trans1[iAtomZ + iaz ][iy][ix]+= potVal;

						BOOST_LOG_TRIVIAL(trace) << boost::format("_trans1[%d ][%d][%d] += %f\n")
							% (iAtomZ+iaz)%((iAtomY+iay+_config->Model.ny) % _config->Model.ny)% ((iAtomX+iax+_config->Model.nx) % _config->Model.nx)%potVal;
						iOffsZ += iOffsStep;
					} // for iaz
				} // if ir < Nr-1

				// advance pointer to next complex potential point:
				//				potPtr += 2;
				// make imaginary part zero for now
				// *potPtr = 0;
				// potPtr++;

				// wrap around when end of y-line is reached:
				//				if (++iay % m_ny == 0)
				//					potPtr -= 2 * m_ny;
				//					iay = 0;
			}
		}
//	} // iaz0+iAtomZ < m_slices
}

void C3DFFTPotential::SliceSetup(){
	C3DPotential::SliceSetup();
	float_tt kmax2;

	if (m_atPot.size() == 0) {
		_nx = 2 * OVERSAMPLING * (int) ceil(_config->Potential.AtomRadiusAngstrom / _config->Model.dx);
		_ny = 2 * OVERSAMPLING * (int) ceil(_config->Potential.AtomRadiusAngstrom / _config->Model.dy);

		// The FFT-resolution in the z-direction must be high enough to avoid
		// artifacts due to premature cutoff of the rec. space scattering factor
		// we will therefore make it roughly the same as the x-resolution
		// However, we will make sure that a single slice contains an integer number
		// of sampling points.
		_nzPerSlice = (int) floor(OVERSAMPLING * _config->Model.sliceThicknessAngstrom / _config->Model.dx);

		// make nzPerSlice odd:
		if (2.0 * (_nzPerSlice >> 1) == _nzPerSlice)	_nzPerSlice += 1;

		// Total number of z-positions should be twice that of atomRadius/sliceThickness
		_nz = (2 * (int) ceil(_config->Potential.AtomRadiusAngstrom / _config->Model.sliceThicknessAngstrom)) * _nzPerSlice;

		m_dkx = 0.5 * OVERSAMPLING / (_nx * _config->Model.dx); // nx*m_dx is roughly 2*_config->Potential.AtomRadiusAngstrom
		m_dkz = _nzPerSlice / (_nz * _config->Model.sliceThicknessAngstrom);
		kmax2 = 0.5 * _nx * m_dkx / OVERSAMPLING;
		// Don't square kmax2 yet!

		scatPar[0][N_SF - 1] = 1.2 * kmax2;
		scatPar[0][N_SF - 2] = 1.1 * kmax2;
		scatPar[0][N_SF - 3] = kmax2;
		// adjust the resolution of the lookup table if necessary
		int ix = 0;
		if (scatPar[0][N_SF - 4] > scatPar[0][N_SF - 3]) {

			if (1) {// TODO: why if(1)? dont use all factors?
				// set additional scattering parameters to zero:
				for (; ix < N_SF - 10; ix++) {
					if (scatPar[0][N_SF - 4 - ix]< scatPar[0][N_SF - 3] - 0.001 * (ix + 1))
						break;
					scatPar[0][N_SF - 4 - ix] = scatPar[0][N_SF - 3] - 0.001 * (ix + 1);
					for (unsigned iy = 1; iy < N_ELEM; iy++) scatPar[iy][N_SF - 4 - ix] = 0;
				}
			} else {
				for (; ix < 20; ix++) {
					if (scatPar[0][N_SF - 4 - ix] < scatPar[0][N_SF - 3])
						break;
					scatPar[0][N_SF - 4 - ix] = scatPar[0][N_SF - 3];
					for (unsigned iy = 1; iy < N_ELEM; iy++)
						scatPar[iy][N_SF - 4 - ix] = scatPar[iy][N_SF - 3];
				}
			}

		}        // end of if (scatPar[0][N_SF-4] > scatPar[0][N_SF-3])

		BOOST_LOG_TRIVIAL(info) << format("Will use %d sampling points per slice, total nz=%d (%d)")
						%	_nzPerSlice % _nz % (_nzPerSlice >> 1);
		BOOST_LOG_TRIVIAL(info) << format("dkx = %g, nx = %d, kmax2 = %g") % m_dkx % _nx % kmax2;
		BOOST_LOG_TRIVIAL(info)<< format("Cutoff scattering angle: kmax=%g (1/A), dk=(%g,%g %g)")
						% kmax2 % m_dkx % m_dkx % m_dkz;
		BOOST_LOG_TRIVIAL(info)<< format("getAtomPotential3D: set resolution of scattering factor to %g/A!")
						% scatPar[0][N_SF - 4 - ix];

		// Now kmax2 is actually kmax**2.
		kmax2 *= kmax2;
	}
}

/********************************************************************************
 * Create Lookup table for 3D potential due to neutral atoms
 ********************************************************************************/
#define PHI_SCALE 47.87658
void C3DFFTPotential::ComputeAtomPotential(std::vector<atom>::iterator &atom){
	float_tt B = m_tds ? 0 : atom->dw;
	float_tt kmax2;
	ComplexArray2D tmp;
	tmp.resize(boost::extents[_nx][_nz]);
	int Znum = atom->Znum;
	std::vector<float_tt> splinb(N_SF), splinc(N_SF), splind(N_SF);
	if (m_atPot.count(Znum) == 0) {
		kmax2 = 0.5 * _nx * m_dkx / OVERSAMPLING;
		kmax2 *= kmax2;

		// setup cubic spline interpolation:
		splinh((scatPar[0]), (scatPar[Znum]), splinb, splinc, splind, N_SF);
		m_atPot[Znum] = ComplexArray2D();
		m_atPot[Znum].resize(boost::extents[_nx/2][_nz/2]);

		float_tt kzmax = m_dkz * _config->Model.nSlices / 2.0;
		// define x-and z-position of atom center:
		// The atom should sit in the top-left corner,
		// however (nzPerSlice+1)/2 above zero in z-direction
		float_tt xPos = -2.0 * PI*0.0; // or m_dx*nx/(OVERSAMPLING), if in center
		int izOffset = (_nzPerSlice - 1) / 2;
		float_tt zPos = -2.0 * PI*(_config->Model.sliceThicknessAngstrom/_nzPerSlice*(izOffset));

		// What this look-up procedure will do is to supply V(r,z) computed from fe(q).
		// Since V(r,z) is rotationally symmetric we might as well compute
		// V(x,y,z) at y=0, i.e. V(x,z).
		// In order to do the proper 3D inverse FT without having to do a complete 3D FFT
		// we will pre-compute the qy-integral for y=0.

		// kzborder = dkz*(nz/(2*OVERSAMP_Z) -1);
		int iz =0, ix=0, iy=0;
		float_tt phase = 0, f=0;
		for (iz = 0; iz < _nz; iz++) {
			float_tt kz = m_dkz * (iz < _nz / 2 ? iz : iz - _nz);
			// We also need to taper off the potential in z-direction in order to avoid cutoff artifacts.
			// zScale = fabs(kz) <= kzborder ? 1.0 : 0.5+0.5*cos(M_PI*(fabs(kz)-kzborder)/(kzmax-kzborder));
			for (ix = 0; ix < _config->Model.nx; ix++) {
				float_tt kx = m_dkx * (ix < _nx / 2 ? ix : ix - _nx);
				float_tt kr = (kx * kx + kz * kz);
				// if this is within the allowed circle:
				if (kr < kmax2) {
					// unsigned ind3d = ix+iz*nx;
					// f = fe3D(Znum,k2,m_tds,1.0,m_scatFactor);
					// multiply scattering factor with Debye-Waller factor:
					// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
					f= seval(scatPar[0], scatPar[Znum], splinb, splinc, splind, N_SF, sqrt(kr))*exp(-kr * B * 0.25);
					// perform the qy-integration for qy <> 0:
					for (iy = 1; iy < _ny; iy++) {
						float_tt s3 = m_dkx * iy;
						s3 = s3 * s3 + kr;
						if (s3 < kmax2) {
							float_tt ss = seval(scatPar[0], scatPar[Znum], splinb, splinc, splind,N_SF, sqrt(s3))* exp(-s3 * B * 0.25);
							f += 2* ss;
						} else
							break;

					}
					f *= m_dkx;
					// note that the factor 2 is missing in the phase (2pi k*r)
					// this places the atoms in the center of the box.
					phase = kx * xPos + kz * zPos;
					//					printf("f = %f\n",f);
					tmp[ix][iz] = std::complex<float_tt>(f * cos(phase), f * sin(phase)); // *zScale
					//					printf("tmp[%d][%d] = (%f,%f)\n",ix,iz,f * cos(phase),f * sin(phase));
					// if ((kx==0) && (ky==0)) printf(" f=%g (%g, [%g, %g])\n",f,f*zScale,atPot[Znum][ind3d][0],atPot[Znum][ind3d][1]);
				}
			}
		} // for iz ...

#if SHOW_SINGLE_POTENTIAL
		// 0 thickness
		imageio = ImageIOPtr(new CImageIO(nz, nx, 0, dkz, dkx));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		sprintf(fileName,"pot_rec_%d.img",Znum);
		imageio->SetThickness(_config->Model.sliceThicknessAngstrom);
		imageio->WriteImage((void**)tmp, fileName);//TODO: modify writeImage
#endif        
		// This converts the 2D kx-kz map of the scattering factor to a 2D real space map.
		time_t time0, time1;
#if FLOAT_PRECISION ==1
		ComplexArray2D out;

		out.resize(boost::extents[_nx][_nz]);
		auto dft = cvmlcpp::DFT<float_tt, 2>(tmp,out,false, FFTW_ESTIMATE,1);
		time(&time0);
		dft.execute();
		time(&time1);
		//    fftwf_complex *ptr=reinterpret_cast<fftwf_complex*>(&temp[0]);
		//    fftwf_plan plan = fftwf_plan_dft_2d(nz,nx,ptr,ptr,FFTW_BACKWARD,FFTW_ESTIMATE);
		//    fftwf_execute(plan);
		//    fftwf_destroy_plan(plan);
#else
		ComplexArray2D out;
		out.resize(boost::extents[m_nx][_nz]);
		auto dft = cvmlcpp::DFT<float_tt, 2>(tmp,out,false, FFTW_ESTIMATE,1);
		dft.execute();

		//		fftw_complex *ptr=reinterpret_cast<fftw_complex*>(&temp[0]);
		//		fftw_plan plan = fftw_plan_dft_2d(nz,nx,ptr,ptr,FFTW_BACKWARD,FFTW_ESTIMATE);
		//		fftw_execute(plan);
		//		fftw_destroy_plan(plan);
#endif
		// We also make sure that the potential touches zero at least somewhere. This will avoid
		// sharp edges that could produce ringing artifacts.
		// It is certainly debatable whether this is a good apprach, or not.
		// printf("Setting up %d x %d potential for Znum=%d, (atom kind %d), Omega=%g\n",nx,nz,Znum,iKind,dkx*dky*dkz);
		// min = atPot[Znum][ny/2+nx/2*ny+nz/2*nx*ny][0];
//		BOOST_LOG_TRIVIAL(trace) << "Before z integration of potential:";
		for (unsigned ix = 0; ix < _nx / 2; ix++){
			for (unsigned iz = 0; iz < _nz / 2; iz++) {
				float_tt zScale = 0;
				// Integrate over nzPerSlice neighboring layers here:::::::::
				for (int iiz = -izOffset; iiz <= izOffset; iiz++) {
					if (iz + izOffset + iiz < _nz / 2){

						float_tt out1 = out[ix][(iz + izOffset + iiz)].real();
//						BOOST_LOG_TRIVIAL(trace) << boost::format("out[%d][%d] = %f") % ix % (iz + izOffset + iiz) % out[ix][(iz + izOffset + iiz)].real();
						zScale += out1*(_nx*_nz);
					}
				}
				if (zScale < 0)
					zScale = 0;

				// assign the iz-th slice the sum of the 3 other slices:
				// and divide by unit cell volume (since this is in 3D):
				// Where does the '1/2' come from??? OVERSAMPLING*OVERSAMP_Y/8 = 1/2
				// if nothing has changed, then OVERSAMPLING=2 OVERSAMP_Z=18.
				// remember, that s=0.5*k;
				// This potential will later again be scaled by lambda*gamma (=0.025*1.39139)
				// Sets real value to this; imaginary value to 0

				m_atPot[Znum][ix][iz] = 47.8658 * m_dkx * m_dkz / (_nz) * zScale;
				//				printf("m_atPot[%d][%d] = 47.8658 * %f * %f / (%d) * %f = %f\n",Znum,ind3d,dkx,dkz,nz,zScale,m_atPot[Znum][ind3d]);
			}
		}
		BOOST_LOG_TRIVIAL(info)<< format("Created 3D (r-z) %d x %d potential array for Z=%d (B=%g, dkx=%g, dky=%g. dkz=%g,sps=%d)")
						% (_nx / 2) % (_nz / 2) % Znum % B % m_dkx %  m_dkx % m_dkz % izOffset;
	}
}

/********************************************************************************
 * Lookup function for 3D potential offset due to charged atoms (ions)
 ********************************************************************************/
void C3DFFTPotential::GetAtomPotentialOffset3D(unsigned Znum, float_tt B,
		unsigned &nzSub, unsigned &Nr, unsigned &Nz_lut, float_tt q,
		ComplexVector &output) {
	/*
	 int ix,iy,iz,iiz,ind3d,iKind,izOffset;
	 double zScale,kzmax,zPos,xPos;
	 fftwf_plan plan;
	 static double f,phase,s2,s3,kmax2,kx,kz,dkx,dky,dkz; // ,dx2,dy2,dz2;
	 static int nx,ny,nz,nzPerSlice;
	 static fftwf_complex **atPot = NULL;
	 static fftwf_complex *temp = NULL;
	 #if SHOW_SINGLE_POTENTIAL == 1
	 ImageIOPtr imageio = ImageIOPtr();
	 static fftwf_complex *ptr = NULL;
	 static char fileName[256];
	 #endif
	 static double *splinb=NULL;
	 static double *splinc=NULL;
	 static double *splind=NULL;
	 */

	// if there is no charge to this atom, return NULL:
	if (q == 0)
		return;
#if !USE_REZ_SFACTS
	BOOST_LOG_TRIVIAL(info)<< format("Using charged atoms only works with scattering factors by Rez et al!") % Znum;
	exit(0);
#endif

	unsigned nx, ny, nz, nzPerSlice;
	float_tt dkx, dky, dkz;
	std::vector<float_tt> splinb(N_SF), splinc(N_SF), splind(N_SF);
	float_tt kmax2;

	ComplexVector temp;

	// scattering factors in:
	// float scatPar[4][30]
	if (m_offsetPot.size() == 0) {

		nx = 2 * OVERSAMPLING * (int) ceil(_config->Potential.AtomRadiusAngstrom / _config->Model.dx);
		ny = 2 * OVERSAMPLING * (int) ceil(_config->Potential.AtomRadiusAngstrom / _config->Model.dy);
		temp.resize(nx * nz);

		// The FFT-resolution in the z-direction must be high enough to avoid
		// artifacts due to premature cutoff of the rec. space scattering factor
		// we will therefore make it roughly the same as the x-resolution
		// However, we will make sure that a single slice contains an integer number
		// of sampling points.
		nzPerSlice = (int) floor(OVERSAMPLING * _config->Model.sliceThicknessAngstrom / _config->Model.dx);
		// make nzPerSlice odd:
		if (2.0 * floor((double) (nzPerSlice >> 1)) == nzPerSlice)
			nzPerSlice += 1;
		// Total number of z-positions should be twice that of atomRadius/sliceThickness
		nz = (2 * (int) ceil(_config->Potential.AtomRadiusAngstrom / _config->Model.sliceThicknessAngstrom)) * nzPerSlice;

		BOOST_LOG_TRIVIAL(info)<< format("Potential offset: will use %d sampling points per slice, total nz=%d (%d)")
						% nzPerSlice % nz % (nzPerSlice >> 1);

		dkx = 0.5 * OVERSAMPLING / (nx * _config->Model.dx);
		dky = 0.5 * OVERSAMPLING / (ny * _config->Model.dy);
		dkz = nzPerSlice / (double) (nz * _config->Model.sliceThicknessAngstrom);
		kmax2 = 0.5 * nx * dkx / (double) OVERSAMPLING; // largest k that we'll admit

		// printf("Cutoff scattering angle:kmax=%g, smax=%g (1/A), dk=(%g,%g %g)\n",kmax2,S_SCALE*kmax2,dkx,dky,dkz);
		scatParOffs[0][N_SF - 1] = 1.2 * kmax2;
		scatParOffs[0][N_SF - 2] = 1.1 * kmax2;
		scatParOffs[0][N_SF - 3] = kmax2;
		if (scatParOffs[0][N_SF - 4] > scatParOffs[0][N_SF - 3]) {
			unsigned ix = 0;
			if (1) {
				// set additional scattering parameters to zero:
				for (ix; ix < N_SF - 10; ix++) {
					if (scatParOffs[0][N_SF - 4 - ix]
					                   < scatParOffs[0][N_SF - 3] - 0.001 * (ix + 1))
						break;
					scatParOffs[0][N_SF - 4 - ix] = scatParOffs[0][N_SF - 3]
					                                               - 0.001 * (ix + 1);
					for (unsigned iy = 1; iy < N_ELEM; iy++)
						scatParOffs[iy][N_SF - 4 - ix] = 0;
				}
			} else {
				for (ix; ix < 20; ix++) {
					if (scatParOffs[0][N_SF - 4 - ix]
					                   < scatParOffs[0][N_SF - 3])
						break;
					scatParOffs[0][N_SF - 4 - ix] = scatParOffs[0][N_SF - 3];
					for (unsigned iy = 1; iy < N_ELEM; iy++)
						scatParOffs[iy][N_SF - 4 - ix] = scatParOffs[iy][N_SF
						                                                 - 3];
				}
			}
			BOOST_LOG_TRIVIAL(info)<< format("getAtomPotentialOffset3D: reduced angular range of scattering factor to %g/A!")
							% scatParOffs[0][N_SF - 4 - ix];
		} // end of if (scatParOffs[0][N_SF-4] > scatParOffs[0][N_SF-3])
		kmax2 *= kmax2;
	}
	// initialize this atom, if it has not been done yet:
	if (m_offsetPot.count(Znum) == 0) {
		// setup cubic spline interpolation:
		splinh(scatParOffs[0], scatParOffs[Znum], splinb, splinc, splind, N_SF);

		m_offsetPot[Znum] = ComplexVector(nx * nz / 4);
		memset((void*) &temp[0], 0, nx * nz * sizeof(complex_tt));

		//complex_tt *temp=(complex_tt *)fftw_malloc(nx*nz*sizeof(complex_tt));
		//memset(temp, 0, nx*nz*sizeof(complex_tt));

		float_tt kzmax = dkz * nz / 2.0;
		// define x-and z-position of atom center:
		// The atom should sit in the top-left corner,
		// however (nzPerSlice+1)/2 above zero in z-direction
		float_tt xPos = -2.0 * PI*0.0; // or m_dx*nx/(OVERSAMPLING), if in center
		unsigned izOffset = (nzPerSlice - 1) / 2;
		float_tt zPos = -2.0 * PI*(_config->Model.sliceThicknessAngstrom/nzPerSlice*(izOffset));

		// kzborder = dkz*(nz/(2*OVERSAMP_Z) -1);
		for (unsigned iz = 0; iz < nz; iz++) {
			float_tt kz = dkz * (iz < nz / 2 ? iz : iz - nz);
			// We also need to taper off the potential in z-direction
			// in order to avoid cutoff artifacts.
			// zScale = fabs(kz) <= kzborder ? 1.0 : 0.5+0.5*cos(M_PI*(fabs(kz)-kzborder)/(kzmax-kzborder));
			// printf("iz=%d, kz=%g, zScale=%g ",iz,kz,zScale);
			for (unsigned ix = 0; ix < nx; ix++) {
				float_tt kx = dkx * (ix < nx / 2 ? ix : ix - nx);
				float_tt s2 = (kx * kx + kz * kz);
				// if this is within the allowed circle:
				if (s2 < kmax2) {
					unsigned ind3d = ix + iz * nx;
					// f = fe3D(Znum,k2,m_tds,1.0,m_scatFactor);
					// multiply scattering factor with Debye-Waller factor:
					// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
					float_tt f = seval(scatParOffs[0], scatParOffs[Znum], splinb, splinc, splind, N_SF, sqrt(s2))
																			* exp(-s2 * B * 0.25);
					// perform the qy-integration for qy <> 0:
					for (unsigned iy = 1; iy < nx; iy++) {
						float_tt s3 = dky * iy;
						s3 = s3 * s3 + s2;
						if (s3 < kmax2) {
							f += 2	* seval(scatPar[0], scatPar[Znum], splinb, splinc, splind,N_SF, sqrt(s3))* exp(-s3 * B * 0.25);
						} else
							break;
					}
					f *= dkx;
					// note that the factor 2 is missing in the phase (2pi k*r)
					// this places the atoms in the center of the box.
					float_tt phase = kx * xPos + kz * zPos;
					temp[ind3d] = complex_tt(f * cos(phase), f * sin(phase)); // *zScale
					//          temp[ind3d][1] = f*sin(phase);        // *zScale
					// if ((kx==0) && (ky==0)) printf(" f=%g (%g, [%g, %g])\n",f,f*zScale,atPot[Znum][ind3d][0],atPot[Znum][ind3d][1]);
				}
			}
		} // for iz....

#if SHOW_SINGLE_POTENTIAL
		imageio = ImageIOPtr(new CImageIO(nz, nx, 0, dkx, dkz, std::vector<double>(),
				"rec. space potential"));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		imageio->SetThickness(_config->Model.sliceThicknessAngstrom);
		imageio->WriteImage((void**)temp, fileName);
#endif        

#if FLOAT_PRECISION ==1
		fftwf_complex *ptr = (fftwf_complex *) &temp[0];
		fftwf_plan plan = fftwf_plan_dft_2d(nz, nx, ptr, ptr, FFTW_BACKWARD,FFTW_ESTIMATE);
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);
#else
		fftw_complex *ptr=(fftw_complex *)&temp[0];
		fftw_plan plan = fftw_plan_dft_2d(nz,nx,ptr,ptr,FFTW_BACKWARD,FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);
#endif
		// We also make sure that the potential touches zero at least somewhere. This will avoid
		// sharp edges that could produce ringing artifacts.
		// It is certainly debatable whether this is a good apprach, or not.
		// printf("Setting up %d x %d potential for Znum=%d, (atom kind %d), Omega=%g\n",nx,nz,Znum,iKind,dkx*dky*dkz);
		// min = atPot[Znum][ny/2+nx/2*ny+nz/2*nx*ny][0];
		for (unsigned ix = 0; ix < nx / 2; ix++)
			for (unsigned iz = 0; iz < nz / 2; iz++) {
				unsigned ind3d = ix + iz * nx / 2;
				float_tt zScale = 0;
				// Integrate over nzPerSlice neighboring layers here:::::::::
				for (int iiz = -izOffset; iiz <= izOffset; iiz++) {
					if (iz + izOffset + iiz < nz / 2)
						zScale += temp[ix + (iz + izOffset + iiz) * nx].real();
				}
				if (zScale < 0)
					zScale = 0;
				// assign the iz-th slice the sum of the 3 other slices:
				// and divide by unit cell volume (since this is in 3D):
				// Where does the '1/2' come from??? OVERSAMPLING*OVERSAMP_Y/8 = 1/2
				// if nothing has changed, then OVERSAMPLING=2 OVERSAMP_Z=18.
				// remember, that s=0.5*k;
				// This potential will later again be scaled by lambda*gamma (=0.025*1.39139)
				//    implicitly sets imaginary part to 0.
				m_offsetPot[Znum][ind3d] = 47.8658 * dkx * dkz / (nz) * zScale;
			}
#if SHOW_SINGLE_POTENTIAL
		imageio = ImageIOPtr(new CImageIO(nz/2, nx/2, 0, m_dx/OVERSAMPLING,
				_config->Model.sliceThicknessAngstrom/nzPerSlice));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		imageio->SetThickness(nz*_config->Model.sliceThicknessAngstrom/nzPerSlice);
		sprintf(fileName,"potentialOffs_rz_%d.img",Znum);
		ptr = atPot[Znum];
		imageio->WriteImage((void**)ptr, fileName);
#endif        
		if (m_printLevel > 1)
			BOOST_LOG_TRIVIAL(info) <<
			format("Created 3D (r-z) %d x %d potential offset array for Z=%d (B=%g, dkx=%g, dky=%g. dkz=%g,sps=%d)")
			%(nx / 2)%( nz / 2)% Znum% B% dkx% dky% dkz% izOffset;
	} // end of creating atom if not exist...
	Nz_lut = nz / 2;
	nzSub = nzPerSlice;
	Nr = nx / 2;
	output = m_offsetPot[Znum];
}
// #undef SHOW_SINGLE_POTENTIAL
// end of fftwf_complex *getAtomPotential3D(...)

}// end namespace QSTEM
