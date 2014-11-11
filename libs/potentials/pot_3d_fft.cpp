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

namespace QSTEM {

C3DFFTPotential::C3DFFTPotential() :C3DPotential() {}

C3DFFTPotential::C3DFFTPotential(const ConfigPtr configReader) :	C3DPotential(configReader), m_dkz(0), m_dkx(0), m_nz(0), m_nzPerSlice(0) {
	for (unsigned i = 0; i < m_nslices; i++) {
		if (m_sliceThicknesses[0] != m_sliceThicknesses[i]) {
			printf("Warning: slice thickness not constant, will give wrong results (iz=%d)!\n",i);
		}
	}
}

void C3DFFTPotential::DisplayParams() {
	CPotential::DisplayParams();
	printf("* Potential calculation: 3D (FFT method)");
}

void C3DFFTPotential::AddAtomNonPeriodic(std::vector<atom>::iterator &atom,	float_tt atomBoxX, unsigned int iAtomX, float_tt atomBoxY, unsigned int iAtomY, float_tt atomZ) {

	unsigned iAtomZ = (int) floor(atomZ / m_sliceThickness + 0.5);
	unsigned nzSub, nRadius, Nz_lut;
	ComplexVector atPotPtr = m_atPot[atom->Znum];

	nRadius = 2 * OVERSAMPLING * (int) ceil(m_atomRadius / m_dx) /2;
	nzSub = (int) floor(OVERSAMPLING * m_sliceThickness / m_dx);
	if (2.0 * (nzSub >> 1) == nzSub)	nzSub += 1;
	Nz_lut = (2 * (int) ceil(m_atomRadius / m_sliceThickness)) * nzSub / 2;

	int xstart = (int) iAtomX - (int) m_iRadX < 0 ? 0 : iAtomX - m_iRadX;
	int xend =(int) iAtomX + (int) m_iRadX >= m_nx ? m_nx - 1 : iAtomX + m_iRadX;
	int ystart = (int) iAtomY - (int) m_iRadY < 0 ? 0 : iAtomY - m_iRadY;
	int yend =	(int) iAtomY + (int) m_iRadY >= m_ny ? m_ny - 1 : iAtomY + m_iRadY;
	int zstart =  (int)iAtomZ -  (int)m_iRadZ < 0 ? -iAtomZ : -m_iRadZ;
	int zend =	 (int)iAtomZ +  (int)m_iRadZ >= m_nslices ?	m_nslices - iAtomZ - 1 : m_iRadZ;

	// if within the potential map range:
	if ((xstart < m_nx) && (xend >= 0) && (ystart < m_ny) && (yend >= 0)) {
		// define range of sampling from atomZ-/+atomRadius
		if (((int) iAtomZ + (int) zstart < (int) m_nslices)	&& ((int) iAtomZ + (int) zend >= 0)) {
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
				float_tt x2 = iax * m_dx - atomBoxX;
				x2 *= x2;
				for (iay = ystart; iay <= yend; iay++) {
					// printf("iax=%d, iay=%d\n",iax,iay);
					float_tt y2 = iay * m_dy - atomBoxY;
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
								- atomZ / m_sliceThickness) * nzSub;
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
									potVal = (1 - ddr)* atPotPtr[ir - iOffsetZ + nRadius].real()+ ddr* atPotPtr[ir + 1 - iOffsetZ + nRadius].real();
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
									potVal =(1 - ddr)* atPotPtr[ir + iOffsetZ].real()+ ddr* atPotPtr[ir + 1+ iOffsetZ].real();
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
	unsigned iAtomX = (int) floor(atomX / m_dx);
	unsigned iAtomY = (int) floor(atomY / m_dy);

	if (m_periodicXY) {
		AddAtomPeriodic(atom, atomX, iAtomX, atomY, iAtomY, atomZ);
	} else {
		AddAtomNonPeriodic(atom, atomX, iAtomX, atomY, iAtomY, atomZ);
	}
}

void C3DFFTPotential::AddAtomPeriodic(std::vector<atom>::iterator &atom,float_tt atomBoxX, unsigned int iAtomX, float_tt atomBoxY,unsigned int iAtomY, float_tt atomZ) {
	unsigned iAtomZ = (int) floor(atomZ / m_sliceThickness + 0.5);
	unsigned iax0 = iAtomX - m_iRadX;
	unsigned iax1 = iAtomX + m_iRadX;
	unsigned iay0 = iAtomY - m_iRadY;
	unsigned iay1 = iAtomY + m_iRadY;

	unsigned nzSub, Nr, Nz_lut;

	Nr = 2 * OVERSAMPLING * (int) ceil(m_atomRadius / m_dx) /2;
	nzSub = (int) floor(OVERSAMPLING * m_sliceThickness / m_dx);
	if (2.0 * (nzSub >> 1) == nzSub)	nzSub += 1;
	Nz_lut = (2 * (int) ceil(m_atomRadius / m_sliceThickness)) * nzSub / 2;

	// define range of sampling from atomZ-/+atomRadius

	int iaz0 = iAtomZ - m_iRadZ < 0 ? -iAtomZ : -m_iRadZ;
	unsigned iaz1 =	iAtomZ + m_iRadZ >= m_nslices ? m_nslices - iAtomZ - 1 : m_iRadZ;

	// if (iatom < 2) printf("iatomZ: %d, %d cz=%g, %g: %d, %d\n",iAtomZ,iaz0,m_sliceThickness,atomZ,(int)(-2.5-atomZ),(int)(atomZ+2.5));
	// printf("%d: iatomZ: %d, %d cz=%g, %g\n",iatom,iAtomZ,iaz0,m_sliceThickness,atomZ);
	bool isAtomZWithinSlices = (iAtomZ + iaz0 < m_nslices) && (iAtomZ + iaz1 >= 0);
	if (isAtomZWithinSlices) {
		// retrieve the pointer for this atom
		ComplexVector atPotPtr = m_atPot[atom->Znum];

#if USE_Q_POT_OFFSETS
		// retrieve the pointer to the array of charge-dependent potential offset
		// This function will return NULL; if the charge of this atom is zero:
		ComplexVector atPotOffsPtr;
		getAtomPotentialOffset3D(atoms[iatom].Znum,muls,m_tds ? 0 : atoms[iatom].dw,&nzSub,&Nr,&Nz_lut,atoms[iatom].q, atPotOffsPtr);
#endif // USE_Q_POT_OFFSETS

		unsigned iOffsLimHi = Nr * (Nz_lut - 1);
		int iOffsLimLo = -Nr * (Nz_lut - 1);
		unsigned iOffsStep = nzSub * Nr;

		// Slices around the slice that this atom is located in must be affected by this atom:
		// iaz must be relative to the first slice of the atom potential box.
		for (unsigned iax = iax0; iax < iax1; iax++) {
			//			unsigned idx = ((iax + 2 * m_nx) % m_nx) * m_ny	+ (iay0 + 2 * m_ny) % m_ny;
			//			complex_tt *potPtr = &(m_trans[iAtomZ + iaz0][idx]);
			float_tt x2 = iax * m_dx - atomBoxX;
			x2 *= x2;

			for (unsigned iay = iay0; iay < iay1; iay++) {
				// printf("iax=%d, iay=%d\n",iax,iay);
				float_tt y2 = iay * m_dy - atomBoxY;
				y2 *= y2;
				float_tt r = sqrt(x2 + y2);
				float_tt ddr = r / m_dr;
				unsigned ir = (int) floor(ddr);
				// add in different slices, once r has been defined
				if (ir < Nr - 1) {
					ddr = ddr - (double) ir;
					//					complex_tt *ptr = potPtr;
					// Include interpolation in z-direction as well (may do it in a very smart way here !):

					float_tt dOffsZ = (iAtomZ + iaz0 - atomZ / m_sliceThickness)* nzSub;
#if Z_INTERPOLATION
					unsigned iOffsZ = (int)dOffsZ;
					float_tt ddz = fabs(dOffsZ - (double)iOffsZ);
#else // Z_INTERPOLATION
					unsigned iOffsZ = (int) (dOffsZ + 0.5);
#endif // Z_INTERPOLATION
					iOffsZ *= Nr;

					for (int iaz = iaz0; iaz <= iaz1; iaz++) {
						float_tt potVal = 0;
						// iOffsZ = (int)(fabs(iAtomZ+iaz-atomZ/m_sliceThickness)*nzSub+0.5);
						if (iOffsZ < 0) {
							if (iOffsZ > iOffsLimLo) {
								// do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
								potVal = (1-ddz)*((1-ddr)*atPotPtr[ir-iOffsZ+Nr][0]+ddr*atPotPtr[ir+1-iOffsZ+Nr][0])+
										ddz *((1-ddr)*atPotPtr[ir-iOffsZ][0]+ddr*atPotPtr[ir+1-iOffsZ][0]);
#if USE_Q_POT_OFFSETS
								// add the charge-dependent potential offset
								if (atPotOffsPtr != NULL) {
									potVal += atoms[iatom].q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir-iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1-iOffsZ+Nr][0])+
											ddz *((1-ddr)*atPotOffsPtr[ir-iOffsZ ][0]+ddr*atPotOffsPtr[ir+1-iOffsZ ][0]));
								}
#endif // USE_Q_POT_OFFSETS
#else // Z_INTERPOLATION
								potVal =(1 - ddr)* atPotPtr[ir - iOffsZ + Nr].real()+ ddr* atPotPtr[ir + 1- iOffsZ + Nr].real();
#if USE_Q_POT_OFFSETS
								// add the charge-dependent potential offset
								if (atPotOffsPtr != NULL) {
									potVal += atoms[iatom].q*((1-ddr)*atPotOffsPtr[ir-iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1-iOffsZ+Nr][0]);
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
								potVal = (1-ddz)*((1-ddr)*atPotPtr[ir+iOffsZ][0]+ddr*atPotPtr[ir+1+iOffsZ][0])+
										ddz *((1-ddr)*atPotPtr[ir+iOffsZ+Nr][0]+ddr*atPotPtr[ir+1+iOffsZ+Nr][0]);
#if USE_Q_POT_OFFSETS
								// add the charge-dependent potential offset
								if (atPotOffsPtr != NULL) {
									potVal += atoms[iatom].q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir+iOffsZ][0]+ddr*atPotOffsPtr[ir+1+iOffsZ][0])+
											ddz *((1-ddr)*atPotOffsPtr[ir+iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1+iOffsZ+Nr][0]));
								}
#endif // USE_Q_POT_OFFSETS
#else // Z_INTERPOLATION
								potVal =(1 - ddr) * atPotPtr[ir + iOffsZ].real()+ ddr* atPotPtr[ir + 1+ iOffsZ].real();
#if USE_Q_POT_OFFSETS
								// add the charge-dependent potential offset
								if (atPotOffsPtr != NULL) {
									potVal += atoms[iatom].q*((1-ddr)*atPotOffsPtr[ir+iOffsZ][0]+ddr*atPotOffsPtr[ir+1+iOffsZ][0]);
								}
#endif // USE_Q_POT_OFFSETS
#endif // Z_INTERPOLATION
							}
						}

						m_trans1[iAtomZ + iaz][(iay+m_ny) % m_ny][(iax+m_nx) % m_nx]+= potVal;
						//						*ptr += potVal;

						//						ptr += m_sliceStep; // advance to the next slice
						// add the remaining potential to the next slice:
						// if (iaz < iaz1) *ptr += (1-ddz)*potVal;
						iOffsZ += iOffsStep;
//						printf("_trans1[%d + %d][%d][%d] += %f\n",iAtomZ,iaz,iay,iax,potVal);
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
	} // iaz0+iAtomZ < m_slices
}

void C3DFFTPotential::SliceSetup(){
	CPotential::SliceSetup();
	float_tt kmax2;

	if (m_atPot.size() == 0) {
		m_nx = 2 * OVERSAMPLING * (int) ceil(m_atomRadius / m_dx);
		m_ny = 2 * OVERSAMPLING * (int) ceil(m_atomRadius / m_dy);

		// The FFT-resolution in the z-direction must be high enough to avoid
		// artifacts due to premature cutoff of the rec. space scattering factor
		// we will therefore make it roughly the same as the x-resolution
		// However, we will make sure that a single slice contains an integer number
		// of sampling points.
		m_nzPerSlice = (int) floor(OVERSAMPLING * m_sliceThickness / m_dx);

		// make nzPerSlice odd:
		if (2.0 * (m_nzPerSlice >> 1) == m_nzPerSlice)	m_nzPerSlice += 1;

		// Total number of z-positions should be twice that of atomRadius/sliceThickness
		m_nz = (2 * (int) ceil(m_atomRadius / m_sliceThickness)) * m_nzPerSlice;



		m_dkx = 0.5 * OVERSAMPLING / (m_nx * m_dx); // nx*m_dx is roughly 2*m_atomRadius
		m_dkz = m_nzPerSlice / (m_nz * m_sliceThickness);
		kmax2 = 0.5 * m_nx * m_dkx / OVERSAMPLING;
		// Don't square kmax2 yet!

		scatPar[0][N_SF - 1] = 1.2 * kmax2;
		scatPar[0][N_SF - 2] = 1.1 * kmax2;
		scatPar[0][N_SF - 3] = kmax2;
		// adjust the resolution of the lookup table if necessary
		unsigned ix = 0;
		if (scatPar[0][N_SF - 4] > scatPar[0][N_SF - 3]) {

			if (1) {
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

		if (m_printLevel > 1){
			printf("Will use %d sampling points per slice, total nz=%d (%d)\n",	m_nzPerSlice, m_nz, m_nzPerSlice >> 1);
			printf("dkx = %g, nx = %d, kmax2 = %g\n", m_dkx, m_nx, kmax2);
			printf("Cutoff scattering angle: kmax=%g (1/A), dk=(%g,%g %g)\n", kmax2, m_dkx, m_dkx, m_dkz);
			printf("getAtomPotential3D: set resolution of scattering factor to %g/A!\n",scatPar[0][N_SF - 4 - ix]);
		}

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

	// largest k that we'll admit

	float_tt kmax2;
	//  ComplexVector temp;
	ComplexArray2D tmp;
	tmp.resize(boost::extents[m_nx][m_nz]);
	int Znum = atom->Znum;


	std::vector<float_tt> splinb(N_SF), splinc(N_SF), splind(N_SF);

	// scattering factors in:
	// float scatPar[4][30]
	// TODO: this section meant to speed things up by skipping initialization.  Can we keep values
	//    as members instead, and move init to another function?



	// initialize this atom, if it has not been done yet:

	if (m_atPot.count(Znum) == 0) {
		kmax2 = 0.5 * m_nx * m_dkx / OVERSAMPLING;
		kmax2 *= kmax2;

		// setup cubic spline interpolation:
		splinh((scatPar[0]), (scatPar[Znum]), splinb, splinc, splind, N_SF);

		// allocate a 3D array:
		m_atPot[Znum] = ComplexVector(m_nx * m_nz / 4);

		float_tt kzmax = m_dkz * m_nz / 2.0;
		// define x-and z-position of atom center:
		// The atom should sit in the top-left corner,
		// however (nzPerSlice+1)/2 above zero in z-direction
		float_tt xPos = -2.0 * PI*0.0; // or m_dx*nx/(OVERSAMPLING), if in center
		int izOffset = (m_nzPerSlice - 1) / 2;
		float_tt zPos = -2.0 * PI*(m_sliceThickness/m_nzPerSlice*(izOffset));

		// What this look-up procedure will do is to supply V(r,z) computed from fe(q).
		// Since V(r,z) is rotationally symmetric we might as well compute
		// V(x,y,z) at y=0, i.e. V(x,z).
		// In order to do the proper 3D inverse FT without having to do a complete 3D FFT
		// we will pre-compute the qy-integral for y=0.

		// kzborder = dkz*(nz/(2*OVERSAMP_Z) -1);
		int iz =0, ix=0, iy=0;
		float_tt phase = 0, f=0;
		for (iz = 0; iz < m_nz; iz++) {
			float_tt kz = m_dkz * (iz < m_nz / 2 ? iz : iz - m_nz);
			// We also need to taper off the potential in z-direction
			// in order to avoid cutoff artifacts.
			// zScale = fabs(kz) <= kzborder ? 1.0 : 0.5+0.5*cos(M_PI*(fabs(kz)-kzborder)/(kzmax-kzborder));
			// printf("iz=%d, kz=%g, zScale=%g ",iz,kz,zScale);
			for (ix = 0; ix < m_nx; ix++) {
				float_tt kx = m_dkx * (ix < m_nx / 2 ? ix : ix - m_nx);
				float_tt s2 = (kx * kx + kz * kz);
				// if this is within the allowed circle:
				if (s2 < kmax2) {
					// unsigned ind3d = ix+iz*nx;
					// f = fe3D(Znum,k2,m_tds,1.0,m_scatFactor);
					// multiply scattering factor with Debye-Waller factor:
					// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
					f= seval(scatPar[0], scatPar[Znum], splinb, splinc, splind, N_SF, sqrt(s2))*exp(-s2 * B * 0.25);
					// perform the qy-integration for qy <> 0:
					for (iy = 1; iy < m_ny; iy++) {
						float_tt s3 = m_dkx * iy;
						s3 = s3 * s3 + s2;
						if (s3 < kmax2) {

							float_tt ss = seval(scatPar[0], scatPar[Znum], splinb, splinc, splind,N_SF, sqrt(s3))* exp(-s3 * B * 0.25);
//							printf("ss = %f\n",ss);
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
		imageio->SetThickness(m_sliceThickness);
		imageio->WriteImage((void**)tmp, fileName);//TODO: modify writeImage
#endif        
		// This converts the 2D kx-kz map of the scattering factor to a 2D real space map.
		time_t time0, time1;
#if FLOAT_PRECISION ==1
		ComplexArray2D out;

		out.resize(boost::extents[m_nx][m_nz]);
		auto dft = cvmlcpp::DFT<float_tt, 2>(tmp,out,false, FFTW_ESTIMATE,1);
		time(&time0);
		dft.execute();
		time(&time1);
//		printf(	"%g sec used for FFT\n",	difftime(time1, time0));
		//    fftwf_complex *ptr=reinterpret_cast<fftwf_complex*>(&temp[0]);
		//    fftwf_plan plan = fftwf_plan_dft_2d(nz,nx,ptr,ptr,FFTW_BACKWARD,FFTW_ESTIMATE);
		//    fftwf_execute(plan);
		//    fftwf_destroy_plan(plan);
#else
		ComplexArray2D out;
		out.resize(boost::extents[m_nx][m_nz]);
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
		for (unsigned ix = 0; ix < m_nx / 2; ix++)
			for (unsigned iz = 0; iz < m_nz / 2; iz++) {
				float_tt zScale = 0;
				unsigned ind3d = ix + iz * m_nx / 2;
				// Integrate over nzPerSlice neighboring layers here:::::::::
				for (int iiz = -izOffset; iiz <= izOffset; iiz++) {
					if (iz + izOffset + iiz < m_nz / 2){

						float_tt out1 = out[ix][(iz + izOffset + iiz)].real();
//						printf("out[%d][%d] = %f\n",ix,iz + izOffset + iiz,out[ix][(iz + izOffset + iiz)].real());
						zScale += out1*(m_nx*m_nz);
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

				m_atPot[Znum][ind3d] = 47.8658 * m_dkx * m_dkz / (m_nz) * zScale;
//				printf("m_atPot[%d][%d] = 47.8658 * %f * %f / (%d) * %f = %f\n",Znum,ind3d,dkx,dkz,nz,zScale,m_atPot[Znum][ind3d]);
			}
		// make sure we don't produce negative potential:
		if (m_printLevel > 1)
			printf("Created 3D (r-z) %d x %d potential array for Z=%d (B=%g, dkx=%g, dky=%g. dkz=%g,sps=%d)\n",m_nx / 2, m_nz / 2, Znum, B, m_dkx, m_dkx, m_dkz, izOffset);
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
	printf(	"Using charged atoms only works with scattering factors by Rez et al!\n",Znum);
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

		nx = 2 * OVERSAMPLING * (int) ceil(m_atomRadius / m_dx);
		ny = 2 * OVERSAMPLING * (int) ceil(m_atomRadius / m_dy);
		temp.resize(nx * nz);

		// The FFT-resolution in the z-direction must be high enough to avoid
		// artifacts due to premature cutoff of the rec. space scattering factor
		// we will therefore make it roughly the same as the x-resolution
		// However, we will make sure that a single slice contains an integer number
		// of sampling points.
		nzPerSlice = (int) floor(OVERSAMPLING * m_sliceThickness / m_dx);
		// make nzPerSlice odd:
		if (2.0 * floor((double) (nzPerSlice >> 1)) == nzPerSlice)
			nzPerSlice += 1;
		// Total number of z-positions should be twice that of atomRadius/sliceThickness
		nz = (2 * (int) ceil(m_atomRadius / m_sliceThickness)) * nzPerSlice;

		if (m_printLevel > 1)
			printf("Potential offset: will use %d sampling points per slice, total nz=%d (%d)\n",nzPerSlice, nz, nzPerSlice >> 1);

		dkx = 0.5 * OVERSAMPLING / (nx * m_dx);
		dky = 0.5 * OVERSAMPLING / (ny * m_dy);
		dkz = nzPerSlice / (double) (nz * m_sliceThickness);
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
			if (m_printLevel > 1)
				printf("getAtomPotentialOffset3D: reduced angular range of scattering factor to %g/A!\n",scatParOffs[0][N_SF - 4 - ix]);
		} // end of if (scatParOffs[0][N_SF-4] > scatParOffs[0][N_SF-3])
		kmax2 *= kmax2;

		//atPot = (fftwf_complex **)malloc((N_ELEM+1)*sizeof(fftwf_complex *));
		//for (unsigned ix=0;ix<=N_ELEM;ix++) atPot[ix] = NULL;
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
		float_tt zPos = -2.0 * PI*(m_sliceThickness/nzPerSlice*(izOffset));

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
		imageio->SetThickness(m_sliceThickness);
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
				m_sliceThickness/nzPerSlice));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		imageio->SetThickness(nz*m_sliceThickness/nzPerSlice);
		sprintf(fileName,"potentialOffs_rz_%d.img",Znum);
		ptr = atPot[Znum];
		imageio->WriteImage((void**)ptr, fileName);
#endif        
		if (m_printLevel > 1)
			printf("Created 3D (r-z) %d x %d potential offset array for Z=%d (B=%g, dkx=%g, dky=%g. dkz=%g,sps=%d)\n",nx / 2, nz / 2, Znum, B, dkx, dky, dkz, izOffset);
	} // end of creating atom if not exist...
	Nz_lut = nz / 2;
	nzSub = nzPerSlice;
	Nr = nx / 2;
	output = m_offsetPot[Znum];
}
// #undef SHOW_SINGLE_POTENTIAL
// end of fftwf_complex *getAtomPotential3D(...)

}// end namespace QSTEM
