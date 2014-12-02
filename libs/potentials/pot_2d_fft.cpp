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

#include "pot_2d_fft.hpp"
#include <boost/format.hpp>
using boost::format;
namespace QSTEM {


C2DFFTPotential::C2DFFTPotential(const ConfigPtr configReader) : CPotential(configReader) {
	m_atPot = std::map<unsigned, ComplexArray2D>();
}


void C2DFFTPotential::DisplayParams() {
	CPotential::DisplayParams();
	BOOST_LOG_TRIVIAL(info)<<format("* Potential calculation: 2D (FFT method)");
}

void C2DFFTPotential::MakeSlices(int nlayer, StructurePtr crystal) {
	CPotential::MakeSlices(nlayer, crystal);
	/* check whether we have constant slice thickness */
	for (unsigned i = 0; i < nlayer; i++) {
		if (m_sliceThicknesses[0] != m_sliceThicknesses[i]) {
			printf( "Warning: slice thickness not constant, will give wrong results (iz=%d)!\n",
					i);
			break;
		}
	}

}

void C2DFFTPotential::AddAtomToSlices(std::vector<atom>::iterator &atom,
		float_tt atomX, float_tt atomY, float_tt atomZ) {
	unsigned iAtomX = (int) floor(atomX / _config->Model.dx);
	unsigned iAtomY = (int) floor(atomY / _config->Model.dy);

	if (_config->Potential.periodicXY) {
		AddAtomPeriodic(atom, atomX, iAtomX, atomY, iAtomY, atomZ);
	} else {
		AddAtomNonPeriodic(atom, atomX, iAtomX, atomY, iAtomY, atomZ);
	}
}
void C2DFFTPotential::CenterAtomZ(std::vector<atom>::iterator &atom, float_tt &z) {

}
void C2DFFTPotential::AddAtomNonPeriodic(std::vector<atom>::iterator &atom,float_tt atomBoxX, int iAtomX, float_tt atomBoxY,int iAtomY, float_tt atomZ) {
	int iAtomZ = (int) floor(atomZ / _config->Model.sliceThicknessAngstrom);
	int iax0, iay0, potentialOffsetX =0, potentialOffsetY=0;
	if(iAtomX - m_iRadX < 0 ){
		iax0 = 0;
		potentialOffsetX = abs(iAtomX - m_iRadX) * OVERSAMPLING;
	} else {
		iax0 = iAtomX - m_iRadX;
	}
	if(iAtomY - m_iRadY < 0 ){
		iay0 =0;
		potentialOffsetY = abs(iAtomY - m_iRadY) * OVERSAMPLING;
	} else {
		iay0 = iAtomY - m_iRadY;
	}
	int iax1 = iAtomX + m_iRadX >= _config->Model.nx ? _config->Model.nx - 1 : iAtomX + m_iRadX;
	int iay1 = iAtomY + m_iRadY >= _config->Model.ny ? _config->Model.ny - 1 : iAtomY + m_iRadY;
	float_tt ddx = (atomBoxX/_config->Model.dx - iAtomX);
	float_tt ddy = (atomBoxY/_config->Model.dy - iAtomY);
	int iOffsX = (int) floor(ddx);
	int iOffsY = (int) floor(ddy);
	ddx -= (double) iOffsX;
	ddy -= (double) iOffsY;
	float_tt s11 = (1 - ddx) * (1 - ddy);
	float_tt s12 = (1 - ddx) * ddy;
	float_tt s21 = ddx * (1 - ddy);
	float_tt s22 = ddx * ddy;
	ComplexArray2D pot = m_atPot[atom->Znum];
	complex_tt added = complex_tt(0,0);
	BOOST_LOG_TRIVIAL(trace)<< format("atom xyz (%-02.3f,%-02.3f,%-02.3f) Ixyz (%-3d,%-3d,%-3d) iax (%-3d .. %-3d) iay (%-3d .. %-3d)")
			% atom->x % atom->y % atom->z % iAtomX % iAtomY % iAtomZ % iax0 %iax1%iay0%iay1;



	for (int iax = iax0; iax < iax1; iax++) {
		int idx = iax * _ny + iay0;
		for (int iay = iay0; iay < iay1; iay++) {
			int xindex = iOffsX + OVERSAMPLING * (iax - iax0) + potentialOffsetX;
			int yindex = iOffsY + OVERSAMPLING * (iay-iay0) + potentialOffsetY;
			float_tt vz = (s11 * pot[xindex][yindex]
						+ s12 * pot[xindex+1][yindex]
						+ s21 * pot[xindex][yindex+1]
						+ s22 * pot[xindex+1][yindex+1]).real();
			m_trans1[iAtomZ][iax][iay] += vz;

			added += vz;
		}
	}
//	BOOST_LOG_TRIVIAL(trace)<< format("added (%g,%g) to potential") % added.real() % added.imag();
	if(added.real()==0 && added.imag()==0)
		added = complex_tt(1,1);
}

void C2DFFTPotential::AddAtomPeriodic(std::vector<atom>::iterator &atom,
		float_tt atomBoxX, int iAtomX, float_tt atomBoxY,
		int iAtomY, float_tt atomZ) {
	unsigned iAtomZ = (int) floor(atomZ / _config->Model.sliceThicknessAngstrom);
	unsigned iax0 = iAtomX - m_iRadX +  _nx;
	unsigned iax1 = iAtomX + m_iRadX +  _nx;
	unsigned iay0 = iAtomY - m_iRadY +  _ny;
	unsigned iay1 = iAtomY + m_iRadY +  _ny;

	float_tt ddx = (-(double) iax0- _nx+ atomBoxX / _config->Model.dx - (double) m_iRadX )* (double) OVERSAMPLING;
	float_tt ddy = (-(double) iay0- _ny+ atomBoxY / _config->Model.dx - (double) m_iRadY )* (double) OVERSAMPLING;
	unsigned iOffsX = (int) floor(ddx);
	unsigned iOffsY = (int) floor(ddy);
	ddx -= (double) iOffsX;
	ddy -= (double) iOffsY;

	float_tt s22 = (1 - ddx) * (1 - ddy);
	float_tt s21 = (1 - ddx) * ddy;
	float_tt s12 = ddx * (1 - ddy);
	float_tt s11 = ddx * ddy;

	ComplexArray2D pot = m_atPot[atom->Znum];

	for (int iax = iax0; iax < iax1; iax++) { // TODO: should use ix += OVERSAMPLING
		for (int iay = iay0; iay < iay1; iay++) {
			int idx = (iax % _nx) * _ny + (iay % _ny);
			int xindex = (iOffsX + OVERSAMPLING * (iax - iax0));
			int yindex = iOffsY + OVERSAMPLING * (iay-iay0);
			m_trans1[iAtomZ][((iAtomX + iax) % _nx)][((iAtomY+iay) % _ny)] +=
					s11 * pot[xindex][yindex]
					                  + s12 * pot[xindex+1][yindex]
					                                        + s21 * pot[xindex][yindex+1]
					                                                            + s22 * pot[xindex+1][yindex+1];
		}
	}
}
void C2DFFTPotential::SliceSetup() {
	CPotential::SliceSetup();
	if (m_atPot.size() == 0) {
		_nx = 2 * OVERSAMPLING * (int) ceil(_config->Potential.AtomRadiusAngstrom / _config->Model.dx);
		_ny = 2 * OVERSAMPLING * (int) ceil(_config->Potential.AtomRadiusAngstrom / _config->Model.dy);
		m_dkx = 0.5 * OVERSAMPLING / ((_nx) * _config->Model.dx);
		m_dky = 0.5 * OVERSAMPLING / ((_ny) * _config->Model.dy);
		m_kmax2 = 0.5 * _nx * m_dkx / (double) OVERSAMPLING; // largest k that we'll admit

		BOOST_LOG_TRIVIAL(info)<< format("Cutoff scattering angle:kmax=%g (1/A)") % m_kmax2;
		scatPar[0][N_SF - 1] = 1.2 * m_kmax2;
		scatPar[0][N_SF - 2] = 1.1 * m_kmax2;
		scatPar[0][N_SF - 3] = m_kmax2;
		if (scatPar[0][N_SF - 4] > scatPar[0][N_SF - 3]) {
			unsigned ix = 0;
			if (1) {
				// set additional scattering parameters to zero:
				for (ix; ix < N_SF - 10; ix++) {
					if (scatPar[0][N_SF - 4 - ix] < scatPar[0][N_SF - 3] - 0.001 * (ix + 1))
						break;
					scatPar[0][N_SF - 4 - ix] = scatPar[0][N_SF - 3] - 0.001 * (ix + 1);
					for (unsigned iy = 1; iy < N_ELEM; iy++)
						scatPar[iy][N_SF - 4 - ix] = 0;
				}
			}

			if (_config->Output.LogLevel < 2)
				printf(
						"getAtomPotential2D: reduced angular range of scattering factor to %g/A!\n",
						scatPar[0][N_SF - 4 - ix]);
		}
		m_kmax2 *= m_kmax2;
	}
}
void  C2DFFTPotential::ComputeAtomPotential(std::vector<atom>::iterator &atom){
	std::vector<float_tt> splinb(N_SF, 0), splinc(N_SF, 0), splind(N_SF, 0);
	float_tt B = _config->Model.UseTDS ? 0 : atom->dw;
	int Znum = atom->Znum;
	if (m_atPot.count(Znum) == 0) {
		// setup cubic spline interpolation:
		splinh(scatPar[0], scatPar[Znum], splinb, splinc, splind,N_SF);
		m_atPot[Znum] = ComplexArray2D();
		m_atPot[Znum].resize(boost::extents[_nx][_ny]);
		for (unsigned ix = 0; ix < _nx; ix++) {
			float_tt kx = m_dkx * (ix < _nx / 2 ? ix : _nx - ix);
			for (unsigned iy = 0; iy < _ny; iy++) {
				float_tt ky = m_dky * (iy < _ny / 2 ? iy : _ny - iy);
				float_tt s2 = (kx * kx + ky * ky);
				// if this is within the allowed circle:
				if (s2 < m_kmax2) {
					// f = fe3D(Znum,k2,m_tds,1.0,m_scatFactor);
					// multiply scattering factor with Debye-Waller factor:
					// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
					float_tt f = seval(scatPar[0], scatPar[Znum], splinb, splinc, splind, N_SF, sqrt(s2))* exp(-s2 * B * 0.25);
					float_tt phase = PI * (kx * _config->Model.dx * _nx + ky * _config->Model.dy * _ny);
					m_atPot[Znum][ix][iy] = std::complex<float_tt>(f * cos(phase),f * sin(phase));
				}
			}
		}
#if SHOW_SINGLE_POTENTIAL == 1
		imageio = ImageIOPtr(new CImageIO(ny, nx, 0, dkx, dky, std::vector<double>(),
				"potential"));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		m_imageIO->SetThickness(m_sliceThickness);
		m_imageIO->WriteImage((void**)atPot[Znum], fileName);
#endif
#if FLOAT_PRECISION == 1
		fftwf_complex *ptr = (fftwf_complex *) (m_atPot[Znum].data());
		fftwf_plan plan = fftwf_plan_dft_2d(_nx, _ny, ptr, ptr, FFTW_BACKWARD,FFTW_ESTIMATE);
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);
#else
		fftw_complex *ptr=(fftw_complex *)(m_atPot[Znum].data());
		fftw_plan plan = fftw_plan_dft_2d(nx,ny,ptr,ptr,FFTW_BACKWARD,FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);
#endif
		for (unsigned ix = 0; ix < _nx; ix++)
			for (unsigned iy = 0; iy < _ny; iy++) {
				m_atPot[Znum][ix][iy] *= m_dkx * m_dky	* (OVERSAMPLING * OVERSAMPLING);
			}
		// make sure we don't produce negative potential:
		// if (min < 0) for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) atPot[Znum][iy+ix*ny][0] -= min;
#if SHOW_SINGLE_POTENTIAL == 1
		imageio = ImageIOPtr(new CImageIO(nx, ny, 0, m_dx/OVERSAMPLING,
				m_dy/OVERSAMPLING, std::vector<double>(), "potential"));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		//imageio->SetThickness(nz*m_sliceThickness/nzPerSlice);
		sprintf(fileName,"potential_%d.img",Znum);
		imageio->WriteImage((void**)atPot[Znum], fileName);
#endif
		BOOST_LOG_TRIVIAL(info)<< format("Created 2D %d x %d potential array for Z=%d (B=%g A^2)")
										% _nx % _ny% Znum% B;
	}

}
#define PHI_SCALE 47.87658
// #define SHOW_SINGLE_POTENTIAL 0
////////////////////////////////////////////////////////////////////////////
// This function should be used yet, because it computes the projected
// potential wrongly, since it doe not yet perform the projection!!!
complex_tt *C2DFFTPotential::GetAtomPotential2D(int Znum, double B) {

}
#undef PHI_SCALE
#undef SHOW_SINGLE_POTENTIAL

} // end namespace QSTEM
