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

#include "wave_convergent.hpp"
#include "boost/log/trivial.hpp"
using boost::format;

namespace QSTEM
{

CConvergentWave::CConvergentWave(const ConfigPtr c) : CBaseWave(c)
{
	// TODO: where does beam current belong?
	//configReader->ReadDoseParameters(m_beamCurrent, m_dwellTime);
	const WaveConfig& w = c->Wave;
	m_Cs = w.Cs;
	m_C5= w.C5;
	m_Cc= w.Cc;
	m_df0= w.Defocus;
	m_Scherzer= "";
	m_astigMag= w.Astigmatism;
	m_a33= w.a_33;
	m_a31= w.a_31;
	m_a44= w.a_44;
	m_a42= w.a_42;
	m_a55= w.a_55;
	m_a53= w.a_53;
	m_a51= w.a_51;
	m_a66= w.a_66;
	m_a64= w.a_64;
	m_a62= w.a_62;
	m_astigAngle= w.AstigmatismAngle;
	m_phi33= w.phi_33;
	m_phi31= w.phi_31;
	m_phi44= w.phi_44;
	m_phi42= w.phi_42;
	m_phi55= w.phi_55;
	m_phi53= w.phi_53;
	m_phi51= w.phi_51;
	m_phi66= w.phi_66;
	m_phi64= w.phi_64;
	m_phi62= w.phi_62;
	m_ismoth= w.Smooth;
	m_gaussScale= w.gaussScale;
	m_gaussFlag= w.Gaussian;
	m_dE_E= w.dE_E;
	m_dI_I= w.dI_I;
	m_dV_V= w.dV_V;
	m_alpha= w.alpha;
	m_aAIS= w.AISaperture;
	m_printLevel = c->Output.LogLevel;
}

/** Copy constructor - used to copy wave just before dispatching multiple threads for STEM simulations */
CConvergentWave::CConvergentWave(const CConvergentWave& other) : CBaseWave(other)
{
	// TODO: need to copy arrays and anything pointed to - anything that needs to be thread-local
}


WavePtr CConvergentWave::Clone()
{
	return WavePtr(new CConvergentWave(*this));
}

void CConvergentWave::DisplayParams()
{
	CBaseWave::DisplayParams();

	printf("* Aperture half angle:  %g mrad\n",m_alpha);
	printf("* AIS aperture:         ");
	if (m_aAIS > 0) printf("%g A\n",m_aAIS);
	else printf("none\n");

	printf("* C_3 (C_s):            %g mm\n",m_Cs*1e-7);
	printf("* C_1 (Defocus):        %g nm%s\n",0.1*m_df0,
			// TODO: what are the available options for defocus?
			(m_Scherzer == "Scherzer") ? " (Scherzer)" : (m_Scherzer=="") ? " (opt.)":"");
	printf("* Astigmatism:          %g nm, %g deg\n",0.1*m_astigMag,RAD2DEG*m_astigAngle);

	// more aberrations:
	if (m_a33 > 0)printf("* a_3,3:                %g nm, phi=%g deg\n",m_a33*1e-1,m_phi33*RAD2DEG);
	if (m_a31 > 0)printf("* a_3,1:                %g nm, phi=%g deg\n",m_a31*1e-1,m_phi31*RAD2DEG);

	if (m_a44 > 0)printf("* a_4,4:                %g um, phi=%g deg\n",m_a44*1e-4,m_phi44*RAD2DEG);
	if (m_a42 > 0)printf("* a_4,2:                %g um, phi=%g deg\n",m_a42*1e-4,m_phi42*RAD2DEG);

	if (m_a55 > 0)printf("* a_5,5:                %g um, phi=%g deg\n",m_a55*1e-4,m_phi55*RAD2DEG);
	if (m_a53 > 0)printf("* a_5,3:                %g um, phi=%g deg\n",m_a53*1e-4,m_phi53*RAD2DEG);
	if (m_a51 > 0)printf("* a_5,1:                %g um, phi=%g deg\n",m_a51*1e-4,m_phi51*RAD2DEG);

	if (m_a66 > 0)printf("* a_6,6:                %g um, phi=%g deg\n",m_a66*1e-7,m_phi66*RAD2DEG);
	if (m_a64 > 0)printf("* a_6,4:                %g um, phi=%g deg\n",m_a64*1e-7,m_phi64*RAD2DEG);
	if (m_a62 > 0)printf("* a_6,2:                %g um, phi=%g deg\n",m_a62*1e-7,m_phi62*RAD2DEG);
	if (m_C5 != 0)printf("* C_5:                  %g mm\n",m_C5*1e-7);

	printf("* C_c:                  %g mm\n",m_Cc*1e-7);
	printf("* Damping dE/E: %g / %g \n",sqrt(m_dE_E*m_dE_E+m_dV_V*m_dV_V+m_dI_I*m_dI_I)*m_v0*1e3,m_v0*1e3);

	/*
    // TODO: where does beam current belong?
  printf("* beam current:         %g pA\n",m_beamCurrent);
  printf("* dwell time:           %g msec (%g electrons)\n",
         m_dwellTime,m_electronScale);
	 */
}

/**********************************************
 * This function creates a incident STEM probe
 * at position (dx,dy)
 * with parameters given in muls
 *
 * The following Abberation functions are being used:
 * 1) ddf = Cc*dE/E + Cc2*(dE/E)^2,
 *    Cc, Cc2 = chrom. Abber. (1st, 2nd order) [1]
 * 2) chi(qx,qy) = (2*pi/lambda)*{0.5*C1*(qx^2+qy^2)+
 *                 0.5*C12a*(qx^2-qy^2)+
 *                 C12b*qx*qy+
 *                 C21a/3*qx*(qx^2+qy^2)+
 *                 ...
 *                 +0.5*C3*(qx^2+qy^2)^2
 *                 +0.125*C5*(qx^2+qy^2)^3
 *                 ... (need to finish)
 *
 *
 *    qx = acos(kx/K), qy = acos(ky/K)
 *
 * References:
 * [1] J. Zach, M. Haider,
 *    "Correction of spherical and Chromatic Abberation
 *     in a low Voltage SEM", Optik 98 (3), 112-118 (1995)
 * [2] O.L. Krivanek, N. Delby, A.R. Lupini,
 *    "Towards sub-Angstroem Electron Beams",
 *    Ultramicroscopy 78, 1-11 (1999)
 *
 *********************************************/
#define SMOOTH_EDGE 5 // make a smooth edge on AIS aperture over +/-SMOOTH_EDGE pixels
void CConvergentWave::FormProbe()
{
	unsigned iy, ixmid, iymid;
	float_tt rmin, rmax, aimin, aimax;
	float_tt k2max, x, y, scale, pixel,alpha;
	float_tt ax = m_nx*m_dx;
	float_tt by = m_ny*m_dy;
	float_tt dx = ax-m_nx/2*m_dx;
	float_tt dy = by-m_ny/2*m_dy;
	float_tt avgRes = sqrt(0.5*(m_dx*m_dx+m_dy*m_dy));
	float_tt edge = SMOOTH_EDGE*avgRes;

	float_tt sum = 0.0;

	/********************************************************
	 * formulas from:
	 * http://cimesg1.epfl.ch/CIOL/asu94/ICT_8.html
	 *
	 * dE_E = dE/E = energy spread of emitted electrons
	 * dV_V = dV/V = acc. voltage fluctuations
	 * dI_I = dI/I = lens current fluctuations
	 * delta defocus in Angstroem (Cc in A)
	 *******************************************************/
	float_tt delta = m_Cc*m_dE_E;
	if (m_printLevel > 2) printf("defocus offset: %g nm (Cc = %g)\n",delta,m_Cc);

	/**********************************************************
	 *  Calculate misc constants
	 *********************************************************/

	float_tt rx = 1.0/ax;
	float_tt rx2 = rx * rx;
	float_tt ry = 1.0/by;
	float_tt ry2 = ry * ry;

	ixmid = m_nx/2;
	iymid = m_ny/2;

	// Start in Fourier space
	ToFourierSpace();

	/* convert convergence angle from mrad to rad */
	alpha = 0.001*m_alpha;
	k2max = sin(alpha)/m_wavlen;  /* = K0*sin(alpha) */
	k2max = k2max * k2max;

	/*   Calculate MTF
       NOTE zero freg is in the bottom left corner and
       expandes into all other corners - not in the center
       this is required for FFT

       PIXEL = diagonal width of pixel squared
       if a pixel is on the apertur boundary give it a weight
       of 1/2 otherwise 1 or 0
	 */
	pixel = ( rx2 + ry2 );
	scale = 1.0/sqrt((double)m_nx*(double)m_ny);

	//#pragma omp parallel for
	for(int iy=0; iy<m_ny; iy++) {
		float_tt ky = (float_tt) iy;
		if( iy > iymid ) ky = (double) (iy-m_ny);
		float_tt ky2 = ky*ky*ry2;
		for(int ix=0; ix<m_nx; ix++) {
			float_tt kx = (double) ix;
			if( ix > ixmid ) kx = (double) (ix-m_nx);
			float_tt k2 = kx*kx*rx2 + ky2;
			float_tt ktheta2 = k2*(m_wavlen*m_wavlen);
			float_tt ktheta = sqrt(ktheta2);
			float_tt phi = atan2(ry*ky,rx*kx);

			// defocus, astigmatism, and shift:
			float_tt chi = ktheta2*(m_df0+delta + m_astigMag*cos(2.0*(phi-m_astigAngle)))/2.0;
			ktheta2 *= ktheta;  // ktheta^3
			if ((m_a33 > 0) || (m_a31 > 0)) {
				//#pragma omp atomic
				chi += ktheta2*(m_a33*cos(3.0*(phi-m_phi33))+m_a31*cos(phi-m_phi31))/3.0;
			}
			ktheta2 *= ktheta;   // ktheta^4
			if ((m_a44 > 0) || (m_a42 > 0) || (m_Cs != 0)) {
				//#pragma omp atomic
				chi += ktheta2*(m_a44*cos(4.0*(phi-m_phi44))+m_a42*cos(2.0*(phi-m_phi42))+m_Cs)/4.0;
			}
			ktheta2 *= ktheta;    // ktheta^5
			if ((m_a55 > 0) || (m_a53 > 0) || (m_a51 > 0)) {
				//#pragma omp atomic
				chi += ktheta2*(m_a55*cos(5.0*(phi-m_phi55))+m_a53*cos(3.0*(phi-m_phi53))+m_a51*cos(phi-m_phi51))/5.0;
			}
			ktheta2 *= ktheta;    // ktheta^6
			if ((m_a66 > 0) || (m_a64 > 0) || (m_a62 = 0) || (m_C5 != 0)) {
				//#pragma omp atomic
				chi += ktheta2*(m_a66*cos(6.0*(phi-m_phi66))+m_a64*cos(4.0*(phi-m_phi64))+m_a62*cos(2.0*(phi-m_phi62))+m_C5)/6.0;
			}
			//#pragma omp atomic
			chi *= 2*M_PI/m_wavlen;
			//#pragma omp atomic
			chi -= 2.0*M_PI*( (dx*kx/ax) + (dy*ky/by) );

			if ( ( m_ismoth != 0) && ( fabs(k2-k2max) <= pixel)) {
				//#pragma omp critical
				//				{
				float x = (float) ( 0.5*scale * cos(chi));
				float y = (float) (-0.5*scale* sin(chi));
				m_wave[ix][iy] = complex_tt(x,y);
				//        m_wave[ix][iy][0]= (float) ( 0.5*scale * cos(chi));
				//        m_wave[ix][iy][1]= (float) (-0.5*scale* sin(chi));
				//				}
			}
			else if ( k2 <= k2max ) {
				//#pragma omp critical
				//				{
				float x = (float) ( scale * cos(chi));
				float y = (float) ( scale* sin(chi));
				m_wave[ix][iy] = complex_tt(x,y);

				//				}
				//        m_wave[ix][iy][0]= (float)  scale * cos(chi);
				//        m_wave[ix][iy][1]= (float) -scale * sin(chi);
			}
			else {
				//#pragma omp critical
				m_wave[ix][iy] =complex_tt(0,0);
			}

		}
	}
	/* Fourier transform into real space */
	ToRealSpace();
	/**********************************************************
	 * display cross section of probe intensity
	 */

	/* multiply with gaussian in Real Space in order to avoid artifacts */
	if (m_gaussFlag) {
#pragma omp parallel for
		for(int ix=0; ix<m_nx; ix++) {
			for(int iy=0; iy<m_ny; iy++) {
				float_tt r = exp(-((ix-m_nx/2)*(ix-m_nx/2)+(iy-m_ny/2)*(iy-m_ny/2))/(m_nx*m_nx*m_gaussScale));
#pragma omp critical
				m_wave[ix][iy] = complex_tt(m_wave[ix][iy].real()*r,m_wave[ix][iy].imag()*r);

				//        m_wave[ix][iy][0] *= (float)r;
				//        m_wave[ix][iy][1] *= (float)r;
			}
		}
	}

	/* Apply AIS aperture in Real Space */
	// printf("center: %g,%g\n",dx,dy);
	if (m_aAIS > 0) {
#pragma omp parallel for
		for(int ix=0; ix<m_nx; ix++) {
			for(int iy=0; iy<m_ny; iy++) {
				x = ix*m_dx-dx;
				y = iy*m_dy-dy;
				float_tt r = sqrt(x*x+y*y);
				delta = r-0.5*m_aAIS+edge;
				if (delta > 0) {
#pragma omp critical
					m_wave[ix][iy] = 0;
				}
				else if (delta >= -edge) {
					scale = 0.5*(1-cos(M_PI*delta/edge));
#pragma omp critical
					m_wave[ix][iy] = complex_tt(scale*m_wave[ix][iy].real(),scale*m_wave[ix][iy].imag());
					//          m_wave[ix][iy][0] = scale*m_wave[ix][iy][0];
					//          m_wave[ix][iy][1] = scale*m_wave[ix][iy][1];
				}
			}
		}
	}

	/*  Normalize probe intensity to unity  */
	for(int ix=0; ix<m_nx; ix++)
		for(int iy=0; iy<m_ny; iy++)
			sum +=  m_wave[ix][iy].real()*m_wave[ix][iy].real()+ m_wave[ix][iy].imag()*m_wave[ix][iy].imag();

	scale = 1.0 / sum;
//	scale = scale * ((double)m_nx) * ((double)m_ny);
	scale = (double) sqrt( scale );

	for(int ix=0; ix<m_nx; ix++)
		for(int iy=0; iy<m_ny; iy++) {
			m_wave[ix][iy] = complex_tt((float) scale*m_wave[ix][iy].real(),(float) scale*m_wave[ix][iy].imag());

//			if(m_wave[ix][iy].real()>1e-5)
//				DLOG(INFO) << format("m_wave[%d][%d] = %2.3f + i %2.3f") % ix % iy % m_wave[ix][iy].real() % m_wave[ix][iy].imag();
		}
	sum = 0;
	for(int ix=0; ix<m_nx; ix++)
		for(int iy=0; iy<m_ny; iy++) {
			sum += m_wave[ix][iy].real()*m_wave[ix][iy].real() + m_wave[ix][iy].imag()*m_wave[ix][iy].imag();
		}
	/*  Output results and find min and max to echo
      remember that complex pix are stored in the file in FORTRAN
      order for compatability
	 */

	rmin = m_wave[0][0].real();
	rmax = rmin;
	aimin = m_wave[0][0].imag();
	aimax = aimin;
#pragma omp parallel for
	for(int iy=0; iy<m_ny; iy++) {
		for(int ix=0; ix<m_nx; ix++) {
#pragma omp critical
			{
				if( m_wave[ix][iy].real() < rmin )
					rmin = m_wave[ix][iy].real();
			}
#pragma omp critical
			{
				if( m_wave[ix][iy].real() > rmax )
					rmax = m_wave[ix][iy].real();
			}
#pragma omp critical
			{
				if( m_wave[ix][iy].imag() < aimin )
					aimin = m_wave[ix][iy].imag();
			}
#pragma omp critical
			{
				if( m_wave[ix][iy].imag() > aimax )
					aimax = m_wave[ix][iy].imag();
			}
		}
	}
	m_rmin = rmin;
	m_rmax = rmax;
	m_aimin = aimin;
	m_aimax = aimax;

	BOOST_LOG_TRIVIAL(trace) << format("wave value range (%f .. %f,i %f ... %f)") % rmin % rmax % aimin % aimax;

	/**********************************************************/

}  /* end probe() */

} // end namespace QSTEM
