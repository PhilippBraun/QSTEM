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

#include "wave_plane.hpp"

static std::string imageFilePrefix="image";
using boost::format;

namespace QSTEM
{

CPlaneWave::CPlaneWave(const ConfigPtr c) : CBaseWave(c)
{
	m_nx = c->Model.nx;
	m_ny = c->Model.ny;
}

CPlaneWave::CPlaneWave() : CBaseWave()
{
}

/** Copy constructor - used to copy wave just before dispatching multiple threads for STEM simulations */
CPlaneWave::CPlaneWave(const CPlaneWave& other) : CBaseWave(other)
{
	// TODO: need to copy arrays and anything pointed to - anything that needs to be thread-local
}

WavePtr CPlaneWave::Clone()
{
	return WavePtr(new CPlaneWave(*this));
}

void CPlaneWave::DisplayParams()
{
	CBaseWave::DisplayParams();

	// TODO: transmit tiltBack status somehow (it's an experiment parameter, not a wave parameter.
	BOOST_LOG_TRIVIAL(info) << format("* Beam tilt:            x=%g deg, y=%g deg")
			% (_config->Wave.tiltX*RAD2DEG) % (_config->Wave.tiltY *RAD2DEG);
}

void CPlaneWave::FormProbe()
{
	float_tt scale = 1/sqrt((float_tt)m_nx*m_ny);
	if ((_config->Wave.tiltX == 0) && (_config->Wave.tiltY == 0)) {
		for (unsigned i=0; i<m_nx; i++)
			for (unsigned j=0; j<m_ny; j++)
			{
				m_wave[i][j]=std::complex<float_tt>(scale,0);

			}
	}
	else {
		TiltBeam();
	}
}

void CPlaneWave::TiltBeam(bool tiltBack)
{
	if ((_config->Wave.tiltX != 0) || (_config->Wave.tiltY != 0))
	{
		float_tt scale = 1/sqrt(m_nx*m_ny);
		int direction = tiltBack ? -1 : 1;

		// produce a tilted wave function (btiltx,btilty):
		float_tt ktx = direction*2.0*M_PI*sin(_config->Wave.tiltX)/GetWavelength();
		float_tt kty = direction*2.0*M_PI*sin(_config->Wave.tiltY)/GetWavelength();
		unsigned px=m_nx*m_ny;
		for (unsigned i=0; i<m_nx; i++)
			for (unsigned j=0; j<m_ny; j++)
			{
				float_tt x = m_dx*(i-m_nx/2);
				float_tt y = m_dy*(j-m_ny/2);
				m_wave[i][j]=std::complex<float_tt>((float_tt)cos(ktx*x+kty*y)*scale,(float_tt)sin(ktx*x+kty*y)*scale);
			}
	}
}

void CPlaneWave::TiltBack()
{
	TiltBeam(true);
}

void CPlaneWave::ReadImage()
{
	m_position.clear();
	m_imageIO->ReadImage(m_image, imageFilePrefix, m_position);
}

void CPlaneWave::WriteImage()
{
	std::map<std::string, double> params;
	std::string comment;
	m_imageIO->WriteImage(m_image, imageFilePrefix, params, comment, m_position);
}

} // end namespace QSTEM
