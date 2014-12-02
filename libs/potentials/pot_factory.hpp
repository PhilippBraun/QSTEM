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

#ifndef POT_FACTORY_H
#define POT_FACTORY_H

#include <map>
#include "pot_interface.hpp"
#include "config_IO/config_reader_factory.hpp"

namespace QSTEM
{

typedef std::map<std::string, CreatePotentialFn> FactoryMap;
// Factory for creating instances of IPotential
class QSTEM_HELPER_DLL_EXPORT CPotFactory
{
protected:
	CPotFactory();
	CPotFactory(const CPotFactory &) { }
	CPotFactory &operator=(const CPotFactory &) { return *this; }
	FactoryMap m_FactoryMap;
	PotPtr GetPotential(const std::string &animalName,const ConfigPtr configReader);
	PotPtr GetPotential(bool _3D, bool fft,const ConfigPtr configReader);

public:
	~CPotFactory() { m_FactoryMap.clear(); }
	static CPotFactory *Get()
	{
		static CPotFactory instance;
		return &instance;
	}
	void Register(const std::string &potentialName, CreatePotentialFn pfnCreate);
	PotPtr GetPotential(const ConfigPtr c);
};

}

#endif
