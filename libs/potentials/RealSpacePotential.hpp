/*
 * RealSpacePotential.hpp
 *
 *  Created on: Dec 2, 2014
 *      Author: philipp
 */

#ifndef REALSPACEPOTENTIAL_HPP_
#define REALSPACEPOTENTIAL_HPP_

#include "pot_base.hpp"
#include <vector>


namespace QSTEM {


class RealSpacePotential : public CPotential
{
public :
	inline RealSpacePotential(const ConfigPtr c) : CPotential(c){};
	virtual ~RealSpacePotential();
protected:
	void AddAtomRealSpace(std::vector<atom>::iterator &atom, float_tt atomX, float_tt atomY, float_tt atomZ);
	virtual void _AddAtomRealSpace(std::vector<atom>::iterator &atom,float_tt atomX, unsigned int ix,
			float_tt atomY, unsigned int iy,float_tt atomZ, unsigned int iatomZ)=0;
};

}
#endif /* REALSPACEPOTENTIAL_HPP_ */
