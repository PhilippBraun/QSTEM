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

#include "experiments.hpp"
#include "experiments/stem.hpp"
#include "experiments/cbed.hpp"
#include "experiments/tem.hpp"

namespace QSTEM
{

ExperimentPtr GetExperiment(QSTEM::Config &config)
{
	switch (config.ExperimentType) {
	case ExperimentType::CBED:
		return ExperimentPtr(new CExperimentCBED(config));
		break;
	case ExperimentType::NBED:
		printf("Unrecognized experiment type: NBED.  Exiting.");
		exit(-1);
		return nullptr;
		break;
	case ExperimentType::STEM:
		return ExperimentPtr(new CExperimentSTEM(config));
		break;
	case ExperimentType::TEM:
		return ExperimentPtr(new CExperimentTEM(config));
		break;
	default:
		printf("Unrecognized experiment type: NONE.  Exiting.");
		exit(-1);
		break;
	}
}

}
