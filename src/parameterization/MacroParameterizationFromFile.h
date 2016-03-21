
/* HEADER MacroParameterizationFromFile ===============================================================
 * author: Paul Cazeaux
 * date: 2014-01-27
 *
 * description:
 *  -> store the parameters of a collection of populations of particles modeling the plasma
 *	-> Particle loading & pushing
 *  
 * ============================================================================= */

#ifndef DEF_PLASMASCALE_MACROPARAMETERIZATIONFROMFILE
#define DEF_PLASMASCALE_MACROPARAMETERIZATIONFROMFILE

#include "parameterization/MacroParameterization.h"
#include "tools/RandomTools.h"
#include "plasma/State.h"
#include <cmath>
#include <cassert>
#include <iostream>
 
#include <eigen3/Eigen/Sparse>

/* Declarations */

class MacroParameterizationFromFile : public MacroParameterization
{
	private:
		/* class members ======================================================================== */

			// Initialization
		std::vector<int>			_group_sizes;
		std::vector<double>			_mean_velocities;
		std::vector<double>			_quiet_start_exponents;
		std::vector<double>			_quiet_mean_thermal_vel;
		std::vector<double>			_random_mean_thermal_vel;

			// Perturbation
		double						_init_occupation;

	public:

		MacroParameterizationFromFile(FILE *& InputDeck);
		
		virtual void 	Load(State & state) const;

		virtual double 	get_initial_thermal_vel(int population_index)	const
			{	return _quiet_mean_thermal_vel.at(population_index) + _random_mean_thermal_vel.at(population_index);	}
};

#endif