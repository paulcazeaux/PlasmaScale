
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
		std::vector<int>			_modes;
		std::vector<double>			_density_amplitudes;
		std::vector<double>			_density_phases;
		std::vector<double>			_velocity_amplitudes;
		std::vector<double>			_velocity_phases;

			// Diagnostics
		std::vector<int>			_bin_numbers;
		std::vector<double>			_upper_velocities;
		std::vector<double>			_lower_velocities;

	public:

		MacroParameterizationFromFile(FILE *& InputDeck);
		
		virtual void 	Load(State & state) const;

		double 			get_initial_thermal_vel(int population_index)	const;

		virtual bool 	HaveVelocityDiagnostics()		const	{return true;	}
		virtual double 	GetBinStart(int index)			const;
		virtual	double	GetBinWidth(int index)			const;
		virtual int		GetNumberOfBins(int index)		const;
};

#endif