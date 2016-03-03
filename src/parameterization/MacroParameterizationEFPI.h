
/* HEADER MacroParameterizationEFPI ===============================================================
 * author: Paul Cazeaux
 * date: 2014-01-27
 *
 * description:
 *  -> store the parameters of a collection of populations of particles modeling the plasma
 *	-> Particle loading & pushing
 *  
 * ============================================================================= */

#ifndef DEF_PLASMASCALE_MACROPARAMETERIZATIONEFPI
#define DEF_PLASMASCALE_MACROPARAMETERIZATIONEFPI

#include "parameterization/MacroParameterization.h"
#include "plasma/State.h"
#include "tools/Tools.h"
#include "tools/CurveDiagnostic.h"
 
#include <iostream>
#include <exception>
#include <stdexcept>
#include <cassert>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

#include <ctime>

/* Declarations */

class MacroParameterizationEFPI : public MacroParameterization
{
	private:
		/* class members ======================================================================== */

		int 									_grid_size;
		int 									_macro_grid_size;
		bool									_record_microsteps;

		/* Spatial quantities */

		typedef std::vector<std::vector<double> > 	spatial_quantity;
		spatial_quantity						_densities;
		spatial_quantity 						_thermal_vel;
		spatial_quantity 						_velocities;

		/* Active variables : assuming that the ion population is the first population */
		/* Data points storing the information used to determine the derivative */
		typedef std::vector<std::vector<double> >		active_variable;
		active_variable 						_stack_ion_density;
		active_variable							_stack_ion_velocity;
		active_variable 						_stack_ion_pressure;

		/* Value for the previous step, used for the leapfrog time integration */
		std::vector<double>	 					_prev_step_ion_density;
		std::vector<double>						_prev_step_ion_velocity;
		std::vector<double>	 					_prev_step_ion_pressure;

		/* Value for the current step, used for the leapfrog time integration */
		std::vector<double>  					_current_step_ion_density;
		std::vector<double> 					_current_step_ion_velocity;
		std::vector<double>  					_current_step_ion_pressure;

		/* Helpers for the quiet start using a discretization of the Maxwellian */
			// TODO = what size ?
		std::vector<double> 					_quiet_start_vel;
		std::vector<double> 					_quiet_start_icdf;

		/* Record arrays for the datapoints from the microsolver */
		std::vector<double>						_record_times;
		active_variable 						_record_ion_density;
		active_variable							_record_ion_velocity;
		active_variable 						_record_ion_pressure;

		/* Parameters for the determination of the passive variables */
		double									_electron_thermal_vel;
		double 									_debye_scaling;

	public:
		/* constructor  ========================================================================= */
		MacroParameterizationEFPI() {}
		MacroParameterizationEFPI(MacroParameterization & parameterization, double electron_thermal_vel);
		virtual ~MacroParameterizationEFPI() {}

		/* move constuctor and assignment ======================================================= */

		MacroParameterizationEFPI(MacroParameterizationEFPI &&parameterization);
		MacroParameterizationEFPI& operator=(MacroParameterizationEFPI &&parameterization);

		/* methods ============================================================================== */

		virtual void Initialize(State & state);
		virtual void Load(State & state) const;

		void RestrictAndPushback(const State & state);
		void ExtrapolateFirstHalfStep();
		void ExtrapolateSecondHalfStep();
		void Lift();

		virtual void Step(State & state);

		virtual void SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics);
		virtual void WriteData(std::fstream & fout);
		
};

#endif