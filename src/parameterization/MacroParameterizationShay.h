
/* HEADER MacroParameterization ===============================================================
 * author: Paul Cazeaux
 * date: 2014-01-27
 *
 * description:
 *  -> store the parameters of a collection of populations of particles modeling the plasma
 *	-> Particle loading & pushing
 *  
 * ============================================================================= */

#ifndef DEF_PLASMASCALE_MACROPARAMETERIZATIONSHAY
#define DEF_PLASMASCALE_MACROPARAMETERIZATIONSHAY

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

/* Declarations */

class MacroParameterizationShay : public MacroParameterization
{
	private:
		/* class members ======================================================================== */

		int 									_grid_size;
		int 									_macro_grid_size;

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

		/* Parameters for the determination of the passive variables */
		double									_electron_thermal_vel;
		double 									_debye_scaling;

	public:
		/* constructor  ========================================================================= */
		MacroParameterizationShay() {}
		MacroParameterizationShay(MacroParameterization & parameterization, double electron_thermal_vel);
		virtual ~MacroParameterizationShay() {}

		/* move constuctor and assignment ======================================================= */

		MacroParameterizationShay(MacroParameterizationShay &&parameterization);
		MacroParameterizationShay& operator=(MacroParameterizationShay &&parameterization);

		/* methods ============================================================================== */

		virtual void Initialize(const State & state);
		virtual void Load(State & state) const;

		void RestrictAndPushback(const State & state);
		void ExtrapolateFirstHalfStep();
		void ExtrapolateSecondHalfStep();
		void Lift();

		virtual void Step(State & state);

		virtual void SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics);
		
};

#endif