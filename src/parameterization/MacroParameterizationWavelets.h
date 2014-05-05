
/* HEADER MacroParameterizationWavelets ===============================================================
 * author: Paul Cazeaux
 * date: 2014-01-27
 *
 * description:
 *  -> store the parameters of a collection of populations of particles modeling the plasma
 *	-> Particle loading & pushing
 *  
 * ============================================================================= */

#ifndef DEF_PLASMASCALE_MACROPARAMETERIZATIONWAVELETS
#define DEF_PLASMASCALE_MACROPARAMETERIZATIONWAVELETS

#include "parameterization/MacroParameterization.h"
#include "plasma/State.h"
#include "tools/Tools.h"
#include "tools/CurveDiagnostic.h"
#include "representation/WaveletRepresentation.h"
#include "representation/WaveletRepresentationP1.h"
#include "representation/MaxwellianRepresentation.h"
#include "representation/MaxwellianRepresentationP1.h"
 
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

class MacroParameterizationWavelets : public MacroParameterization
{
	private:
		/* class members ======================================================================== */

		int 										_grid_size;
		int 										_macro_grid_size;
		int 										_depth;
		double 										_ion_vmax;
		bool										_record_microsteps;

		/* Fine grid wavelet coefficients by population */
		std::vector<std::unique_ptr<Representation> > 	_distributions;

		/* Active variables : assuming that the ion population is the first population */
		/* Data points storing the information used to determine the derivative */
		typedef WaveletRepresentation				active_variable;
		std::vector<active_variable> 				_stack_ion_distribution;

		/* Value for the previous step, used for the leapfrog time integration */
		active_variable	 							_prev_step_ion_distribution;

		/* Value for the current step, used for the leapfrog time integration */
		active_variable								_current_step_ion_distribution;

		/* Record arrays for the datapoints from the microsolver */
		std::vector<double>							_record_times;
		std::vector<active_variable> 				_record_ion_distribution;

		/* Parameters for the determination of the passive variables */
		double										_electron_thermal_vel;
		double 										_debye_scaling;

		/* Arrays for the diagnostics */
		std::vector<double> 						_ion_density;
		std::vector<double> 						_ion_velocity;
		std::vector<double> 						_ion_pressure;


	public:
		/* constructor  ========================================================================= */
		MacroParameterizationWavelets() {}
		MacroParameterizationWavelets(MacroParameterization & parameterization, double electron_thermal_vel, double ion_vmax, int depth);
		virtual ~MacroParameterizationWavelets() {}

		/* move constuctor and assignment ======================================================= */

		MacroParameterizationWavelets(MacroParameterizationWavelets &&parameterization);
		MacroParameterizationWavelets& operator=(MacroParameterizationWavelets &&parameterization);

		/* methods ============================================================================== */

		virtual void Initialize(const State & state);
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