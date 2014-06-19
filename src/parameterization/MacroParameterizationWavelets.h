
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

typedef WaveletRepresentationP1 ActiveWaveletRepresentation;
typedef MaxwellianRepresentationP1 ActiveMaxwellianRepresentation;

class MacroParameterizationWavelets : public MacroParameterization
{
	private:
		/* class members ======================================================================== */

		int 										_grid_size;
		int 										_macro_grid_size;
		int 										_depth;
		int 										_cutoff;
		double 										_ion_vmax;
		bool										_record_microsteps;

		/* Fine grid wavelet coefficients by population */
		std::vector<std::unique_ptr<Representation> > 	_distributions;

		/* Active variables : assuming that the ion population is the first population */
		/* Data points storing the information used to determine the derivative */
		std::vector<ActiveWaveletRepresentation> 	_stack_ion_distribution;
		std::vector<ActiveMaxwellianRepresentation> _stack_electron_distribution;
		int 										_stack_index;

		/* Value for the previous step, used for the leapfrog time integration */
		ActiveWaveletRepresentation					_prev_step_ion_distribution;
		ActiveMaxwellianRepresentation				_prev_step_electron_distribution;

		/* Value for the current step, used for the leapfrog time integration */
		ActiveWaveletRepresentation					_current_step_ion_distribution;
		ActiveMaxwellianRepresentation				_current_step_electron_distribution;

		/* Record arrays for the datapoints from the microsolver */
		std::vector<double>							_record_times;
		std::vector<ActiveWaveletRepresentation> 	_record_ion_distribution;
		std::vector<ActiveMaxwellianRepresentation> _record_electron_distribution;

		/* Arrays for the diagnostics */
		std::vector<double> 						_ion_density;
		std::vector<double> 						_ion_velocity;
		std::vector<double> 						_ion_pressure;
		std::vector<double> 						_electron_density;
		std::vector<double> 						_electron_velocity;
		std::vector<double> 						_electron_pressure;


	public:
		/* constructor  ========================================================================= */
		MacroParameterizationWavelets() {}
		MacroParameterizationWavelets(MacroParameterization & parameterization, double ion_vmax);
		virtual ~MacroParameterizationWavelets() {}

		/* move constuctor and assignment ======================================================= */

		MacroParameterizationWavelets(MacroParameterizationWavelets &&parameterization);
		MacroParameterizationWavelets& operator=(MacroParameterizationWavelets &&parameterization);

		/* methods ============================================================================== */

		virtual void Initialize(State & state);
		virtual void Load(State & state) const;

		void RestrictAndPushback(const State & state, const double delay);
		void Extrapolate(const double ratio);
		void Lift();

		virtual void Step(State & state);

		virtual void SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics);
		virtual void WriteData(std::fstream & fout);
		
};

#endif