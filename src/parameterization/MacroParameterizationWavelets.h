
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
#include <iomanip>
#include <exception>
#include <stdexcept>
#include <cassert>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
 
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/UmfPackSupport>

#include <ctime>
#include <chrono>

/* Declarations */

typedef std::chrono::milliseconds timeunit;
typedef WaveletRepresentationP1 ActiveWaveletRepresentation;
typedef MaxwellianRepresentationP1 ActiveMaxwellianRepresentation;

class MacroParameterizationWavelets : public MacroParameterization
{
	private:
		/* class members ======================================================================== */

		int 														_grid_end;
		int 														_macro_grid_end;
		int 														_depth;
		int 														_cutoff;
		bool														_record_microsteps;

		/* Fine grid wavelet coefficients by population */
		std::vector<std::unique_ptr<Representation> > 				_distributions;

		/* Active variables : assuming that the ion population is the first population */
		/* Data points storing the information used to determine the derivative */
		std::vector<ActiveWaveletRepresentation> 					_stack_ion_distribution;
		std::vector<double> 											_stack_ion_number;
		std::vector<double> 											_stack_electron_number;
		std::vector<double> 										_stack_electron_energy;
		int 														_stack_index;

		/* Record arrays for the datapoints from the microsolver */
		std::vector<double>											_record_times;
		std::vector<ActiveWaveletRepresentation> 					_record_ion_distribution;

		/* Parameters for the determination of the passive variables */
		double														_total_moment;
		Eigen::SparseMatrix<double> 								_J;

		/* Arrays for the diagnostics */
		std::vector<double> 										_ion_density;
		std::vector<double> 										_ion_velocity;
		std::vector<double> 										_ion_pressure;

		/* Vector to store an acceleration profile for the approximation of characteristics */
		std::vector<double> 										_accfield;

	public:
		/* constructor  ========================================================================= */
		MacroParameterizationWavelets() {}
		MacroParameterizationWavelets(MacroParameterization & parameterization);
		virtual ~MacroParameterizationWavelets() {}

		/* move constuctor and assignment ======================================================= */

		MacroParameterizationWavelets(MacroParameterizationWavelets &&parameterization);
		MacroParameterizationWavelets& operator=(MacroParameterizationWavelets &&parameterization);

		/* methods ============================================================================== */

		virtual void Initialize(State & state);
		virtual void Load(State & state) const;
		virtual void SetAccField(State & state);
		void RestrictAndPushback(const State & state, const double delay);
		void Lift();

		virtual void Step(State & state);

		virtual void SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics);
		virtual void WriteData(State & state, std::fstream & fout);

		void CalculateTotalMoment(const State & state);
		void CalculateIonMoment(const State & state, double & ion_particle_moment, double & ion_distr_moment);
		void CalculateElectronMoment(const State & state, double & electron_particle_moment);
};

#endif