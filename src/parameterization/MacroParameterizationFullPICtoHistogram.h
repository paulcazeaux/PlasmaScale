
/* HEADER MacroParameterizationFullPICtoHistogram ===============================================================
 * author: Paul Cazeaux
 * date: 2014-01-27
 *
 * description:
 *  -> store the parameters of a collection of populations of particles modeling the plasma
 *	-> Particle loading & pushing
 *  
 * ============================================================================= */

#ifndef DEF_PLASMASCALE_MACROPARAMETERIZATIONFULLPICTOHISTOGRAM
#define DEF_PLASMASCALE_MACROPARAMETERIZATIONFULLPICTOHISTOGRAM

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

class MacroParameterizationFullPICtoHistogram : public MacroParameterization
{
	private:
		/* class members ======================================================================== */

		int 										_grid_size;
		int 										_number_of_bins;
		int 										_min_depth;
		int 										_max_depth;
		double 										_ion_vmax;

		/* Spatial quantities */

		std::vector<std::vector<double>	>			_histogram;

	public:
		/* constructor  ========================================================================= */
		MacroParameterizationFullPICtoHistogram() {}
		MacroParameterizationFullPICtoHistogram(MacroParameterization & parameterization, double ion_vmax);
		virtual ~MacroParameterizationFullPICtoHistogram() {}

		/* move constuctor and assignment ======================================================= */

		MacroParameterizationFullPICtoHistogram(MacroParameterizationFullPICtoHistogram &&parameterization);
		MacroParameterizationFullPICtoHistogram& operator=(MacroParameterizationFullPICtoHistogram &&parameterization);

		/* methods ============================================================================== */

		virtual void Initialize(State & state);
		virtual void Load(State & state) const {}

		void ComputeVariables(const State & state);
		void Coarsen();
		virtual void Step(State & state);

		virtual void SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics);
		virtual void WriteData(std::fstream & fout);
};

#endif