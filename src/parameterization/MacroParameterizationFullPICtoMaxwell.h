
/* HEADER MacroParameterizationFullPICtoMaxwell ===============================================================
 * author: Paul Cazeaux
 * date: 2014-01-27
 *
 * description:
 *  -> store the parameters of a collection of populations of particles modeling the plasma
 *	-> Particle loading & pushing
 *  
 * ============================================================================= */

#ifndef DEF_PLASMASCALE_MACROPARAMETERIZATIONFULLPICTOMAXWELL
#define DEF_PLASMASCALE_MACROPARAMETERIZATIONFULLPICTOMAXWELL

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

class MacroParameterizationFullPICtoMaxwell : public MacroParameterization
{
	private:
		/* class members ======================================================================== */

		int 						_grid_size;

		/* Spatial quantities */

		std::vector<double>			_ion_density;
		std::vector<double> 		_ion_pressure;
		std::vector<double> 		_ion_thermal_velocity;
		std::vector<double> 		_ion_velocity;

	public:
		/* constructor  ========================================================================= */
		MacroParameterizationFullPICtoMaxwell() {}
		MacroParameterizationFullPICtoMaxwell(MacroParameterization & parameterization);
		virtual ~MacroParameterizationFullPICtoMaxwell() {}

		/* move constuctor and assignment ======================================================= */

		MacroParameterizationFullPICtoMaxwell(MacroParameterizationFullPICtoMaxwell &&parameterization);
		MacroParameterizationFullPICtoMaxwell& operator=(MacroParameterizationFullPICtoMaxwell &&parameterization);

		/* methods ============================================================================== */

		virtual void Initialize(State & state);
		virtual void Load(State & state) const {}

		void ComputeVariables(const State & state);
		virtual void Step(State & state);

		virtual void SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics);
		virtual void WriteData(std::fstream & fout);
};

#endif