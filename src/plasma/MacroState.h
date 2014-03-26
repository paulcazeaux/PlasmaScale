/* HEADER  MacroState ====================================================================================
 * author: Paul Cazeaux
 * date: 2014-01-27
 * 
 * description:
 *  -> store the macroscopic state of the simulation
 * 
 * ============================================================================================== */


#ifndef DEF_PLASMASCALE_MACROSTATE
#define DEF_PLASMASCALE_MACROSTATE

/* includes ===================================================================================== */
#include <iostream>
#include <vector>
#include <memory>
#include <exception>
#include <iostream>
#include <stdexcept>

#include "plasma/Plasma.h"
#include "plasma/State.h"
#include "parameterization/MacroParameterization.h"
#include "parameterization/MacroParameterizationFromFile.h"
#include "parameterization/MacroParameterizationEFPI.h"
 #include "parameterization/MacroParameterizationFullPIC.h"

class MacroState
{
	friend class History;

	private:
		/* class members ======================================================================== */

		/* pointer on object */
		std::shared_ptr<const Plasma>				_plasma;
		std::unique_ptr<MacroParameterization>		_parameterization;
		std::unique_ptr<State>						_micro_state;

		/* members */
		std::shared_ptr<double>						_simulation_time;
		std::shared_ptr<int>						_macro_iteration;

		double 										_macro_dt;

		/* singleton ============================================================================ */
		MacroState(const MacroState&);
		MacroState& operator=(const MacroState&);
		MacroState& operator=(MacroState);

	public:
		/* constuctor and destructor ============================================================ */
		MacroState(	FILE *& InputDeck );
		~MacroState() {}

		/* getter */

		std::shared_ptr<double> 		get_simulation_time()	const 	{ return _simulation_time;	}
		std::shared_ptr<const Plasma> 	get_plasma() 			const	{ return _plasma;	}

		/* operators ============================================================================= */
		friend std::ostream& operator<<( std::ostream& os, const MacroState& state);

		/* methods =============================================================================== */
		void Step();
		void SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic>	> &diagnostics)
		{
			_micro_state->SetupDiagnostics(diagnostics);
			_parameterization->SetupDiagnostics(diagnostics);
		}
		void WriteData(std::fstream & fout);
};


 #endif