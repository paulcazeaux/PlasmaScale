/* HEADER Simulation ===============================================================
 * author: Paul Cazeaux
 * date: 2013-11-28
 *
 * description:
 *  -> transparent container for the simulation: plasma library, state, and history
 *  -> User-friendly interface
 *  
 * ============================================================================= */

#ifndef DEF_PLASMASCALE_SIMULATION
#define DEF_PLASMASCALE_SIMULATION

#include <memory>

#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "plasma/Plasma.h"
#include "plasma/State.h"
#include "plasma/History.h"

class Simulation
{
	private:
		/* class members ======================================================= */
		bool							_is_initialized;

	public:

		/* pointer on object */

		std::shared_ptr<const Plasma> 	_plasma;
		std::unique_ptr<State> 			_state;
		std::unique_ptr<History>		_history;

		std::vector<std::unique_ptr<Diagnostic>	>	_diagnostics;


		/* Constructor, destructor ============================================= */
		Simulation() {_is_initialized = false;}
		~Simulation() {}

		/* test    ============================================================= */
		bool is_initialized()			const { return _is_initialized;	}

		/* methods ============================================================= */
		void Setup(int argc, char ** argv);
		void Step();
		void InitWindows();
		double * get_simulation_time_ptr()		const { return _state->get_simulation_time().get();	}
};

 #endif