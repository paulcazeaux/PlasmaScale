/* HEADER  History ==============================================================================
 * author: Paul Cazeaux
 * date: 2013-11-29
 * 
 * description:
 *  -> store the history of the simulation (energy, moment...)
 * 
 * ============================================================================================== */

#ifndef DEF_PLASMASCALE_HISTORY
#define DEF_PLASMASCALE_HISTORY

/* includes ===================================================================================== */
#include <iostream>
#include <vector>
#include <cassert>

#include "plasma/State.h"
#include "plasma/Plasma.h"
#include "fields/PlasmaFields.h"
#include "tools/CurveDiagnostic.h"

class History
{
	private:
		/* class members ======================================================================== */
		/* pointer on object */
		std::shared_ptr<const Plasma>			_plasma;

		/* size */
		const int 								_max_size;
		int 									_size;

		/* time */
		std::vector<double>						_time_array;
		
		const int 								_velocity_accumulation_interval;

		/* accumulation interval */

		int 									_interval;

		/* aggregated energy */
		std::vector<double>						_electrostatic_energy;
		std::vector<double> 					_kinetic_energy;
		std::vector<double> 					_total_energy;
		std::vector<double>						_moment;

		/* electrostatic energy by mode */
		const int 								_max_mode;
		std::vector<std::vector<double>	>		_electrostatic_energy_by_mode;

		/* moment and energy by population */
		const int 								_number_of_populations;
		std::vector<std::vector<double>	>		_kinetic_energy_by_population;
		std::vector<std::vector<double>	>		_moment_by_population;

		/* singleton ============================================================================ */
		History(const History&);
		History& operator=(const History&);
		History& operator=(History);

	public:
		/* constuctor and destructor ============================================================ */
		History(std::shared_ptr<const Plasma> plasma);
		~History() {}

		/* methods ============================================================================== */
		void	Resize(int size);
		void	Compute(const State& state);
		void	Compute(const MacroState& state);
		void	Comb();
		void 	SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics);
};


 #endif