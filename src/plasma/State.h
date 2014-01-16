/* HEADER  State ====================================================================================
 * author: Paul Cazeaux
 * date: 2013-11-29
 * 
 * description:
 *  -> store the state of the simulation
 * 
 * ============================================================================================== */

#ifndef DEF_PLASMASCALE_STATE
#define DEF_PLASMASCALE_STATE

/* includes ===================================================================================== */
#include <iostream>
#include <vector>
#include <memory>

#include "plasma/Plasma.h"
#include "particles/CollectionOfPopulations.h"
#include "fields/PlasmaFields.h"
#include "tools/CurveDiagnostic.h"
#include "tools/ScatterDiagnostic.h"

class State
{
	friend class History;

	private:
		/* class members ======================================================================== */

		/* pointer on object */
		std::shared_ptr<const Plasma>				_plasma;
		std::unique_ptr<PlasmaFields>				_fields;
		std::unique_ptr<CollectionOfPopulations>	_populations;

		/* members */
		std::shared_ptr<double>						_simulation_time;
		std::shared_ptr<int>						_iteration;

		/* singleton ============================================================================ */
		State(const State&);
		State& operator=(const State&);
		State& operator=(State);

	public:
		/* constuctor and destructor ============================================================ */
		State(	std::shared_ptr<const Plasma> plasma );
		State(	std::shared_ptr<const Plasma> plasma, FILE *& InputDeck );
		State(	std::shared_ptr<const Plasma> plasma, std::unique_ptr<CollectionOfPopulations> populations );
		~State() {}

		/* getter */

		std::shared_ptr<double> get_simulation_time()	const { return _simulation_time;	}

		/* operator ============================================================================= */
		friend std::ostream& operator<<( std::ostream& os, const State& state);

		/* method =============================================================================== */
		void Step();
		void SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic>	> &diagnostics);
};


 #endif