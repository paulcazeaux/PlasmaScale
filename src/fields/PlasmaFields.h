/* HEADER PlasmaFields ===============================================================
 * author: Paul Cazeaux
 * date: 2013-11-28
 *
 * description:
 *  -> store & compute the fields (potential, energy, charge)
 *  
 * ============================================================================= */

#ifndef DEF_PLASMASCALE_PLASMAFIELDS
#define DEF_PLASMASCALE_PLASMAFIELDS

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "plasma/Plasma.h"

/* Declarations */

class PlasmaFields
{
	friend class State;
	friend class PopulationOfParticles;

	private:
		/* class members ======================================================================== */
		/* pointer on object */
		std::shared_ptr<const Plasma>				_plasma;

		/* Simulation parameters */
		const int 									_grid_end;

		/* Time */
		std::shared_ptr<double>						_simulation_time;
		std::shared_ptr<int>						_iteration;

		/* Fields */
		std::vector<double> 						_charge;
		std::vector<double> 						_electrical_field;
		std::vector<double> 						_potential;

		/* Computed energy */
		double 										_electrostatic_energy;

	public:
		/* constuctor and destructor ============================================================ */
		PlasmaFields(std::shared_ptr<const Plasma> plasma,
					std::shared_ptr<double> simulation_time,
					std::shared_ptr<int> iteration);
		~PlasmaFields() {}

		/* getter =========================================================================*/
		double * get_density_ptr()						{return _charge.data();				}
		double * get_electrical_field_ptr()				{return _electrical_field.data();	}
		double * get_potential_ptr()					{return _potential.data();		}

		/* methods ======================================================================= */
		void PushBackElectrostaticEnergy(std::vector<double>& electrostatic_energy) const;
		void Compute();
		void ResetCharge()
		{
			for (double & rho: _charge) rho = 0.;
		}
};

#endif