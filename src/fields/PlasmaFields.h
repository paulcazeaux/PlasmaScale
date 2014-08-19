/* HEADER PlasmaFields ===============================================================
 * author: Paul Cazeaux
 * date: 2013-11-28
 *
 * description:
 *  -> store & compute the fields (potential, energy, charge)
 *	-> Poisson solver
 *  
 * ============================================================================= */

#ifndef DEF_PLASMASCALE_PLASMAFIELDS
#define DEF_PLASMASCALE_PLASMAFIELDS

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "plasma/Plasma.h"
#include "fields/Field.h"

#define AKPERP               2.405

/* Declarations */

class PlasmaFields
{
	private:
		/* class members ======================================================================== */
		/* pointer on object */
		std::shared_ptr<const Plasma>				_plasma;

		/* Physical parameters */
		const double								_length;
		const double								_epsilon;
		const double								_E0;
		const double								_w0;

		/* Simulation parameters */
		const int 									_grid_size;
		const double								_dx;

		/* Modal filter & Fourier Poisson kernel */
		std::unique_ptr<const std::vector<double> > _filtered_kernel;
		std::unique_ptr<const std::vector<double> > _filter;
		const double								_filter_parameter_1;
		const double								_filter_parameter_2;

		/* Time */
		std::shared_ptr<double>						_simulation_time;
		std::shared_ptr<int>						_iteration;

		/* Fields */
		Field 										_charge;
		Field 										_electrical_field;
		Field 										_potential;

		/* Computed energy */
		const int 									_max_mode;
		double 										_electrostatic_energy_total;
		std::vector<double>							_electrostatic_energy_by_mode;

	public:
		/* constuctor and destructor ============================================================ */
		PlasmaFields(std::shared_ptr<const Plasma> plasma,
					std::shared_ptr<double> simulation_time,
					std::shared_ptr<int> iteration);
		~PlasmaFields() {}

		/* getter =========================================================================*/
		double * get_density_ptr()						const 	{return _charge.get_ptr();				}
		double * get_electrical_field_ptr()				const	{return _electrical_field.get_ptr();	}
		double * get_potential_ptr()					const 	{return _potential.get_ptr();		}

		/* methods ======================================================================= */
		void PushBackElectrostaticEnergy(std::vector<double>& electrostatic_energy, std::vector<std::vector<double>	>& electrostatic_energy_by_mode) const;
		void ComputeAndFilter();
		void GetEField(std::vector<double> & Efield);
		void ResetCharge()
		{
			_charge.Reset();
		}
		void SubstractMeanDensity(const double rho0)
		{
			for (int i=0; i<_grid_size; i++)
			{
				_charge._values[i] -= rho0;
			}
		}

		void WeighParticle(double position, double charge);
		
};

#endif