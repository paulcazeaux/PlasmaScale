/* HEADER PopulationOfParticles ===============================================================
 * author: Paul Cazeaux
 * date: 2013-11-28
 *
 * description:
 *  -> store the parameters of a population of particles (mass, charge...)
 *  -> store the variables of a population of particles (positions, velocities, weights)
 *	-> Particle loading & pushing
 *  
 * ============================================================================= */

#ifndef DEF_PLASMASCALE_POPULATIONOFPARTICLES
#define DEF_PLASMASCALE_POPULATIONOFPARTICLES

#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <memory>


#include "fields/PlasmaFields.h"
#include "plasma/Plasma.h"
#include "parameterization/MacroParameterization.h"
#include "tools/Tools.h"
#include "tools/RandomTools.h"

/* Declarations */

class PopulationOfParticles
{
	friend class State;
	friend class PlasmaFields;

	private:
		/* class members ======================================================================== */
		/* pointer on object ============================================================= */
		std::shared_ptr<const Plasma>			_plasma;
		std::shared_ptr<int>					_iteration;

		/* Parameters ==================================================================== */
		int 				 					_population_size;
		const double							_unit_mass;
		double									_total_mass;
		double 									_total_weight;
		double 									_mean_density;
		const double							_unit_charge;


		/* Variables ===================================================================== */
		std::vector<double>						_position;
		std::vector<double>						_velocity;
		std::vector<double>						_weights;
		
		/* Energy ======================================================================== */
		double 									_moment;
		double 									_kinetic_energy;


		/* Velocity bins     ============================================================= */

		int 			 						_number_of_bins;
		double									_bin_width;
		double									_bin_start;

		/* velocity profiles ============================================================ */
		bool									_profiling_active;
		int 									_velocity_accumulation_interval;
		double 									_accumulation_weight;
		int 									_count;

		/* Pointers required for diagnostics =========================================== */
		std::unique_ptr<std::vector<double> >	_mid_bin_array;
		std::unique_ptr<std::vector<double> >	_partial_velocity_profile;
		std::unique_ptr<std::vector<double> >	_accumulated_velocity_profile;

	public:
		/* constuctor and destructor =================================================== */
		PopulationOfParticles(std::shared_ptr<const Plasma> plasma,
			std::shared_ptr<int> iteration,
			double population_size, double unit_mass, double unit_charge);
		PopulationOfParticles(const MacroParameterization & parameterization, const int index, std::shared_ptr<int> iteration);

		~PopulationOfParticles() {}

		/* getters ===============================================================================*/
		double	get_unit_charge()							const { return _unit_charge;}
		double	get_unit_mass()								const { return _unit_mass;	}
		int 	get_size()									const { return _population_size;	}

		double	get_kinetic_energy()						const { return _kinetic_energy;	}
		double	get_moment()								const { return _moment;	}

		bool	is_velocity_profiling_active()				const { return _profiling_active;	}
		int    	get_number_of_bins()						const { return _number_of_bins;	}

		/* setters ===============================================================================*/
		void 	set_new_number_of_particles(const int np);

		/* methods ============================================================================== */
		void 	Reset();
		bool	CheckParameters(const MacroParameterization & parameterization,const int index);
		void 	ComputeAggregateParameters();
		void 	Accelerate(const PlasmaFields &fields, double factor = 1.0);
		void	Move();
		void	Weigh(PlasmaFields	&fields);

		void 	Prepare(const PlasmaFields &fields, const bool toggle_half_step = true);

		void 	SetupVelocityDiagnostics(int nbins, int velocity_accumulation_interval, double vupper, double vlower, double v0, double vt1, double vt2);
		void 	SetupVelocityDiagnostics(const MacroParameterization & parameterization, const int index);
		void	ComputeVelocityProfile();
		void 	ComputeVelocityProfile(std::vector<double> & velocity_profile);
};

#endif

