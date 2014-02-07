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

#include "fields/PlasmaFields.h"
#include "plasma/Plasma.h"
#include "parameterization/MacroParameterization.h"
#include "tools/RandomTools.h"

/* Declarations */

class PopulationOfParticles
{
	friend class State;

	private:
		/* class members ======================================================================== */
		/* pointer on object */
		std::shared_ptr<const Plasma>			_plasma;
		std::shared_ptr<int>					_iteration;

		/* Parameters */
		std::unique_ptr<int> 					_population_size;
		const double							_unit_mass;
		double									_total_mass;
		double 									_total_weight;
		double 									_mean_density;
		const double							_unit_charge;
		const double							_cyclotronic_rotation_parameter;
		const bool 								_magnetized;


		/* Variables */
		std::vector<double>						_position;
		std::vector<double>						_velocity_x;
		std::vector<double>						_velocity_y;
		std::vector<double>						_weights;

		/* Energy */
		double 									_moment;
		double 									_kinetic_energy;


		/* Velocity bins     ==================================================================== */

		std::unique_ptr<int> 					_number_of_bins;
		double									_bin_width;
		double									_bin_start;
		std::unique_ptr<std::vector<double> >	_mid_bin_array;

		/* Density profiles  ==================================================================== */

		std::unique_ptr<std::vector<double> >	_conditional_velocity_density_function;


		/* velocity profiles ==================================================================== */
		bool									_profiling_active;

		std::unique_ptr<std::vector<double> >	_partial_velocity_profile;
		std::unique_ptr<std::vector<double> >	_accumulated_velocity_profile;
		int 									_velocity_accumulation_interval;
		double 									_accumulation_weight;
		int 									_count;

	public:
		/* constuctor and destructor ============================================================ */
		PopulationOfParticles(std::shared_ptr<const Plasma> plasma,
			std::shared_ptr<int> iteration,
			double population_size, double unit_mass, double unit_charge, double cyclotronic_rotation_parameter);
		PopulationOfParticles(const MacroParameterization & parameterization, const int index, std::shared_ptr<int> iteration);

		~PopulationOfParticles() {}

		/* getters */
		double	get_unit_charge()							const { return _unit_charge;}
		double	get_unit_mass()								const { return _unit_mass;	}
		double	get_size()									const { return *_population_size;	}
		double	get_cyclotronic_parameter()					const { return _cyclotronic_rotation_parameter;	}

		double	get_kinetic_energy()						const { return _kinetic_energy;	}
		double	get_moment()								const { return _moment;	}

		bool	is_velocity_profiling_active()				const { return _profiling_active;	}
		int    	get_number_of_bins()						const { return *_number_of_bins;	}

		/* methods ============================================================================== */
		void 	Reset();
		bool	CheckParameters(const int size, const double unit_mass, const double unit_charge, const double cyclotronic_rotation_parameter);
		void 	ComputeAggregateParameters();
		void 	Accelerate(const PlasmaFields &fields, double factor = 1.0);
		void	Move();
		void	Weigh(PlasmaFields	&fields);

		void	Prepare(const PlasmaFields	&fields);

		void 	SetupVelocityDiagnostics(int nbins, int velocity_accumulation_interval, double vupper, double vlower, double v0, double vt1, double vt2);
		void 	SetupVelocityDiagnostics(const MacroParameterization & parameterization, const int index);
		void	ComputeVelocityProfile();
		
};

#endif

