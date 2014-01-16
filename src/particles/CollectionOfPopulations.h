
/* HEADER PopulationOfParticles ===============================================================
 * author: Paul Cazeaux
 * date: 2013-11-28
 *
 * description:
 *  -> store the populations of particles modeling the plasma
 *	-> Particle loading & pushing
 *  
 * ============================================================================= */

#ifndef DEF_PLASMASCALE_COLLECTIONOFPOPULATIONS
#define DEF_PLASMASCALE_COLLECTIONOFPOPULATIONS

#include <iostream>
#include <string>
#include <vector>

#include "plasma/Plasma.h"
#include "particles/PopulationOfParticles.h"
#include "fields/PlasmaFields.h"

/* Declarations */

class CollectionOfPopulations
{
	private:
		/* class members ======================================================================== */
		/* pointer to objects */
		std::shared_ptr<const Plasma>							_plasma;

		/* Populations */
		const int												_number_of_populations;
		std::vector<std::unique_ptr<PopulationOfParticles> >	_populations_of_particles;

	public:
		/* constuctor and destructor ============================================================ */
		CollectionOfPopulations( std::shared_ptr<const Plasma>	plasma);
		CollectionOfPopulations( std::shared_ptr<const Plasma>	plasma,
			std::shared_ptr<int> iteration, std::shared_ptr<double> simulation_time, FILE *& InputDeck);
		~CollectionOfPopulations() {}

		/* getters for the diagnostics ========================================================== */
		std::vector<std::vector<double> * >	get_vector_of_position_arrays();
		std::vector<std::vector<double> * >	get_vector_of_x_velocity_arrays();
		std::vector<std::vector<double> * >	get_vector_of_y_velocity_arrays();
		std::vector<bool> 					get_vector_of_magnetizations();
		std::vector<int *>					get_vector_of_sizes();

		std::vector<std::vector<double> * > get_vector_of_bin_arrays();
		std::vector<std::vector<double> * > get_vector_of_velocity_profiles();
		std::vector<int *>					get_vector_of_number_of_bins();

		/* methods ============================================================================== */
		void 	ChangeTime(std::shared_ptr<double> simulation_time, std::shared_ptr<int> iteration)
		{
			for (int i=0; i<_number_of_populations; i++)
			{
				_populations_of_particles.at(i)->ChangeTime(simulation_time, iteration);
			}
		}

		void	Prepare(PlasmaFields& fields)
		{
			fields.ComputeAndFilter();
			for (int i=0; i<_number_of_populations; i++)
			{
				_populations_of_particles.at(i)->Prepare(fields);
			}

		}

		void	Weigh(PlasmaFields& fields)
		{
			fields.ResetCharge();
			for (int i=0; i<_number_of_populations; i++)
			{
				_populations_of_particles.at(i)->Weigh(fields);
			}
		}

		void	Accelerate(const PlasmaFields& fields, double factor = 1.0)
		{
			for (int i=0; i<_number_of_populations; i++)
			{
				_populations_of_particles.at(i)->Accelerate(fields, factor);
			}
		}

		void	Move(PlasmaFields& fields)
		{
			for (int i=0; i<_number_of_populations; i++)
			{
				_populations_of_particles.at(i)->Move();
			}
			this->Weigh(fields);
		}

		void	ComputeVelocityProfile();

		void	PushBackKineticEnergy(std::vector<double> &kinetic_energy, std::vector<std::vector<double> >	&kinetic_energy_by_population);
		void	PushBackMoment(std::vector<double> &moment, std::vector<std::vector<double> >	&moment_by_population);
		
};

#endif