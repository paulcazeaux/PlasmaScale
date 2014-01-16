#include "CollectionOfPopulations.h"

/* constuctor and destructor ============================================================ */
CollectionOfPopulations::CollectionOfPopulations( std::shared_ptr<const Plasma>	plasma) :
		_plasma(plasma),
		_number_of_populations(plasma->get_number_of_populations())
{
	_populations_of_particles = std::vector<std::unique_ptr<PopulationOfParticles> > (_number_of_populations);
}

CollectionOfPopulations::CollectionOfPopulations( std::shared_ptr<const Plasma>	plasma,
			std::shared_ptr<int> iteration, std::shared_ptr<double> simulation_time, FILE *& InputDeck) :
		_plasma(plasma), 
		_number_of_populations(plasma->get_number_of_populations())
{
	_populations_of_particles = std::vector<std::unique_ptr<PopulationOfParticles> > (_number_of_populations);

	for (int population_index = 0; population_index < _number_of_populations; population_index++)
	{

		int n, nv2, nlg, mode, nbins;
		double wp, wc, qm, vt1, vt2, v0;
		double x1, v1, thetax, thetav;
		double vlower, vupper;

		char a_char[80];

		printf("\n");
		while (fscanf(InputDeck, "%d %d %d %d ", &n, &nv2, &nlg, &mode) < 4)  /*  added nbins  */
		{
			fscanf(InputDeck,"%s",a_char);
		}
		while (fscanf(InputDeck, "%lg %lg %lg %lg %lg %lg", &wp, &wc, &qm, &vt1, &vt2, &v0) < 6)
		{
			fscanf(InputDeck,"%s",a_char);
		}
		while (fscanf(InputDeck,"%lg %lg %lg %lg", &x1, &v1, &thetax, &thetav) < 4)
		{
			fscanf(InputDeck,"%s",a_char);
		}
		while(fscanf(InputDeck," %d %lg %lg ",&nbins,&vlower,&vupper) < 3)
		{
			fscanf(InputDeck,"%s",a_char);
		}
		printf("\n");
		printf(" n = %4d  nv2 = %4d  nlg = %4d  mode = %4d \n", n, nv2, nlg, mode);
		printf(" wp = %6.3f   wc = %6.3f   qm = %6.3f \n", wp, wc, qm);
		printf(" vt1 = %6.3f  vt2 = %6.3f  v0 = %6.3f \n", vt1, vt2, v0);
		printf(" x1 = %6.3f   v1 = %6.3f   thetax = %6.3f   thetav = %6.3f \n", x1, v1, thetax, thetav);
		printf(" nbins = %4d   vlower = %6.3f  vupper = %6.3f ", nbins, vlower, vupper);
		printf("\n");

		double cyclotronic_rotation_parameter = tan(-0.5*wc*plasma->get_dt());
		double unit_charge 	= plasma->get_length()*wp*wp/(plasma->get_epsilon()*n*qm);
		double unit_mass 	= unit_charge/qm;
		int population_size = n;
		int group_size = n / nlg;
		int quiet_start_exponent = nv2;

		_populations_of_particles[population_index] = std::unique_ptr<PopulationOfParticles> 
				(new PopulationOfParticles(plasma, iteration, simulation_time,
						population_size, unit_mass, unit_charge, cyclotronic_rotation_parameter));
		_populations_of_particles.at(population_index)->Load(group_size, v0, vt1, vt2, 
													mode, x1, thetax, v1, thetav, quiet_start_exponent);
		_populations_of_particles.at(population_index)->SetupVelocityDiagnostics(nbins, plasma->get_velocity_accumulation_interval(), vupper, vlower, v0, vt1, vt2);
		
	}

}

/* getters for the diagnostics ========================================================== */
std::vector<std::vector<double> * >	CollectionOfPopulations::get_vector_of_position_arrays()
{
	std::vector<std::vector<double> * > positions = std::vector<std::vector<double> * > (_number_of_populations);
	for (int i = 0; i < _number_of_populations; i++)
	{
		positions[i] = (_populations_of_particles.at(i)->_position).get();
	}
	return positions;
}

std::vector<std::vector<double> * >	CollectionOfPopulations::get_vector_of_x_velocity_arrays()
{
	std::vector<std::vector<double> * > velocities = std::vector<std::vector<double> * > (_number_of_populations);
	for (int i = 0; i < _number_of_populations; i++)
	{
		velocities[i] = (_populations_of_particles.at(i)->_velocity_x).get();
	}
	return velocities;
}

std::vector<std::vector<double> * >	CollectionOfPopulations::get_vector_of_y_velocity_arrays()
{
	std::vector<std::vector<double> * > velocities = std::vector<std::vector<double> * > (_number_of_populations);
	for (int i = 0; i < _number_of_populations; i++)
	{
		velocities[i] = (_populations_of_particles.at(i)->_velocity_y).get();
	}
	return velocities;
}

std::vector<bool> CollectionOfPopulations::get_vector_of_magnetizations()
{
	std::vector<bool> magnetizations = std::vector<bool>(_number_of_populations);
	for (int i = 0; i < _number_of_populations; i++)
	{
		magnetizations[i] = _populations_of_particles.at(i)->_magnetized;
	}
	return magnetizations;
}

std::vector<int *>	CollectionOfPopulations::get_vector_of_sizes()
{
	std::vector<int *> sizes = std::vector<int*>(_number_of_populations);
	for (int i = 0; i < _number_of_populations; i++)
	{
		sizes[i] = _populations_of_particles.at(i)->_population_size.get();
	}
	return sizes;
}

std::vector<std::vector<double> * >	CollectionOfPopulations::get_vector_of_bin_arrays()
{
	std::vector<std::vector<double> * > bin_arrays = std::vector<std::vector<double> * > (_number_of_populations);
	for (int i = 0; i < _number_of_populations; i++)
	{
		bin_arrays[i] = (_populations_of_particles.at(i)->_mid_bin_array).get();
	}
	return bin_arrays;
}

std::vector<std::vector<double> * >	CollectionOfPopulations::get_vector_of_velocity_profiles()
{
	std::vector<std::vector<double> * > velocity_profiles = std::vector<std::vector<double> * > (_number_of_populations);
	for (int i = 0; i < _number_of_populations; i++)
	{
		velocity_profiles[i] = (_populations_of_particles.at(i)->_accumulated_velocity_profile).get();
	}
	return velocity_profiles;
}

std::vector<int *>	CollectionOfPopulations::get_vector_of_number_of_bins()
{
	std::vector<int *> number_of_bins = std::vector<int *>(_number_of_populations);
	for (int i = 0; i < _number_of_populations; i++)
	{
		number_of_bins[i] = _populations_of_particles.at(i)->_number_of_bins.get();
	}
	return number_of_bins;
}


/* methods ============================================================================== */

void 	CollectionOfPopulations::ComputeVelocityProfile()
{
	for (int i=0; i<_number_of_populations; i++)
	{
		_populations_of_particles.at(i)->ComputeVelocityProfile();
	}
}


void	CollectionOfPopulations::PushBackKineticEnergy(std::vector<double> &kinetic_energy, std::vector<std::vector<double> >	&kinetic_energy_by_population)
{
	double ke_total = 0.0;
	for (int i=0; i<_number_of_populations; i++)
	{
		double ke_pop = _populations_of_particles.at(i)->get_kinetic_energy();
		kinetic_energy_by_population.at(i).push_back(ke_pop + 1e-30);
		ke_total += ke_pop;
	}
	kinetic_energy.push_back(ke_total);
}

void	CollectionOfPopulations::PushBackMoment(std::vector<double> &moment, std::vector<std::vector<double> >	&moment_by_population)
{
	double p_total = 0.0;
	for (int i=0; i<_number_of_populations; i++)
	{
		double p_pop = _populations_of_particles.at(i)->get_moment();
		moment_by_population.at(i).push_back(p_pop);
		p_total += p_pop;
	}
	moment.push_back(p_total);
}
