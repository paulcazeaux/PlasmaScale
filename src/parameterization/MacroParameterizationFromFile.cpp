#include "parameterization/MacroParameterizationFromFile.h"

/* constuctor and destructor ============================================================ */

MacroParameterizationFromFile::MacroParameterizationFromFile(FILE *& InputDeck)
{
	/* Parameters for the plasma */
	double length, dt, epsilon, la, e0, w0;
	int number_of_populations, grid_size, macro_grid_size, number_of_microsteps, macro_to_micro_dt_ratio, velocity_accumulation_interval, max_mode;
	double filter_parameter_1, filter_parameter_2;

	char a_char[80];

	/* read lines until we get to numbers */
	while (std::fscanf(InputDeck,"%d %lg %d %d %d", &number_of_populations, &dt, &number_of_microsteps, &macro_to_micro_dt_ratio, &velocity_accumulation_interval) <5)
	{
		std::fscanf(InputDeck, "%s", a_char);
	}
	while (std::fscanf(InputDeck,"%lg %d %d %d %lg %lg", &length, &grid_size, &macro_grid_size, &max_mode, &filter_parameter_1, &filter_parameter_2) <6)
	{
		std::fscanf(InputDeck, "%s", a_char);
	}
	while (std::fscanf(InputDeck,"%lg %lg %lg %lg", &epsilon, &la, &e0, &w0) < 4)
	{
		std::fscanf(InputDeck, "%s", a_char);
	}
		/* note: la is l/a */

	if(velocity_accumulation_interval<0) 
	{ 
		std::printf("\nError:  accum can't be negative! \n"); exit(1);
	}

	_plasma = std::make_shared<const Plasma>(length, dt, number_of_microsteps, macro_to_micro_dt_ratio,
				epsilon, la, 0.0, e0, w0,
				number_of_populations, grid_size, macro_grid_size, velocity_accumulation_interval, max_mode,
				filter_parameter_1, filter_parameter_2);
	std::cout << *_plasma << std::endl;

	_number_of_populations = number_of_populations;
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


		_cyclotronic_rotation_parameters.push_back(std::tan(-0.5*wc*_plasma->get_dt()));
		_unit_charges.push_back(_plasma->get_length()*wp*wp/(_plasma->get_epsilon()*n*qm));
		_unit_masses.push_back(_unit_charges.at(population_index)/qm);
		_population_sizes.push_back(n);

		_group_sizes.push_back(n / nlg);
		_mean_velocities.push_back(v0);
		_quiet_start_exponents.push_back(nv2);
		_random_mean_thermal_vel.push_back(vt1);
		_quiet_mean_thermal_vel.push_back(vt2);

		_modes.push_back(mode);
		_density_amplitudes.push_back(x1);
		_density_phases.push_back(thetax);
		_velocity_amplitudes.push_back(v1);
		_velocity_phases.push_back(thetav);

		_bin_numbers.push_back( (nbins>2 ? nbins : 2));
		_upper_velocities.push_back(vupper);
		_lower_velocities.push_back(vlower);
	}
}

void MacroParameterizationFromFile::Load(State & state) const
/* Fill the particle arrays to initialize the microscopic state */
{
	std::vector<std::vector<double> * > positions 		= state.get_vector_of_position_arrays();
	std::vector<std::vector<double> * > x_velocities 	= state.get_vector_of_x_velocity_arrays();
	std::vector<std::vector<double> * > y_velocities 	= state.get_vector_of_y_velocity_arrays();
	std::vector<std::vector<double> * > weights 		= state.get_vector_of_y_velocity_arrays();

	assert(positions.size() ==_number_of_populations);
	assert(x_velocities.size() ==_number_of_populations);
	assert(y_velocities.size() ==_number_of_populations);
	assert(weights.size() ==_number_of_populations);

	state.Reset();

	for (int population_index=0; population_index < _number_of_populations; population_index++)
	{
		int group_size = _group_sizes.at(population_index);
		int population_size = _population_sizes.at(population_index);

		std::vector<double> * position 		= positions.at(population_index);
		std::vector<double> * velocity_x 	= x_velocities.at(population_index);
		std::vector<double> * velocity_y 	= y_velocities.at(population_index);
		std::vector<double> * weight 		= weights.at(population_index);

		assert(position->size() == population_size);
		assert(velocity_x->size() == population_size);
		assert(velocity_y->size() == population_size);
		assert(weight->size() == population_size);

		double v0 = _mean_velocities.at(population_index);
		double vt1 = _random_mean_thermal_vel.at(population_index);
		double vt2 = _quiet_mean_thermal_vel.at(population_index);
		double quiet_start_exponent = _quiet_start_exponents.at(population_index);

		double population_density = _plasma->get_length() / static_cast<double>(population_size);
		double group_density = _plasma->get_length() * static_cast<double>(group_size) / static_cast<double>(population_size);
		static std::vector<double> helper = std::vector<double> (population_size);

		for (int i=0; i < group_size; i++)
		{
			position->at(i) 	= ( static_cast<double>(i) + 0.5 ) * population_density;
			velocity_x->at(i) 	= v0;
		}

		if (vt2 != 0.0) 
		{
			double vmax = 5.0*vt2;
			double dv = 2.0*vmax/(static_cast<double>(population_size - 1));
			double vvnv2 = 1.0;
			helper.front() = 0.0;

			for (int i=1; i < helper.size(); i++)
			{
				double vv = ((i - 0.5)*dv - vmax)/vt2;
				if (quiet_start_exponent != 0)
				{
					vvnv2 = std::pow(vv,static_cast<double>(quiet_start_exponent));
				}
				helper.at(i) = helper.at(i-1) + vvnv2*std::exp(-0.5*vv*vv);
			}

			double df = helper.back()/group_size;

			auto it_helper = helper.begin();
			for (int i=0; i<group_size; i++)
			{
				double fv = (i + 0.5)*df;
				while (fv >= *(it_helper+1)) 
				{
					it_helper++;
					if (it_helper == helper.end()) 
						perror("distribution function error");
				}
				velocity_x->at(i) += dv* (std::distance(helper.begin(), it_helper)
										 + (fv - *it_helper)/(*(it_helper+1) - *it_helper)) - vmax;
			}

			double xs = 0.0;
			for (int i=0; i<group_size; i++)
			{
				position->at(i) = xs*group_density + 0.5*population_density;
				/* bit-reversed scrambling to reduce the correlation with the positions */
				double xsi = 0.5;
				xs -= 0.5;
				while (xs >= 0.0)
				{
					xsi *= 0.5;
					xs -= xsi;
				} 
				xs += 2.0*xsi;
			}
		}
		bool magnetized = (_cyclotronic_rotation_parameters.at(population_index) != 0.) ;
		if (magnetized) 
		{
			for (int i=0; i < group_size; i++) 
			{
				double v = velocity_x->at(i);
				double theta = RandomTools::Generate_randomly_uniform(0., 2.*M_PI);
				velocity_x->at(i) = v*std::cos(theta);
				velocity_y->at(i) = v*std::sin(theta);
			}
		}
		if (group_size < population_size)
		{
			double xs = 0.0;
			for (int k=group_size; k< population_size; k+=group_size)
			{
				xs += group_density;
				for (int j=0; j < group_size; j++) 
				{
					 position->at(j+k) = position->at(j) + xs;
					 velocity_x->at(j+k) = velocity_x->at(j);
					 if (magnetized) 
					 {
						velocity_x->at(j+k) = velocity_y->at(j);
					 }
				}
			}
		}

		if (vt1 != 0.0)
		{
			auto it_vel_y = velocity_y->begin();
			for (auto & vel_x : *velocity_x)
			{
				vel_x 	+= vt1*RandomTools::Generate_randomly_normal(0.0, 1.0);
				if (magnetized)
				{
					*it_vel_y++ 	+= vt1*RandomTools::Generate_randomly_normal(0.0, 1.0);
				}
			}
		}

		/* Add the perturbation */

		double w = 2*M_PI*_modes.at(population_index)/_plasma->get_length();
		double x1 = _density_amplitudes.at(population_index);
		double v1 = _velocity_amplitudes.at(population_index);
		double thetax = _density_phases.at(population_index);
		double thetav = _velocity_phases.at(population_index);

		auto it_vel_x 	= velocity_x->begin();
		auto it_weights = weight->begin();

		for (auto & x : *position)
		{
			double theta = w * x;
			*it_weights++ = 1. + x1 * std::sin(theta + thetax);
			*it_vel_x++  += v1 * std::sin(theta + thetav);
		}
	}
}

double MacroParameterizationFromFile::get_initial_thermal_vel(int population_index) const
{
	return _quiet_mean_thermal_vel.at(population_index) + _random_mean_thermal_vel.at(population_index);
}

double MacroParameterizationFromFile::GetBinStart(int population_index) const
{
	double vupper = _upper_velocities.at(population_index);
	double vlower = _lower_velocities.at(population_index);
	double vt1 = _random_mean_thermal_vel.at(population_index);
	double vt2 = _quiet_mean_thermal_vel.at(population_index);
	double v0 = _mean_velocities.at(population_index);

	if(vupper-vlower < 0.0)
	{
		printf("\nInitialization error: vupper must be > vlower! \n");
		exit(1);
	}
	if(vt1<0 ||vt2<0)
	{ 
		printf("\nInitialization error: can't have negative thermal voltages!\n");
		exit(1);
	}
	if(vupper-vlower > 0.0)
	{
		return vlower;
	}
	else if(vt1 + vt2 > 0.0) 
	{
		return v0 - 5.0*(vt1 + vt2);
	}
	else if (v0!=0)
	{
		return (v0 < 0 ? 2.0*v0 : 0);
	}
	else
	{
		return 0.0;
	}

}

double MacroParameterizationFromFile::GetBinWidth(int population_index)	const
{
	double vupper = _upper_velocities.at(population_index);
	double vlower = _lower_velocities.at(population_index);
	double vt1 = _random_mean_thermal_vel.at(population_index);
	double vt2 = _quiet_mean_thermal_vel.at(population_index);
	double v0 = _mean_velocities.at(population_index);
	int nbins = _bin_numbers.at(population_index);

	if(vupper-vlower < 0.0)
	{
		printf("\nInitialization error: vupper must be > vlower! \n");
		exit(1);
	}
	if(vt1<0 ||vt2<0)
	{ 
		printf("\nInitialization error: can't have negative thermal voltages!\n");
		exit(1);
	}
	if(vupper-vlower > 0.0)
	{
		return (vupper-vlower)/static_cast<double>(nbins);
	}
	else if(vt1 + vt2 > 0.0) 
	{
		return 10.0 * (vt1 + vt2)/static_cast<double>(nbins);  /*  so that the distribution goes from
							v0-5*vt to v0+5*vt  */
	}
	else if (v0!=0)
	{
		return 2.0*std::abs(v0)/static_cast<double>(nbins);
	}
	else
	{
		return 1.0/static_cast<double>(nbins);
	}

}


int MacroParameterizationFromFile::GetNumberOfBins(int population_index) const
{
	return _bin_numbers.at(population_index);
}
