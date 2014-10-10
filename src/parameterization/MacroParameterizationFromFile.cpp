#include "parameterization/MacroParameterizationFromFile.h"

/* constuctor and destructor ============================================================ */

MacroParameterizationFromFile::MacroParameterizationFromFile(FILE *& InputDeck)
{
	/* Parameters for the plasma */
	double length, dt, epsilon, la, e0, w0;
	int number_of_populations, grid_size, macro_grid_size, number_of_microsteps, macro_to_micro_dt_ratio, velocity_accumulation_interval, depth, cutoff, max_mode;
	double filter_parameter_1, filter_parameter_2, intensity;

	int max_size_history, use_full_PIC, record_microsteps;
	char a_char[80];

	/* read lines until we get to numbers */
	while (std::fscanf(InputDeck,"%d %lg %d %d %d %d %d %lg", &number_of_populations, &dt, &number_of_microsteps, &macro_to_micro_dt_ratio, &velocity_accumulation_interval, &depth, &cutoff, &intensity) <8)
	{
		std::fscanf(InputDeck, "%s", a_char);
	}
	while (std::fscanf(InputDeck,"%lg %d %d %d %lg %lg %d %d %d", &length, &grid_size, &macro_grid_size, &max_mode, &filter_parameter_1, &filter_parameter_2, &max_size_history, &use_full_PIC, &record_microsteps) <9)
	{
		std::fscanf(InputDeck, "%s", a_char);
	}
	while (std::fscanf(InputDeck,"%lg %lg %lg %lg", &epsilon, &la, &e0, &w0) < 4)
	{
		std::fscanf(InputDeck, "%s", a_char);
	}
		/* note: la is l/a */

	/* Parameters for the initialization */
	while (std::fscanf(InputDeck,"%lg", &_init_occupation) < 1)
	{
		std::fscanf(InputDeck, "%s", a_char);
	}


	if(velocity_accumulation_interval<0) 
	{ 
		std::printf("\nError:  accum can't be negative! \n"); exit(1);
	}

	_plasma = std::make_shared<const Plasma>(length, dt, number_of_microsteps, macro_to_micro_dt_ratio,
				epsilon, la, 0.0, e0, w0,
				number_of_populations, grid_size, macro_grid_size, velocity_accumulation_interval, max_mode, depth, cutoff, intensity,
				filter_parameter_1, filter_parameter_2, max_size_history, static_cast<bool>(use_full_PIC), static_cast<bool>(record_microsteps));
	std::cout << *_plasma << std::endl;
	std::cout << "Initialization of a plasma expansion :" << std::endl;
	std::cout << "Initial size: " << _init_occupation * _plasma->get_length() << " out of a total length: " << _plasma->get_length() << std::endl;

	_number_of_populations = number_of_populations;
	if (_number_of_populations !=2)
	{
		std::cout << "Wrong number of particle populations! We cannot proceed with the test." << std::endl;
		exit(1);
	}
	for (int population_index = 0; population_index < _number_of_populations; population_index++)
	{
		int n, nv2, nlg, nbins;
		double wp, wc, qm, vt1, vt2, v0;
		double vlower, vupper;

		char a_char[80];

		printf("\n");
		while (fscanf(InputDeck, "%d %d %d ", &n, &nv2, &nlg) < 3)
		{
			fscanf(InputDeck,"%s",a_char);
		}
		while (fscanf(InputDeck, "%lg %lg %lg %lg %lg %lg", &wp, &wc, &qm, &vt1, &vt2, &v0) < 6)
		{
			fscanf(InputDeck,"%s",a_char);
		}
		while(fscanf(InputDeck," %d %lg %lg ",&nbins,&vlower,&vupper) < 3)
		{
			fscanf(InputDeck,"%s",a_char);
		}
		std::cout << "Parameters for population " << (population_index+1) << std::endl << "---------------------------" << std::endl;
		std::cout << "Number of particles:\t" << n << "\tNumber of init. groups:\t" << nlg <<  "\tQuiet start exponent:\t" << nv2 << std::endl;
		std::cout << "Plasma pulsation:\t" << wp << "\tCycl. pulsation:\t" << wc << "\tq/m ratio:\t\t" << qm << std::endl;
		std::cout << "Random therm. veloc.:\t" << vt1 << "\tQuiet therm. veloc.:\t" << vt2 << "\tMean velocity:\t\t" << v0 << std::endl;
		std::cout << "Velocity init.:\t# bins:\t" << nbins << "\tLower bound:\t\t" << vlower << "\tUpper bound:\t\t" << vupper << std::endl;
		std::cout << "---------------------------" << std::endl;


		_cyclotronic_rotation_parameters.push_back(std::tan(-0.5*wc*_plasma->get_dt()));
		_unit_charges.push_back(_plasma->get_length()*wp*wp/(_plasma->get_epsilon()*n*qm));
		_unit_masses.push_back(_unit_charges.at(population_index)/qm);
		_plasma_pulsations.push_back(wp);
		_population_sizes.push_back(n);

		_group_sizes.push_back(n / nlg);
		_mean_velocities.push_back(v0);
		_quiet_start_exponents.push_back(nv2);
		_random_mean_thermal_vel.push_back(vt1);
		_quiet_mean_thermal_vel.push_back(vt2);

		_bin_numbers.push_back( (nbins>2 ? nbins : 2));
		_upper_velocities.push_back(vupper);
		_lower_velocities.push_back(vlower);
	}

	_debye_scaling = std::pow((_random_mean_thermal_vel.back()+_quiet_mean_thermal_vel.back())/_plasma_pulsations.back(), 2.);

}

MacroParameterizationFromFile::MacroParameterizationFromFile(FILE *& InputDeck, int number_of_microsteps)
{
	/* Parameters for the plasma */
	double length, dt, epsilon, la, e0, w0;
	int number_of_populations, grid_size, macro_grid_size, unused_number_of_microsteps, macro_to_micro_dt_ratio, velocity_accumulation_interval, depth, cutoff, max_mode;
	double filter_parameter_1, filter_parameter_2, intensity;

	int max_size_history, use_full_PIC, record_microsteps;
	char a_char[80];

	/* read lines until we get to numbers */
	while (std::fscanf(InputDeck,"%d %lg %d %d %d %d %d %lg", &number_of_populations, &dt, &unused_number_of_microsteps, &macro_to_micro_dt_ratio, &velocity_accumulation_interval, &depth, &cutoff, &intensity) <8)
	{
		std::fscanf(InputDeck, "%s", a_char);
	}
	while (std::fscanf(InputDeck,"%lg %d %d %d %lg %lg %d %d %d", &length, &grid_size, &macro_grid_size, &max_mode, &filter_parameter_1, &filter_parameter_2, &max_size_history, &use_full_PIC, &record_microsteps) <9)
	{
		std::fscanf(InputDeck, "%s", a_char);
	}
	while (std::fscanf(InputDeck,"%lg %lg %lg %lg", &epsilon, &la, &e0, &w0) < 4)
	{
		std::fscanf(InputDeck, "%s", a_char);
	}
		/* note: la is l/a */

	/* Parameters for the initialization */
	while (std::fscanf(InputDeck,"%lg", &_init_occupation) < 1)
	{
		std::fscanf(InputDeck, "%s", a_char);
	}


	if(velocity_accumulation_interval<0) 
	{ 
		std::printf("\nError:  accum can't be negative! \n"); exit(1);
	}

	_plasma = std::make_shared<const Plasma>(length, dt, number_of_microsteps, macro_to_micro_dt_ratio,
				epsilon, la, 0.0, e0, w0,
				number_of_populations, grid_size, macro_grid_size, velocity_accumulation_interval, max_mode, depth, cutoff, intensity,
				filter_parameter_1, filter_parameter_2, max_size_history, static_cast<bool>(use_full_PIC), static_cast<bool>(record_microsteps));
	std::cout << *_plasma << std::endl;
	std::cout << "Initialization of a plasma expansion :" << std::endl;
	std::cout << "Initial size: " << _init_occupation * _plasma->get_length() << " out of a total length: " << _plasma->get_length() << std::endl;

	_number_of_populations = number_of_populations;
	if (_number_of_populations !=2)
	{
		std::cout << "Wrong number of particle populations! We cannot proceed with the test." << std::endl;
		exit(1);
	}
	for (int population_index = 0; population_index < _number_of_populations; population_index++)
	{
		int n, nv2, nlg, nbins;
		double wp, wc, qm, vt1, vt2, v0;
		double vlower, vupper;

		char a_char[80];

		printf("\n");
		while (fscanf(InputDeck, "%d %d %d ", &n, &nv2, &nlg) < 3)
		{
			fscanf(InputDeck,"%s",a_char);
		}
		while (fscanf(InputDeck, "%lg %lg %lg %lg %lg %lg", &wp, &wc, &qm, &vt1, &vt2, &v0) < 6)
		{
			fscanf(InputDeck,"%s",a_char);
		}
		while(fscanf(InputDeck," %d %lg %lg ",&nbins,&vlower,&vupper) < 3)
		{
			fscanf(InputDeck,"%s",a_char);
		}
		std::cout << "Parameters for population " << (population_index+1) << std::endl << "---------------------------" << std::endl;
		std::cout << "Number of particles:\t" << n << "\tNumber of init. groups:\t" << nlg <<  "\tQuiet start exponent:\t" << nv2 << std::endl;
		std::cout << "Plasma pulsation:\t" << wp << "\tCycl. pulsation:\t" << wc << "\tq/m ratio:\t\t" << qm << std::endl;
		std::cout << "Random therm. veloc.:\t" << vt1 << "\tQuiet therm. veloc.:\t" << vt2 << "\tMean velocity:\t\t" << v0 << std::endl;
		std::cout << "Velocity init.:\t# bins:\t" << nbins << "\tLower bound:\t\t" << vlower << "\tUpper bound:\t\t" << vupper << std::endl;
		std::cout << "---------------------------" << std::endl;


		_cyclotronic_rotation_parameters.push_back(std::tan(-0.5*wc*_plasma->get_dt()));
		_unit_charges.push_back(_plasma->get_length()*wp*wp/(_plasma->get_epsilon()*n*qm));
		_unit_masses.push_back(_unit_charges.at(population_index)/qm);
		_plasma_pulsations.push_back(wp);
		_population_sizes.push_back(n);

		_group_sizes.push_back(n / nlg);
		_mean_velocities.push_back(v0);
		_quiet_start_exponents.push_back(nv2);
		_random_mean_thermal_vel.push_back(vt1);
		_quiet_mean_thermal_vel.push_back(vt2);

		_bin_numbers.push_back( (nbins>2 ? nbins : 2));
		_upper_velocities.push_back(vupper);
		_lower_velocities.push_back(vlower);
	}

	_debye_scaling = std::pow((_random_mean_thermal_vel.back()+_quiet_mean_thermal_vel.back())/_plasma_pulsations.back(), 2.);

}

void MacroParameterizationFromFile::Load(State & state) const
/* Fill the particle arrays to initialize the microscopic state */
{
	std::vector<std::vector<double> * > positions 		= state.get_vector_of_position_arrays();
	std::vector<std::vector<double> * > x_velocities 	= state.get_vector_of_x_velocity_arrays();
	std::vector<std::vector<double> * > y_velocities 	= state.get_vector_of_y_velocity_arrays();
	std::vector<std::vector<double> * > weights 		= state.get_vector_of_weight_arrays();

	assert(positions.size() ==_number_of_populations);
	assert(x_velocities.size() ==_number_of_populations);
	assert(y_velocities.size() ==_number_of_populations);
	assert(weights.size() ==_number_of_populations);

	state.Reset();

	for (int population_index=0; population_index < _number_of_populations; population_index++)
	{
		int group_size = _group_sizes.at(population_index);
		int population_size = _population_sizes.at(population_index);
		int nbins = _bin_numbers.at(population_index);

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

		static std::vector<double> helper;
		helper.resize(nbins);

		for (int i=0; i < group_size; i++)
		{
			position->at(i) 	= ( static_cast<double>(i) + 0.5 ) * population_density;
			velocity_x->at(i) 	= v0;
		}

		if (vt2 != 0.0) 
		{
			double vmax = 5.0*vt2;
			double dv = 2.0*vmax/(static_cast<double>(nbins - 1));
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

	}

	/* Implement the population expansion */
	int size = _plasma->get_grid_size();
	double tol2 = 1e-20, iter_max = 100;
	double offset = (1.-_init_occupation)*_plasma->get_length()/2.;
	double population_density = static_cast<double>(size)/static_cast<double>(_population_sizes.front());

	std::vector<double> 	ion_density = std::vector<double>(size),
							potential = std::vector<double>(size),
						   	exp_potential = std::vector<double>(size), 
						   	update = std::vector<double>(size);
	std::vector<double> residual = std::vector<double>(size),
							conj_dir = std::vector<double>(size),
							dir = std::vector<double>(size);

	for (int i=0; i<_population_sizes.front(); i++)
	{
		double x 			= _init_occupation*positions.front()->at(i) + offset;
		int bin 			= _plasma->find_index_on_grid(x);
		double cellpos		= _plasma->find_position_in_cell(x);

		positions.front()->at(i) = x;
		weights.front()->at(i) = 1.;

		ion_density[bin] += (1. - cellpos);
		if (bin+1< size)
			ion_density[bin+1] += cellpos;
		else
			ion_density[0] += cellpos;
	}

	/* Determine the electron density from the Boltzmann approximation */
	/* Initialization */

	for (int i=0; i<size; i++)
	{
		ion_density.at(i) = population_density*ion_density.at(i);
		potential.at(i) = std::log(ion_density.at(i)+1e-10);
	}


	/* Newton's method loop */
	double scaling = - _debye_scaling/std::pow(_plasma->get_dx(), 2.);
	for (int count=0; count<iter_max; count++)
	{
		/* Inner solve : CG algorithm */

		for (int i=0; i<size; i++)
		{
			exp_potential.at(i) = std::exp(potential.at(i));
				/* Initial value for the CG algorithm : zero. */
			residual.at(i)	=  exp_potential.at(i)  +  scaling * (
											  potential.at((i>0?i-1:size-1)) 
											+ potential.at((i<size-1?i+1:0)) 
											- 2*potential.at(i)
												    			 )
									- ion_density.at(i);
			conj_dir.at(i) = residual.at(i);
		}
		std::fill(update.begin(), update.end(), 0.);


		/* Inner loop */

		double rnorm2 = 0, delta, alpha, beta;
		for (int i=0; i<size; i++)
		{
			rnorm2 += residual.at(i)*residual.at(i);
		}

		for (int count_cg=0; count_cg<iter_max; count_cg++)
		{
			/* Compute the matrix-vector product */
			delta = 0;
			for (int i=0; i<size; i++)
			{
				dir.at(i) = (exp_potential.at(i)-2*scaling)*conj_dir.at(i) 
										+ scaling * (	conj_dir.at((i>0?i-1:size-1)) 
													  + conj_dir.at((i<size-1?i+1:0))
													);
				delta += dir.at(i) * conj_dir.at(i);
			}

			beta = 1./rnorm2;
			alpha = rnorm2/delta;
			rnorm2 = 0;
			for (int i=0; i<size; i++)
			{
				update.at(i) += alpha * conj_dir.at(i);
				residual.at(i) -= alpha* dir.at(i);
				rnorm2 += residual.at(i)*residual.at(i);
			}
			if (rnorm2 < tol2/2.)
				break;
			beta *= rnorm2;
			for (int i=0; i<size; i++)
			{
				conj_dir.at(i) *= beta;
				conj_dir.at(i) += residual.at(i);
			}
		}

		double err=0;
		for (int i=0; i<size; i++)
		{
			err += update.at(i) * update.at(i);
			potential.at(i) -= update.at(i);
		}
		if (err < tol2)
			{
				std::cout << "Newton error: " << err << std::endl;
				break;
			}
	}
	for (int i=0; i<size; i++)
		exp_potential.at(i) = std::exp(potential.at(i));

	auto it_weights = weights.at(1)->begin();
	for (auto & x : *positions.back())
	{
		int bin = _plasma->find_index_on_grid(x);
		double cellpos = _plasma->find_position_in_cell(x);
		*it_weights++ = Tools::EvaluateP1Function(exp_potential, bin, cellpos);
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

	if(vupper - vlower < 0.0)
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

double MacroParameterizationFromFile::GetBinEnd(int population_index) const
{
	double vupper = _upper_velocities.at(population_index);
	double vlower = _lower_velocities.at(population_index);
	double vt1 = _random_mean_thermal_vel.at(population_index);
	double vt2 = _quiet_mean_thermal_vel.at(population_index);
	double v0 = _mean_velocities.at(population_index);

	if(vupper - vlower < 0.0)
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
		return vupper;
	}
	else if(vt1 + vt2 > 0.0) 
	{
		return v0 + 5.0*(vt1 + vt2);
	}
	else if (v0!=0)
	{
		return (v0 > 0 ? 2.0*v0 : 0);
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
