#include "parameterization/MacroParameterizationFromFile.h"

/* constuctor and destructor ============================================================ */

MacroParameterizationFromFile::MacroParameterizationFromFile(FILE *& InputDeck)
{
	/* Parameters for the plasma */
	double length, dt, epsilon, e0, w0;
	int number_of_populations, grid_end, macro_grid_end, number_of_microsteps, macro_to_micro_dt_ratio, velocity_accumulation_interval, depth, cutoff;
	double intensity;

	int max_size_history, use_full_PIC, record_microsteps;
	char a_char[80];

	std::cout << std::scientific;

	/* read lines until we get to numbers */
	while (std::fscanf(InputDeck,"%d %lg %d %d %d %d %d %lg", &number_of_populations, &dt, &number_of_microsteps, &macro_to_micro_dt_ratio, &velocity_accumulation_interval, &depth, &cutoff, &intensity) <8)
	{
		std::fscanf(InputDeck, "%s", a_char);
	}
	while (std::fscanf(InputDeck,"%lg %d %d %d %d %d", &length, &grid_end, &macro_grid_end, &max_size_history, &use_full_PIC, &record_microsteps) <6)
	{
		std::fscanf(InputDeck, "%s", a_char);
	}
	while (std::fscanf(InputDeck,"%lg %lg %lg", &epsilon, &e0, &w0) < 3)
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

	if (use_full_PIC == 1)
	{
		number_of_microsteps = 0;
	}
	_plasma = std::make_shared<const Plasma>(length, dt, number_of_microsteps, macro_to_micro_dt_ratio,
				epsilon, e0, w0,
				number_of_populations, grid_end, macro_grid_end, velocity_accumulation_interval, depth, cutoff, intensity,
				max_size_history, static_cast<bool>(use_full_PIC), static_cast<bool>(record_microsteps));
	std::cout << *_plasma << std::endl;
	std::cout << "Initialization of a plasma expansion :" << std::endl;
	std::cout << "Initial size: " << _init_occupation * _plasma->_length << " out of a total length: " << _plasma->_length << std::endl;

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
		while (fscanf(InputDeck, "%lg %lg %lg %lg %lg", &wp, &qm, &vt1, &vt2, &v0) < 5)
		{
			fscanf(InputDeck,"%s",a_char);
		}
		while(fscanf(InputDeck," %d %lg %lg ",&nbins,&vlower,&vupper) < 3)
		{
			fscanf(InputDeck,"%s",a_char);
		}
		std::cout << "Parameters for population " << (population_index+1) << std::endl << "---------------------------" << std::endl;
		std::cout << "Number of particles:\t" << n << "\t\tNumber of init. groups:\t" << nlg <<  "\tQuiet start exponent:\t" << nv2 << std::endl;
		std::cout << "Plasma pulsation:\t" << wp << "\tq/m ratio:\t" << qm << std::endl;
		std::cout << "Random therm. veloc.:\t" << vt1 << "\tQuiet therm. veloc.:\t" << vt2 << "\tMean velocity:\t\t" << v0 << std::endl;
		std::cout << "Velocity init.:\t# bins:\t" << nbins << "\t\tLower bound:\t\t" << vlower << "\tUpper bound:\t\t" << vupper << std::endl;
		std::cout << "---------------------------" << std::endl;

		if(vupper-vlower < 0.0)
		{
			printf("\nInitialization error: vupper must be > vlower! \n");
			exit(1);
		}
		if(vt1<0 || vt2<0)
		{ 
			printf("\nInitialization error: can't have negative thermal voltages!\n");
			exit(1);
		}

		_reference_densities.push_back(M_2_SQRTPI*n/(_init_occupation*_plasma->_length));
		_unit_charges.push_back(_plasma->_epsilon*wp*wp/(_reference_densities.at(population_index)*qm));
		_unit_masses.push_back(_unit_charges.at(population_index)/qm);
		_plasma_pulsations.push_back(wp);
		_population_sizes.push_back(n);

		_group_sizes.push_back(n / nlg);
		_mean_velocities.push_back(v0);
		_quiet_start_exponents.push_back(nv2);
		_random_mean_thermal_vel.push_back(vt1);
		_quiet_mean_thermal_vel.push_back(vt2);
		_thermal_velocities.push_back(vt1+vt2);

		_bin_numbers.push_back( (nbins>2 ? nbins : 2));
		_upper_velocities.push_back(vupper);
		_lower_velocities.push_back(vlower);
	}
}

void MacroParameterizationFromFile::Load(State & state) const
/* Fill the particle arrays to initialize the microscopic state */
{
	std::vector<std::vector<double> * > positions 	= state.get_vector_of_position_arrays();
	std::vector<std::vector<double> * > velocities 	= state.get_vector_of_velocity_arrays();
	std::vector<std::vector<double> * > weights 	= state.get_vector_of_weight_arrays();

	assert(positions.size() ==_number_of_populations);
	assert(velocities.size() ==_number_of_populations);
	assert(weights.size() ==_number_of_populations);

	for (int population_index=0; population_index < _number_of_populations; population_index++)
	{
		int population_size = _population_sizes.at(population_index);

		assert(positions.at(population_index)->size() == population_size);
		assert(velocities.at(population_index)->size() == population_size);
		assert(weights.at(population_index)->size() == population_size);
	}


	state.Reset();

	/**************************************/
	/* Now we implement a Gaussian plasma */
	/**************************************/

	/* Characteristic radius of the spatial Gaussian profile */
	double r0 = _init_occupation*_plasma->_length;

	double xs = 0.0;
	double x0 = .5/_population_sizes.front();
	for (int i=0; i<_population_sizes.front(); i++)
	{
		double x = _plasma->_length*(xs + x0);
		positions.front()->at(i)	= x;
		x /= r0;
		weights.front()->at(i) 		= std::exp(-x*x);

		/* bit-reversed scrambling to reduce the correlation between positions and velocities */
		double xsi = 0.5;
		xs -= 0.5;
		while (xs >= 0.0)
		{
			xsi *= 0.5;
			xs -= xsi;
		} 
		xs += 2.0*xsi;
	}

	double w0 = _population_sizes.front() / std::accumulate(weights.front()->begin(), weights.front()->end(), 0.);
	for (double & w: *weights.front()) w *= w0;

	/*----------------------------------------------------------------*/
	/* Newton's method solve for the electrostatic potential assuming */
	/* 			a Boltzmann distribution for the electrons 			  */
	/*----------------------------------------------------------------*/
		/* Initialization */

	double tol2 = 1e-30, iter_max = 100;
	double debye_length = (_random_mean_thermal_vel.back()+_quiet_mean_thermal_vel.back())/_plasma_pulsations.back();
	double scaling = std::pow(debye_length/_plasma->_dx, 2.);
	int grid_end = _plasma->_grid_end;

	std::vector<double> 	ion_density = std::vector<double>(grid_end+1),
						   	exp_potential = std::vector<double>(grid_end+1);
	Eigen::Map<Eigen::VectorXd> Ion_density(ion_density.data(), grid_end+1);
	Eigen::Map<Eigen::VectorXd> Exp_potential(exp_potential.data(), grid_end+1);
	Eigen::VectorXd Residual(grid_end+1), Potential(grid_end+1), Update(grid_end+1);

	std::vector<Eigen::Triplet<double> > J_T;
	Eigen::SparseMatrix<double> J;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > solver;

	/* Neumann condition at x = 0 */
	J_T.push_back(Eigen::Triplet<double>(0,0, scaling));
	J_T.push_back(Eigen::Triplet<double>(0,1, -scaling));
	for (int i=1; i<grid_end; i++)
	{
		J_T.push_back(Eigen::Triplet<double>(i, i-1, -scaling));
		J_T.push_back(Eigen::Triplet<double>(i, i,  2*scaling));
		J_T.push_back(Eigen::Triplet<double>(i, i+1, -scaling));
	}
	/* Neumann condition at x = _length */
	J_T.push_back(Eigen::Triplet<double>(grid_end, grid_end, scaling));
	J_T.push_back(Eigen::Triplet<double>(grid_end, grid_end-1, -scaling));

	J = Eigen::SparseMatrix<double>(grid_end+1, grid_end+1);
	J.setFromTriplets(J_T.begin(), J_T.end());

	/* Measure the ion density */
	for (int i=0; i<_population_sizes.front(); i++)
	{
		double x 			= positions.front()->at(i);
		double w 			= weights.front()->at(i);
		int bin 			= _plasma->find_index_on_grid(x);

		double s = _plasma->find_position_in_cell(x);
		ion_density[bin  ] += w*(1. - s);
		ion_density[bin+1] += w*s;
	}
	ion_density[0] *= 2;
	ion_density[grid_end] *= 2;
	double n0 = 1./(_reference_densities.front()*_plasma->_dx);
	for (double & rho: ion_density) 
		rho *= n0;

	Potential = (1e-6 + Ion_density.array()).log();
	Exp_potential = Potential.array().exp();
	Residual = Exp_potential - Ion_density ;
	Residual(0) *= .5;  Residual(grid_end) *= .5;
	Residual += J*Potential;

	/* Newton's method loop */
	for (int count=0; count<iter_max; count++)
	{
		Eigen::SparseMatrix<double> Je = J;
		Je.coeffRef(0,0) += .5* Exp_potential(0); 
		for (int i=1; i<grid_end; i++) Je.coeffRef(i,i) += Exp_potential(i); 
		Je.coeffRef(grid_end, grid_end) += .5* Exp_potential(grid_end);
		
		solver.compute(Je);

		Update = solver.solve(Residual);

		Potential -= Update;
		Exp_potential = Potential.array().exp();
		Residual = Exp_potential - Ion_density ;
		Residual(0) *= .5;  Residual(grid_end) *= .5;
		Residual += J*Potential;
	
		if (Residual.squaredNorm() < tol2)
			break;
	}

	xs = 0.0;
	x0 = .5/_population_sizes.back();

	for (int i=0; i<_population_sizes.back(); i++)
	{
		double x 		= _plasma->_length*(xs + x0);
		int bin 		= _plasma->find_index_on_grid(x);
		double cellpos 	= _plasma->find_position_in_cell(x);

		positions.back()->at(i) = x;
		weights.back()->at(i) 	= Tools::EvaluateP1Function(exp_potential, bin, cellpos);

		/* bit-reversed scrambling to reduce the correlation between positions and velocities */
		double xsi = 0.5;
		xs -= 0.5;
		while (xs >= 0.0)
		{
			xsi *= 0.5;
			xs -= xsi;
		} 
		xs += 2.0*xsi;
	}

	w0 = _population_sizes.back() / std::accumulate(weights.back()->begin(), weights.back()->end(), 0.);
	for (double & w: *weights.back()) w *= w0;

	/***************************************************************************************/
	/* Finally we implement spatially constant temperature particle velocity distributions */
	/***************************************************************************************/
	for (int population_index=0; population_index < _number_of_populations; population_index++)
	{
		std::vector<double> * velocity 	= velocities.at(population_index);
		int population_size = _population_sizes.at(population_index);
		int nbins 			= _bin_numbers.at(population_index);

		double v0 = _mean_velocities.at(population_index);
		double vt1 = _random_mean_thermal_vel.at(population_index);
		double vt2 = _quiet_mean_thermal_vel.at(population_index);
		double quiet_start_exponent = _quiet_start_exponents.at(population_index);

		static std::vector<double> helper;
		helper.resize(nbins);

		for (double & v: *velocity)
			v = v0;

		if (vt2 != 0.0) 
		{
			double vmax = 5.0*vt2;
			double dv = 2.0*vmax/(static_cast<double>(nbins - 1));
			helper.front() = 0.0;

			double vv = -vmax/vt2;
			for (int i=1; i < helper.size(); i++)
			{
				if (quiet_start_exponent == 0)
				{
					helper.at(i) = helper.at(i-1) + std::exp(-0.5*vv*vv);
					vv += dv/vt2;
					helper.at(i) +=  std::exp(-0.5*vv*vv);
				}
				else
				{					
					double vvnv2 = std::pow(vv, static_cast<double>(quiet_start_exponent));
					helper.at(i) = helper.at(i-1) + vvnv2*std::exp(-0.5*vv*vv);
					vv += dv/vt2;
					vvnv2 = std::pow(vv, static_cast<double>(quiet_start_exponent));
					helper.at(i) +=  vvnv2*std::exp(-0.5*vv*vv);
				}
			}

			double df = helper.back()/population_size;

			auto it_helper = helper.begin();
			for (int i=0; i<population_size; i++)
			{
				double fv = (i + 0.5)*df;
				while (fv >= *(it_helper+1)) 
				{
					it_helper++;
					if (it_helper == helper.end()) 
						perror("distribution function error");
				}
				velocity->at(i) += dv* (std::distance(helper.begin(), it_helper)
										 + (fv - *it_helper)/(*(it_helper+1) - *it_helper)) - vmax;
			}
		}

		if (vt1 != 0.0)
			for (auto & vel_x : *velocity)
				vel_x 	+= vt1*RandomTools::Generate_randomly_normal(0.0, 1.0);
	}

}