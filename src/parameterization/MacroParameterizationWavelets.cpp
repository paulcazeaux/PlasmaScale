#include "parameterization/MacroParameterizationWavelets.h"


MacroParameterizationWavelets::MacroParameterizationWavelets(MacroParameterizationWavelets &&parameterization) :
			MacroParameterization(std::move(parameterization)),
			_grid_end(parameterization._grid_end),
			_macro_grid_end(parameterization._macro_grid_end),
			_depth(parameterization._depth),
			_cutoff(parameterization._cutoff),
			_record_microsteps(parameterization._record_microsteps),
			_distributions(std::move(parameterization._distributions)),
			_stack_ion_distribution(std::move(parameterization._stack_ion_distribution)),
			_stack_ion_number(std::move(parameterization._stack_ion_number)),
			_stack_electron_number(std::move(parameterization._stack_electron_number)),
			_stack_electron_energy(std::move(parameterization._stack_electron_energy)),
            _stack_index(parameterization._stack_index),
			_record_times(std::move(parameterization._record_times)),
			_record_ion_distribution(std::move(parameterization._record_ion_distribution)),
			_total_moment(parameterization._total_moment),
			_J(parameterization._J)
		{}

MacroParameterizationWavelets& MacroParameterizationWavelets::operator=(MacroParameterizationWavelets &&parameterization)
{
	if (this == &parameterization)
		return *this;
	
	MacroParameterization::operator=(std::move(parameterization));
	_grid_end = parameterization._plasma->_grid_end;
   	_macro_grid_end = parameterization._plasma->_macro_grid_end;
   	_depth = parameterization._depth;
   	_cutoff = parameterization._cutoff;
   	_record_microsteps = parameterization._plasma->_record_microsteps;

	_distributions = std::move(parameterization._distributions);
	_stack_ion_distribution = std::move(parameterization._stack_ion_distribution);
	_stack_ion_number = std::move(parameterization._stack_ion_number);
	_stack_electron_number = std::move(parameterization._stack_electron_number);
	_stack_electron_energy = std::move(parameterization._stack_electron_energy);
    _stack_index = parameterization._stack_index;

	_record_times = std::move(parameterization._record_times);
	_record_ion_distribution = std::move(parameterization._record_ion_distribution);

	_total_moment = parameterization._total_moment;
	_J = parameterization._J;
	return *this;
}

MacroParameterizationWavelets::MacroParameterizationWavelets(MacroParameterization & parameterization) :
	MacroParameterization(std::move(parameterization))
{
	if (_number_of_populations != 2)
		throw std::runtime_error("There are " + std::to_string(_number_of_populations) + " and not 2 populations as required.\n");

	_grid_end = _plasma->_grid_end;
	_macro_grid_end = _plasma->_macro_grid_end;
	_depth = _plasma->_depth;
	_cutoff = _plasma->_cutoff;
	_record_microsteps = _plasma->_record_microsteps;

	_distributions.clear();
	_distributions.emplace_back(new ActiveWaveletRepresentation(_plasma, _lower_velocities.at(0), _upper_velocities.at(0), _depth, _reference_densities.at(0)));
	_distributions.emplace_back(new ActiveMaxwellianRepresentation(_plasma, _reference_densities.at(1)));
	MaxwellianRepresentation::InitializeQuietStartArrays(std::pow(2, _depth));

	int number_of_microsteps = _plasma->_number_of_microsteps;
	_stack_ion_distribution.reserve(number_of_microsteps+1);
	_stack_ion_number.reserve(number_of_microsteps+1);
	_stack_electron_number.reserve(number_of_microsteps+1);
	_stack_electron_energy.reserve(number_of_microsteps+1);
    _stack_index = 0;

	if (_record_microsteps)
	{
        _record_ion_distribution.clear();
		_record_times.reserve(2*number_of_microsteps+2);
		for (int i = 0; i < 2*number_of_microsteps+2; i++)
		{
			_record_ion_distribution.emplace_back(_plasma, _lower_velocities.at(0), _upper_velocities.at(0), _depth, _reference_densities.at(0));
		}
	}

	/* Create the Jacobian matrix used for the self-consistent electron initialization */
	double debye_length = _thermal_velocities.back()/_plasma_pulsations.back();
	double scaling = std::pow(debye_length/_plasma->_dx, 2.);

	std::vector<Eigen::Triplet<double> > J_T;

		/* Neumann condition at x = 0 */
	J_T.push_back(Eigen::Triplet<double>(0,0, scaling));
	J_T.push_back(Eigen::Triplet<double>(0,1, -scaling));
		/* Bulk values */
	for (int i=1; i<_grid_end; i++)
	{
		J_T.push_back(Eigen::Triplet<double>(i, i-1, -scaling));
		J_T.push_back(Eigen::Triplet<double>(i, i,  2*scaling));
		J_T.push_back(Eigen::Triplet<double>(i, i+1, -scaling));
	}
		/* (inhomogeneous) Neumann condition at x = _length */
	J_T.push_back(Eigen::Triplet<double>(_grid_end, _grid_end, scaling));
	J_T.push_back(Eigen::Triplet<double>(_grid_end, _grid_end-1, -scaling));

	_J = Eigen::SparseMatrix<double>(_grid_end+1, _grid_end+1);
	_J.setFromTriplets(J_T.begin(), J_T.end());

}



void MacroParameterizationWavelets::Initialize(State & state)
{
	_stack_index = 0;
	this->RestrictAndPushback(state, 0.);
	_stack_ion_distribution.front().GetDensityVelocityPressure(_ion_density, _ion_velocity, _ion_pressure);
	this->CalculateTotalMoment(state);
}

void MacroParameterizationWavelets::Load(State & state) const
/* Fill the particle arrays to initialize the microscopic state */
{	
	/* Initialize the particle arrays */
	std::vector<std::vector<double> * > positions 		= state.get_vector_of_position_arrays();
	std::vector<std::vector<double> * > velocities	 	= state.get_vector_of_velocity_arrays();
	std::vector<std::vector<double> * > weights 		= state.get_vector_of_weight_arrays();

	assert(positions.size() == 2);
	assert(velocities.size() == 2);
	assert(weights.size() == 2);

	for (int population_index=0; population_index < 2; population_index++)
	{
		std::vector<double>::iterator position 		= positions.at(population_index)->begin();
		std::vector<double>::iterator velocity 		= velocities.at(population_index)->begin();
		std::vector<double>::iterator weight 		= weights.at(population_index)->begin();
		
		int population_size = _population_sizes.at(population_index);
		state.set_new_number_of_particles(population_index, population_size);
		_distributions.at(population_index)->Load(population_size, position, velocity, weight);
	}
	state.Prepare(true);
}

void MacroParameterizationWavelets::SetAccField(State & state)
{
	state.GetEField(_accfield);
	double e_to_acc_factor = _unit_charges.front()/_unit_masses.front() *_plasma->_dt*_plasma->_dt;
	for (auto & val : _accfield) 	val *= e_to_acc_factor;
}


void MacroParameterizationWavelets::RestrictAndPushback(const State & state, const double delay)
{
    if (_stack_index >= _stack_ion_distribution.size())
        _stack_ion_distribution.emplace_back(_plasma, _lower_velocities.at(0), _upper_velocities.at(0), _depth, _reference_densities.at(0));
	if (_stack_index >= _stack_ion_number.size())
		_stack_ion_number.emplace_back(0.);
	if (_stack_index >= _stack_electron_number.size())
		_stack_electron_number.emplace_back(0.);
	if (_stack_index >= _stack_electron_energy.size())
		_stack_electron_energy.emplace_back(0.);

	/************************************************************************/
	/* Now we weigh the particles and compute the moments on the macro grid */
	/************************************************************************/
    /* First, we initialize the arrays */
	std::vector<double>::iterator 	ion_position 		= state.get_vector_of_position_arrays().front()->begin();
	std::vector<double>::iterator	ion_velocity 		= state.get_vector_of_velocity_arrays().front()->begin();
	std::vector<double>::iterator 	ion_weight 			= state.get_vector_of_weight_arrays().front()->begin();
	int ion_population_size = state.get_number_of_particles(0);
    
    _stack_ion_distribution.at(_stack_index).set_grid_end(_macro_grid_end);
	if (delay != 0)
		_stack_ion_distribution.at(_stack_index).Weigh(ion_population_size, ion_position, ion_velocity, ion_weight, delay, _accfield);
	else
		_stack_ion_distribution.at(_stack_index).Weigh(ion_population_size, ion_position, ion_velocity, ion_weight);

	_stack_ion_number.at(_stack_index) = state.get_total_weight(0);
	_stack_electron_number.at(_stack_index) = state.get_total_weight(1);
	_stack_electron_energy.at(_stack_index) = state.get_kinetic_energy(1);	

	if (_record_microsteps)
	{
		_record_times.push_back(*state.get_simulation_time());
		int m = _record_times.size()-1;
		_record_ion_distribution.at(m) = _stack_ion_distribution.at(_stack_index);
	}
	_stack_index++;
}

void MacroParameterizationWavelets::Lift()
{
	int size = _grid_end+1;
	*dynamic_cast<ActiveWaveletRepresentation*>(_distributions.front().get()) = _stack_ion_distribution.front();

	/*-------------------------------------------------------------*/
	/*	  We lift all these quantities to the fine grid    */
	/*-------------------------------------------------------------*/

	while (_distributions.front()->get_grid_end() < _grid_end)	{_distributions.front()->Refine();}

	/*----------------------------------------------------------------*/
	/* 		Initialize the electron density on the fine grid: 	  */
	/* Newton's method solve for the electrostatic potential assuming */
	/* 			a Boltzmann distribution for the electrons 			  */
	/*----------------------------------------------------------------*/

	static std::vector<double>	ion_density = std::vector<double>(size),
							  	ion_velocity = std::vector<double>(size),
							 	exp_potential = std::vector<double>(size);
	static 	Eigen::UmfPackLU<Eigen::SparseMatrix<double> > PoissonSolver;

	Eigen::Map<Eigen::VectorXd> Ion_density(ion_density.data(), size);
	Eigen::Map<Eigen::VectorXd> Exp_potential(exp_potential.data(), size);
	Eigen::VectorXd Residual(size), Potential(size), Update(size);

	_distributions.front()->GetDensityVelocity(ion_density, ion_velocity);

	/* Calculate thermal velocity for the electrons from the projected kinetic energy */
	double vt2 =  2*_stack_electron_energy.front()/(_unit_masses.at(1) * _population_sizes.at(1));

	/* Approximate the correction from mean velocity effect using the ion density and velocity */
	double N = std::accumulate(ion_density.begin(), ion_density.end(), 0.);
	for (int i=0; i<_grid_end; i++)
	{
		vt2 -= ion_density.at(i)*std::pow(ion_velocity.at(i), 2.0)/N;
	}
	double thermal_velocity = std::sqrt(vt2);

	/* Calculate derivative at x = _length due to escaped particles */
	double end_deriv = _stack_ion_number.front()/_reference_densities.at(0) - _stack_electron_number.front()/_reference_densities.at(1);
	end_deriv /= _plasma->_dx;
	
	/* Initialize solver variables */
	Eigen::SparseMatrix<double> Jr = std::pow(thermal_velocity/_thermal_velocities.at(1), 2.0) * _J;
	Potential = (Ion_density.array()+1e-6).log();
	Exp_potential = Potential.array().exp();
	Residual = Exp_potential - Ion_density;
	Residual(0) *= .5;  Residual(size-1) *= .5;
	Residual(size-1) += end_deriv;
	Residual += Jr*Potential;

	/* Newton's method loop */
	double tol2 = 1e-20, iter_max = 100;
	for (int count=0; count<iter_max; count++)
	{
		Eigen::SparseMatrix<double> Je = Jr;
		Je.coeffRef(0,0) += .5* Exp_potential(0); 
		for (int i=1; i<size-1; i++) Je.coeffRef(i,i) += Exp_potential(i); 
		Je.coeffRef(size-1, size-1) += .5* Exp_potential(size-1);
		
		PoissonSolver.compute(Je);
		Update = PoissonSolver.solve(Residual);

		Potential -= Update;
		Exp_potential = Potential.array().exp();
		Residual = Exp_potential - Ion_density;
		Residual(0) *= .5;  Residual(size-1) *= .5;
		Residual(size-1) += end_deriv;
		Residual += Jr*Potential;
	
		if (count > 90)
			std::cout << "Newton iteration count: " << count << "\t|\tResidual norm:\t" << Residual.norm() << "\t|\tStep norm:\t" << Update.norm() << std::endl;
		if (Update.squaredNorm() < tol2)
			break;
	}

	/*-------------------------------------------------------------*/
	/* Deduce the electron density from the Boltzmann distribution */
	/*-------------------------------------------------------------*/
	_distributions.back()->SetAdiabaticValues(exp_potential, ion_velocity, thermal_velocity);

	/* Set the new sizes of each population */
	_population_sizes.at(0) = static_cast<int>(_stack_ion_number.front()+.5);
	_population_sizes.at(1) = static_cast<int>(_stack_electron_number.front()+.5);
}

void MacroParameterizationWavelets::Step(State & state)
{

	static int count_steps = 0;
	timeunit field_time(0), PIC_time(0), algebra_time(0), load_time(0), lift_time(0), pushback_time(0), cutoff_time(0);
	std::chrono::high_resolution_clock::time_point start, stop;

	int number_of_microsteps = _plasma->_number_of_microsteps;
	if (number_of_microsteps == 0)
	{
		/***************************/
		/* Case of full PIC method */
		/***************************/
		std::shared_ptr<double> simulation_time = state.get_simulation_time();
		double current_time = *simulation_time;

			start = std::chrono::high_resolution_clock::now();
		for (int i=0; i<_plasma->_macro_to_micro_dt_ratio; i++)	{state.Step();}

			stop = std::chrono::high_resolution_clock::now();
			PIC_time += std::chrono::duration_cast<timeunit>(stop - start);
		*simulation_time = current_time + _plasma->_macro_to_micro_dt_ratio*_plasma->_dt;
	}
	else if (number_of_microsteps == _plasma->_macro_to_micro_dt_ratio)
	{
		/********************************************************************************/
		/* Case of periodic electronic reset to a self-consistent Boltzmann equilibrium */
		/********************************************************************************/
		std::shared_ptr<double> simulation_time = state.get_simulation_time();
		double current_time = *simulation_time;
			start = std::chrono::high_resolution_clock::now();
		this->SetAccField(state);
			stop = std::chrono::high_resolution_clock::now();
			field_time += std::chrono::duration_cast<timeunit>(stop - start);

		if (_record_microsteps)
			_record_times.clear();
		_stack_index = 0;
			start = std::chrono::high_resolution_clock::now();
		for (int i=0; i<number_of_microsteps; i++)	{ state.Step(); }
			stop = std::chrono::high_resolution_clock::now();
			PIC_time += std::chrono::duration_cast<timeunit>(stop - start);

			start = std::chrono::high_resolution_clock::now();
		this->RestrictAndPushback(state, 0.);
			stop = std::chrono::high_resolution_clock::now();
			pushback_time += std::chrono::duration_cast<timeunit>(stop - start);

		_stack_index = 1;
			start = std::chrono::high_resolution_clock::now();
		this->Lift();
			stop = std::chrono::high_resolution_clock::now();
			lift_time += std::chrono::duration_cast<timeunit>(stop - start);

		/*		Here we load the electrons only, without using 		*/
		/*	the normal MacroParameterizationWavelets::Load function	*/

			start = std::chrono::high_resolution_clock::now();
		std::vector<double>::iterator position 		= state.get_vector_of_position_arrays().at(1)->begin();
		std::vector<double>::iterator velocity 		= state.get_vector_of_velocity_arrays().at(1)->begin();
		std::vector<double>::iterator weight 		= state.get_vector_of_weight_arrays().at(1)->begin();

		int population_size = static_cast<int>(_stack_electron_number.front()+.5);
		_population_sizes.at(1) = population_size;
		state.set_new_number_of_particles(1, population_size);

		_distributions.at(1)->Load(population_size, position, velocity, weight);
		state.Prepare(1, false);
			stop = std::chrono::high_resolution_clock::now();
			load_time += std::chrono::duration_cast<timeunit>(stop - start);

		*simulation_time = current_time + _plasma->_macro_to_micro_dt_ratio*_plasma->_dt;
	    if (_record_microsteps)
		{
			_record_ion_distribution.at(_record_times.size()) = _stack_ion_distribution.front();
			_record_times.push_back(*simulation_time);
		}

	}
	else
	{

		/***********************************************************/
		/* Case of Lagrangian equation-free projective integration */
		/***********************************************************/
		std::shared_ptr<double> simulation_time = state.get_simulation_time();
		double current_time = *simulation_time;
	    double ratio = static_cast<double>(_plasma->_macro_to_micro_dt_ratio - number_of_microsteps);
		static ActiveWaveletRepresentation slope;
		double slope_ion_number, slope_electron_number, slope_energy;
		/* Explicit Euler integration */
		/* Obtain the ion acceleration field for approximation of the characteristics */
			start = std::chrono::high_resolution_clock::now();
		this->SetAccField(state);
			stop = std::chrono::high_resolution_clock::now();
			field_time += std::chrono::duration_cast<timeunit>(stop - start);
		if (_record_microsteps)
			_record_times.clear();
		/* Now we advance while measuring the coarse data for a set number of microsteps */
	    _stack_index = 0;
				start = std::chrono::high_resolution_clock::now();
		this->RestrictAndPushback(state, ratio + number_of_microsteps);
				stop = std::chrono::high_resolution_clock::now();
				pushback_time += std::chrono::duration_cast<timeunit>(stop - start);
		for (int i=0; i<number_of_microsteps; i++)
		{
				start = std::chrono::high_resolution_clock::now();
			state.Step();
				stop = std::chrono::high_resolution_clock::now();
				PIC_time += std::chrono::duration_cast<timeunit>(stop - start);
				start = std::chrono::high_resolution_clock::now();
			this->RestrictAndPushback(state, ratio + number_of_microsteps-i-1);
				stop = std::chrono::high_resolution_clock::now();
				pushback_time += std::chrono::duration_cast<timeunit>(stop - start);
		}

		/**************************/
		/* Extrapolate using EFPI */
		/**************************/

		/* First we calculate the slope using the mean-square estimate */
			start = std::chrono::high_resolution_clock::now();
		double slope_factor = static_cast<double>(6*ratio)/static_cast<double>(number_of_microsteps*(number_of_microsteps+1)*(number_of_microsteps+2));
		slope 					= slope_factor*number_of_microsteps*(_stack_ion_distribution.at(number_of_microsteps)-_stack_ion_distribution.front());

		slope_ion_number 		= slope_factor*number_of_microsteps*(_stack_ion_number.at(number_of_microsteps)-_stack_ion_number.front());
		slope_electron_number 	= slope_factor*number_of_microsteps*(_stack_electron_number.at(number_of_microsteps)-_stack_electron_number.front());
		slope_energy 			= slope_factor*number_of_microsteps*(_stack_electron_energy.at(number_of_microsteps)-_stack_electron_energy.front());

		for (int i=1; 2*i<number_of_microsteps; i++)
		{
			slope 					+= slope_factor*(number_of_microsteps-2*i)*(_stack_ion_distribution.at(number_of_microsteps-i)-_stack_ion_distribution.at(i));

			slope_ion_number 		+= slope_factor*(number_of_microsteps-2*i)*(_stack_ion_number.at(number_of_microsteps-i)-_stack_ion_number.at(i));
			slope_electron_number 	+= slope_factor*(number_of_microsteps-2*i)*(_stack_electron_number.at(number_of_microsteps-i)-_stack_electron_number.at(i));
			slope_energy 			+= slope_factor*(number_of_microsteps-2*i)*(_stack_electron_energy.at(number_of_microsteps-i)-_stack_electron_energy.at(i));
		}

		/* We take a projective Euler step following this slope */
	    std::swap(_stack_ion_distribution.front(), _stack_ion_distribution.at(number_of_microsteps));
		_stack_ion_distribution.front() += slope;

		_stack_ion_number.front() 		= _stack_ion_number.at(number_of_microsteps) + slope_ion_number;
		_stack_electron_number.front() 	= _stack_electron_number.at(number_of_microsteps) + slope_electron_number;
		_stack_electron_energy.front() 	= _stack_electron_energy.at(number_of_microsteps) + slope_energy;
			stop = std::chrono::high_resolution_clock::now();
			algebra_time += std::chrono::duration_cast<timeunit>(stop - start);

		/* We denoise the resulting distribution using a fixed wavelet resolution cutoff */
			start = std::chrono::high_resolution_clock::now();
	    _stack_ion_distribution.front().Cutoff(_cutoff);
			stop = std::chrono::high_resolution_clock::now();
			cutoff_time += std::chrono::duration_cast<timeunit>(stop - start);

		/* Finally we lift the distribution to the fine grid and we load the particles */
		_stack_index = 1;
			start = std::chrono::high_resolution_clock::now();
		this->Lift();
			stop = std::chrono::high_resolution_clock::now();
			lift_time += std::chrono::duration_cast<timeunit>(stop - start);

			start = std::chrono::high_resolution_clock::now();
		this->Load(state);
			stop = std::chrono::high_resolution_clock::now();
			load_time += std::chrono::duration_cast<timeunit>(stop - start);

		*simulation_time = current_time + _plasma->_macro_to_micro_dt_ratio*_plasma->_dt;
	    if (_record_microsteps)
		{
			_record_ion_distribution.at(_record_times.size()) = _stack_ion_distribution.front();
			_record_times.push_back(*simulation_time);
		}
	}

	_stack_index = 0;
	this->RestrictAndPushback(state, 0.);
	_stack_ion_distribution.front().GetDensityVelocityPressure(_ion_density, _ion_velocity, _ion_pressure);
	/* Calculate thermal velocity for the electrons from the projected kinetic energy */
	double vt2 =  2*_stack_electron_energy.front()/(_unit_masses.at(1) * _population_sizes.at(1));

	/* Approximate the correction from mean velocity effect using the ion density and velocity */
	double N = 0.5*(_ion_density.front() + _ion_density.back()) + std::accumulate(_ion_density.begin()+1, _ion_density.end()-1, 0.);
	vt2 -= 0.5/N*(_ion_density.front()*std::pow(_ion_velocity.front(), 2.0) + _ion_density.back()*std::pow(_ion_velocity.back(), 2.0));
	for (int i=1; i<_macro_grid_end; i++)
	{
		vt2 -= _ion_density.at(i)*std::pow(_ion_velocity.at(i), 2.0)/N;
	}

	/* Diagnostics: output timing info */
	if ((count_steps % 100) == 0)
		std::cout << "|\tStep\t|\t Total\t|\t   PIC\t|\t Field\t|"
						<< "\tPushback\t|\tAlgebra\t|\tCutoff\t|"
						<< "\t  Load\t|\t  Lift\t||\t Thermal vel\t||"
						<< "\tGlob. charge\t|" << std::endl;
	std::cout << "|\t" 	<< std::setw(4) << ++count_steps  
		<< "\t|\t" 		<< std::setw(6) << (PIC_time + field_time + pushback_time + algebra_time + cutoff_time + load_time + lift_time).count() 
		<< "\t|\t" 		<< std::setw(6) << PIC_time.count() 
		<< "\t|\t" 		<< std::setw(6) << field_time.count() 
		<< "\t|\t" 		<< std::setw(8) << pushback_time.count() 
		<< "\t|\t" 		<< std::setw(7) << algebra_time.count() 
		<< "\t|\t" 		<< std::setw(6) << cutoff_time.count()  
		<< "\t|\t" 		<< std::setw(6) << load_time.count() 
		<< "\t|\t" 		<< std::setw(6) << lift_time.count() 
		<< "\t||\t" 	<< std::sqrt(vt2) 
		<< "\t||\t" 	<< _stack_ion_number.front()/_reference_densities.at(0) - _stack_electron_number.front()/_reference_densities.at(1) 
		<< "\t|" 		<< std::endl;
}


void MacroParameterizationWavelets::WriteData(State & state, std::fstream & fout)
{
	const int downsampling = 0;

	_stack_index = 0;
	this->RestrictAndPushback(state, 0.);
	_stack_ion_distribution.front().GetDensityVelocityPressure(_ion_density, _ion_velocity, _ion_pressure);

		/* Diagnostic: extract electron parameters */
	std::vector<double>::iterator 	electron_position 		= state.get_vector_of_position_arrays().back()->begin();
	std::vector<double>::iterator	electron_velocity 		= state.get_vector_of_velocity_arrays().back()->begin();
	std::vector<double>::iterator 	electron_weight 		= state.get_vector_of_weight_arrays().back()->begin();
	int electron_population_size = state.get_number_of_particles(1);
    _distributions.back()->set_grid_end(_macro_grid_end);
	_distributions.back()->Weigh(electron_population_size, electron_position, electron_velocity, electron_weight);

	if (!_record_microsteps)
	{
		_stack_ion_distribution.front().DWT(downsampling);
		fout << _stack_ion_distribution.front(); 
		_stack_ion_distribution.front().iDWT();
		fout << *_distributions.back();
	}
	else
	{
		for (int i=0; i<_record_times.size(); i++)
		{
			fout << "t = " << _record_times.at(i) << std::endl;
			_record_ion_distribution.at(i).DWT(downsampling);
			fout << _record_ion_distribution.at(i);
		}
	}
}


void MacroParameterizationWavelets::SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics)
{
	double * x_array 	= _plasma->get_x_grid_ptr();
	int * grid_end 	= _plasma->get_grid_end_ptr();

	diagnostics.emplace_back(new CurveDiagnostic(
				"linlin", "X", "Ion density", 0, 340));
	diagnostics.back()->AddData(x_array, _ion_density.data(), grid_end, 2);

	diagnostics.emplace_back(new CurveDiagnostic(
				"linlin", "X", "Ion velocity", 410, 340));
	diagnostics.back()->AddData(x_array, _ion_velocity.data(), grid_end, 2);

	diagnostics.emplace_back(new CurveDiagnostic(
				"linlin", "X", "Ion pressure", 820, 340));
	diagnostics.back()->AddData(x_array, _ion_pressure.data(), grid_end, 4);
}


void MacroParameterizationWavelets::CalculateTotalMoment(const State & state)
{
	_total_moment = 0;
	double dt = _plasma->_dt;
	std::vector<std::vector<double> * > velocities	 	= state.get_vector_of_velocity_arrays();
	std::vector<std::vector<double> * > weights 		= state.get_vector_of_weight_arrays();

	assert(velocities.size() == 2);
	assert(weights.size() == 2);

	for (int population_index=0; population_index < 2; population_index++)
	{
		std::vector<double>::iterator velocity 		= velocities.at(population_index)->begin();
		std::vector<double>::iterator weight 		= weights.at(population_index)->begin();
		int population_size = state.get_number_of_particles(population_index);
		double m = _unit_masses.at(population_index)/dt;

		for (int i=0; i < population_size; i++) 
		{
			double v = *(velocity+i);
			double w = *(weight+i);
			_total_moment += m*w*v;
		}
	}
}

void MacroParameterizationWavelets::CalculateIonMoment(const State & state, double & ion_particle_moment, double & ion_distr_moment)
{
	ion_particle_moment = 0;
	std::vector<std::vector<double> * > velocities	 	= state.get_vector_of_velocity_arrays();
	std::vector<std::vector<double> * > weights 		= state.get_vector_of_weight_arrays();

	assert(velocities.size() == 2);
	assert(weights.size() == 2);

	std::vector<double>::iterator velocity 		= velocities.front()->begin();
	std::vector<double>::iterator weight 		= weights.front()->begin();
	int population_size = state.get_number_of_particles(0);

	for (int i=0; i < population_size; i++) 
	{
		double v = *(velocity+i);
		double w = *(weight+i);
		ion_particle_moment += w*v;
	}
	ion_particle_moment *= this->get_unit_mass(0)/_plasma->_dt;
	_stack_ion_distribution.at(_stack_index-1).GetVelocityMoment(ion_distr_moment);
	ion_distr_moment *= this->get_unit_mass(0)*population_size;
}

void MacroParameterizationWavelets::CalculateElectronMoment(const State & state, double & electron_particle_moment)
{
	electron_particle_moment = 0;
	std::vector<std::vector<double> * > velocities 	= state.get_vector_of_velocity_arrays();
	std::vector<std::vector<double> * > weights 		= state.get_vector_of_weight_arrays();

	assert(velocities.size() == 2);
	assert(weights.size() == 2);

	std::vector<double>::iterator velocity 		= velocities.back()->begin();
	std::vector<double>::iterator weight 		= weights.back()->begin();
	int population_size = state.get_number_of_particles(1);
	for (int i=0; i < population_size; i++) 
	{
		double v = *(velocity+i);
		double w = *(weight+i);
		electron_particle_moment += w*v;
	}
	electron_particle_moment *= this->get_unit_mass(1)/_plasma->_dt;
}