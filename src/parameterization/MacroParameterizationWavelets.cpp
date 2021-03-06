#include "parameterization/MacroParameterizationWavelets.h"


MacroParameterizationWavelets::MacroParameterizationWavelets(MacroParameterizationWavelets &&parameterization) :
			MacroParameterization(std::move(parameterization)),
			_grid_size(parameterization._grid_size),
			_macro_grid_size(parameterization._macro_grid_size),
			_depth(parameterization._depth),
			_cutoff(parameterization._cutoff),
			_ion_vmax(parameterization._ion_vmax),
			_record_microsteps(parameterization._record_microsteps),
			_distributions(std::move(parameterization._distributions)),
			_stack_ion_distribution(std::move(parameterization._stack_ion_distribution)),
            _stack_index(parameterization._stack_index),
			_current_step_ion_distribution(std::move(parameterization._current_step_ion_distribution)),
			_record_times(std::move(parameterization._record_times)),
			_record_ion_distribution(std::move(parameterization._record_ion_distribution)),
			_electron_thermal_vel(parameterization._electron_thermal_vel),
			_debye_scaling(parameterization._debye_scaling),
			_total_moment(parameterization._total_moment)
		{}

MacroParameterizationWavelets& MacroParameterizationWavelets::operator=(MacroParameterizationWavelets &&parameterization)
{
	if (this == &parameterization)
		return *this;
	
	MacroParameterization::operator=(std::move(parameterization));
	_grid_size = parameterization._plasma->get_grid_size();
   	_macro_grid_size = parameterization._plasma->get_macro_grid_size();
   	_depth = parameterization._depth;
   	_cutoff = parameterization._cutoff;
   	_ion_vmax = parameterization._ion_vmax;
   	_record_microsteps = parameterization._plasma->get_record_microsteps();

	_distributions = std::move(parameterization._distributions);
	_stack_ion_distribution = std::move(parameterization._stack_ion_distribution);
    _stack_index = parameterization._stack_index;
	_current_step_ion_distribution = std::move(parameterization._current_step_ion_distribution);

	_record_times = std::move(parameterization._record_times);
	_record_ion_distribution = std::move(parameterization._record_ion_distribution);


	_electron_thermal_vel = parameterization._electron_thermal_vel;
	_debye_scaling = parameterization._debye_scaling;
	_total_moment = parameterization._total_moment;
	return *this;
}

MacroParameterizationWavelets::MacroParameterizationWavelets(MacroParameterization & parameterization,
				double electron_thermal_vel, double ion_vmax) :
	MacroParameterization(std::move(parameterization)), _ion_vmax(ion_vmax), _electron_thermal_vel(electron_thermal_vel)
{
	if (_number_of_populations != 2)
		throw std::runtime_error("There are " + std::to_string(_number_of_populations) + " and not 2 populations as required.\n");

	_grid_size = _plasma->get_grid_size();
	_macro_grid_size = _plasma->get_macro_grid_size();
	_depth = _plasma->get_wavelet_depth();
	_cutoff = _plasma->get_wavelet_cutoff();
	_record_microsteps = _plasma->get_record_microsteps();

	_current_step_ion_distribution	= ActiveWaveletRepresentation(_plasma, _ion_vmax, _depth, _macro_grid_size);

	_distributions.clear();
	_distributions.emplace_back(new ActiveWaveletRepresentation(_plasma, _ion_vmax, _depth, _grid_size));
	_distributions.emplace_back(new ActiveMaxwellianRepresentation(_plasma, _grid_size));

	int number_of_microsteps = _plasma->get_number_of_microsteps();
	_stack_ion_distribution.reserve(number_of_microsteps+1);
    _stack_index = 0;

	if (_record_microsteps)
	{
        _record_ion_distribution.clear();
		_record_times.reserve(2*number_of_microsteps+4);
		for (int i = 0; i < 2*number_of_microsteps+4; i++)
		{
			_record_ion_distribution.emplace_back(_plasma, _ion_vmax, _depth, _macro_grid_size);
		}
	}

	MaxwellianRepresentation::InitializeQuietStartArrays(std::pow(2, _depth));
	_debye_scaling = std::pow(_electron_thermal_vel/_plasma_pulsations.back(), 2.);
	_total_moment = 0;
}

void MacroParameterizationWavelets::Initialize(State & state)
{
	this->CalculateTotalMoment(state);
	this->SetAccField(state);
	this->RestrictAndPushback(state, 0.);
	this->Lift();
	this->Load(state);
    _stack_ion_distribution.front().GetDensityVelocityPressure(_ion_density, _ion_velocity, _ion_pressure);
}


void MacroParameterizationWavelets::Load(State & state) const
/* Fill the particle arrays to initialize the microscopic state */
{	
	/* Initialize the particle arrays */
	std::vector<std::vector<double> * > positions 		= state.get_vector_of_position_arrays();
	std::vector<std::vector<double> * > x_velocities 	= state.get_vector_of_x_velocity_arrays();
	std::vector<std::vector<double> * > y_velocities 	= state.get_vector_of_y_velocity_arrays();
	std::vector<std::vector<double> * > weights 		= state.get_vector_of_weight_arrays();

	assert(positions.size() == 2);
	assert(x_velocities.size() == 2);
	assert(y_velocities.size() == 2);
	assert(weights.size() == 2);

	for (int population_index=0; population_index < 2; population_index++)
	{
		std::vector<double>::iterator position 		= positions.at(population_index)->begin();
		std::vector<double>::iterator velocity_x 	= x_velocities.at(population_index)->begin();
		std::vector<double>::iterator weight 		= weights.at(population_index)->begin();
		int population_size = this->get_population_size(population_index);
		_distributions.at(population_index)->Load(population_size, position, velocity_x, weight);

		if (_cyclotronic_rotation_parameters.at(population_index) != 0.) 
		{
			std::vector<double>::iterator velocity_y 	= y_velocities.at(population_index)->begin();
			for (int i=0; i < population_size; i++) 
			{
				double v = *(velocity_x+i);
				double theta 	= RandomTools::Generate_randomly_uniform(0., 2.*M_PI);
				velocity_x[i] = v*std::cos(theta);
				velocity_y[i] = v*std::sin(theta);
			}
		}
	}
	state.Prepare(true);
}

void MacroParameterizationWavelets::SetAccField(State & state)
{
	state.GetEField(_accfield);
	double e_to_acc_factor = _unit_charges.front()/_unit_masses.front() *_plasma->get_dt()*_plasma->get_dt();
	for (auto & val : _accfield)
	{
		val *= e_to_acc_factor;
	}
}

void MacroParameterizationWavelets::CalculateTotalMoment(const State & state)
{
	_total_moment = 0;
	double dt = _plasma->get_dt();
	std::vector<std::vector<double> * > x_velocities 	= state.get_vector_of_x_velocity_arrays();
	std::vector<std::vector<double> * > weights 		= state.get_vector_of_weight_arrays();

	assert(x_velocities.size() == 2);
	assert(weights.size() == 2);

	for (int population_index=0; population_index < 2; population_index++)
	{
		std::vector<double>::iterator velocity_x 	= x_velocities.at(population_index)->begin();
		std::vector<double>::iterator weight 		= weights.at(population_index)->begin();
		int population_size = this->get_population_size(population_index);
		double m = _unit_masses.at(population_index)/dt;

		for (int i=0; i < population_size; i++) 
			{
				double v = *(velocity_x+i);
				double w = *(weight+i);
				_total_moment += m*w*v;
			}
	}
}

void MacroParameterizationWavelets::CalculateIonMoment(const State & state, double & ion_particle_moment, double & ion_distr_moment)
{
	ion_particle_moment = 0;
	std::vector<std::vector<double> * > x_velocities 	= state.get_vector_of_x_velocity_arrays();
	std::vector<std::vector<double> * > weights 		= state.get_vector_of_weight_arrays();

	assert(x_velocities.size() == 2);
	assert(weights.size() == 2);

	std::vector<double>::iterator velocity_x 	= x_velocities.front()->begin();
	std::vector<double>::iterator weight 		= weights.front()->begin();
	int population_size = this->get_population_size(0);

	for (int i=0; i < population_size; i++) 
	{
		double v = *(velocity_x+i);
		double w = *(weight+i);
		ion_particle_moment += w*v;
	}
	ion_particle_moment *= this->get_unit_mass(0)/_plasma->get_dt();
	_stack_ion_distribution.at(_stack_index-1).GetVelocityMoment(ion_distr_moment);
	ion_distr_moment *= this->get_unit_mass(0)*population_size;
}

void MacroParameterizationWavelets::CalculateElectronMoment(const State & state, double & electron_particle_moment)
{
	electron_particle_moment = 0;
	std::vector<std::vector<double> * > x_velocities 	= state.get_vector_of_x_velocity_arrays();
	std::vector<std::vector<double> * > weights 		= state.get_vector_of_weight_arrays();

	assert(x_velocities.size() == 2);
	assert(weights.size() == 2);

	std::vector<double>::iterator velocity_x 	= x_velocities.back()->begin();
	std::vector<double>::iterator weight 		= weights.back()->begin();
	int population_size = this->get_population_size(1);
	for (int i=0; i < population_size; i++) 
	{
		double v = *(velocity_x+i);
		double w = *(weight+i);
		electron_particle_moment += w*v;
	}
	electron_particle_moment *= this->get_unit_mass(1)/_plasma->get_dt();
}

void MacroParameterizationWavelets::RestrictAndPushback(const State & state, const double delay)
{
	/***********************************************************************/
	/* Now we weigh the particles and compute the moments on the fine grid */
	/***********************************************************************/
    /* First, we initialize the arrays */
	std::vector<double>::iterator 	ion_position 		= state.get_vector_of_position_arrays().front()->begin();
	std::vector<double>::iterator	ion_velocity 		= state.get_vector_of_x_velocity_arrays().front()->begin();
	std::vector<double>::iterator 	ion_weight 			= state.get_vector_of_weight_arrays().front()->begin();
	int ion_population_size = state.get_vector_of_position_arrays().front()->size();
    if (_stack_index >= _stack_ion_distribution.size())
        _stack_ion_distribution.emplace_back(_plasma, _ion_vmax, _depth, _grid_size);
    
    int size = _stack_ion_distribution.at(_stack_index).get_grid_size();
    
    if (size != _grid_size)
    {
        _stack_ion_distribution.at(_stack_index).set_grid_size(_grid_size);
        size = _grid_size;
    }
    
	/* Then we weigh the particles and restrict the values to the macroscopic grid using a linear smoothing. */
	_stack_ion_distribution.at(_stack_index).Weigh(ion_population_size, ion_position, ion_velocity, ion_weight, delay, _accfield);

	while (size > _macro_grid_size)
	{
		_stack_ion_distribution.at(_stack_index).Coarsen();
		size = _stack_ion_distribution.at(_stack_index).get_grid_size();
	}
	assert(size == _macro_grid_size);


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
	int size = _macro_grid_size;
	*dynamic_cast<ActiveWaveletRepresentation*>(_distributions.front().get()) = _stack_ion_distribution.front();


	/*----------------------------------------------------------------*/
	/* 		Initialize the electron density on the coarse grid 		  */
	/*----------------------------------------------------------------*/

	static std::vector<double>	ion_density = std::vector<double>(size),
							  	ion_velocity = std::vector<double>(size),
							 	potential = std::vector<double>(size),
							 	exp_potential = std::vector<double>(size), 
							 	update = std::vector<double>(size),
							 	electron_velocity = std::vector<double>(size);
	static std::vector<double> residual = std::vector<double>(size),
								conj_dir = std::vector<double>(size),
								dir = std::vector<double>(size),
							 	rhs = std::vector<double>(size);

	/*----------------------------------------------------------------*/
	/* Newton's method solve for the electrostatic potential assuming */
	/* 			a Boltzmann distribution for the electrons 			  */
	/*----------------------------------------------------------------*/

	double tol2 = 1e-12, iter_max = 20;
	_distributions.front()->GetDensityVelocity(ion_density, ion_velocity);

	/* Initialization */
	for (int i=0; i<size; i++)
	{
		potential.at(i) = std::log(ion_density.at(i));
	}
	/* Newton's method loop */
	double scaling = _debye_scaling/std::pow(_plasma->get_macro_dx(), 2.);
	for (int count=0; count<iter_max; count++)
	{
		/* Inner solve : CG algorithm */

		for (int i=0; i<size; i++)
		{
			exp_potential.at(i) = std::exp(potential.at(i));
				/* Initial value for the CG algorithm : zero. */
			residual.at(i)	=  
				  2./3. * exp_potential.at(i) + 1./6. * (exp_potential.at((i>0?i-1:size-1)) + exp_potential.at((i<size-1?i+1:0)))
				+ scaling * ( 2*potential.at(i) - potential.at((i>0?i-1:size-1)) - potential.at((i<size-1?i+1:0)) )
				- 2./3. * ion_density.at(i) - 1./6. * (ion_density.at((i>0?i-1:size-1)) + ion_density.at((i<size-1?i+1:0)));
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
				dir.at(i) =			(2./3.*exp_potential.at(i)+2*scaling)*conj_dir.at(i) 
										+ (1/6.*exp_potential.at((i>0?i-1:size-1)) - scaling) * conj_dir.at((i>0?i-1:size-1)) 
										+ (1/6.*exp_potential.at((i<size-1?i+1:0)) - scaling) * conj_dir.at((i<size-1?i+1:0));
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
			break;
	}

	/*---------------------------------------*/
	/* Now compute the electron density flux */
	/*---------------------------------------*/

	static const double total_moment = _total_moment;
	double mi = this->get_unit_mass(0)*this->get_population_size(0), me = this->get_unit_mass(1)*this->get_population_size(1);
	double ion_moment = 0;
	for (int i=0; i<size; i++)
	{
		ion_moment += ion_density.at(i)*ion_velocity.at(i);
	}
	ion_moment /= size;

	/* Assemble right-hand side */
	for (int i=0; i<size; i++)
	{
		rhs.at(i) = (ion_density.at(i)*ion_velocity.at(i) + total_moment/me - (1 + mi/me)*ion_moment);
	}
	/* Assemble matrix diagonals */
	static std::vector<double> D = std::vector<double>(size), D1 = std::vector<double>(size);
	for (int i=0; i<size; i++)
	{
		D.at(i) =  2./3. + scaling * ( std::exp(-.5*(potential.at(i)+potential.at((i>0?i-1:size-1)))) 
								   	  + std::exp(-.5*(potential.at(i)+potential.at((i<size-1?i+1:0)))));
		D1.at(i) = 1./6. - scaling * std::exp(-.5*(potential.at(i)+potential.at((i<size-1?i+1:0))));
	}

	double rnorm2 = 0, delta, alpha, beta;
	for (int i=0; i<size; i++)
	{
		exp_potential.at(i) = std::exp(potential.at(i));
			/* Initial value for the CG algorithm : zero. */
		residual.at(i)	= 2./3.*rhs.at(i) + 1./6.*(rhs.at((i>0?i-1:size-1)) + rhs.at((i<size-1?i+1:0)));
		conj_dir.at(i) = residual.at(i);
		rnorm2 += residual.at(i)*residual.at(i);
	}
	std::fill(electron_velocity.begin(), electron_velocity.end(), 0.);


	/* Start conjugate gradient loop */
	for (int count_cg=0; count_cg<iter_max; count_cg++)
	{
		/* Compute the matrix-vector product */
		delta = 0;
		for (int i=0; i<size; i++)
		{
			dir.at(i) = D.at(i)*conj_dir.at(i)	+ D1.at(i)*conj_dir.at((i<size-1?i+1:0))
												+ D1.at((i>0?i-1:size-1))*conj_dir.at((i>0?i-1:size-1));
			delta += dir.at(i) * conj_dir.at(i);
		}

		beta = 1./rnorm2;
		alpha = rnorm2/delta;
		rnorm2 = 0;
		for (int i=0; i<size; i++)
		{
			electron_velocity.at(i) += alpha * conj_dir.at(i);
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

	/*-------------------------------------------------------------*/
	/* Deduce the electron density from the Boltzmann distribution */
	/*-------------------------------------------------------------*/
	static std::vector<double>	electron_thermal_vel = std::vector<double>(size);

	for (int i=0; i<size; i++)
	{
		exp_potential.at(i) = std::exp(potential.at(i));
		electron_velocity.at(i) = electron_velocity.at(i)/exp_potential.at(i);
	}

	_distributions.back()->SetAdiabaticValues(exp_potential, electron_velocity, _electron_thermal_vel);

	/*-------------------------------------------------------------*/
	/*	  Finally we lift all these quantities to the fine grid    */
	/*-------------------------------------------------------------*/

	while (size < _grid_size)
	{
		_distributions.front()->Refine();
		_distributions.back()->Refine();
		size *= 2;
	}
	assert(size == _grid_size);
}

void MacroParameterizationWavelets::Step(State & state)
{
	std::cout.precision(8);

	static int count_steps = 0;
	timeunit field_time(0), PIC_time(0), algebra_time(0), load_time(0), lift_time(0), pushback_time(0), cutoff_time(0);
	std::chrono::high_resolution_clock::time_point start, stop;

	int number_of_microsteps = _plasma->get_number_of_microsteps();
	if (number_of_microsteps == 0)
	{		
		std::shared_ptr<double> simulation_time = state.get_simulation_time();
		double current_time = *simulation_time;

			start = std::chrono::high_resolution_clock::now();
		this->SetAccField(state);
			stop = std::chrono::high_resolution_clock::now();
			field_time += std::chrono::duration_cast<timeunit>(stop - start);

			start = std::chrono::high_resolution_clock::now();
		for (int i=0; i<_plasma->get_macro_to_micro_dt_ratio(); i++)
			state.Step();

			stop = std::chrono::high_resolution_clock::now();
			PIC_time += std::chrono::duration_cast<timeunit>(stop - start);
			start = std::chrono::high_resolution_clock::now();
		_stack_index = 0;
		this->RestrictAndPushback(state, 0.);
			stop = std::chrono::high_resolution_clock::now();
			pushback_time += std::chrono::duration_cast<timeunit>(stop - start);
			start = std::chrono::high_resolution_clock::now();
	    _stack_ion_distribution.front().Cutoff(_cutoff);
			stop = std::chrono::high_resolution_clock::now();
			cutoff_time += std::chrono::duration_cast<timeunit>(stop - start);
			start = std::chrono::high_resolution_clock::now();
		this->Lift();
			stop = std::chrono::high_resolution_clock::now();
			lift_time += std::chrono::duration_cast<timeunit>(stop - start);
		*simulation_time = current_time + _plasma->get_macro_to_micro_dt_ratio()*_plasma->get_dt();
	}
	else if (number_of_microsteps == _plasma->get_macro_to_micro_dt_ratio())
	{

		std::shared_ptr<double> simulation_time = state.get_simulation_time();
		double current_time = *simulation_time;
			start = std::chrono::high_resolution_clock::now();
		this->SetAccField(state);
			stop = std::chrono::high_resolution_clock::now();
			field_time += std::chrono::duration_cast<timeunit>(stop - start);

		if (_record_microsteps)
			_record_times.clear();
		_stack_index = 0;
		for (int i=0; i<number_of_microsteps; i++)
		{
				start = std::chrono::high_resolution_clock::now();
			state.Step();
				stop = std::chrono::high_resolution_clock::now();
				PIC_time += std::chrono::duration_cast<timeunit>(stop - start);
		}
			start = std::chrono::high_resolution_clock::now();
		this->RestrictAndPushback(state, 0.);
			stop = std::chrono::high_resolution_clock::now();
			pushback_time += std::chrono::duration_cast<timeunit>(stop - start);
			start = std::chrono::high_resolution_clock::now();
	    _stack_ion_distribution.front().Cutoff(_cutoff);
			stop = std::chrono::high_resolution_clock::now();
			cutoff_time += std::chrono::duration_cast<timeunit>(stop - start);
		_stack_index = 1;
			start = std::chrono::high_resolution_clock::now();
		this->Lift();
			stop = std::chrono::high_resolution_clock::now();
			lift_time += std::chrono::duration_cast<timeunit>(stop - start);

		/* Here we load the electrons only, without using the normal MacroParameterizationWavelets::Load function */
			start = std::chrono::high_resolution_clock::now();
		std::vector<double>::iterator position 		= state.get_vector_of_position_arrays().at(1)->begin();
		std::vector<double>::iterator velocity_x 	= state.get_vector_of_x_velocity_arrays().at(1)->begin();
		std::vector<double>::iterator velocity_y 	= state.get_vector_of_y_velocity_arrays().at(1)->begin();
		std::vector<double>::iterator weight 		= state.get_vector_of_weight_arrays().at(1)->begin();
		int population_size = this->get_population_size(1);

		_distributions.at(1)->Load(population_size, position, velocity_x, weight);
		if (_cyclotronic_rotation_parameters.at(1) != 0.) 
		{
			for (int i=0; i < population_size; i++) 
			{
				double v = *(velocity_x+i);
				double theta = RandomTools::Generate_randomly_uniform(0., 2.*M_PI);
				velocity_x[i] = v*std::cos(theta);
				velocity_y[i] = v*std::sin(theta);
			}
		}
		state.Prepare(1, false);
			stop = std::chrono::high_resolution_clock::now();
			load_time += std::chrono::duration_cast<timeunit>(stop - start);

		*simulation_time = current_time + _plasma->get_macro_to_micro_dt_ratio()*_plasma->get_dt();
	    if (_record_microsteps)
		{
			_record_ion_distribution.at(_record_times.size()) = _stack_ion_distribution.front();
			_record_times.push_back(*simulation_time);
		}
	}
	else
	{
		std::shared_ptr<double> simulation_time = state.get_simulation_time();
		double current_time = *simulation_time;
	    double ratio = static_cast<double>(_plasma->get_macro_to_micro_dt_ratio() - number_of_microsteps);

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
		static ActiveWaveletRepresentation first_slope = ActiveWaveletRepresentation(_plasma, _ion_vmax, _depth, _macro_grid_size);
		double slope_factor = static_cast<double>(6*ratio)/static_cast<double>(number_of_microsteps*(number_of_microsteps+1)*(number_of_microsteps+2));
		first_slope = slope_factor*number_of_microsteps*(_stack_ion_distribution.at(number_of_microsteps)-_stack_ion_distribution.front());
		for (int i=1; 2*i<number_of_microsteps; i++)
			first_slope += slope_factor*(number_of_microsteps-2*i)*(_stack_ion_distribution.at(number_of_microsteps-i)-_stack_ion_distribution.at(i));

		/* We take a projective Euler step following this slope */
	    std::swap(_stack_ion_distribution.front(), _stack_ion_distribution.at(number_of_microsteps));
		_stack_ion_distribution.front() += first_slope;
			stop = std::chrono::high_resolution_clock::now();
			algebra_time += std::chrono::duration_cast<timeunit>(stop - start);

		/* We denoise the resulting distribution using a fixed wavelet resolution cutoff */
			start = std::chrono::high_resolution_clock::now();
	    _stack_ion_distribution.front().Cutoff(_cutoff);
			stop = std::chrono::high_resolution_clock::now();
			cutoff_time += std::chrono::duration_cast<timeunit>(stop - start);

		/* Finally we lift the distribution to the fine grid and we load the particles */
		_stack_index = 1;
			double ion_particle_moment, ion_distr_moment, electron_moment; 
			this->CalculateIonMoment(state, ion_particle_moment, ion_distr_moment);
			this->CalculateElectronMoment(state, electron_moment);
			std::cout << ++count_steps <<  " | Final particle ion moment: " << ion_particle_moment << "; Distribution ion moment: " << ion_distr_moment  << "; Electron moment: " << electron_moment << "; Total moment: " << ion_particle_moment + electron_moment << std::endl;


			start = std::chrono::high_resolution_clock::now();
		this->Lift();
			stop = std::chrono::high_resolution_clock::now();
			lift_time += std::chrono::duration_cast<timeunit>(stop - start);
			start = std::chrono::high_resolution_clock::now();
		this->Load(state);
			stop = std::chrono::high_resolution_clock::now();
			load_time += std::chrono::duration_cast<timeunit>(stop - start);

		*simulation_time = current_time + _plasma->get_macro_to_micro_dt_ratio()*_plasma->get_dt();
	    if (_record_microsteps)
		{
			_record_ion_distribution.at(_record_times.size()) = _stack_ion_distribution.front();
			_record_times.push_back(*simulation_time);
		}
	}


	_distributions.front()->GetDensityVelocityPressure(_ion_density, _ion_velocity, _ion_pressure);
	std::swap(_current_step_ion_distribution, _stack_ion_distribution.front());

	/* Diagnostics: output timing info */
	// std::cout << ++count_steps <<  " | PIC time: " << PIC_time.count() << ", Field time: " << field_time.count() << ", Pushback time: " << pushback_time.count() << ", Algebra time: " << algebra_time.count() << ", Cutoff time: " << cutoff_time.count()  << ", Load time: "<< load_time.count() << ", Lift time: " << lift_time.count() << std::endl;

	/* Diagnostic: extract electron parameters */
	std::vector<double>::iterator 	electron_position 		= state.get_vector_of_position_arrays().back()->begin();
	std::vector<double>::iterator	electron_velocity 		= state.get_vector_of_x_velocity_arrays().back()->begin();
	std::vector<double>::iterator 	electron_weight 			= state.get_vector_of_weight_arrays().back()->begin();
	int electron_population_size = state.get_vector_of_position_arrays().back()->size();

	_distributions.back()->Weigh(electron_population_size, electron_position, electron_velocity, electron_weight);


	// /* Diagnostic: extract electron velocity distribution */
	// std::vector<std::vector<double> > profile_by_population;
	// state.ComputeVelocityProfile(profile_by_population);
	// std::cout << "t = " << count_steps <<  " | Velocity profile: " << std::endl;
	// for (auto val: profile_by_population.back())
	// 	std::cout << val << "\t";
	// std::cout << std::endl;
}

void MacroParameterizationWavelets::SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics)
{
	double * x_array 	= _plasma->get_x_grid_ptr();
	int * grid_size 	= _plasma->get_grid_size_ptr();

	diagnostics.emplace_back(new CurveDiagnostic(
				"linlin", "X", "Ion density", 0, 340));
	diagnostics.back()->AddData(x_array, _ion_density.data(), grid_size, 2);

	diagnostics.emplace_back(new CurveDiagnostic(
				"linlin", "X", "Ion velocity", 410, 340));
	diagnostics.back()->AddData(x_array, _ion_velocity.data(), grid_size, 2);

	diagnostics.emplace_back(new CurveDiagnostic(
				"linlin", "X", "Ion pressure", 820, 340));
	diagnostics.back()->AddData(x_array, _ion_pressure.data(), grid_size, 4);
}

void MacroParameterizationWavelets::WriteData(std::fstream & fout)
{
	const int downsampling = 2;

	if (!_record_microsteps)
	{
		_current_step_ion_distribution.DWT(downsampling);
		fout << _current_step_ion_distribution; 
		fout << *_distributions.back();
	}
	else
	{
		if (_record_times.size()!=_plasma->get_number_of_microsteps()*2+4)
			return;

		for (int i=0; i<_record_times.size(); i++)
		{
			fout << "t = " << _record_times.at(i) << std::endl;
			_record_ion_distribution.at(i).DWT(downsampling);
			fout << _record_ion_distribution.at(i);
		}

	}
}