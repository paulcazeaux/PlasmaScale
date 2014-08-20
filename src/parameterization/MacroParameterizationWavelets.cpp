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
			_prev_step_ion_distribution(std::move(parameterization._prev_step_ion_distribution)),
			_current_step_ion_distribution(std::move(parameterization._current_step_ion_distribution)),
			_record_times(std::move(parameterization._record_times)),
			_record_ion_distribution(std::move(parameterization._record_ion_distribution)),
			_electron_thermal_vel(parameterization._electron_thermal_vel),
			_debye_scaling(parameterization._debye_scaling)
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
	_prev_step_ion_distribution = std::move(parameterization._prev_step_ion_distribution);
	_current_step_ion_distribution = std::move(parameterization._current_step_ion_distribution);

	_record_times = std::move(parameterization._record_times);
	_record_ion_distribution = std::move(parameterization._record_ion_distribution);


	_electron_thermal_vel = parameterization._electron_thermal_vel;
	_debye_scaling = parameterization._debye_scaling;
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

	_prev_step_ion_distribution 	= ActiveWaveletRepresentation(_plasma, _ion_vmax, _depth, _macro_grid_size);
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
		_record_times.reserve(2*number_of_microsteps+5);
		for (int i = 0; i < 2*number_of_microsteps+5; i++)
		{
			_record_ion_distribution.emplace_back(_plasma, _ion_vmax, _depth, _macro_grid_size);
		}
	}

	MaxwellianRepresentation::InitializeQuietStartArrays(std::pow(2, _depth));
	_debye_scaling = std::pow(_electron_thermal_vel/_plasma_pulsations.at(1), 2.);
}

void MacroParameterizationWavelets::Initialize(State & state)
{
    _stack_index = 0;
	std::shared_ptr<double> simulation_time = state.get_simulation_time();
	double current_time = *simulation_time;
    
	if (_record_microsteps)
	{
		_record_times.clear();
	}
	this->RestrictAndPushback(state);
	_current_step_ion_distribution = _stack_ion_distribution.front();
	_prev_step_ion_distribution = _current_step_ion_distribution;
    
    /* Take an Euler step backwards to initialize the distribution at the previous step */
    
    int number_of_microsteps = _plasma->get_number_of_microsteps();
	if (_record_microsteps)
	{
		_record_times.clear();
	}
	/* Leapfrog integration : using two-stage integration */
    /* Stage 1 */
    // Step 1: Run the fine solver
	for (int i=0; i<number_of_microsteps; i++)
	{
		state.Step();
		this->RestrictAndPushback(state);
	}
	/* Compute the derivative by least-squares and step backwards */
    _prev_step_ion_distribution = _current_step_ion_distribution;
    _prev_step_ion_distribution -= _plasma->get_macro_to_micro_dt_ratio() * Tools::EvaluateSlope<ActiveWaveletRepresentation>(_stack_ion_distribution, _stack_index);
	_stack_index = 1;
    
	this->Lift();
	state.Load(*this);
    _distributions.front()->GetDensityVelocityPressure(_ion_density, _ion_velocity, _ion_pressure);
    *simulation_time = current_time;
    if (_record_microsteps)
	{
		_record_times.clear();
	}
    
}

void MacroParameterizationWavelets::Load(State & state) const
/* Fill the particle arrays to initialize the microscopic state */
{	
	state.Reset();
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
		std::vector<double>::iterator velocity_y 	= y_velocities.at(population_index)->begin();
		std::vector<double>::iterator weight 		= weights.at(population_index)->begin();
		int population_size = positions.at(population_index)->size();

		_distributions.at(population_index)->Load(population_size, position, velocity_x, weight);

		if (_cyclotronic_rotation_parameters.at(population_index) != 0.) 
		{
			for (int i=0; i < population_size; i++) 
			{
				double v = *(velocity_x+i);
				double theta 	= RandomTools::Generate_randomly_uniform(0., 2.*M_PI);
				velocity_x[i] = v*std::cos(theta);
				velocity_y[i] = v*std::sin(theta);
			}
		}
	}
}

void MacroParameterizationWavelets::RestrictAndPushback(const State & state)
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
	_stack_ion_distribution.at(_stack_index).Weigh(ion_population_size, ion_position, ion_velocity, ion_weight);
    
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

void MacroParameterizationWavelets::ExtrapolateFirstHalfStep(const double ratio)
{
	/* First, compute the derivative by least-squares and step forward */
	_stack_ion_distribution.front() = Tools::EvaluateSlope<ActiveWaveletRepresentation>(_stack_ion_distribution, _stack_index);
    _stack_ion_distribution.front() *= ratio;
	_stack_ion_distribution.front() += 0.5*(_prev_step_ion_distribution + _current_step_ion_distribution);
    _stack_ion_distribution.front().Cutoff(_cutoff);
	_stack_index = 1;
    
    if (_record_microsteps)
	{
		_record_times.push_back(_record_times.front() + ratio/2.*_plasma->get_dt());
		int m = _record_times.size()-1;
		_record_ion_distribution.at(m) = _stack_ion_distribution.front();
	}
}

void MacroParameterizationWavelets::ExtrapolateSecondHalfStep(const double ratio)
{
	/* First, compute the derivative by least-squares */
	_stack_ion_distribution.front() = Tools::EvaluateSlope<ActiveWaveletRepresentation>(_stack_ion_distribution, _stack_index);
    _stack_ion_distribution.front() *= ratio;
	_stack_ion_distribution.front()  += _current_step_ion_distribution;
    _stack_ion_distribution.front().Cutoff(_cutoff);
	_stack_index = 1;
    
    if (_record_microsteps)
	{
		_record_times.push_back(_record_times.front()+ratio*_plasma->get_dt());
		int m = _record_times.size()-1;
		_record_ion_distribution.at(m) = _stack_ion_distribution.front();
	}
}

void MacroParameterizationWavelets::Lift()
{
	int size = _macro_grid_size;
	*dynamic_cast<WaveletRepresentation*>(_distributions.front().get()) = _stack_ion_distribution.front();

	while (size < _grid_size)
	{
		_distributions.front()->Refine();
		size = _distributions.front()->get_grid_size();
	}
	assert(size == _grid_size);


	/* Initialize the electron density on the coarse grid */

	static std::vector<double> ion_density = std::vector<double>(size),
							   ion_velocity = std::vector<double>(size),
								potential = std::vector<double>(size),
							   	exp_potential = std::vector<double>(size), 
							   	update = std::vector<double>(size);
	static std::vector<double> residual = std::vector<double>(size),
								conj_dir = std::vector<double>(size),
								dir = std::vector<double>(size);

	/* Newton's method solve for the electrostatic potential assuming a Boltzmann distribution for the electrons */
	double tol2 = 1e-20, iter_max = 10;
	_distributions.front()->GetDensityVelocity(ion_density, ion_velocity);

	/* Initialization */
	for (int i=0; i<size; i++)
	{
		potential.at(i) = std::log(ion_density.at(i));
	}

	/* Newton's method loop */
	double scaling = - 0.25 * _debye_scaling/std::pow(_plasma->get_macro_dx(), 2.);
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
			break;
	}
	for (int i=0; i<size; i++)
	{
		exp_potential.at(i) = std::exp(potential.at(i));
	}
	
	/* Deduce the electron density from the Boltzmann distribution */
	_distributions.at(1)->SetAdiabaticValues(exp_potential, ion_velocity, _electron_thermal_vel);
}

void MacroParameterizationWavelets::Step(State & state)
{	
	int number_of_microsteps = _plasma->get_number_of_microsteps();
    double ratio = static_cast<double>(_plasma->get_macro_to_micro_dt_ratio());
	std::shared_ptr<double> simulation_time = state.get_simulation_time();
	double current_time = *simulation_time;
	if (_record_microsteps)
	{
		_record_times.clear();
	}
	/* Leapfrog integration : using two-stage integration */
		/* Stage 1 */
			// Step 1: Run the fine solver
    _stack_index = 0;
	this->RestrictAndPushback(state);
	std::swap(_prev_step_ion_distribution, _current_step_ion_distribution);
	_current_step_ion_distribution = _stack_ion_distribution.front();
    
	for (int i=0; i<number_of_microsteps; i++)
	{
		state.Step();
		this->RestrictAndPushback(state);
	}
			// Step 2: Extrapolate using EFPI
	this->ExtrapolateFirstHalfStep(ratio);
	this->Lift();
	state.Load(*this);
	*simulation_time = current_time + .5*ratio*_plasma->get_dt();

    /* And repeat for stage 2 */
    _stack_index = 0;
	this->RestrictAndPushback(state);
    
	for (int i=0; i<number_of_microsteps; i++)
	{
		state.Step();
		this->RestrictAndPushback(state);
	}
	this->ExtrapolateSecondHalfStep(ratio);
	this->Lift();
	state.Load(*this);
	_distributions.front()->GetDensityVelocityPressure(_ion_density, _ion_velocity, _ion_pressure);
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
	if (!_record_microsteps)
	{
		fout << _current_step_ion_distribution; 
	}
	else
	{
		if (_record_times.size()!=_plasma->get_number_of_microsteps()*2+4)
			return;

		for (int i=0; i<_record_times.size(); i++)
		{
			fout << "t = " << _record_times.at(i) << std::endl;
			fout << _record_ion_distribution.at(i);
		}

	}
}





