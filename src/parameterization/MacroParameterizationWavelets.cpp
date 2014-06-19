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
			_stack_electron_distribution(std::move(parameterization._stack_electron_distribution)),
            _stack_index(parameterization._stack_index),
			_prev_step_ion_distribution(std::move(parameterization._prev_step_ion_distribution)),
			_prev_step_electron_distribution(std::move(parameterization._prev_step_electron_distribution)),
			_current_step_ion_distribution(std::move(parameterization._current_step_ion_distribution)),
			_current_step_electron_distribution(std::move(parameterization._current_step_electron_distribution)),
			_record_times(std::move(parameterization._record_times)),
			_record_ion_distribution(std::move(parameterization._record_ion_distribution)),
			_record_electron_distribution(std::move(parameterization._record_electron_distribution))
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
	_stack_electron_distribution = std::move(parameterization._stack_electron_distribution);
    _stack_index = parameterization._stack_index;
	_prev_step_ion_distribution = std::move(parameterization._prev_step_ion_distribution);
	_prev_step_electron_distribution = std::move(parameterization._prev_step_electron_distribution);
	_current_step_ion_distribution = std::move(parameterization._current_step_ion_distribution);
	_current_step_electron_distribution = std::move(parameterization._current_step_electron_distribution);

	_record_times = std::move(parameterization._record_times);
	_record_ion_distribution = std::move(parameterization._record_ion_distribution);
	_record_electron_distribution = std::move(parameterization._record_electron_distribution);

	return *this;
}

MacroParameterizationWavelets::MacroParameterizationWavelets(MacroParameterization & parameterization,
				double ion_vmax) :
	MacroParameterization(std::move(parameterization)), _ion_vmax(ion_vmax)
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
	_prev_step_electron_distribution 	= ActiveMaxwellianRepresentation(_plasma, _macro_grid_size);
	_current_step_electron_distribution	= ActiveMaxwellianRepresentation(_plasma, _macro_grid_size);


	_distributions.clear();
	_distributions.emplace_back(new ActiveWaveletRepresentation(_plasma, _ion_vmax, _depth, _grid_size));
	_distributions.emplace_back(new ActiveMaxwellianRepresentation(_plasma, _grid_size));

	int number_of_microsteps = _plasma->get_number_of_microsteps();
	_stack_ion_distribution.reserve(number_of_microsteps+1);
    _stack_index = 0;

	if (_record_microsteps)
	{
        _record_ion_distribution.clear();
		_record_times.reserve(2*number_of_microsteps+2);
		for (int i = 0; i < 2*number_of_microsteps+2; i++)
		{
			_record_ion_distribution.emplace_back(_plasma, _ion_vmax, _depth, _macro_grid_size);
			_record_electron_distribution.emplace_back(_plasma, _macro_grid_size);
		}
	}
	MaxwellianRepresentation::InitializeQuietStartArrays(std::pow(2, _depth));
}

void MacroParameterizationWavelets::Initialize(State & state)
{
    double ratio = static_cast<double>(_plasma->get_macro_to_micro_dt_ratio());
	std::shared_ptr<double> simulation_time = state.get_simulation_time();
	double current_time = *simulation_time;
    
	if (_record_microsteps)
	{
		_record_times.clear();
	}
    _stack_index = 0;
	this->RestrictAndPushback(state, ratio);
	std::swap(_current_step_ion_distribution, _stack_ion_distribution.front());
	std::swap(_current_step_electron_distribution, _stack_electron_distribution.front());
    _stack_index = 0;
	this->RestrictAndPushback(state, 0.);
    
    /* Take an Euler step backwards to initialize the distribution at the previous step */
    
    int number_of_microsteps = _plasma->get_number_of_microsteps();
	if (_record_microsteps)
	{
		_record_times.clear();
	}
	for (int i=0; i<number_of_microsteps; i++)
	{
		state.Step();
		this->RestrictAndPushback(state, ratio-i-1);
	}
	/* Compute the derivative by least-squares and step backwards */
    _prev_step_ion_distribution = _current_step_ion_distribution;
    _prev_step_ion_distribution -= ratio * Tools::EvaluateSlope<ActiveWaveletRepresentation>(_stack_ion_distribution, _stack_index);
    _prev_step_electron_distribution = _current_step_electron_distribution;
    _prev_step_electron_distribution -= ratio * Tools::EvaluateSlope<ActiveMaxwellianRepresentation>(_stack_electron_distribution, _stack_index);

	_stack_index = 1;
	this->Lift();
	state.Load(*this);
    _distributions.front()->GetDensityVelocityPressure(_ion_density, _ion_velocity, _ion_pressure);
    _distributions.back()->GetDensityVelocityPressure(_electron_density, _electron_velocity, _electron_pressure);
    *simulation_time = current_time;
    if (_record_microsteps)
	{
		_record_times.clear();
	}
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

void MacroParameterizationWavelets::RestrictAndPushback(const State & state, const double delay)
{
	/***********************************************************************/
	/* Now we weigh the particles and compute the moments on the fine grid */
	/***********************************************************************/
    /* First, we initialize the arrays */
	/* We start with the ions */
	std::vector<double>::iterator 	position 		= state.get_vector_of_position_arrays().front()->begin();
	std::vector<double>::iterator	velocity 		= state.get_vector_of_x_velocity_arrays().front()->begin();
	std::vector<double>::iterator 	weight 			= state.get_vector_of_weight_arrays().front()->begin();
	int population_size = state.get_vector_of_position_arrays().front()->size();
    if (_stack_index >= _stack_ion_distribution.size())
        _stack_ion_distribution.emplace_back(_plasma, _ion_vmax, _depth, _grid_size);
    
    int size = _stack_ion_distribution.at(_stack_index).get_grid_size();
    
    if (size != _grid_size)
    {
        _stack_ion_distribution.at(_stack_index).set_grid_size(_grid_size);
        size = _grid_size;
    }
	/* We weigh the ions */
	_stack_ion_distribution.at(_stack_index).Weigh(population_size, position, velocity, weight, delay);
   

	/* Repeat for the electrons */
	position 		= state.get_vector_of_position_arrays().back()->begin();
	velocity 		= state.get_vector_of_x_velocity_arrays().back()->begin();
	weight 			= state.get_vector_of_weight_arrays().back()->begin();
	population_size = state.get_vector_of_position_arrays().back()->size();
    if (_stack_index >= _stack_electron_distribution.size())
        _stack_electron_distribution.emplace_back(_plasma, _grid_size);
    
    size = _stack_electron_distribution.at(_stack_index).get_grid_size();
    
    if (size != _grid_size)
    {
        _stack_electron_distribution.at(_stack_index).set_grid_size(_grid_size);
        size = _grid_size;
    }
    /* We weigh the electrons */
	_stack_electron_distribution.at(_stack_index).Weigh(population_size, position, velocity, weight, delay);
    
	/* Then we restrict the values to the macroscopic grid using a linear smoothing. */
	while (size > _macro_grid_size)
	{
		_stack_ion_distribution.at(_stack_index).Coarsen();
		_stack_electron_distribution.at(_stack_index).Coarsen();
		size /= 2;
	}
	assert(size == _macro_grid_size);

	if (_record_microsteps)
	{
		_record_times.push_back(*state.get_simulation_time());
		int m = _record_times.size()-1;
		_record_ion_distribution.at(m) = _stack_ion_distribution.at(_stack_index);
		_record_electron_distribution.at(m) = _stack_electron_distribution.at(_stack_index);
	}

	_stack_index++;
}

void MacroParameterizationWavelets::Extrapolate(const double ratio)
{
	/* First, compute the derivative by least-squares and step forward */
	_stack_ion_distribution.front() += ratio * Tools::EvaluateSlope<ActiveWaveletRepresentation>(_stack_ion_distribution, _stack_index);
	_stack_electron_distribution.front() += ratio * Tools::EvaluateSlope<ActiveMaxwellianRepresentation>(_stack_electron_distribution, _stack_index);
    _stack_ion_distribution.front().Cutoff(_cutoff);
	_stack_index = 1;
}

void MacroParameterizationWavelets::Lift()
{
	int size = _macro_grid_size;
	*dynamic_cast<ActiveWaveletRepresentation*>(_distributions.front().get()) = _stack_ion_distribution.front();
	*dynamic_cast<ActiveMaxwellianRepresentation*>(_distributions.back().get()) = _stack_electron_distribution.front();

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
	int number_of_microsteps = _plasma->get_number_of_microsteps();
    double ratio = static_cast<double>(_plasma->get_macro_to_micro_dt_ratio());
	std::shared_ptr<double> simulation_time = state.get_simulation_time();
	double current_time = *simulation_time;

	/* Prepare the next step for leap-frog integration */
	_current_step_ion_distribution = _prev_step_ion_distribution;
	_current_step_electron_distribution = _prev_step_electron_distribution;
	if (_record_microsteps)
		_record_times.clear();

    _stack_index = 0;
	this->RestrictAndPushback(state, 2.*ratio);
	_prev_step_ion_distribution = _stack_ion_distribution.front();
	_prev_step_electron_distribution = _stack_electron_distribution.front();

	/* Leapfrog integration : using two-stage integration */
		/* Stage 1 */
    _stack_index = 0;
	this->RestrictAndPushback(state, ratio);
	_current_step_ion_distribution += _stack_ion_distribution.front();
	_current_step_electron_distribution += _stack_electron_distribution.front();
	_current_step_ion_distribution *= 0.5;
	_current_step_electron_distribution *= 0.5;

    if (_record_microsteps)
	{
		_record_times.clear();
		_record_times.push_back(current_time - .5*ratio*_plasma->get_dt());
		_record_ion_distribution.front() = _current_step_ion_distribution;
		_record_electron_distribution.front() = _current_step_electron_distribution;
	}

	std::swap(_current_step_ion_distribution, _stack_ion_distribution.front());
	std::swap(_current_step_electron_distribution, _stack_electron_distribution.front());

	/* Initial value for this half-step in stored in _stack_ion_distribution.front() */
	for (int i=0; i<number_of_microsteps; i++)
	{
		state.Step();
		this->RestrictAndPushback(state, ratio-i-1);
	}
			// Step 2: Extrapolate using EFPI
	this->Extrapolate(ratio);
	this->Lift();
	state.Load(*this);
	*simulation_time = current_time + .5*ratio*_plasma->get_dt();

    /* And repeat for stage 2 */

    _stack_index = 1;
	_stack_ion_distribution.front() = _current_step_ion_distribution;
	_stack_electron_distribution.front() = _current_step_electron_distribution;

    if (_record_microsteps)
	{
		_record_times.push_back(current_time);
		int m = _record_times.size()-1;
		_record_ion_distribution.at(m) = _current_step_ion_distribution;
		_record_electron_distribution.at(m) = _current_step_electron_distribution;
	}

    /* Initial value for this half-step in stored in _stack_ion_distribution.front() */
	for (int i=0; i<number_of_microsteps; i++)
	{
		state.Step();
		this->RestrictAndPushback(state, .5*ratio-i-1);
	}
	this->Extrapolate(ratio);
	this->Lift();
	state.Load(*this);
	_distributions.front()->GetDensityVelocityPressure(_ion_density, _ion_velocity, _ion_pressure);
	_distributions.back()->GetDensityVelocityPressure(_electron_density, _electron_velocity, _electron_pressure);
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


	diagnostics.emplace_back(new CurveDiagnostic(
				"linlin", "X", "Electron density", 0, 340));
	diagnostics.back()->AddData(x_array, _electron_density.data(), grid_size, 2);

	diagnostics.emplace_back(new CurveDiagnostic(
				"linlin", "X", "Electron velocity", 410, 340));
	diagnostics.back()->AddData(x_array, _electron_velocity.data(), grid_size, 2);

	diagnostics.emplace_back(new CurveDiagnostic(
				"linlin", "X", "Electron pressure", 820, 340));
	diagnostics.back()->AddData(x_array, _electron_pressure.data(), grid_size, 4);
}

void MacroParameterizationWavelets::WriteData(std::fstream & fout)
{
	if (!_record_microsteps)
	{
		fout << _current_step_ion_distribution; 
	}
	else
	{
		if (_record_times.size()!=_plasma->get_number_of_microsteps()*2+2)
			return;

		for (int i=0; i<_record_times.size(); i++)
		{
			fout << "t = " << _record_times.at(i) << std::endl;
			fout << _record_electron_distribution.at(i);
		}

	}
}





