#include "parameterization/MacroParameterizationEPFI.h"


MacroParameterizationEPFI::MacroParameterizationEPFI(MacroParameterizationEPFI &&parameterization) :
			MacroParameterization(std::move(parameterization)),
			_grid_size(parameterization._grid_size),
			_macro_grid_size(parameterization._macro_grid_size),
			_record_microsteps(parameterization._record_microsteps),
			_densities(std::move(parameterization._densities)),
			_thermal_vel(std::move(parameterization._thermal_vel)),
			_velocities(std::move(parameterization._velocities)),
			_stack_ion_density(std::move(parameterization._stack_ion_density)),
			_stack_ion_velocity(std::move(parameterization._stack_ion_velocity)),
			_stack_ion_pressure(std::move(parameterization._stack_ion_pressure)),
			_prev_step_ion_density(std::move(parameterization._prev_step_ion_density)),
			_prev_step_ion_velocity(std::move(parameterization._prev_step_ion_velocity)),
			_prev_step_ion_pressure(std::move(parameterization._prev_step_ion_pressure)),
			_current_step_ion_density(std::move(parameterization._current_step_ion_density)),
			_current_step_ion_velocity(std::move(parameterization._current_step_ion_velocity)),
			_current_step_ion_pressure(std::move(parameterization._current_step_ion_pressure)),
			_quiet_start_vel(parameterization._quiet_start_vel),
			_quiet_start_icdf(parameterization._quiet_start_icdf),
			_record_times(std::move(parameterization._record_times)),
			_record_ion_density(std::move(parameterization._record_ion_density)),
			_record_ion_velocity(std::move(parameterization._record_ion_velocity)),
			_record_ion_pressure(std::move(parameterization._record_ion_pressure)),
			_electron_thermal_vel(parameterization._electron_thermal_vel),
			_debye_scaling(parameterization._debye_scaling)
		{}

MacroParameterizationEPFI& MacroParameterizationEPFI::operator=(MacroParameterizationEPFI &&parameterization)
{
	if (this == &parameterization)
		return *this;
	
	MacroParameterization::operator=(std::move(parameterization));
	_grid_size = parameterization._plasma->get_grid_size();
   	_macro_grid_size = parameterization._plasma->get_macro_grid_size();
   	_record_microsteps = parameterization._plasma->get_record_microsteps();

	_densities = std::move(parameterization._densities);
	_thermal_vel = std::move(parameterization._thermal_vel);
	_velocities = std::move(parameterization._velocities);


	_stack_ion_density = std::move(parameterization._stack_ion_density);
	_stack_ion_velocity = std::move(parameterization._stack_ion_velocity);
	_stack_ion_pressure = std::move(parameterization._stack_ion_pressure);

	_prev_step_ion_density = std::move(parameterization._prev_step_ion_density);
	_prev_step_ion_velocity = std::move(parameterization._prev_step_ion_velocity);
	_prev_step_ion_pressure = std::move(parameterization._prev_step_ion_pressure);

	_current_step_ion_density = std::move(parameterization._current_step_ion_density);
	_current_step_ion_velocity = std::move(parameterization._current_step_ion_velocity);
	_current_step_ion_pressure = std::move(parameterization._current_step_ion_pressure);

	_quiet_start_vel = parameterization._quiet_start_vel;
	_quiet_start_icdf = parameterization._quiet_start_icdf;

	_record_times = std::move(parameterization._record_times);
	_record_ion_density = std::move(parameterization._stack_ion_density);
	_record_ion_velocity = std::move(parameterization._stack_ion_velocity);
	_record_ion_pressure = std::move(parameterization._stack_ion_pressure);


	_electron_thermal_vel = parameterization._electron_thermal_vel;
	_debye_scaling = parameterization._debye_scaling;
	return *this;
}

MacroParameterizationEPFI::MacroParameterizationEPFI(MacroParameterization & parameterization, double electron_thermal_vel) :
	MacroParameterization(std::move(parameterization)), _electron_thermal_vel(electron_thermal_vel)
{
	if (_number_of_populations != 2)
		throw std::runtime_error("There are " + std::to_string(_number_of_populations) + " and not 2 populations as required.\n");

	_grid_size = _plasma->get_grid_size();
	_macro_grid_size = _plasma->get_macro_grid_size();
	_record_microsteps = _plasma->get_record_microsteps();

	_densities	 	= std::vector<std::vector<double> >(_number_of_populations);
	_thermal_vel	= std::vector<std::vector<double> >(_number_of_populations);
	_velocities 	= std::vector<std::vector<double> >(_number_of_populations);

	_stack_ion_density 	= std::vector<std::vector<double> >(_macro_grid_size);
	_stack_ion_velocity = std::vector<std::vector<double> >(_macro_grid_size);
	_stack_ion_pressure = std::vector<std::vector<double> >(_macro_grid_size);

	_prev_step_ion_density 	= std::vector<double>(_macro_grid_size);
	_prev_step_ion_velocity = std::vector<double>(_macro_grid_size);
	_prev_step_ion_pressure = std::vector<double>(_macro_grid_size);

	_current_step_ion_density 	= std::vector<double>(_macro_grid_size);
	_current_step_ion_velocity 	= std::vector<double>(_macro_grid_size);
	_current_step_ion_pressure 	= std::vector<double>(_macro_grid_size);

	for (int i = 0; i < _number_of_populations; i++)
	{
		_densities.at(i).resize(_grid_size);
		_thermal_vel.at(i).resize(_grid_size);
		_velocities.at(i).resize(_grid_size);
	}

	int number_of_microsteps = _plasma->get_number_of_microsteps();
	for (int i = 0; i < _macro_grid_size; i++)
	{
		_stack_ion_density.at(i).reserve(number_of_microsteps+1);
		_stack_ion_velocity.at(i).reserve(number_of_microsteps+1);
		_stack_ion_pressure.at(i).reserve(number_of_microsteps+1);
	}

	if (_record_microsteps)
	{
		_record_times.reserve(2*number_of_microsteps+2);
		_record_ion_density  = std::vector<std::vector<double> >(2*number_of_microsteps+2);
		_record_ion_velocity = std::vector<std::vector<double> >(2*number_of_microsteps+2);
		_record_ion_pressure = std::vector<std::vector<double> >(2*number_of_microsteps+2);
		for (int i = 0; i < 2*number_of_microsteps+2; i++)
		{
			_record_ion_density.at(i).resize(_macro_grid_size);
			_record_ion_velocity.at(i).resize(_macro_grid_size);
			_record_ion_pressure.at(i).resize(_macro_grid_size);
		}
	}

	int helper_size = _grid_size;		 						// TODO : change this
	_quiet_start_vel = std::vector<double>(helper_size);
	_quiet_start_icdf = std::vector<double>(helper_size);

	double vmax = 5.0;
	double dv = 2.0*vmax/(static_cast<double>(helper_size-1));
	_quiet_start_vel.front() = -vmax;
	_quiet_start_icdf.front() = 0.0;
	for (int i=1; i < helper_size; i++)
	{
		double vv = ((static_cast<double>(i) - 0.5)*dv - vmax);
		_quiet_start_vel.at(i) 	= static_cast<double>(i)*dv - vmax;
		_quiet_start_icdf.at(i) = _quiet_start_icdf.at(i-1) + std::exp(-0.5*vv*vv);
	}
	double norm = 1./_quiet_start_icdf.back();
	for (int i=1; i<helper_size; i++)
	{
		_quiet_start_icdf.at(i) *= norm;
	}

	_debye_scaling = -std::pow(_electron_thermal_vel, 2.) *  _unit_masses.front() * _plasma->get_length()
			/ (_plasma->get_epsilon() * static_cast<double>(_population_sizes.front()) * std::pow(_unit_charges.front(), 2.));
}

void MacroParameterizationEPFI::Initialize(const State & state)
{
	for (int bin=0; bin<_macro_grid_size; bin++)
	{	
		_stack_ion_density.at(bin).resize(0);
		_stack_ion_velocity.at(bin).resize(0);
		_stack_ion_pressure.at(bin).resize(0);
	}
	if (_record_microsteps)
	{
		_record_times.resize(0);
	}
	this->RestrictAndPushback(state);
	this->Lift();

	for (int bin=0; bin<_macro_grid_size; bin++)
	{	
		_current_step_ion_density.at(bin) = _stack_ion_density.at(bin).front();
		_current_step_ion_velocity.at(bin) = _stack_ion_velocity.at(bin).front();
		_current_step_ion_pressure.at(bin) = _stack_ion_pressure.at(bin).front();
		_stack_ion_density.at(bin).resize(0);
		_stack_ion_velocity.at(bin).resize(0);
		_stack_ion_pressure.at(bin).resize(0);
	}

	_prev_step_ion_density 	= _current_step_ion_density;
	_prev_step_ion_velocity = _current_step_ion_velocity;
	_prev_step_ion_pressure = _current_step_ion_pressure;

}
void MacroParameterizationEPFI::Load(State & state) const
/* Fill the particle arrays to initialize the microscopic state */
{	
	state.Reset();
	double dx = _plasma->get_dx();

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
		std::vector<double> * position 		= positions.at(population_index);
		std::vector<double> * velocity_x 	= x_velocities.at(population_index);
		std::vector<double> * velocity_y 	= y_velocities.at(population_index);
		std::vector<double> * weight 		= weights.at(population_index);

		int population_size = position->size();
		double mean_bin_size 	= static_cast<double>(population_size) / static_cast<double>(_grid_size) ;

		int bin_start_index = 0;
		for (int bin = 0; bin < _grid_size; bin++)
		{
			int bin_end_index = std::ceil(static_cast<double>(bin+1)*mean_bin_size-0.5);
			int bin_size = bin_end_index - bin_start_index;
			assert(bin_size > 1);

			auto it_vel = _quiet_start_vel.begin();
			auto it_icdf = _quiet_start_icdf.begin();
			double dn = 1./static_cast<double>(bin_size);
			double dv = *(it_vel+1) - *it_vel;
			
			double xs = 0.;
			for (int i=0; i<bin_size; i++)
			{
				double fv = (static_cast<double>(i) + 0.5)*dn;
				while (fv >= *(it_icdf+1)) 
				{
					it_vel++;
					it_icdf++;
				}
								/* bit-reversed scrambling to reduce the correlation with the positions */
				double xsi = 0.5;
				xs -= 0.5;
				while (xs >= 0.0)
				{
					xsi *= 0.5;
					xs -= xsi;
				} 
				xs += 2.0*xsi;
				double cellpos = xs + 0.5/static_cast<double>(bin_size);

				position->at(i+bin_start_index) = (static_cast<double>(bin)+cellpos) * dx;
				velocity_x->at(i+bin_start_index) = Tools::EvaluateP1Function(_thermal_vel.at(population_index), bin, cellpos) * (*it_vel + dv*(fv - *it_icdf)/(*(it_icdf+1) - *it_icdf));

						/* Add the perturbation */
				weight->at(i+bin_start_index) 		 = Tools::EvaluateP1Function(_densities.at(population_index), bin, cellpos);
				velocity_x->at(i+bin_start_index) 	+= Tools::EvaluateP1Function(_velocities.at(population_index), bin, cellpos);

			}
			bin_start_index = bin_end_index;
		}

		if (_cyclotronic_rotation_parameters.at(population_index) != 0.) 
		{
			for (int i=0; i < population_size; i++) 
			{
				double v = velocity_x->at(i);
				double theta = RandomTools::Generate_randomly_uniform(0., 2.*M_PI);
				velocity_x->at(i) = v*std::cos(theta);
				velocity_y->at(i) = v*std::sin(theta);
			}
		}
	}
}

void MacroParameterizationEPFI::RestrictAndPushback(const State & state)
{
	std::vector<double> * ion_position 		= state.get_vector_of_position_arrays().front();
	std::vector<double> * ion_velocity 		= state.get_vector_of_x_velocity_arrays().front();
	std::vector<double> * ion_weight 		= state.get_vector_of_weight_arrays().front();
	int ion_population_size = ion_position->size();
	double ion_population_density = static_cast<double>(_grid_size)/static_cast<double>(ion_population_size);

	static std::vector<int> bins = std::vector<int>(ion_population_size);
	static std::vector<double> cellpos = 		std::vector<double>(ion_population_size);
	static std::vector<double> left_weight = 	std::vector<double>(ion_population_size);
	static std::vector<double> right_weight = 	std::vector<double>(ion_population_size);

	// First pass
	for (int i=0; i<ion_population_size; i++)
	{
		double pos = ion_position->at(i);
		double weight = ion_weight->at(i);

		bins.at(i) = _plasma->find_index_on_grid(pos);
		cellpos.at(i) = _plasma->find_position_in_cell(pos);
		right_weight.at(i) = cellpos.at(i) * weight;
		left_weight.at(i) = weight - right_weight.at(i);
	}

	static std::vector<double> working_ion_density = std::vector<double>(_grid_size);
	static std::vector<double> working_ion_pressure = std::vector<double>(_grid_size);
	static std::vector<double> working_ion_velocity = std::vector<double>(_grid_size);

	// Now we weigh the particles and compute the moments on the fine grid

	std::fill(working_ion_density.begin(), working_ion_density.end(), 0.);
	std::fill(working_ion_velocity.begin(), working_ion_velocity.end(), 0.);

	double dt = _plasma->get_dt();
	for (int i=0; i<ion_population_size; i++)
	{
		int bin = bins.at(i);
		double velocity = ion_velocity->at(i) / dt;

		working_ion_density.at(bin) += left_weight.at(i);
		working_ion_velocity.at(bin) += left_weight.at(i) * velocity;
		if (bin < _grid_size-1)
		{
			working_ion_density.at(bin+1) += right_weight.at(i);
			working_ion_velocity.at(bin+1) += right_weight.at(i) * velocity;
		}
		else
		{
			working_ion_density.at(0) += right_weight.at(i);
			working_ion_velocity.at(0) += right_weight.at(i) * velocity;
		}
	}


	for (int bin=0; bin<_grid_size; bin++)
	{
		working_ion_velocity.at(bin) /= working_ion_density.at(bin);
		working_ion_density.at(bin) *= ion_population_density;
	}

	std::fill(working_ion_pressure.begin(), working_ion_pressure.end(), 0.);
	for (int i=0; i<ion_population_size; i++)
	{
		int bin = bins.at(i);
		double velocity = ion_velocity->at(i) / dt;

		if (bin < _grid_size-1)
		{
			double velsquare = std::pow( velocity - Tools::EvaluateP1Function(working_ion_velocity, bin, cellpos.at(i)), 2.0);
			working_ion_pressure.at(bin) += left_weight.at(i) * velsquare;
			working_ion_pressure.at(bin+1) += right_weight.at(i) * velsquare;
		}
		else
		{
			double velsquare = std::pow( velocity - Tools::EvaluateP1Function(working_ion_velocity, bin, cellpos.at(i)), 2.0);
			working_ion_pressure.at(bin) += left_weight.at(i) * velsquare;
			working_ion_pressure.at(0) += right_weight.at(i) * velsquare;
		}			
	}

	for (int bin=0; bin<_grid_size; bin++)
	{
		working_ion_pressure.at(bin) *= ion_population_density;
	}

	// Then we restrict the values to the macroscopic grid using a linear smoothing.
	int size = _grid_size;

	while (size > _macro_grid_size)
	{
		size /= 2;
		double dens0 = working_ion_density.at(0), vel0 = working_ion_velocity.at(0), pres0 = working_ion_pressure.at(0);
		for (int i=0; i<size-1; i++)
		{
			working_ion_density.at(i) = 0.25 * working_ion_density.at(2*i) + 0.5 * working_ion_density.at(2*i+1) + 0.25 * working_ion_density.at(2*i+2);
			working_ion_velocity.at(i) = 0.25 * working_ion_velocity.at(2*i) + 0.5 * working_ion_velocity.at(2*i+1) + 0.25 * working_ion_velocity.at(2*i+2);
			working_ion_pressure.at(i) = 0.25 * working_ion_pressure.at(2*i) + 0.5 * working_ion_pressure.at(2*i+1) + 0.25 * working_ion_pressure.at(2*i+2);
		}
		{
			int i=size-1;
			working_ion_density.at(i) = 0.25 * working_ion_density.at(2*i) + 0.5 * working_ion_density.at(2*i+1) + 0.25 * dens0;
			working_ion_velocity.at(i) = 0.25 * working_ion_velocity.at(2*i) + 0.5 * working_ion_velocity.at(2*i+1) + 0.25 * vel0;
			working_ion_pressure.at(i) = 0.25 * working_ion_pressure.at(2*i) + 0.5 * working_ion_pressure.at(2*i+1) + 0.25 * pres0;
		}
	}

	// Finally, we pushback into the data points.
	assert(size == _macro_grid_size);
	for (int i=0; i<size; i++)
	{
		_stack_ion_density.at(i).push_back(working_ion_density.at(i));
		_stack_ion_velocity.at(i).push_back(working_ion_velocity.at(i));
		_stack_ion_pressure.at(i).push_back(working_ion_pressure.at(i));
	}
	if (_record_microsteps)
	{
		_record_times.push_back(*state.get_simulation_time());
		int m = _record_times.size()-1;
		for (int i=0; i<size; i++)
		{
			_record_ion_density.at(m).at(i) = working_ion_density.at(i);
			_record_ion_velocity.at(m).at(i) = working_ion_velocity.at(i);
			_record_ion_pressure.at(m).at(i) = working_ion_pressure.at(i);
		}
	}
}

void MacroParameterizationEPFI::ExtrapolateFirstHalfStep()
{
	/* First, compute the derivative by least-squares and step forward */
	double macro_to_micro_dt_ratio = static_cast<double>(_plasma->get_macro_to_micro_dt_ratio());

	for (int bin=0; bin<_macro_grid_size; bin++)
	{
		_stack_ion_density.at(bin).front()  = macro_to_micro_dt_ratio * Tools::EvaluateSlope(_stack_ion_density.at(bin));
		_stack_ion_velocity.at(bin).front() = macro_to_micro_dt_ratio * Tools::EvaluateSlope(_stack_ion_velocity.at(bin));
		_stack_ion_pressure.at(bin).front() = macro_to_micro_dt_ratio * Tools::EvaluateSlope(_stack_ion_pressure.at(bin));

		_stack_ion_density.at(bin).resize(1);
		_stack_ion_velocity.at(bin).resize(1);
		_stack_ion_pressure.at(bin).resize(1);

		_stack_ion_density.at(bin).front()  += 0.5*(_prev_step_ion_density.at(bin) + _current_step_ion_density.at(bin));
		_stack_ion_velocity.at(bin).front() += 0.5*(_prev_step_ion_velocity.at(bin) + _current_step_ion_velocity.at(bin));
		_stack_ion_pressure.at(bin).front() += 0.5*(_prev_step_ion_pressure.at(bin) + _current_step_ion_pressure.at(bin));
	}
}

void MacroParameterizationEPFI::ExtrapolateSecondHalfStep()
{
	/* First, compute the derivative by least-squares */
	double macro_to_micro_dt_ratio = static_cast<double>(_plasma->get_macro_to_micro_dt_ratio());

	for (int bin=0; bin<_macro_grid_size; bin++)
	{
		_stack_ion_density.at(bin).front()  = macro_to_micro_dt_ratio * Tools::EvaluateSlope(_stack_ion_density.at(bin));
		_stack_ion_velocity.at(bin).front() = macro_to_micro_dt_ratio * Tools::EvaluateSlope(_stack_ion_velocity.at(bin));
		_stack_ion_pressure.at(bin).front() = macro_to_micro_dt_ratio * Tools::EvaluateSlope(_stack_ion_pressure.at(bin));

		_stack_ion_density.at(bin).resize(1);
		_stack_ion_velocity.at(bin).resize(1);
		_stack_ion_pressure.at(bin).resize(1);

		_stack_ion_density.at(bin).front()  += _current_step_ion_density.at(bin); 
		_stack_ion_velocity.at(bin).front() += _current_step_ion_velocity.at(bin);
		_stack_ion_pressure.at(bin).front() += _current_step_ion_pressure.at(bin);
	}

	/* Initialize the next step */

	std::swap(_prev_step_ion_density, _current_step_ion_density);
	std::swap(_prev_step_ion_velocity, _current_step_ion_velocity);
	std::swap(_prev_step_ion_pressure, _current_step_ion_pressure);

	for (int bin=0; bin<_macro_grid_size; bin++)
	{
		_current_step_ion_density.at(bin) = _stack_ion_density.at(bin).front();
		_current_step_ion_velocity.at(bin) = _stack_ion_velocity.at(bin).front();
		_current_step_ion_pressure.at(bin) = _stack_ion_pressure.at(bin).front();

	}
}

void MacroParameterizationEPFI::Lift()
{
	/* First, lift the ion density, velocity and thermal velocity and the electron density to the fine grid */
	int size = _macro_grid_size;
	for (int i=0; i<size; i++)
	{
		_densities.front().at(i) = _stack_ion_density.at(i).front();
		_velocities.front().at(i) = _stack_ion_velocity.at(i).front();
		_thermal_vel.front().at(i) = std::sqrt(_stack_ion_pressure.at(i).front() / _stack_ion_density.at(i).front());
		_stack_ion_density.at(i).resize(0);
		_stack_ion_velocity.at(i).resize(0);
		_stack_ion_pressure.at(i).resize(0);
	}



	// Another possibility : lift the pressure to the fine grid and then compute the thermal velocity

	while (size < _grid_size)
	{
		for (int i=size-1; i>=0; i--)
		{
			_densities.front().at(2*i+1) 	= _densities.front().at(i);
			_velocities.front().at(2*i+1) 	= _velocities.front().at(i);
			_thermal_vel.front().at(2*i+1) 	= _thermal_vel.front().at(i);
		}
		{
			_densities.front().front() 		= 0.5*(_densities.front().at(1) + _densities.front().at(2*size-1));
			_velocities.front().front() 	= 0.5*(_velocities.front().at(1) + _velocities.front().at(2*size-1));
			_thermal_vel.front().front() 	= 0.5*(_thermal_vel.front().at(1) + _thermal_vel.front().at(2*size-1));
		}
		for (int i=1; i<size; i++)
		{
			_densities.front().at(2*i) 		= 0.5*(_densities.front().at(2*i-1) + _densities.front().at(2*i+1));
			_velocities.front().at(2*i) 	= 0.5*(_velocities.front().at(2*i-1) + _velocities.front().at(2*i+1));
			_thermal_vel.front().at(2*i) 	= 0.5*(_thermal_vel.front().at(2*i-1) + _thermal_vel.front().at(2*i+1));
		}
		size *= 2;
	}
	assert(size == _grid_size);

	/* Initialize the electron density on the fine grid */
	static std::vector<double> log_ion_density = std::vector<double>(size);

	double scaling = _debye_scaling/_plasma->get_macro_dx();
	for (int i=0; i<size; i++)
		log_ion_density.at(i) = std::log(_densities.front().at(i));
	for (int i=0; i<size; i++)
	{
		_densities.at(1).at(i)  = _densities.front().at(i) 
						+ scaling * (
								  0.5 *log_ion_density.at(i)
								- 0.25*log_ion_density.at((i>0 ? i-1: size-1)) 
								- 0.25*log_ion_density.at((i+1<size ? i+1: 0)) 
									);
	}

	/* Finally, fill the electron velocity and thermal velocity on the fine grid */
	for (int i=0; i<size; i++)
	{
		_velocities.at(1).at(i) = _velocities.front().at(i);
	}
	std::fill(_thermal_vel.at(1).begin(), _thermal_vel.at(1).end(), _electron_thermal_vel);

}

void MacroParameterizationEPFI::Step(State & state)
{	
	int number_of_microsteps = _plasma->get_number_of_microsteps();
	std::shared_ptr<double> simulation_time = state.get_simulation_time();
	double current_time = *simulation_time;
	if (_record_microsteps)
	{
		_record_times.resize(0);
	}
	/* Leapfrog integration : using two-stage integration */
		/* Stage 1 */
			// Step 1: Run the fine solver
	this->RestrictAndPushback(state);
	for (int i=0; i<number_of_microsteps; i++)
	{
		state.Step();
		this->RestrictAndPushback(state);
	}
			// Step 2: Extrapolate using EPFI
	this->ExtrapolateFirstHalfStep();
	this->Lift();
	state.Load(*this);
	*simulation_time = current_time + _plasma->get_macro_to_micro_dt_ratio()/2.*_plasma->get_dt();

		/* And reapeat for the stage 2 */
	this->RestrictAndPushback(state);
	for (int i=0; i<number_of_microsteps; i++)
	{
		state.Step();
		this->RestrictAndPushback(state);
	}
	this->ExtrapolateSecondHalfStep();
	this->Lift();
	state.Load(*this);
}

void MacroParameterizationEPFI::SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics)
{
	double * x_array 	= _plasma->get_x_grid_ptr();
	int * grid_size 	= _plasma->get_grid_size_ptr(); 

	double dt = _plasma->get_dt();
	double length = _plasma->get_length();

	diagnostics.emplace_back(new CurveDiagnostic(
				"linlin", "X", "Ion density", 0, 340));
	diagnostics.back()->AddData(x_array, _densities.front().data(), grid_size, 2);

	diagnostics.emplace_back(new CurveDiagnostic(
				"linlin", "X", "Ion velocity", 410, 340));
	diagnostics.back()->AddData(x_array, _velocities.front().data(), grid_size, 2);

	diagnostics.emplace_back(new CurveDiagnostic(
				"linlin", "X", "Ion thermal velocity", 820, 340));
	diagnostics.back()->AddData(x_array, _thermal_vel.front().data(), grid_size, 4);
}

void MacroParameterizationEPFI::WriteData(std::fstream & fout)
{
	if (!_record_microsteps)
	{
		fout << "Density:" << std::endl;
		for (double & density : _current_step_ion_density)
			fout << density << "\t";
		fout << std::endl << "Velocity:" << std::endl;
		for (double & velocity : _current_step_ion_velocity)
			fout << velocity << "\t";
		fout << std::endl << "Pressure:" << std::endl;
		for (double & pressure : _current_step_ion_pressure)
			fout << pressure << "\t";
		fout << std::endl;
	}
	else
	{
		if (_record_times.size()!=_plasma->get_number_of_microsteps()*2+2)
			return;

		for (int i=0; i<_record_times.size(); i++)
		{
			fout << "t = " << _record_times.at(i) << std::endl;
			fout << "Density:" << std::endl;
			for (double & density : _record_ion_density.at(i))
				fout << density << "\t";
			fout << std::endl << "Velocity:" << std::endl;
			for (double & velocity : _record_ion_velocity.at(i))
				fout << velocity << "\t";
			fout << std::endl << "Pressure:" << std::endl;
			for (double & pressure : _record_ion_pressure.at(i))
				fout << pressure << "\t";
			fout << std::endl;
		}

	}
}




