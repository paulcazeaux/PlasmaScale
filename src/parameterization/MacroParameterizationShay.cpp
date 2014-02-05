#include "parameterization/MacroParameterizationShay.h"


MacroParameterizationShay::MacroParameterizationShay(MacroParameterizationShay &&parameterization) :
			MacroParameterization(std::move(parameterization)),
			_grid_size(parameterization._grid_size),
			_macro_grid_size(parameterization._macro_grid_size),
			_densities(std::move(parameterization._densities)),
			_temperatures(std::move(parameterization._temperatures)),
			_velocities(std::move(parameterization._velocities)),
			_ion_density(std::move(parameterization._ion_density)),
			_ion_velocity(std::move(parameterization._ion_velocity)),
			_ion_pressure(std::move(parameterization._ion_pressure)),
			_quiet_start_vel(parameterization._quiet_start_vel),
			_quiet_start_icdf(parameterization._quiet_start_icdf),
			_electron_temperature(parameterization._electron_temperature),
			_debye_scaling(parameterization._debye_scaling)
		{}

MacroParameterizationShay& MacroParameterizationShay::operator=(MacroParameterizationShay &&parameterization)
{
	if (this == &parameterization)
		return *this;
	
	MacroParameterization::operator=(std::move(parameterization));
	_grid_size = parameterization._plasma->get_grid_size();

	_densities = std::move(parameterization._densities);
	_temperatures = std::move(parameterization._temperatures);
	_velocities = std::move(parameterization._velocities);


	_ion_density = std::move(parameterization._ion_density);
	_ion_velocity = std::move(parameterization._ion_velocity);
	_ion_pressure = std::move(parameterization._ion_pressure);

	_quiet_start_vel = parameterization._quiet_start_vel;
	_quiet_start_icdf = parameterization._quiet_start_icdf;


	_electron_temperature = parameterization._electron_temperature;
	_debye_scaling = parameterization._debye_scaling;
	return *this;
}

MacroParameterizationShay::MacroParameterizationShay(MacroParameterization & parameterization, double electron_temperature, int number_of_microsteps) :
	MacroParameterization(std::move(parameterization)), _electron_temperature(electron_temperature)
{
	if (_number_of_populations != 2)
		throw std::runtime_error("There are " + std::to_string(_number_of_populations) + " and not 2 populations as required.\n");

	_grid_size = _plasma->get_grid_size();
	_macro_grid_size = _plasma->get_macro_grid_size();

	_densities	 	= std::unique_ptr<std::vector<std::vector<double> > >(new std::vector<std::vector<double> >(_number_of_populations));
	_temperatures	= std::unique_ptr<std::vector<std::vector<double> > >(new std::vector<std::vector<double> >(_number_of_populations));
	_velocities 	= std::unique_ptr<std::vector<std::vector<double> > >(new std::vector<std::vector<double> >(_number_of_populations));

	_ion_density 	= std::unique_ptr<std::vector<std::vector<double> > >(new std::vector<std::vector<double> >(_macro_grid_size));
	_ion_velocity 	= std::unique_ptr<std::vector<std::vector<double> > >(new std::vector<std::vector<double> >(_macro_grid_size));
	_ion_pressure 	= std::unique_ptr<std::vector<std::vector<double> > >(new std::vector<std::vector<double> >(_macro_grid_size));

	for (int i = 0; i < _number_of_populations; i++)
	{
		_densities->at(i).resize(_grid_size);
		_temperatures->at(i).resize(_grid_size);
		_velocities->at(i).resize(_grid_size);
	}
	for (int i = 0; i < _macro_grid_size; i++)
	{
		_ion_density->at(i).reserve(number_of_microsteps);
		_ion_velocity->at(i).reserve(number_of_microsteps);
		_ion_pressure->at(i).reserve(number_of_microsteps);
	}

	int helper_size = _macro_grid_size;		 						// To be changed
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
	double df = 1./_quiet_start_icdf.back();
	for (int i=1; i<helper_size; i++)
	{
		_quiet_start_icdf.at(i) *= df;
	}

	_debye_scaling = _electron_temperature *  _unit_masses->front() * _plasma->get_length()
			/ (_plasma->get_epsilon() * static_cast<double>(_population_sizes->front()) * std::pow(_unit_charges->front(), 2.));
}

void MacroParameterizationShay::Load(State * state) const
/* Fill the particle arrays to initialize the microscopic state */
{	

	state->Reset();
	double dx = _plasma->get_dx();

	/* Initialize the particle arrays */

	std::vector<std::vector<double> * > positions 		= state->get_vector_of_position_arrays();
	std::vector<std::vector<double> * > x_velocities 	= state->get_vector_of_x_velocity_arrays();
	std::vector<std::vector<double> * > y_velocities 	= state->get_vector_of_y_velocity_arrays();
	std::vector<std::vector<double> * > weights 		= state->get_vector_of_y_velocity_arrays();

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

			double vt = _temperatures->at(population_index).at(bin);	// P1 interpolation ? TODO

			for (int i=0; i<bin_size; i++)
			{
				double fv = (static_cast<double>(i) + 0.5)*dn;
				while (fv >= *(it_icdf+1)) 
				{
					it_vel++;
					it_icdf++;
				}
				double cellpos = RandomTools::Generate_randomly_uniform(0.,1.);
				position->at(i+bin_start_index) = (static_cast<double>(bin)+cellpos) * dx;
				velocity_x->at(i+bin_start_index) = vt * (*it_vel + (fv - *it_icdf)/(*(it_icdf+1) - *it_icdf));
			}
			bin_start_index = bin_end_index;
		}

		if (_cyclotronic_rotation_parameters->at(population_index) != 0.) 
		{
			for (int i=0; i < population_size; i++) 
			{
				double v = velocity_x->at(i);
				double theta = RandomTools::Generate_randomly_uniform(0., 2.*M_PI);
				velocity_x->at(i) = v*std::cos(theta);
				velocity_y->at(i) = v*std::sin(theta);
			}
		}

		/* Add the perturbation */

		bin_start_index = 0;
		for (int bin = 0; bin < _grid_size; bin++)
		{
			int bin_end_index = std::ceil(static_cast<double>(bin+1)*mean_bin_size-0.5);
			
			double d = _densities->at(population_index).at(bin);
			double v = _velocities->at(population_index).at(bin);		// P1 interpolation ? TODO
			for (int i=0; i<bin_end_index-bin_start_index; i++)
			{
				weight->at(i+bin_start_index) 		 = d;
				velocity_x->at(i+bin_start_index) 	+= v;
			}
			bin_start_index = bin_end_index;
		}
	}
}

void MacroParameterizationShay::RestrictAndPushback(State * state)
{
	std::vector<double> * ion_position 		= state->get_vector_of_position_arrays().front();
	std::vector<double> * ion_velocity 		= state->get_vector_of_x_velocity_arrays().front();
	std::vector<double> * ion_weight 		= state->get_vector_of_y_velocity_arrays().front();
	int ion_population_size = ion_position->size();
	double ion_population_density = static_cast<double>(_grid_size)/static_cast<double>(ion_population_size);

	static std::vector<int> bins = std::vector<int>(ion_population_size);
	static std::vector<double> cellpos = std::vector<double>(ion_population_size);
	static std::vector<double> left_weight = std::vector<double>(ion_population_size);
	static std::vector<double> right_weight = std::vector<double>(ion_population_size);

	// First pass
	for (int i=0; i<ion_population_size; i++)
	{
		double pos = ion_position->at(i);
		double weight = ion_weight->at(i);

		bins.at(i) = _plasma->find_index_on_grid(pos);
		cellpos.at(i) = _plasma->find_position_in_cell(pos);
		right_weight.at(i) = _plasma->find_position_in_cell(pos) * weight;
		left_weight.at(i) = weight - right_weight.at(i);
	}

	static std::vector<double> working_ion_density = std::vector<double>(_grid_size);
	static std::vector<double> working_ion_pressure = std::vector<double>(_grid_size);
	static std::vector<double> working_ion_velocity = std::vector<double>(_grid_size);

	// Now we weigh the particles and compute the moments on the fine grid

	std::fill(working_ion_density.begin(), working_ion_density.end(), 0.);
	std::fill(working_ion_velocity.begin(), working_ion_velocity.end(), 0.);

	for (int i=0; i<ion_population_size; i++)
	{
		int bin = bins.at(i);
		double velocity = ion_velocity->at(i);

		working_ion_density.at(bin) += left_weight.at(i);
		working_ion_velocity.at(bin) += left_weight.at(i) * velocity;
		if (bin < ion_population_size - 1)
		{
			working_ion_density.at(bin+1) += right_weight.at(i);
			working_ion_velocity.at(bin+1) += right_weight.at(i) * velocity;
		}
		else
		{
			working_ion_density.front() += right_weight.at(i);
			working_ion_velocity.front() += right_weight.at(i) * velocity;
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
		double velocity = ion_velocity->at(i);

		if (bin < ion_population_size - 1)
		{
			double velsquare = std::pow( velocity - (1.-cellpos.at(i)) * working_ion_velocity.at(bin) 
											  			-cellpos.at(i) * working_ion_velocity.at(bin+1), 2.0);
			working_ion_pressure.at(bin) += left_weight.at(i) * velsquare;
			working_ion_pressure.at(bin+1) += right_weight.at(i) * velsquare;
		}
		else
		{
			double velsquare = std::pow( velocity - (1.-cellpos.at(i)) * working_ion_velocity.at(bin) 
											  		   -cellpos.at(i)  * working_ion_velocity.front(), 2.0);
			working_ion_pressure.at(bin) += left_weight.at(i) * velsquare;
			working_ion_pressure.front() += right_weight.at(i) * velsquare;
		}			
	}

	// Then we restrict the values to the macroscopic grid using a linear smoothing.
	int size = _grid_size;
	while (size > _macro_grid_size)
	{
		for (int i=0; i<size; i++)
		{
			working_ion_density.at(i) = 0.25 * working_ion_density.at(2*i) + 0.5 * working_ion_density.at(2*i+1) + 0.25 * working_ion_density.at(2*i+2);
			working_ion_velocity.at(i) = 0.25 * working_ion_velocity.at(2*i) + 0.5 * working_ion_velocity.at(2*i+1) + 0.25 * working_ion_velocity.at(2*i+2);
			working_ion_pressure.at(i) = 0.25 * working_ion_pressure.at(2*i) + 0.5 * working_ion_pressure.at(2*i+1) + 0.25 * working_ion_pressure.at(2*i+2);
		}
		size /= 2;
	}

	// Finally, we pushback into the data points.
	assert(size == macro_grid_size);
	for (int i=0; i<size; i++)
	{
		_ion_density->at(i).push_back(working_ion_density.at(i));
		_ion_velocity->at(i).push_back(working_ion_velocity.at(i));
		_ion_pressure->at(i).push_back(working_ion_pressure.at(i));
	}
}

void MacroParameterizationShay::ExtrapolateAndLift(int macro_to_micro_dt_ratio)
{
	/* First, compute the derivative by least-squares */

	for (int bin=0; bin<_macro_grid_size; bin++)
	{
		_ion_density->at(bin).front()  += macro_to_micro_dt_ratio * Tools::EvaluateSlope(_ion_density->at(bin));
		_ion_velocity->at(bin).front() += macro_to_micro_dt_ratio * Tools::EvaluateSlope(_ion_velocity->at(bin));
		_ion_pressure->at(bin).front() += macro_to_micro_dt_ratio * Tools::EvaluateSlope(_ion_pressure->at(bin));

		_ion_density->at(bin).resize(1);
		_ion_velocity->at(bin).resize(1);
		_ion_pressure->at(bin).resize(1);
	}


	/* Next, lift the ion density, velocity and temperature to the fine grid */
	int size = _macro_grid_size;
	for (int i=0; i<size; i++)
	{
		_densities->front().at(i) = _ion_density->at(i).front();
		_velocities->front().at(i) = _ion_velocity->at(i).front();
		_temperatures->front().at(i) = _ion_pressure->at(i).front() / _ion_density->at(i).front();
	}

	while (size < _grid_size)
	{
		for (int i=0; i<size; i++)
		{
			_densities->front().at(2*i+1) = _densities->front().at(i);
			_velocities->front().at(2*i+1) = _velocities->front().at(i);
			_temperatures->front().at(2*i+1) = _temperatures->front().at(i);
		}
		_densities->front().front() = 0.5*(_densities->front().at(1) + _densities->front().at(size-1));
		_velocities->front().front() = 0.5*(_velocities->front().at(1) + _velocities->front().at(size-1));
		_temperatures->front().front() = 0.5*(_temperatures->front().at(1) + _temperatures->front().at(size-1));
		for (int i=1; i<size; i++)
		{
			_densities->front().at(2*i) = 0.5*(_densities->front().at(2*i-1) + _densities->front().at(2*i+1));
			_velocities->front().at(2*i) = 0.5*(_velocities->front().at(2*i-1) + _velocities->front().at(2*i+1));
			_temperatures->front().at(2*i) = 0.5*(_temperatures->front().at(2*i-1) + _temperatures->front().at(2*i+1));
		}
		size *= 2;
	}
	assert(size == _grid_size);

	/* Finally, compute the electron density, velocity and temperature on the fine grid */
	static std::vector<double> log_ion_density = std::vector<double>(_grid_size);

	for (int i=0; i<_grid_size; i++)
		log_ion_density.at(i) = std::log(_densities->front().at(i));

	double scaling = _debye_scaling/_plasma->get_dx();
	for (int i=0; i<_grid_size; i++)
	{
		_densities->at(1).at(i)  = _densities->front().at(i) 
						+ scaling * (
								  0.5 *log_ion_density.at(i)
								- 0.25*log_ion_density.at((i>0 ? i-1: _grid_size-1)) 
								- 0.25*log_ion_density.at((i+1<_grid_size? i+1: 0)) 
									);
		_velocities->at(1).at(i) = _velocities->front().at(i);
	}
	std::fill(_temperatures->at(1).begin(), _temperatures->at(1).end(), _electron_temperature);
}
