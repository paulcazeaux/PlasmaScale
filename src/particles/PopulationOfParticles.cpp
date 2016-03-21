#include "PopulationOfParticles.h"

/* Constructors ========================================================================================= */
PopulationOfParticles::PopulationOfParticles(std::shared_ptr<const Plasma> plasma,
			std::shared_ptr<int> iteration,
			double population_size, double unit_mass, double unit_charge) :
		_plasma(plasma), 
		_iteration(iteration),
		_population_size(population_size),
		_unit_mass(unit_mass),
		_total_mass(unit_mass * population_size),
		_total_weight(population_size),
		_mean_density(population_size*unit_charge/plasma->_length),
		_unit_charge(unit_charge)
{
	_weights = std::vector<double> (population_size, 1.0);
	_position = std::vector<double> (population_size);
	_velocity = std::vector<double> (population_size);

	_profiling_active = false;
}

PopulationOfParticles::PopulationOfParticles(const MacroParameterization & parameterization, const int index,
					std::shared_ptr<int> iteration) :
		_plasma(parameterization.get_plasma()),
		_iteration(iteration),
		_population_size(parameterization.get_population_size(index)),
		_unit_mass(parameterization.get_unit_mass(index)),
		_unit_charge(parameterization.get_unit_charge(index))

{
	_weights = std::vector<double> (_population_size, 1.0);
	_position = std::vector<double> (_population_size);
	_velocity = std::vector<double> (_population_size);

	_total_mass = _unit_mass * _population_size;
	_total_weight = _population_size;
	_mean_density = _population_size*_unit_charge/_plasma->_length;

	this->SetupVelocityDiagnostics(parameterization, index);
}

void PopulationOfParticles::set_new_number_of_particles(const int np)
{
	_population_size = np;
	_weights.resize(np);
	_position.resize(np);
	_velocity.resize(np);

	_total_mass = _unit_mass * _population_size;
	_total_weight = _population_size;
	_mean_density = _population_size*_unit_charge/_plasma->_length;
}

void PopulationOfParticles::Reset()
{
	std::fill(_weights.begin(), 	_weights.end(), 	1.);
	std::fill(_position.begin(), 	_position.end(), 	0.);
	std::fill(_velocity.begin(), 	_velocity.end(), 	0.);
}

bool PopulationOfParticles::CheckParameters(const MacroParameterization & parameterization, const int index)
{
	return ((parameterization.get_population_size(index)==_population_size)
			&& (parameterization.get_unit_charge(index)==_unit_charge)
			&& (parameterization.get_unit_mass(index)==_unit_mass));
}

void PopulationOfParticles::ComputeAggregateParameters()
{
	double total_weight = 0.;
	for (auto & weight : _weights )
	{
		total_weight += weight;
	}
	_total_weight = total_weight;
	_total_mass = total_weight * _unit_mass;
	_mean_density = _total_weight * _unit_charge/_plasma->_length;
}

void PopulationOfParticles::Move()
{	
	int i=0, N = _population_size;
	while (i<N)
 	{
		_position[i] += _velocity[i];

		if (_position[i] < 0.0)
		{
			/* Reflected particle */
			_position[i] = - _position[i];
			_velocity[i] = -_velocity[i];
			i++;
		}
		else if (_position[i]  >= _plasma->_length) 
		{
			/* Particle escapes */
			_total_weight  -= _weights[i];
			_total_mass 	= _total_weight * _unit_mass;
			_mean_density 	= _total_weight * _unit_charge/_plasma->_length;

			_position[i] = _position[N-1];
			_velocity[i] = _velocity[N-1];
			_weights[i]  = _weights[N-1];
			N--;
		}
		else
		{
			i++;
		}
	}
	_population_size = N;
}

void PopulationOfParticles::Accelerate(const PlasmaFields& fields, double factor /* = 1.0 */ )
{
	/* Compute the acceleration from the Lorentz forces */
	static std::vector<double> acceleration;
	std::vector<double>::iterator it_acc, it_vel, it_weights;

	double e_to_acc_factor = _unit_charge/_unit_mass * std::pow(_plasma->_dt, 2) * factor;
	double sum_v = 0.0, sum_v_square = 0.0;

	acceleration.resize(_plasma->_grid_end+1);
 	for (int i=0; i<=_plasma->_grid_end; i++)
 		acceleration.at(i) = e_to_acc_factor * fields._electrical_field.at(i);  // do full steps

	int grid_end = _plasma->_grid_end;
	for (int i=0; i<_population_size; i++)
	{
		double	vold =  _velocity.at(i), vnew = vold;
		int 	bin = _plasma->find_index_on_grid(_position.at(i));
		double	s = _plasma->find_position_in_cell(_position.at(i));		// Belongs to [0, 1)

		if (bin >= 0 && bin < grid_end)
		{
			vnew += (1. - s) * acceleration.at(bin);
			vnew += s * acceleration.at(bin+1);
		}

		sum_v 			+= _weights.at(i) * vnew;
		sum_v_square 	+= _weights.at(i) * vold*vnew;
		_velocity.at(i) = vnew;
	}
	_moment  = 		_unit_mass * sum_v / _plasma->_dt;
	_kinetic_energy = 0.5 * _unit_mass * sum_v_square / std::pow(_plasma->_dt, 2.0);
}

void PopulationOfParticles::Weigh(PlasmaFields& fields)
{
	double qdx = _unit_charge / _plasma->_dx;
	int grid_end = _plasma->_grid_end;

	for (int i=0; i<_population_size; i++)
	{
		double x = _position[i];
		int bin 	= _plasma->find_index_on_grid(x);
		if (bin>=0 && bin<grid_end)
		{
			double s = _plasma->find_position_in_cell(x);
			double w = _weights[i]*qdx;

			fields._charge[bin]   += (1.-s)*w;
			fields._charge[bin+1] += s*w;
		}
	}
}

void PopulationOfParticles::Prepare(const PlasmaFields &fields, const bool toggle_half_step /* = true */)
{
 	/* Rescaling */
	double dt = _plasma->_dt;
	for (auto & v : _velocity) 
  		v *= dt;
	  	
	/* Prepare the particles velocities by taking one-half step back with the acceleration induced by the fields */
	if (toggle_half_step)
		this->Accelerate(fields, -0.5);
	else
	{
		double sum_v = 0.0, sum_v_square = 0.0;
		for (int i=0; i<_population_size; i++)
		{
			double	v =  _velocity.at(i);
			sum_v 			+= _weights.at(i) * v;
			sum_v_square 	+= _weights.at(i) * v * v;
		}
		_moment  = 		_unit_mass * sum_v / _plasma->_dt;
		_kinetic_energy = 0.5 * _unit_mass * sum_v_square / std::pow(_plasma->_dt, 2.0);
	}
	this->ComputeVelocityProfile();
}


void PopulationOfParticles::SetupVelocityDiagnostics(int nbins, int velocity_accumulation_interval, double vupper, double vlower, double v0, double vt1, double vt2)
{
	if (nbins==0 || velocity_accumulation_interval==0)
	{
		_profiling_active = false;
		return;
	}

	_profiling_active = true;
	_number_of_bins = (nbins>2 ? nbins : 2);
	_velocity_accumulation_interval = velocity_accumulation_interval;
	_count 							= 0;
	_accumulation_weight			= 1./static_cast<double>(velocity_accumulation_interval);
	_mid_bin_array 					= std::make_unique<std::vector<double> > (_number_of_bins);
	_partial_velocity_profile 		= std::make_unique<std::vector<double> > (_number_of_bins);
	_accumulated_velocity_profile	= std::make_unique<std::vector<double> > (_number_of_bins);

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
		_bin_start = vlower;
		_bin_width = (vupper-vlower)/static_cast<double>(_number_of_bins);
	}
	else if(vt1 + vt2 > 0.0) 
	{
		_bin_start = v0 - 5.0*(vt1 + vt2);  
		_bin_width = 10.0 * (vt1 + vt2)/static_cast<double>(_number_of_bins);  /*  so that the distribution goes from
							v0-5*vt to v0+5*vt  */
	}
	else if (v0!=0)
	{
		_bin_start = (v0 < 0 ? 2.0*v0 : 0);
		_bin_width = 2.0*std::abs(v0)/static_cast<double>(_number_of_bins);
	}
	else
	{
		_bin_start = 0.0;
		_bin_width = 1.0/static_cast<double>(_number_of_bins);
	}

		/*  setup _mid_bin_array for this population  */
	for(int i=0; i<_number_of_bins; i++)
	{
		_mid_bin_array->at(i) = (_bin_start + _bin_width*(static_cast<double>(i)+0.5));
	}

		/* Rescaling */
	_bin_width *= _plasma->_dt;
	_bin_start *= _plasma->_dt;

}

void PopulationOfParticles::SetupVelocityDiagnostics(const MacroParameterization & parameterization, const int index)
{
	if (!parameterization.HaveVelocityDiagnostics())
	{
		_profiling_active = false;
		return;
	}

	_profiling_active = true;
	_number_of_bins = parameterization.GetNumberOfBins(index);
	_count 							= 0;
	_accumulation_weight			= 1./static_cast<double>(_plasma->_velocity_accumulation_interval);
	_mid_bin_array 					= std::make_unique<std::vector<double> > (_number_of_bins);
	_partial_velocity_profile 		= std::make_unique<std::vector<double> > (_number_of_bins);
	_accumulated_velocity_profile	= std::make_unique<std::vector<double> > (_number_of_bins);

	_bin_start = parameterization.GetBinStart(index);
	_bin_width = parameterization.GetBinWidth(index);

		/*  setup _mid_bin_array for this population  */
	for(int i=0; i<_number_of_bins; i++)
	{
		_mid_bin_array->at(i) = (_bin_start + _bin_width*(static_cast<double>(i)+0.5));
	}

		/* Rescaling */
	_bin_width *= _plasma->_dt;
	_bin_start *= _plasma->_dt;
}

void PopulationOfParticles::ComputeVelocityProfile()
{
	if (!_profiling_active)
		return;

	/* Add the current velocity distribution to the partial accumulated distribution */
	if (_count < _velocity_accumulation_interval)
	{
		for (int i=0; i<_population_size; i++)
		{
			int v_index = static_cast<int>(std::floor((_velocity.at(i)-_bin_start)/_bin_width));
			if(v_index>=0 && v_index<_number_of_bins)
			{
				_partial_velocity_profile->at(v_index) += _accumulation_weight * _weights.at(i);
			}
		}
		_count++;
	}
	/* Have we reached the accumulation interval ? */
	if (_count == _velocity_accumulation_interval)
	{
	/* Yes ! Now we update the accumulated profile */
		for (int i=0; i<_number_of_bins; i++)
			_accumulated_velocity_profile->at(i) =  _partial_velocity_profile->at(i) * _unit_mass;

			/* Reset the partial accumulated distribution to the current velocity distribution */
		std::fill(_partial_velocity_profile->begin(), _partial_velocity_profile->end(), 0.);
		for (int i=0; i<_population_size; i++)
		{
			int v_index = static_cast<int>(std::floor((_velocity.at(i)-_bin_start)/_bin_width));
			if(v_index>=0 && v_index<_number_of_bins)
			{
				_partial_velocity_profile->at(v_index) += _accumulation_weight * _weights.at(i);
			}
		}
	/* Reset the accumulation counter */
		_count = 1;
	}
}

void PopulationOfParticles::ComputeVelocityProfile(std::vector<double> & velocity_profile)
{
	velocity_profile.resize(_number_of_bins);
	std::fill(velocity_profile.begin(), velocity_profile.end(), 0.);

	for (int i=0; i<_population_size; i++)
	{
		int v_index = static_cast<int>(std::floor((_velocity.at(i)-_bin_start)/_bin_width));
		if(v_index>=0 && v_index<_number_of_bins)
		{
			velocity_profile.at(v_index) += _unit_mass * _weights.at(i);
		}
	}
}
