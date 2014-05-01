#include "PopulationOfParticles.h"

/* Constructors ========================================================================================= */
PopulationOfParticles::PopulationOfParticles(std::shared_ptr<const Plasma> plasma,
			std::shared_ptr<int> iteration,
			double population_size, double unit_mass, double unit_charge, double cyclotronic_rotation_parameter) :
		_plasma(plasma), 
		_iteration(iteration),
		_unit_mass(unit_mass),
		_total_mass(unit_mass * population_size),
		_total_weight(population_size),
		_mean_density(population_size*unit_charge/plasma->get_length()),
		_unit_charge(unit_charge),
		_cyclotronic_rotation_parameter(cyclotronic_rotation_parameter),
		_magnetized(_cyclotronic_rotation_parameter!=0.)
{
	_population_size = std::unique_ptr<int> (new int(population_size));
	_weights = std::vector<double> (population_size, 1.0);
	_position = std::vector<double> (population_size);
	_velocity_x = std::vector<double> (population_size);
	_velocity_y = std::vector<double> (population_size);

	_profiling_active = false;
}

PopulationOfParticles::PopulationOfParticles(const MacroParameterization & parameterization, const int index,
					std::shared_ptr<int> iteration) :
		_plasma(parameterization.get_plasma()),
		_iteration(iteration),
		_unit_mass(parameterization.get_unit_mass(index)),
		_unit_charge(parameterization.get_unit_charge(index)),
		_cyclotronic_rotation_parameter(parameterization.get_cyclotronic_rotation_parameter(index)),
		_magnetized(_cyclotronic_rotation_parameter!=0.)

{
	_population_size = std::unique_ptr<int> (new int(parameterization.get_population_size(index)));
	_weights = std::vector<double> (*_population_size, 1.0);
	_position = std::vector<double> (*_population_size);
	_velocity_x = std::vector<double> (*_population_size);
	_velocity_y = std::vector<double> (*_population_size);

	_total_mass = _unit_mass * (*_population_size);
	_total_weight = *_population_size;
	_mean_density = (*_population_size)*_unit_charge/_plasma->get_length();

	this->SetupVelocityDiagnostics(parameterization, index);
}

void PopulationOfParticles::Reset()
{
	std::fill(_weights.begin(), 	_weights.end(), 		1.);
	std::fill(_position.begin(), 	_position.end(), 		0.);
	std::fill(_velocity_x.begin(), _velocity_x.end(), 	0.);
	std::fill(_velocity_y.begin(), _velocity_y.end(), 	0.);
}

bool PopulationOfParticles::CheckParameters(const MacroParameterization & parameterization, const int index)
{
	return ((parameterization.get_population_size(index)==*_population_size)
			&& (parameterization.get_unit_charge(index)==_unit_charge)
			&& (parameterization.get_unit_mass(index)==_unit_mass) 
			&& (parameterization.get_cyclotronic_rotation_parameter(index)==_cyclotronic_rotation_parameter));
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
	_mean_density = _total_weight * _unit_charge/_plasma->get_length();
}

void PopulationOfParticles::Move()
{
	static double plasma_length = static_cast<double>(_plasma->get_length());
	std::vector<double>::iterator it_vel_x = _velocity_x.begin();
	
	for (auto & position : _position) 
 	{
		position += *(it_vel_x++);
		while (position < 0.0)
		{
			position += plasma_length;
		}
		while (position  >= plasma_length)
		{
			position -= plasma_length;
		}
 	}
}

void PopulationOfParticles::Accelerate(const PlasmaFields& fields, double factor /* = 1.0 */ )
{
	/* Compute the acceleration from the Lorentz forces */
	static double dt = _plasma->get_dt();
 	static int 	grid_size = _plasma->get_grid_size();
	static std::vector<double> acceleration(grid_size);
	static double previous_e_to_acc_factor = 1.;
	static int previous_iteration = 0;
	std::vector<double>::iterator it_acc, it_vel_x, it_vel_y, it_weights;

	double e_to_acc_factor = _unit_charge/_unit_mass *dt*dt* factor;
	double sum_v = 0.0, sum_v_square = 0.0;

 	if (!_magnetized)
   	{
   		if (e_to_acc_factor != previous_e_to_acc_factor || (*_iteration) != previous_iteration) /* check if it is necessary to recompute the array */
		{ 
			double * e_field_array = fields.get_electrical_field_ptr();

			int i=0;
		 	for (auto & acc : acceleration)
		 	{
		 		assert((i>=0)&&(i<grid_size));
		 		acc = e_to_acc_factor * e_field_array[i++];  // do full steps
		 	}
		 	previous_e_to_acc_factor = e_to_acc_factor;
 			previous_iteration = *_iteration;
		}

 		it_vel_x 	= _velocity_x.begin();
 		it_weights	= _weights.begin();
 		for (auto & position : _position)
 		{
			double	vold = *it_vel_x, vnew = *it_vel_x;
 			int 	bin = _plasma->find_index_on_grid(position);				// Guaranteed to belong to [[0 ... grid_size -1]]
 			double	position_in_cell = _plasma->find_position_in_cell(position);		// Belongs to [0, 1)
		
			if (bin < grid_size - 1)
			{
				vnew += (1. - position_in_cell) * acceleration.at(bin);
				vnew += position_in_cell * acceleration.at(bin+1);
			}
			else
			{
				vnew += (1. - position_in_cell) * acceleration.at(bin);
				vnew += position_in_cell * acceleration.at(bin+1-grid_size);  
			}

			sum_v 			+= (*it_weights) * vnew;
			sum_v_square 	+= std::pow(*it_weights, 2.0)*vold*vnew;
			*(it_vel_x++) = vnew;
			it_weights++;
		}
		_moment  = 		_unit_mass * sum_v 		 /dt;
 		_kinetic_energy = 0.5 * _unit_mass * sum_v_square / std::pow(dt, 2.0);
	}
	else
   	{
   		if (e_to_acc_factor != previous_e_to_acc_factor || (*_iteration) != previous_iteration) /* check if it is necessary to recompute the array */
		{ 
			double * e_field_array = fields.get_electrical_field_ptr();
			int i=0;
		 	for (auto & acc : acceleration)
		 	{
		 		assert((i>=0)&&(i<grid_size));
		 		acc = 0.5 * e_to_acc_factor * e_field_array[i++];  // do half steps
		 	}
		 	previous_e_to_acc_factor = e_to_acc_factor;
 			previous_iteration = *_iteration;
		}
 		double s = 2.0*_cyclotronic_rotation_parameter/(1.0 + std::pow(_cyclotronic_rotation_parameter, 2.0));

 		it_vel_x 	= _velocity_x.begin();
 		it_vel_y 	= _velocity_y.begin();
 		it_weights	= _weights.begin();
 		for (auto & position : _position)
 		{			
 			double	vx, vy, acc;
 			int 	bin = _plasma->find_index_on_grid(position);				// Guaranteed to belong to [[0 ... grid_size -1]]
 			double	position_in_cell = _plasma->find_position_in_cell(position);		// Belongs to [0, 1)
			if (bin < grid_size - 1)
			{
				acc = (1. - position_in_cell) * acceleration.at(bin) + position_in_cell * acceleration.at(bin+1);  
			}
			else
			{
				acc = (1. - position_in_cell) * acceleration.at(bin) + position_in_cell * acceleration.at(bin+1-grid_size);  
			}
			vy	= *it_vel_y;
			vx	= *it_vel_x - _cyclotronic_rotation_parameter * vy + acc;
			vy += s * vx;
			vx -= _cyclotronic_rotation_parameter * vy;

			sum_v_square += std::pow(*it_weights++, 2.0) * (vx*vx + vy*vy);

			*it_vel_x++ = vx + acc;
			*it_vel_y++ = vy;
		}

 		_kinetic_energy = 0.5 * _unit_mass * sum_v_square / std::pow(dt, 2.0);
	}
}

void PopulationOfParticles::Weigh(PlasmaFields& fields)
{
	auto it_weights = _weights.begin();
	for (auto & position : _position)
	{
		fields.WeighParticle( position, _unit_charge * (*it_weights++));
	}
	fields.SubstractMeanDensity(_mean_density);
}

void PopulationOfParticles::Prepare(const PlasmaFields &fields)
{
	static double dt = _plasma->get_dt();
 
 	/* Rescaling */
	if (!_magnetized)
	{
		for (auto & it : _velocity_x) 
	  		it *= dt;
	}
	else
	{
		double cos_cyclotron = 1.0 /std::sqrt(1.0 + std::pow(_cyclotronic_rotation_parameter, 2.0));
		double sin_cyclotron = cos_cyclotron * _cyclotronic_rotation_parameter;

		std::vector<double>::iterator it_vel_x = _velocity_x.begin();
		std::vector<double>::iterator it_vel_y = _velocity_y.begin();
		for (; it_vel_x!= _velocity_x.end(); it_vel_x++, it_vel_y++) 
		{
			double vx = *it_vel_x, vy = *it_vel_y;
			*it_vel_x =  dt * ( sin_cyclotron * vy + cos_cyclotron * vx);
			*it_vel_y =  dt * ( cos_cyclotron * vy - sin_cyclotron * vx);
		}
	}

	/* Prepare the particles velocities by taking one-half step back with the acceleration induced by the fields */
	this->Accelerate(fields, -0.5);
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
	_number_of_bins = std::unique_ptr<int>(new int( nbins>2 ? nbins : 2));
	_velocity_accumulation_interval = velocity_accumulation_interval;
	_count 							= 0;
	_accumulation_weight			= 1./static_cast<double>(velocity_accumulation_interval);
	_mid_bin_array 					= std::unique_ptr<std::vector<double> > (new std::vector<double>(*_number_of_bins));
	_partial_velocity_profile 		= std::unique_ptr<std::vector<double> > (new std::vector<double>(*_number_of_bins, 0.));
	_accumulated_velocity_profile	= std::unique_ptr<std::vector<double> > (new std::vector<double>(*_number_of_bins, 0.));

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
		_bin_width = (vupper-vlower)/static_cast<double>(*_number_of_bins);
	}
	else if(vt1 + vt2 > 0.0) 
	{
		_bin_start = v0 - 5.0*(vt1 + vt2);  
		_bin_width = 10.0 * (vt1 + vt2)/static_cast<double>(*_number_of_bins);  /*  so that the distribution goes from
							v0-5*vt to v0+5*vt  */
	}
	else if (v0!=0)
	{
		_bin_start = (v0 < 0 ? 2.0*v0 : 0);
		_bin_width = 2.0*std::abs(v0)/static_cast<double>(*_number_of_bins);
	}
	else
	{
		_bin_start = 0.0;
		_bin_width = 1.0/static_cast<double>(*_number_of_bins);
	}

		/*  setup _mid_bin_array for this population  */
	for(int i=0; i<*_number_of_bins; i++)
	{
		_mid_bin_array->at(i) = (_bin_start + _bin_width*(static_cast<double>(i)+0.5));
	}

		/* Rescaling */
	_bin_width *= _plasma->get_dt();
	_bin_start *= _plasma->get_dt();

}

void PopulationOfParticles::SetupVelocityDiagnostics(const MacroParameterization & parameterization, const int index)
{
	if (!parameterization.HaveVelocityDiagnostics())
	{
		_profiling_active = false;
		return;
	}

	_profiling_active = true;
	_number_of_bins = std::unique_ptr<int>(new int( parameterization.GetNumberOfBins(index)));
	_velocity_accumulation_interval = _plasma->get_velocity_accumulation_interval();
	_count 							= 0;
	_accumulation_weight			= 1./static_cast<double>(_velocity_accumulation_interval);
	_mid_bin_array 					= std::unique_ptr<std::vector<double> > (new std::vector<double>(*_number_of_bins));
	_partial_velocity_profile 		= std::unique_ptr<std::vector<double> > (new std::vector<double>(*_number_of_bins, 0.));
	_accumulated_velocity_profile	= std::unique_ptr<std::vector<double> > (new std::vector<double>(*_number_of_bins, 0.));

	_bin_start = parameterization.GetBinStart(index);
	_bin_width = parameterization.GetBinWidth(index);

		/*  setup _mid_bin_array for this population  */
	for(int i=0; i<*_number_of_bins; i++)
	{
		_mid_bin_array->at(i) = (_bin_start + _bin_width*(static_cast<double>(i)+0.5));
	}

		/* Rescaling */
	_bin_width *= _plasma->get_dt();
	_bin_start *= _plasma->get_dt();
}

void PopulationOfParticles::ComputeVelocityProfile()
{
	if (!_profiling_active)
		return;
	/* Add the current velocity distribution to the partial accumulated distribution */
	if (_count < _velocity_accumulation_interval)
	{
		auto it_weights = _weights.begin();
		for (auto & v : _velocity_x)
		{
			int v_index = static_cast<int>(std::floor((v-_bin_start)/_bin_width));
			if(v_index>=0 && v_index<*_number_of_bins)
			{
				_partial_velocity_profile->at(v_index) += _accumulation_weight * (*it_weights++);
			}
		}
		_count++;
	}
	/* Have we reached the accumulation interval ? */
	if (_count == _velocity_accumulation_interval)
	{
	/* Yes ! Now we update the accumulated profile */
		auto it_avp = _accumulated_velocity_profile->begin();
		for (auto & v : *_partial_velocity_profile)
		{
			*(it_avp++) =  v * _unit_mass;
		}

			/* Reset the partial accumulated distribution to the current velocity distribution */
		std::fill(_partial_velocity_profile->begin(), _partial_velocity_profile->end(), 0.);
		auto it_weights = _weights.begin();
		for (auto & v : _velocity_x)
		{
			int v_index = static_cast<int>(std::floor((v-_bin_start)/_bin_width));
			if(v_index>=0 && v_index<*_number_of_bins)
			{
				_partial_velocity_profile->at(v_index) += _accumulation_weight * (*it_weights++);
			}
		}
	/* Reset the accumulation counter */
		_count = 1;
	}

}
