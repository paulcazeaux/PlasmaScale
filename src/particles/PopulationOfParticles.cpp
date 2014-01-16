#include "PopulationOfParticles.h"

/* Constructors ========================================================================================= */
PopulationOfParticles::PopulationOfParticles(std::shared_ptr<const Plasma> plasma,
			std::shared_ptr<int> iteration, std::shared_ptr<double> simulation_time,
			double population_size, double unit_mass, double unit_charge, double cyclotronic_rotation_parameter) :
		_plasma(plasma), 
		_iteration(iteration),
		_simulation_time(simulation_time),
		_unit_mass(unit_mass),
		_total_mass(unit_mass * population_size),
		_total_weight(population_size),
		_mean_density(population_size*unit_charge/plasma->get_length()),
		_unit_charge(unit_charge),
		_cyclotronic_rotation_parameter(cyclotronic_rotation_parameter),
		_magnetized(_cyclotronic_rotation_parameter!=0.)
{
	_population_size = std::unique_ptr<int> (new int(population_size));
	_weights = std::unique_ptr<std::vector<double> > (new std::vector<double> (population_size, 1.0));
	_position = std::unique_ptr<std::vector<double> > (new std::vector<double> (population_size));
	_velocity_x = std::unique_ptr<std::vector<double> > (new std::vector<double> (population_size));
	_velocity_y = std::unique_ptr<std::vector<double> > (new std::vector<double> (population_size));

	_moment = 0.;
	_kinetic_energy = 0.;

	_profiling_active = false;
}

void PopulationOfParticles::Reset()
{
	std::fill(_weights->begin(), 	_weights->end(), 		1.);
	std::fill(_position->begin(), 	_position->end(), 		0.);
	std::fill(_velocity_x->begin(), _velocity_x->end(), 	0.);
	std::fill(_velocity_y->begin(), _velocity_y->end(), 	0.);
}

void PopulationOfParticles::ChangeTime(std::shared_ptr<double> simulation_time, std::shared_ptr<int> iteration)
{
	_iteration = iteration;
	_simulation_time = simulation_time;
}

void PopulationOfParticles::Load(int group_size, double v0, double vt1, double vt2, 
					int mode, double x1, double thetax, double v1, double thetav, int quiet_start_exponent)
{
	double population_density = _plasma->get_length() / static_cast<double>(*_population_size);
	double group_density = _plasma->get_length() * static_cast<double>(group_size) / static_cast<double>(*_population_size);
	static std::vector<double> helper = std::vector<double> (*_population_size);

	this->Reset();

	for (int i=0; i < group_size; i++)
	{
		_position->at(i) 	= ( static_cast<double>(i) + 0.5 ) * population_density;
		_velocity_x->at(i) 	= v0;
	}

	if (vt2 != 0.0) 
	{
		double vmax = 5.0*vt2;
		double dv = 2.0*vmax/(static_cast<double>(*_population_size - 1));
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
			_velocity_x->at(i) += dv* (std::distance(helper.begin(), it_helper)
									 + (fv - *it_helper)/(*(it_helper+1) - *it_helper)) - vmax;
		}

		double xs = 0.0;
		for (int i=0; i<group_size; i++)
		{
			_position->at(i) = xs*group_density + 0.5*population_density;
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

	if (_magnetized) 
	{
		for (int i=0; i < group_size; i++) 
		{
			double v = _velocity_x->at(i);
			double theta = RandomTools::Generate_randomly_uniform(0., 2.*M_PI);
			_velocity_x->at(i) = v*std::cos(theta);
			_velocity_y->at(i) = v*std::sin(theta);
		}
	}
	if (group_size < *_population_size)
	{
		double xs = 0.0;
		for (int k=group_size; k< *_population_size; k+=group_size)
		{
			xs += group_density;
			for (int j=0; j < group_size; j++) 
			{
				 _position->at(j+k) = _position->at(j) + xs;
				 _velocity_x->at(j+k) = _velocity_x->at(j);
				 if (_magnetized)
				 {
					_velocity_x->at(j+k) = _velocity_y->at(j);
				 }
			}
		}
	}

	if (vt1 != 0.0)
	{
		auto it_vel_y = _velocity_y->begin();
		for (auto & vel_x : *_velocity_x)
		{
			vel_x 	+= vt1*RandomTools::Generate_randomly_normal(0.0, 1.0);
			if (_magnetized)
			{
				*it_vel_y++ 	+= vt1*RandomTools::Generate_randomly_normal(0.0, 1.0);
			}
		}
	}

	/* Add the perturbation */
	this->Perturbate(mode, x1, v1, thetax, thetav);
}

void PopulationOfParticles::Perturbate(const int mode, const double x1, const double v1, const double thetax /* = 0. */, const double thetav /* = 0. */)
{
	double plasma_length = _plasma->get_length();
	auto it_vel_x 	= _velocity_x->begin();
	auto it_weights = _weights->begin();
	_total_weight = 0;

	for (auto & position : *_position)
	{
		double theta = (2*M_PI*mode/plasma_length) * position;
		*it_weights 	+= x1 * std::sin(theta + thetax);
		_total_weight 	+= *it_weights++; 
		*it_vel_x++ 	+= v1 * std::sin(theta + thetav);
	}
	_total_mass 	= _unit_mass * _total_weight;
	_mean_density 	= _total_weight * _unit_charge/_plasma->get_length();
}


void PopulationOfParticles::Move()
{
	static double plasma_length = static_cast<double>(_plasma->get_length());
	std::vector<double>::iterator it_vel_x = _velocity_x->begin();
	
	for (auto & position : *_position) 
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
	static double previous_e_to_acc_factor = 1., previous_iteration = -1.;
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
		}

 		it_vel_x 	= _velocity_x->begin();
 		it_weights	= _weights->begin();
 		for (auto & position : *_position)
 		{			
			double	vold, vnew;
 			int 	bin = _plasma->find_index_on_grid(position);				// Guaranteed to belong to [[0 ... grid_size -1]]
 			double	position_in_cell = _plasma->find_position_in_cell(position);		// Belongs to [0, 1)
		
			vold = *it_vel_x;
			if (bin < grid_size - 1)
			{
				vnew = vold + (1. - position_in_cell) * acceleration.at(bin) + position_in_cell * acceleration.at(bin+1);
			}
			else
			{
				vnew = vold + (1. - position_in_cell) * acceleration.at(bin) + position_in_cell * acceleration.at(bin+1-grid_size);  
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
		}
 		double s = 2.0*_cyclotronic_rotation_parameter/(1.0 + std::pow(_cyclotronic_rotation_parameter, 2.0));

 		it_vel_x 	= _velocity_x->begin();
 		it_vel_y 	= _velocity_y->begin();
 		it_weights	= _weights->begin();
 		for (auto & position : *_position)
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
	auto it_weights = _weights->begin();
	for (auto & position : *_position) 
	{
		fields.WeighParticle( position, _unit_charge * (*it_weights++));
	}
	fields.SubstractMeanDensity(_mean_density);
}

void PopulationOfParticles::Prepare(const PlasmaFields &fields)
{
	static double dt = _plasma->get_dt();
 
	if (!_magnetized)
	{
		for (auto & it : *_velocity_x) 
	  		it *= dt;
	}
	else
	{
		double cos_cyclotron = 1.0 /std::sqrt(1.0 + std::pow(_cyclotronic_rotation_parameter, 2.0));
		double sin_cyclotron = cos_cyclotron * _cyclotronic_rotation_parameter;

		std::vector<double>::iterator it_vel_x = _velocity_x->begin();
		std::vector<double>::iterator it_vel_y = _velocity_y->begin();
		for (; it_vel_x!= _velocity_x->end(); it_vel_x++, it_vel_y++) 
		{
			double vx = *it_vel_x, vy = *it_vel_y;
			*it_vel_x =  dt * ( sin_cyclotron * vy + cos_cyclotron * vx);
			*it_vel_y =  dt * ( cos_cyclotron * vy - sin_cyclotron * vx);
		}
	}
	this->Accelerate(fields, -0.5);

	if (_profiling_active)
	{
		_bin_width *= dt;
		_bin_start *= dt;

			/*  This code does the velocity-distribution stuff : Set up the initial 
					velocity distribution diagnostic to reflect the startup values.  */

		auto it_weights = _weights->begin();
		for (auto & v : *_velocity_x)
		{
			int v_index = static_cast<int>(std::floor((v-_bin_start)/_bin_width));
			if(v_index>=0 && v_index<*_number_of_bins)
			{
				_partial_velocity_profile->at(v_index) += _accumulation_weight * (*it_weights++);
			}
		}

		auto it_avp = _accumulated_velocity_profile->begin();
		for (auto & v : *_partial_velocity_profile)
		{
			*(it_avp++) =  v * _unit_mass * static_cast<double>(_velocity_accumulation_interval);
		}
		_count++;
	}
}


void PopulationOfParticles::SetupVelocityDiagnostics(int nbins, int velocity_accumulation_interval, double vupper, double vlower, double v0, double vt1, double vt2)
{

	// if(nbins>NVBINMAX)
	// {
	// 	nbins=NVBINMAX;
	// }
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

}

/*							LEGACY XES1 COMMENT 					 */
/**********************************************************************
	velocity()
	This code maintains the velocity distribution code.  
	It works thusly:  I wanted to accumulate velocities in bins for
	ten iterations, then slap up a velocity distribution to be seen
	by the users of xes1.

	To do that I had to play the following tricks: 
	I do all the work in vbin_inst, and I'm always displaying vbin.
	When I've accumulated velocities for accum timesteps
	I copy the current accumulated vbin_inst into vbin.
***********************************************************************/

void PopulationOfParticles::ComputeVelocityProfile()
{
	if (!_profiling_active)
		return;

	/* Have we reached the accumulation interval ? */
	if (_count == _velocity_accumulation_interval)
	{
	/* Yes ! Now we update the accumulated profile */
		auto it_avp = _accumulated_velocity_profile->begin();
		for (auto & v : *_partial_velocity_profile)
		{
			*(it_avp++) =  v * _unit_mass;
		}

	/* Reset the accumulation counter */
		_count = 0;
		std::fill(_partial_velocity_profile->begin(), _partial_velocity_profile->end(), 0.);
	}

	/* Add the current velocity distribution to the partial accumulated distribution */
	auto it_weights = _weights->begin();
	for (auto & v : *_velocity_x)
	{
		int v_index = static_cast<int>(std::floor((v-_bin_start)/_bin_width));
		if(v_index>=0 && v_index<*_number_of_bins)
		{
			_partial_velocity_profile->at(v_index) += _accumulation_weight * (*it_weights++);
		}
	}
	_count++;
}
