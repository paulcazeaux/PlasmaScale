#include "State.h"



/* constuctor and destructor ============================================================ */
State::State(	std::shared_ptr<const Plasma> plasma ) :
	_plasma(plasma), _number_of_populations(plasma->get_number_of_populations())
{
	_simulation_time 	= std::make_shared<double>(0.0);
	_iteration			= std::make_shared<int>(0);
	_fields				= std::unique_ptr<PlasmaFields> (new PlasmaFields(_plasma, _simulation_time, _iteration));
}

State::State(const MacroParameterization & parameterization) :
	_plasma(parameterization.get_plasma()), 
	_number_of_populations(parameterization.get_number_of_populations())
{
	_simulation_time 	= std::make_shared<double>(0.0);
	_iteration			= std::make_shared<int>(0);
	_fields				= std::unique_ptr<PlasmaFields> (new PlasmaFields(_plasma, _simulation_time, _iteration));
	_populations 		= std::vector<std::unique_ptr<PopulationOfParticles> > ();

	for (int index = 0; index<_number_of_populations; index++)
	{
		_populations.push_back(std::unique_ptr<PopulationOfParticles> 
				(new PopulationOfParticles(parameterization, index, _iteration)));		
	}
	this->Load(parameterization);
}

/* operator << */
std::ostream& operator<<( std::ostream& os, const State& state)
{
	os << "            STATE:                 " << std::endl;
	os << "===================================" << std::endl << std::endl;

	os << "Number of iterations:	"	<< *(state._iteration) 			<< std::endl;
	os << "Current simulation time:	"	<< *(state._simulation_time) 	<< std::endl;

	os << "===================================" << std::endl << std::endl;

	os << state._plasma;

	os << std::endl;
	os << "===================================" << std::endl << std::endl;

	return os;
}

void State::Load(const MacroParameterization & parameterization)
{
	/* First we check that the size and parameters are right */
	for (int index=0; index<_number_of_populations; index++)
	{
		if (!_populations.at(index)->CheckParameters(parameterization, index))
		{
            std::cout << "Attention ! Wrong parameters for the population. Creating a new one, this could break diagnostics pointers !!" << std::endl;
			_populations.at(index).reset(new PopulationOfParticles(parameterization, index, _iteration));
		}
	}
	/* Then we let the parameterization fill the particle arrays */
	parameterization.Load(*this);
	for (auto & population : _populations)
	{
		population->ComputeAggregateParameters();
	}
	this->Prepare();
	this->Weigh();
	_fields->ComputeAndFilter();
}

void State::Load(const MacroParameterization & parameterization, const bool toggle_half_step)
{
	/* First we check that the size and parameters are right */
	for (int index=0; index<_number_of_populations; index++)
	{
		if (!_populations.at(index)->CheckParameters(parameterization, index))
		{
            std::cout << "Attention ! Wrong parameters for the population. Creating a new one, this could break diagnostics pointers !!" << std::endl;
			_populations.at(index).reset(new PopulationOfParticles(parameterization, index, _iteration));
		}
	}
	/* Then we let the parameterization fill the particle arrays */
	parameterization.Load(*this);
	this->Weigh();
	_fields->ComputeAndFilter();
	for (auto & population : _populations)
	{
		population->ComputeAggregateParameters();
	}
	if (toggle_half_step)
	{
		*_iteration = -1;	// Forcing the recomputation of the acceleration field
		for (auto & population : _populations)
		{
			population->Prepare(*_fields);
		}
		*_iteration = 0;
	}
	else
	{
		_fields->ComputeAndFilter();
		*_iteration = -1;	// Forcing the recomputation of the acceleration field
		_populations.front()->Prepare(*_fields, false);
		_populations.back()->Prepare(*_fields, true);
		*_iteration = 0;
	}
}

void State::Step()
{
	this->Accelerate();
	this->Move();
	_fields->ComputeAndFilter();

	(*_iteration)++;
	*_simulation_time += _plasma->get_dt();
}

void State::GetEField(std::vector<double> & Efield)
{
	_fields->ComputeAndFilter();
	_fields->GetEField(Efield);
}

void State::SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics)
{
	char buffer[25];
	double dt = _plasma->get_dt();
	double length = _plasma->get_length();


	std::vector<std::vector<double> * > positions 		= this->get_vector_of_position_arrays();
	std::vector<std::vector<double> * > x_velocities 	= this->get_vector_of_x_velocity_arrays();
	std::vector<std::vector<double> * > y_velocities 	= this->get_vector_of_y_velocity_arrays();
	std::vector<bool> magnetizations					= this->get_vector_of_magnetizations();
	std::vector<int *> sizes 							= this->get_vector_of_sizes();

	/*  set up  windows for the phase space distribution of particles.  */
	diagnostics.emplace_back(new ScatterDiagnostic(
				"linlin", "X", "Vx-X Phase Space",
				 410, 0, 
				 1.0, 1.0/dt,
				 false, false, 0.0, length, -3.0, 3.0, std::string("open")));
	for (int i = 0; i < _number_of_populations; i++)
	{
		diagnostics.back()->AddData(positions.at(i), x_velocities.at(i), sizes.at(i), i);
	}


	/*  set up  windows for velocity distributions.  */

	for (int i = 0; i<_number_of_populations; i++)
	{
		if (magnetizations.at(i))
		{
			std::sprintf(buffer, "Vy vs Vx species %d", i+1);
			diagnostics.emplace_back(new ScatterDiagnostic("linlin", "Velocity", buffer, 410, 335));
			diagnostics.back()->AddData(x_velocities.at(i), y_velocities.at(i), sizes.at(i), i);
		}
	}

	if (_plasma->is_velocity_profiling_active())
	{
		std::vector<std::vector<double> *> mid_bin_arrays 		= this->get_vector_of_bin_arrays();
		std::vector<std::vector<double> *> velocity_profiles 	= this->get_vector_of_velocity_profiles();
		std::vector<int * >	number_of_bins						= this->get_vector_of_number_of_bins();
		
		/*  set up  windows for velocity distributions.  */
		for (int i = 0; i<_number_of_populations; i++)
		{
			std::sprintf(buffer, "Species %d f(Vx)", i+1);
			diagnostics.emplace_back(new CurveDiagnostic("linlin", "Velocity", buffer, 820, 0));
			diagnostics.back()->AddData(mid_bin_arrays.at(i), velocity_profiles.at(i), number_of_bins.at(i), i);
		}

		/********************************************/
		/*  this graph puts up a curve of ALL the velocity distributions. */
		diagnostics.emplace_back(new CurveDiagnostic("linlin","Velocity","f(v) ALL", 820, 0));
		for (int i = 0; i<_number_of_populations; i++)
		{
			diagnostics.back()->AddData(mid_bin_arrays.at(i), velocity_profiles.at(i), number_of_bins.at(i), i);
		}
	}


	double * x_array 	= _plasma->get_x_grid_ptr();
	int * grid_size 	= _plasma->get_grid_size_ptr(); 
	double * rho		= _fields->get_density_ptr();
	double * e 			= _fields->get_electrical_field_ptr();
	double * phi 		= _fields->get_potential_ptr();

	diagnostics.emplace_back(new CurveDiagnostic("linlin", "X", "rho(x)", 0, 670, 1.0, 1.0, false, true, 0.0, length, 0.0, 0.0));
	diagnostics.back()->AddData(x_array, rho, grid_size, 2);

	diagnostics.emplace_back(new CurveDiagnostic("linlin", "X", "E field(x)", 410, 670, 1.0, 1.0, false, true, 0.0, length, 0.0, 0.0));
	diagnostics.back()->AddData(x_array, e, grid_size, 3);

	diagnostics.emplace_back(new CurveDiagnostic("linlin", "X", "Potential(x)", 820, 670, 1.0, 1.0, false, true, 0.0, length, 0.0, 0.0));
	diagnostics.back()->AddData(x_array, phi, grid_size, 4);

	/* SKIP THE GRAPH OF Potential phi(k) : IT IS FALSE IN THE ORIGINAL XES1 */
}


/* getters for the diagnostics ========================================================== */
std::vector<std::vector<double> * >	State::get_vector_of_position_arrays() const
{
	std::vector<std::vector<double> * > positions = std::vector<std::vector<double> * > (_populations.size());
	for (int i = 0; i < _number_of_populations; i++)
	{
		positions.at(i) = & (_populations.at(i)->_position);
	}
	return positions;
}

std::vector<std::vector<double> * >	State::get_vector_of_x_velocity_arrays() const
{
	std::vector<std::vector<double> * > velocities = std::vector<std::vector<double> * > (_populations.size());
	for (int i = 0; i < _number_of_populations; i++)
	{
		velocities.at(i) = & (_populations.at(i)->_velocity_x);
	}
	return velocities;
}

std::vector<std::vector<double> * >	State::get_vector_of_y_velocity_arrays() const
{
	std::vector<std::vector<double> * > velocities = std::vector<std::vector<double> * > (_populations.size());
	for (int i = 0; i < _number_of_populations; i++)
	{
		velocities.at(i) = & (_populations.at(i)->_velocity_y);
	}
	return velocities;
}

std::vector<std::vector<double> * >	State::get_vector_of_weight_arrays() const
{
	std::vector<std::vector<double> * > weights = std::vector<std::vector<double> * > (_populations.size());
	for (int i = 0; i < _number_of_populations; i++)
	{
		weights.at(i) = & (_populations.at(i)->_weights);
	}
	return weights;
}

std::vector<bool> State::get_vector_of_magnetizations() const
{
	std::vector<bool> magnetizations = std::vector<bool>(_populations.size());
	for (int i = 0; i < _number_of_populations; i++)
	{
		magnetizations.at(i) = _populations.at(i)->_magnetized;
	}
	return magnetizations;
}

std::vector<int *>	State::get_vector_of_sizes() const
{
	std::vector<int *> sizes = std::vector<int *>(_number_of_populations);
	for (int i = 0; i < _number_of_populations; i++)
	{
		sizes.at(i) = (_populations.at(i)->_population_size).get();
	}
	return sizes;
}

std::vector<std::vector<double> * >	State::get_vector_of_bin_arrays() const
{
	std::vector<std::vector<double> * > bin_arrays = std::vector<std::vector<double> * > (_number_of_populations);
	for (int i = 0; i < _number_of_populations; i++)
	{
		bin_arrays.at(i) = (_populations.at(i)->_mid_bin_array).get();
	}
	return bin_arrays;
}

std::vector<std::vector<double> * >	State::get_vector_of_velocity_profiles() const 
{
	std::vector<std::vector<double> * > velocity_profiles = std::vector<std::vector<double> * > (_number_of_populations);
	for (int i = 0; i < _number_of_populations; i++)
	{
		velocity_profiles.at(i) = (_populations.at(i)->_accumulated_velocity_profile).get();
	}
	return velocity_profiles;
}

std::vector<int *>	State::get_vector_of_number_of_bins() const
{
	std::vector<int *> number_of_bins = std::vector<int *>(_number_of_populations);
	for (int i = 0; i < _number_of_populations; i++)
	{
		number_of_bins.at(i) = (_populations.at(i)->_number_of_bins).get();
	}
	return number_of_bins;
}


/* methods ============================================================================== */

void 	State::ComputeVelocityProfile()
{
	for (auto & population : _populations)
	{
		population->ComputeVelocityProfile();
	}
}

void	State::PushBackElectrostaticEnergy(std::vector<double>& electrostatic_energy, std::vector<std::vector<double>	>& electrostatic_energy_by_mode) const
{
	_fields->PushBackElectrostaticEnergy(electrostatic_energy, electrostatic_energy_by_mode);
}


void	State::PushBackKineticEnergy(std::vector<double> &kinetic_energy, std::vector<std::vector<double> >	&kinetic_energy_by_population) const
{
	double ke_total = 0.0;
	for (int i=0; i<_number_of_populations; i++)
	{
		double ke_pop = _populations.at(i)->get_kinetic_energy();
		kinetic_energy_by_population.at(i).push_back(ke_pop + 1e-30);
		ke_total += ke_pop;
	}
	kinetic_energy.push_back(ke_total);
}

void	State::PushBackMoment(std::vector<double> &moment, std::vector<std::vector<double> >	&moment_by_population) const
{
	double p_total = 0.0;
	for (int i=0; i<_number_of_populations; i++)
	{
		double p_pop = _populations.at(i)->get_moment();
		moment_by_population.at(i).push_back(p_pop);
		p_total += p_pop;
	}
	moment.push_back(p_total);
}




