#include "Plasma.h"

/* Constructor */

Plasma::Plasma( double length, double dt, int number_of_microsteps, int macro_to_micro_dt_ratio,
				double epsilon, double la, double rho0, double E0, double w0, 
				int number_of_populations, int grid_size, int macro_grid_size, int velocity_accumulation_interval, int max_mode,
				double filter_parameter_1, double filter_parameter_2, int max_size_history )
				: _length(length), 
					_dt(dt),
					_number_of_microsteps(number_of_microsteps),
					_macro_to_micro_dt_ratio(macro_to_micro_dt_ratio),
					_epsilon(epsilon),
					_rho0(rho0),
					_E0(E0), 
					_w0(w0),
					_number_of_populations(number_of_populations), 
					_velocity_accumulation_interval(velocity_accumulation_interval),
					_max_mode(max_mode),
					_filter_parameter_1(filter_parameter_1),
					_filter_parameter_2(filter_parameter_2),
					_max_size_history(max_size_history)
{
	_grid_size = std::unique_ptr<int>(new int(grid_size));
	_macro_grid_size = std::unique_ptr<int>(new int(macro_grid_size));
	_dx = _length / static_cast<double>(*_grid_size);
	_highest_mode = (*_grid_size) / 2; 
	_inverse_of_particle_radius = la / _length;
	_profiling_active = (_velocity_accumulation_interval>0);

	std::vector<double> x_grid = std::vector<double>(*_grid_size+1);
	std::vector<double> k_grid = std::vector<double>(_highest_mode);
	for (int i=1; i < *_grid_size+1; i++)
	{
		x_grid[i] = _dx * i;
	}
	for (int k = 1; k < _highest_mode; k++)
	{
		k_grid[k] = 2 * M_PI /_length * k;
	}
	_x_grid = std::unique_ptr<std::vector<double> > (new std::vector<double>(x_grid));
	_k_grid = std::unique_ptr<std::vector<double> > (new std::vector<double>(k_grid));
}

/* operator << */
std::ostream& operator<<( std::ostream& os, const Plasma& plasma)
{
	os << "PLASMA:" << std::endl;
	os << "=======" << std::endl << std::endl;

	os << std::endl;
	os << "Steps : " << std::endl;
	os << "-----   " << std::endl;
	os << "\t Timestep: \t" 	<< plasma._dt 						<< std::endl;
	os << "\t Grid step: \t"	<< plasma._dx 						<< std::endl;

	os << std::endl;
	os << "Spatial dimensions: " 	<< std::endl;
	os << "------------------	" 	<< std::endl;
	os << "\t System length: \t" 	<< plasma._length 						<< std::endl;
	os << "\t Particle radius: \t"	<< 1/plasma._inverse_of_particle_radius	<< std::endl;

	os << std::endl;
	os << "Grid parameters: "	<< std::endl;
	os << "---------------	"	<< std::endl;
	os << "\t Grid size: \t"	<< *(plasma._grid_size)				<< std::endl;
	os << "\t Grid length: \t"	<< plasma._x_grid->back()			<< std::endl;
	os << "\t Highest mode: \t"	<< plasma._highest_mode 			<< std::endl;
	os << "\t Highest frequency: \t"	<< plasma._k_grid->back()			<< std::endl;

	os << std::endl;
	os << "Plasma parameters:	"	<< std::endl;
	os << "-----------------	"	<< std::endl;
	os << "\t Epsilon: \t"	<< plasma._epsilon	<< std::endl;
	os << "\t Mean density: \t"	<< plasma._rho0 	<< std::endl;
	os << "\t External E field intensity: \t" 		<< plasma._E0 		<< std::endl;
	os << "\t External E field frequency: \t" 		<< plasma._w0 		<< std::endl;
	os << std::endl;
	os << "\t Number of particle populations: \t" << plasma._number_of_populations << std::endl;

	os << std::endl;
	os << "Filtering:	" 	<< std::endl;
	os << "---------	"	<< std::endl;
	os << "\t Parameter 1: \t"	<< plasma._filter_parameter_1 		<< std::endl;
	os << "\t Parameter 2: \t"	<< plasma._filter_parameter_2 		<< std::endl;

	os << std::endl;
	os << "Diagnostics parameters:	"	<< std::endl;
	os << "----------------------	"	<< std::endl;
	os << "\t History stack size: \t"	<< plasma._max_size_history  				<< std::endl;
	os << "\t Number of modes: \t"	<< plasma._max_mode							<< std::endl;
	os << "\t Velocity accumulation interval: \t"	<< plasma._velocity_accumulation_interval 	<< std::endl;

	return os;
}