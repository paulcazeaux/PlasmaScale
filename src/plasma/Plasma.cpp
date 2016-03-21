#include "Plasma.h"

/* Constructor */

Plasma::Plasma( double length, double dt, int number_of_microsteps, int macro_to_micro_dt_ratio,
				double epsilon, double E0, double w0, 
				int number_of_populations, int grid_end, int macro_grid_end, int velocity_accumulation_interval, int depth, int cutoff, double intensity,
				int max_size_history, bool use_full_PIC, bool record_microsteps)
				: 
					_dt(dt),
					_dx(length / static_cast<double>(grid_end)),
					_macro_dx(length / static_cast<double>(macro_grid_end)),
					_number_of_microsteps(number_of_microsteps),
					_macro_to_micro_dt_ratio(macro_to_micro_dt_ratio),
					_record_microsteps(record_microsteps),
					_use_full_PIC(use_full_PIC),
                    _length(length),
                    _grid_end(grid_end),
                    _macro_grid_end(macro_grid_end),
                    _depth(depth),
                    _cutoff(cutoff),
                    _intensity(intensity),
					_epsilon(epsilon),
					_E0(E0), 
					_w0(w0),
					_number_of_populations(number_of_populations),
					_max_size_history(max_size_history),
					_velocity_accumulation_interval(velocity_accumulation_interval),
					_profiling_active(_velocity_accumulation_interval>0)
{
	_grid_end_ptr = std::make_unique<int>(_grid_end);
	_macro_grid_end_ptr = std::make_unique<int>(_macro_grid_end);
	_x_grid = std::make_unique<std::vector<double> >(_grid_end+1);
	_macro_x_grid = std::make_unique<std::vector<double> >(_macro_grid_end+1);

	for (int i=1; i < _grid_end+1; i++)
	{
		_x_grid->at(i) = _dx * i;
	}
	for (int i=1; i < _macro_grid_end+1; i++)
	{
		_macro_x_grid->at(i) = _macro_dx * i;
	}

	std::vector<Eigen::Triplet<double>> M;

	/* Dirichlet condition at x = 0 */
	double scaling = _epsilon/(_dx*_dx);
	for (int i=0; i<grid_end-1; i++)
	{
		if (i>0) M.push_back(Eigen::Triplet<double>(i, i-1, -scaling));
		M.push_back(Eigen::Triplet<double>(i, i,  2*scaling));
		M.push_back(Eigen::Triplet<double>(i, i+1, -scaling));
	}
	/* Neumann condition at x = _length */
	M.push_back(Eigen::Triplet<double>(grid_end-1, grid_end-1, scaling));
	M.push_back(Eigen::Triplet<double>(grid_end-1, grid_end-2, -scaling));

	_M = SpMat(grid_end, grid_end);
	_M.setFromTriplets(M.begin(), M.end());
	_solver.compute(_M);
}

/* operator << */
std::ostream& operator<<( std::ostream& os, const Plasma& plasma)
{
	os << "PLASMA PARAMETERS:" << std::endl;
	os << "==================" << std::endl;

	os << "Steps : " << std::endl;
	os << "-----   " << std::endl;
	os << "\t Timestep: \t" 	<< plasma._dt 						<< "\t\t\t";
	os << "\t Grid step: \t"	<< plasma._dx 						<< std::endl;

	os << std::endl;
	os << "EFPI parameters : " << std::endl;
	os << "-----   " << std::endl;
	os << "\t Microscopic steps computed for each projection: \t" 	<< plasma._number_of_microsteps 	<< std::endl;
	os << "\t Ratio of microscopic to macroscopic timesteps: \t"	<< plasma._macro_to_micro_dt_ratio	<< std::endl;
	os << "\t Using full PIC timestepping without extrapolation: \t"<< plasma._use_full_PIC 			<< std::endl;
	os << "\t Recording microstepping: \t\t\t\t"					<< plasma._record_microsteps 		<< std::endl;

	os << std::endl;
	os << "Spatial dimensions: " 	<< std::endl;
	os << "------------------	" 	<< std::endl;
	os << "\t System length: \t" 	<< plasma._length;

	os << std::endl;
	os << "Grid parameters: "	<< std::endl;
	os << "---------------	"	<< std::endl;
	os << "\t Grid size: \t"	<< plasma._grid_end;
	os << "\t\t Grid length: \t"	<< plasma._x_grid->back();
	os << "\t\t Wavelet tree depth: \t" << plasma._depth 			<< std::endl;
    os << "\t\t Wavelet tree cutoff: \t" << plasma._cutoff 			<< std::endl;
    os << "\t\t PURE denoising intensity: \t" << plasma._intensity  << std::endl;

	os << std::endl;
	os << "Plasma parameters:	"	<< std::endl;
	os << "-----------------	"	<< std::endl;
	os << "\t " << plasma._number_of_populations << " population of particles;" << std::endl;
	os << "\t Epsilon: \t"	<< plasma._epsilon;
	os << "\t External E field intensity: \t" 		<< plasma._E0;
	os << "\t\t External E field frequency: \t" 		<< plasma._w0 		<< std::endl;

	os << std::endl;
	os << "Diagnostics parameters:	"	<< std::endl;
	os << "----------------------	"	<< std::endl;
	os << "\t History stack size: \t"	<< plasma._max_size_history ;
	os << "\t Velocity accumulation interval: \t"	<< plasma._velocity_accumulation_interval 	<< std::endl;

	return os;
}

void 	Plasma::Poisson_Solver(std::vector<double> & charge, std::vector<double> & potential) const
{
	assert(charge.size() == _grid_end+1 && potential.size() == _grid_end+1);

	/* Account for BCs */
	potential.at(0) = 0;
	charge.at(_grid_end) *= .5;

	/* Use Eigen Cholesky solver */
	Eigen::Map<Eigen::VectorXd> charge_map(charge.data()+1, _grid_end);
	Eigen::Map<Eigen::VectorXd> potential_map(potential.data()+1, _grid_end);
	potential_map = _solver.solve(charge_map);

	/* Reset to the original charge */
	charge.at(_grid_end) *= 2.;
}
