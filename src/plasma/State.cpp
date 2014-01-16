#include "State.h"



/* constuctor and destructor ============================================================ */
State::State(	std::shared_ptr<const Plasma> plasma ) :
	_plasma(plasma)
{
	_simulation_time 	= std::make_shared<double>(0.0);
	_iteration			= std::make_shared<int>(0);
	_fields				= std::unique_ptr<PlasmaFields> (new PlasmaFields(_plasma, _simulation_time, _iteration));
}

State::State(	std::shared_ptr<const Plasma> plasma, FILE *& InputDeck ) :
	_plasma(plasma)
{
	_simulation_time 	= std::make_shared<double>(0.0);
	_iteration			= std::make_shared<int>(0);
	_fields				= std::unique_ptr<PlasmaFields> (new PlasmaFields(_plasma, _simulation_time, _iteration));
	_populations 		= std::unique_ptr<CollectionOfPopulations> (new CollectionOfPopulations(_plasma, _iteration, _simulation_time, InputDeck));

	_populations->Prepare(*_fields);
	_populations->Weigh(*_fields);
	_fields->ComputeAndFilter();
}

State::State(	std::shared_ptr<const Plasma> plasma, std::unique_ptr<CollectionOfPopulations> populations ) :
	_plasma(plasma),
	_populations(std::move(populations))
{
	_simulation_time 	= std::make_shared<double>(0.0);
	_iteration			= std::make_shared<int>(0);
	_fields				= std::unique_ptr<PlasmaFields> (new PlasmaFields(_plasma, _simulation_time, _iteration));

	_populations->ChangeTime(_simulation_time, _iteration);
	_populations->Prepare(*_fields);
	_populations->Weigh(*_fields);
	_fields->ComputeAndFilter();
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


void State::Step()
{
	_populations->Accelerate(*_fields);
	_populations->Move(*_fields);
	_fields->ComputeAndFilter();

	(*_iteration)++;
	*_simulation_time += _plasma->get_dt();
}

void State::SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics)
{
	char buffer[25];
	double dt = _plasma->get_dt();
	double length = _plasma->get_length();
	int number_of_populations = _plasma->get_number_of_populations();


	std::vector<std::vector<double> * > positions 		= _populations->get_vector_of_position_arrays();
	std::vector<std::vector<double> * > x_velocities 	= _populations->get_vector_of_x_velocity_arrays();
	std::vector<std::vector<double> * > y_velocities 	= _populations->get_vector_of_y_velocity_arrays();
	std::vector<bool> magnetizations					= _populations->get_vector_of_magnetizations();
	std::vector<int *> sizes 							= _populations->get_vector_of_sizes();

	/*  set up  windows for velocity distributions.  */
	diagnostics.emplace_back(new ScatterDiagnostic(
				"linlin", "X", "Vx-X Phase Space",
				 300, 10, 
				 1.0, 1.0/dt,
				 false, false, 0.0, length, -3.0, 3.0));
	for (int i = 0; i < number_of_populations; i++)
	{
		diagnostics.back()->AddData(positions[i], x_velocities[i], sizes[i], i);
	}


	/*  set up  windows for velocity distributions.  */

	for (int i = 0; i<number_of_populations; i++)
	{
		if (magnetizations[i])
		{
			std::sprintf(buffer, "Vy vs Vx species %d", i+1);
			diagnostics.emplace_back(new ScatterDiagnostic("linlin", "Velocity", buffer, 700, 400));
			diagnostics.back()->AddData(x_velocities[i], y_velocities[i], sizes[i], i);
		}
	}

	if (_plasma->is_velocity_profiling_active())
	{
		std::vector<std::vector<double> *> mid_bin_arrays 		= _populations->get_vector_of_bin_arrays();
		std::vector<std::vector<double> *> velocity_profiles 	= _populations->get_vector_of_velocity_profiles();
		std::vector<int * >	number_of_bins						= _populations->get_vector_of_number_of_bins();
		
		/*  set up  windows for velocity distributions.  */
		for (int i = 0; i<number_of_populations; i++)
		{
			std::sprintf(buffer, "Species %d f(Vx)", i+1);
			diagnostics.emplace_back(new CurveDiagnostic("linlin", "Velocity", buffer, 700, 400));
			diagnostics.back()->AddData(mid_bin_arrays[i], velocity_profiles[i], number_of_bins[i], i);
		}

		/********************************************/
		/*  this graph puts up a curve of ALL the velocity distributions. */
		diagnostics.emplace_back(new CurveDiagnostic("linlin","Velocity","f(v) ALL species", 700, 400));		
		for (int i = 0; i<number_of_populations; i++)
		{
			diagnostics.back()->AddData(mid_bin_arrays[i], velocity_profiles[i], number_of_bins[i], i);
		}
	}

		/* SKIP THE GRAPH OF f(v) TOTAL : IT IS FALSE IN THE ORIGINAL XES1 */

	double * x_array 	= _plasma->get_x_grid_ptr();
	int * grid_size 	= _plasma->get_grid_size_ptr(); 
	double * rho		= _fields->get_density_ptr();
	double * e 			= _fields->get_electrical_field_ptr();
	double * phi 		= _fields->get_potential_ptr();

	diagnostics.emplace_back(new CurveDiagnostic("linlin", "X", "rho(x)", 400, 100, 1.0, 1.0, false, true, 0.0, length, 0.0, 0.0));
	diagnostics.back()->AddData(x_array, rho, grid_size, 2);

	diagnostics.emplace_back(new CurveDiagnostic("linlin", "X", "E field(x)", 10, 500, 1.0, 1.0, false, true, 0.0, length, 0.0, 0.0));
	diagnostics.back()->AddData(x_array, e, grid_size, 3);

	diagnostics.emplace_back(new CurveDiagnostic("linlin", "X", "Potential(x)", 400, 300, 1.0, 1.0, false, true, 0.0, length, 0.0, 0.0));
	diagnostics.back()->AddData(x_array, phi, grid_size, 4);

	/* SKIP THE GRAPH OF Potential phi(k) : IT IS FALSE IN THE ORIGINAL XES1 */
}