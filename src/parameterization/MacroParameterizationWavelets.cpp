#include "parameterization/MacroParameterizationWavelets.h"


MacroParameterizationWavelets::MacroParameterizationWavelets(MacroParameterizationWavelets &&parameterization) :
			MacroParameterization(std::move(parameterization)),
			_grid_size(parameterization._grid_size),
			_macro_grid_size(parameterization._macro_grid_size),
			_record_microsteps(parameterization._record_microsteps),
			_distributions(std::move(parameterization._distributions)),
			_stack_ion_distribution(std::move(parameterization._stack_ion_distribution)),
			_prev_step_distribution(std::move(parameterization._prev_step_distribution)),
			_current_step_distribution(std::move(parameterization._current_step_distribution)),
			_record_times(std::move(parameterization._record_times)),
			_record_wave_expansion(std::move(parameterization._record_wave_expansion)),
			_electron_thermal_vel(parameterization._electron_thermal_vel),
			_debye_scaling(parameterization._debye_scaling)
		{}

MacroParameterizationWavelets& MacroParameterizationWavelets::operator=(MacroParameterizationWavelets &&parameterization)
{
	if (this == &parameterization)
		return *this;
	
	MacroParameterization::operator=(std::move(parameterization));
	_grid_size = parameterization._plasma->get_grid_size();
   	_macro_grid_size = parameterization._plasma->get_macro_grid_size();
   	_record_microsteps = parameterization._plasma->get_record_microsteps();

	_distributions = std::move(parameterization._distributions);
	_stack_ion_distribution = std::move(parameterization._stack_ion_distribution);
	_prev_step_distribution = std::move(parameterization._prev_step_distribution);
	_current_step_distribution = std::move(parameterization._current_step_distribution);

	_record_times = std::move(parameterization._record_times);
	_record_wave_expansion = std::move(parameterization._record_wave_expansion);


	_electron_thermal_vel = parameterization._electron_thermal_vel;
	_debye_scaling = parameterization._debye_scaling;
	return *this;
}

MacroParameterizationWavelets::MacroParameterizationWavelets(MacroParameterization & parameterization, double electron_thermal_vel, std::vector<double> vmax, int depth) :
	MacroParameterization(std::move(parameterization)), _electron_thermal_vel(electron_thermal_vel)
{
	if (_number_of_populations != 2)
		throw std::runtime_error("There are " + std::to_string(_number_of_populations) + " and not 2 populations as required.\n");

	_grid_size = _plasma->get_grid_size();
	_macro_grid_size = _plasma->get_macro_grid_size();
	_record_microsteps = _plasma->get_record_microsteps();

	_distributions	 			= std::vector<std::vector<std::unique_ptr<Representation> > >(_number_of_populations);
	_stack_ion_distribution 	= std::vector<std::vector<WaveletRepresentation> >(_macro_grid_size);
	_prev_step_distribution 	= std::vector<WaveletRepresentation>(_macro_grid_size);
	_current_step_distribution	= std::vector<WaveletRepresentation>(_macro_grid_size);


	for (int j = 0; j < _grid_size; j++)
	{
		_distributions.at(0).emplace_back(new WaveletRepresentation(_plasma, vmax.at(0), depth));
	}	
	for (int j = 0; j < _grid_size; j++)
	{
		_distributions.at(1).emplace_back(new MaxwellianRepresentation(_plasma));
	}

	int number_of_microsteps = _plasma->get_number_of_microsteps();
	for (int i = 0; i < _macro_grid_size; i++)
	{
		_stack_ion_distribution.at(i).reserve(number_of_microsteps+1);
	}

	if (_record_microsteps)
	{
		_record_times.reserve(2*number_of_microsteps+2);
		_record_ion_distribution  = std::vector<std::vector<std::vector<double> > >(2*number_of_microsteps+2);
		for (int i = 0; i < 2*number_of_microsteps+2; i++)
		{
			_record_ion_distribution.at(i).resize(_macro_grid_size);
		}
	}

	MaxwellianRepresentation::InitializeQuietStartArrays(std::pow(2, depth));
	_debye_scaling = std::pow(_electron_thermal_vel/_plasma_pulsations.at(1), 2.);
}

void MacroParameterizationWavelets::Initialize(const State & state)
{
	for (int bin=0; bin<_macro_grid_size; bin++)
	{	
		_stack_ion_distribution.at(bin).resize(0);
	}
	if (_record_microsteps)
	{
		_record_times.resize(0);
	}
	this->RestrictAndPushback(state);
	this->Lift();

	for (int bin=0; bin<_macro_grid_size; bin++)
	{	
		_current_step_ion_distribution.at(bin) = _stack_ion_distribution.at(bin).front();
		_stack_ion_distribution.at(bin).resize(0);
	}
	_prev_step_ion_distribution = _current_step_ion_distribution;
}

void MacroParameterizationWavelets::Load(State & state) const
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
		std::vector<double>::iterator position 		= positions.at(population_index)->begin();
		std::vector<double>::iterator velocity_x 	= x_velocities.at(population_index)->begin();
		std::vector<double>::iterator velocity_y 	= y_velocities.at(population_index)->begin();
		std::vector<double>::iterator weight 		= weights.at(population_index)->begin();
		int population_size = positions.at(population_index)->size();

		distribution.Load(population_size, position, velocity_x, weight);

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

void MacroParameterizationWavelets::RestrictAndPushback(const State & state)
{
	std::vector<double>::iterator 	ion_position 		= state.get_vector_of_position_arrays().front()->begin();
	std::vector<double>::iterator	ion_velocity 		= state.get_vector_of_x_velocity_arrays().front()->begin();
	std::vector<double>::iterator 	ion_weight 			= state.get_vector_of_weight_arrays().front()->begin();
	int ion_population_size = state.get_vector_of_position_arrays().front()->size();

	static std::vector<WaveletRepresentation> working_ion_distribution = std::vector<double>(_grid_size);

	/***********************************************************************/
	/* Now we weigh the particles and compute the moments on the fine grid */
	/***********************************************************************/

	working_ion_distribution.Weigh(ion_population_size, ion_position->begin(), ion_velocity->begin(), ion_weight->begin());
	/* Then we restrict the values to the macroscopic grid using a linear smoothing. */
	int size = working_ion_distribution.get_grid_size();
	assert(size == _grid_size);

	while (size > _macro_grid_size)
	{
		working_ion_distribution.Coarsen();
		int size = working_ion_distribution.get_grid_size();
	}

	/* Finally, we pushback into the data points. */
	assert(size == _macro_grid_size);
	_stack_ion_distribution.push_back(working_ion_distribution);

	if (_record_microsteps)
	{
		_record_times.push_back(*state.get_simulation_time());
		int m = _record_times.size()-1;
		std::swap(_record_ion_distribution.at(m), working_ion_distribution);
	}
}

void MacroParameterizationWavelets::ExtrapolateFirstHalfStep()
{
	/* First, compute the derivative by least-squares and step forward */
	double macro_to_micro_dt_ratio = static_cast<double>(_plasma->get_macro_to_micro_dt_ratio());


	_stack_ion_distribution.front()  = macro_to_micro_dt_ratio * Tools::EvaluateSlope(_stack_ion_distribution);
	_stack_ion_distribution.resize(1);
	_stack_ion_distribution.front()  += 0.5*(_prev_step_ion_distribution.at(bin) + _current_step_ion_distribution);

}

void MacroParameterizationWavelets::ExtrapolateSecondHalfStep()
{
	/* First, compute the derivative by least-squares */
	double macro_to_micro_dt_ratio = static_cast<double>(_plasma->get_macro_to_micro_dt_ratio());

	_stack_ion_distribution.front()  = macro_to_micro_dt_ratio * Tools::EvaluateSlope(_stack_ion_distribution);
	_stack_ion_distribution.resize(1);
	_stack_ion_distribution.front()  += _current_step_ion_distribution;

	/* Initialize the next step */

	std::swap(_prev_step_ion_distribution, _current_step_ion_distribution);
	_current_step_ion_distribution = _stack_ion_distribution.front();
}

void MacroParameterizationWavelets::Lift()
{
	/* First, lift the ion density, velocity and thermal velocity and the electron density to the fine grid */
	int size = _macro_grid_size;
	*_distributions.front() = _stack_ion_distribution.front();
	_stack_ion_distribution.resize(0);
	
	// Another possibility : lift the pressure to the fine grid and then compute the thermal velocity

	while (size < _grid_size)
	{
		_distributions.front()->Refine();
	}
	assert(size == _grid_size);


	/* Initialize the electron density on the coarse grid */

	static std::vector<double> potential = std::vector<double>(size),
							   	exp_potential = std::vector<double>(size), 
							   	update = std::vector<double>(size);
	static std::vector<double> residual = std::vector<double>(size),
								conj_dir = std::vector<double>(size),
								dir = std::vector<double>(size);

	/* Newton's method solve for the electrostatic potential assuming a Boltzmann distribution for the electrons */
	double tol2 = 1e-20, iter_max = 10;

	/* Initialization */
	for (int i=0; i<size; i++)
	{
		potential.at(i) = std::log(_densities.front().at(i));
	}

	/* Newton's method loop */
	double scaling = -0.25 * _debye_scaling/std::pow(_plasma->get_dx(), 2.);
	for (int count=0; count<iter_max; count++)
	{
		/* Inner solve : CG algorithm */

		for (int i=0; i<size; i++)
		{
			exp_potential.at(i) = std::exp(potential.at(i));
				/* Initial value for the CG algorithm : zero. */
			residual.at(i)	=  exp_potential.at(i)  +  scaling * (	
														  potential.at((i>0?i-1:size-1)) 
														+ potential.at((i<size-1?i+1:0)) 
														- 2*potential.at(i)
															     )
									- _densities.front().at(i);
			conj_dir.at(i) = residual.at(i);
		}
		std::fill(update.begin(), update.end(), 0.);


		/* Inner loop */

		double rnorm2 = 0, delta, alpha, beta;
		for (int i=0; i<size; i++)
		{
			rnorm2 += residual.at(i)*residual.at(i);
		}

		for (int count_cg=0; count_cg<iter_max; count_cg++)
		{
			/* Compute the matrix-vector product */
			delta = 0;
			for (int i=0; i<size; i++)
			{
				dir.at(i) = (exp_potential.at(i)-2*scaling)*conj_dir.at(i) 
										+ scaling * (	conj_dir.at((i>0?i-1:size-1)) 
													  + conj_dir.at((i<size-1?i+1:0)) 
														
															     			);
				delta += dir.at(i) * conj_dir.at(i);
			}

			beta = 1./rnorm2;
			alpha = rnorm2/delta;
			rnorm2 = 0;
			for (int i=0; i<size; i++)
			{
				update.at(i) += alpha * conj_dir.at(i);
				residual.at(i) -= alpha* dir.at(i);
				rnorm2 += residual.at(i)*residual.at(i);
			}
			if (rnorm2 < tol2/2.)
				break;
			beta *= rnorm2;
			for (int i=0; i<size; i++)
			{
				conj_dir.at(i) *= beta;
				conj_dir.at(i) += residual.at(i);
			}
		}

		double err=0;
		for (int i=0; i<size; i++)
		{
			err += update.at(i) * update.at(i);
			potential.at(i) -= update.at(i);
		}
		if (err < tol2)
			break;
	}
	
	/* Deduce the electron density from the Boltzmann distribution */
	for (int i=0; i<size; i++)
	{
		_densities.at(1).at(i) = std::exp(potential.at(i));
	}

	/* Finally, fill the electron velocity and thermal velocity on the fine grid */
	for (int i=0; i<size; i++)
	{
		_velocities.at(1).at(i) = _velocities.front().at(i);
	}
	std::fill(_thermal_vel.at(1).begin(), _thermal_vel.at(1).end(), _electron_thermal_vel);
}

void MacroParameterizationWavelets::Step(State & state)
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
			// Step 2: Extrapolate using EFPI
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

void MacroParameterizationWavelets::SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics)
{
	double * x_array 	= _plasma->get_x_grid_ptr();
	int * grid_size 	= _plasma->get_grid_size_ptr(); 

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

void MacroParameterizationEFPI::WriteData(std::fstream & fout)
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





