#include "parameterization/MacroParameterizationFullPIC.h"


MacroParameterizationFullPIC::MacroParameterizationFullPIC(MacroParameterizationFullPIC &&parameterization) :
			MacroParameterization(std::move(parameterization)),
			_grid_size(parameterization._grid_size),
			_ion_density(std::move(parameterization._ion_density)),
			_ion_thermal_velocity(std::move(parameterization._ion_thermal_velocity)),
			_ion_velocity(std::move(parameterization._ion_velocity))
		{}

MacroParameterizationFullPIC& MacroParameterizationFullPIC::operator=(MacroParameterizationFullPIC &&parameterization)
{
	if (this == &parameterization)
		return *this;
	
	MacroParameterization::operator=(std::move(parameterization));
	_grid_size = parameterization._plasma->get_grid_size();

	_ion_density = std::move(parameterization._ion_density);
	_ion_thermal_velocity = std::move(parameterization._ion_thermal_velocity);
	_ion_velocity = std::move(parameterization._ion_velocity);

	return *this;
}

MacroParameterizationFullPIC::MacroParameterizationFullPIC(MacroParameterization & parameterization) :
	MacroParameterization(std::move(parameterization))
{
	_grid_size = _plasma->get_grid_size();

	_ion_density	 		= std::vector<double>(_grid_size);
	_ion_thermal_velocity	= std::vector<double>(_grid_size);
	_ion_velocity 			= std::vector<double>(_grid_size);
}

void MacroParameterizationFullPIC::Initialize(const State & state)
{
	this->ComputeVariables(state);
}

void MacroParameterizationFullPIC::ComputeVariables(const State & state)
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

	// Now we weigh the particles and compute the moments on the fine grid

	std::fill(_ion_density.begin(), _ion_density.end(), 0.);
	std::fill(_ion_velocity.begin(), _ion_velocity.end(), 0.);

	double dt = _plasma->get_dt();
	for (int i=0; i<ion_population_size-1; i++)
	{
		int bin = bins.at(i);
		double velocity = ion_velocity->at(i) / dt;

		_ion_density.at(bin) += left_weight.at(i);
		_ion_velocity.at(bin) += left_weight.at(i) * velocity;
		if (bin < _grid_size-1)
		{
			_ion_density.at(bin+1) += right_weight.at(i);
			_ion_velocity.at(bin+1) += right_weight.at(i) * velocity;
		}
		else
		{
			_ion_density.at(0) += right_weight.at(i);
			_ion_velocity.at(0) += right_weight.at(i) * velocity;
		}
	}


	for (int bin=0; bin<_grid_size; bin++)
	{
		_ion_velocity.at(bin) /= _ion_density.at(bin);
	}

	std::fill(_ion_thermal_velocity.begin(), _ion_thermal_velocity.end(), 0.);
	for (int i=0; i<ion_population_size; i++)
	{
		int bin = bins.at(i);
		double velocity = ion_velocity->at(i) / dt;

		if (bin < _grid_size-1)
		{
			double velsquare 				= std::pow( velocity - (1.-cellpos.at(i)) * _ion_velocity.at(bin) 
											  			-cellpos.at(i) * _ion_velocity.at(bin+1), 2.0);
			_ion_thermal_velocity.at(bin  ) += left_weight.at(i) * velsquare;
			_ion_thermal_velocity.at(bin+1) += right_weight.at(i) * velsquare;
		}
		else
		{
			double velsquare 				= std::pow( velocity - (1.-cellpos.at(i)) * _ion_velocity.at(bin) 
											  		   -cellpos.at(i)  * _ion_velocity.at(0), 2.0);
			_ion_thermal_velocity.at(bin) 	+= left_weight.at(i) * velsquare;
			_ion_thermal_velocity.at(0  ) 	+= right_weight.at(i) * velsquare;
		}			
	}

	for (int bin=0; bin<_grid_size; bin++)
	{
		_ion_thermal_velocity.at(bin) = std::sqrt(_ion_thermal_velocity.at(bin) / _ion_density.at(bin));
		_ion_density.at(bin) *= ion_population_density;
	}
}


void MacroParameterizationFullPIC::Step(State & state)
{	
	for (int i=0; i< _plasma->get_macro_to_micro_dt_ratio(); i++)
	{
		state.Step();
	}
	this->ComputeVariables(state);
}

void MacroParameterizationFullPIC::SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics)
{
	double * x_array 	= _plasma->get_x_grid_ptr();
	int * grid_size 	= _plasma->get_grid_size_ptr(); 

	double dt = _plasma->get_dt();
	double length = _plasma->get_length();

	diagnostics.emplace_back(new CurveDiagnostic(
				"linlin", "X", "Ion density", 0, 340));
	diagnostics.back()->AddData(x_array, _ion_density.data(), grid_size, 2);

	diagnostics.emplace_back(new CurveDiagnostic(
				"linlin", "X", "Ion velocity", 410, 340));
	diagnostics.back()->AddData(x_array, _ion_velocity.data(), grid_size, 2);

	diagnostics.emplace_back(new CurveDiagnostic(
				"linlin", "X", "Ion thermal velocity", 820, 340));
	diagnostics.back()->AddData(x_array, _ion_thermal_velocity.data(), grid_size, 4);
}





