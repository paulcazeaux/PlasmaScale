#include "parameterization/MacroParameterization.h"



MacroParameterization::MacroParameterization(MacroParameterization &&parameterization) :
	_plasma(std::move(parameterization._plasma)),
	_number_of_populations(parameterization._number_of_populations),
	_population_sizes(std::move(parameterization._population_sizes)),
	_unit_charges(std::move(parameterization._unit_charges)),
	_unit_masses(std::move(parameterization._unit_masses)),
	_plasma_pulsations(std::move(parameterization._plasma_pulsations)),
	_cyclotronic_rotation_parameters(std::move(parameterization._cyclotronic_rotation_parameters))
{}



MacroParameterization& MacroParameterization::operator=(MacroParameterization &&parameterization)
{
	_plasma = std::move(parameterization._plasma);
	_number_of_populations = parameterization._number_of_populations;

	_population_sizes = std::move(parameterization._population_sizes);
	_unit_charges = std::move(parameterization._unit_charges);
	_unit_masses = std::move(parameterization._unit_masses);
	_plasma_pulsations = std::move(parameterization._plasma_pulsations);
	_cyclotronic_rotation_parameters = std::move(parameterization._cyclotronic_rotation_parameters);
	return *this;
}