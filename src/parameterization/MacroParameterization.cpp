#include "parameterization/MacroParameterization.h"



MacroParameterization::MacroParameterization(MacroParameterization &&parameterization) :
	_plasma(std::move(parameterization._plasma)),
	_number_of_populations(parameterization._number_of_populations),
	_population_sizes(std::move(parameterization._population_sizes)),
	_reference_densities(std::move(parameterization._reference_densities)),
	_unit_charges(std::move(parameterization._unit_charges)),
	_unit_masses(std::move(parameterization._unit_masses)),
	_plasma_pulsations(std::move(parameterization._plasma_pulsations)),
	_thermal_velocities(std::move(parameterization._thermal_velocities)),
	_bin_numbers(std::move(parameterization._bin_numbers)),
	_upper_velocities(std::move(parameterization._upper_velocities)),
	_lower_velocities(std::move(parameterization._lower_velocities))
{}



MacroParameterization& MacroParameterization::operator=(MacroParameterization &&parameterization)
{
	_plasma = std::move(parameterization._plasma);
	_number_of_populations = parameterization._number_of_populations;

	_population_sizes = std::move(parameterization._population_sizes);
	_reference_densities = std::move(parameterization._reference_densities);
	_unit_charges = std::move(parameterization._unit_charges);
	_unit_masses = std::move(parameterization._unit_masses);
	_plasma_pulsations = std::move(parameterization._plasma_pulsations);
	_thermal_velocities = std::move(parameterization._thermal_velocities);
	
	_bin_numbers = std::move(parameterization._bin_numbers);
	_upper_velocities = std::move(parameterization._upper_velocities);
	_lower_velocities = std::move(parameterization._lower_velocities);
	return *this;
}

bool 	MacroParameterization::HaveVelocityDiagnostics() const
 {
 	if ((_upper_velocities.at(0)>_lower_velocities.at(0) || _thermal_velocities.at(0)>0.0)
 		&& (_upper_velocities.at(1)>_lower_velocities.at(1) || _thermal_velocities.at(1)>0.0))
 		return true;
	else
 		return false;  
 }


double MacroParameterization::GetBinStart(int population_index) const
{
	double vupper = _upper_velocities.at(population_index);
	double vlower = _lower_velocities.at(population_index);
	double vt = _thermal_velocities.at(population_index);

	if(vupper > vlower)
	{
		return vlower;
	}
	else if(vt > 0.0) 
	{
		return -5.0*vt;
	}
	else
	{
		return 0.0;
	}

}

double MacroParameterization::GetBinEnd(int population_index) const
{
	double vupper = _upper_velocities.at(population_index);
	double vlower = _lower_velocities.at(population_index);
	double vt = _thermal_velocities.at(population_index);

	if(vupper > vlower)
	{
		return vupper;
	}
	else if(vt > 0.0) 
	{
		return 5.0*vt;
	}
	else
	{
		return 0.0;
	}

}

double MacroParameterization::GetBinWidth(int population_index)	const
{
	double vupper = _upper_velocities.at(population_index);
	double vlower = _lower_velocities.at(population_index);
	double vt = _thermal_velocities.at(population_index);
	int nbins = _bin_numbers.at(population_index);

	if(vupper > vlower)
	{
		return (vupper-vlower)/static_cast<double>(nbins);
	}
	else if(vt > 0.0) 
	{
		return 10.0 * vt/static_cast<double>(nbins);  
		/*  so that the distribution goes from -5*vt to +5*vt  */
	}
	else
	{
		return 0.;
	}

}


int MacroParameterization::GetNumberOfBins(int population_index) const
{
	return _bin_numbers.at(population_index);
}