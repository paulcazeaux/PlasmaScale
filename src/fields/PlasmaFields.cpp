#include "PlasmaFields.h"



/* Constructor =============================================================== */

PlasmaFields::PlasmaFields(std::shared_ptr<const Plasma> plasma,
					std::shared_ptr<double> simulation_time,
					std::shared_ptr<int> iteration)
					: 	_plasma(plasma),
						_grid_end  (plasma->_grid_end),
					  	_simulation_time(simulation_time),
					  	_iteration(iteration),
						_electrostatic_energy(0.)
{
	_charge = std::vector<double>(_grid_end+1);
	_potential = std::vector<double>(_grid_end+1);
	_electrical_field = std::vector<double>(_grid_end+1);
}

void PlasmaFields::Compute()
{	
	double dx = _plasma->_dx;

	/* Poisson solve */
	_plasma->Poisson_Solver(_charge, _potential);

	/* Calculate the electrostatic energy */
	_electrostatic_energy = 0;
	for (int i=0; i<_grid_end; i++)
	{
		_electrostatic_energy += (_charge.at(i)*_potential.at(i) + _charge.at(i+1)*_potential.at(i+1));
	}
	_electrostatic_energy *= 0.25*dx;

	/* Compute the electrical field */
	double external_e_field = _plasma->_E0*cos(_plasma->_w0*(*_simulation_time));

	for (int i=1; i<_grid_end; i++)
	{
		_electrical_field[i] = 0.5/dx*(_potential[i-1] - _potential[i+1]) + external_e_field;
	}
	_electrical_field[0] = (_potential[0]-_potential[1])/dx + (0.5*dx/_plasma->_epsilon)*_charge[0] + external_e_field;
	_electrical_field[_grid_end] = (_potential[_grid_end-1]-_potential[_grid_end])/dx - (0.5*dx/_plasma->_epsilon)*_charge[_grid_end] + external_e_field;
}

void PlasmaFields::PushBackElectrostaticEnergy(std::vector<double>& electrostatic_energy) const
{
	electrostatic_energy.push_back(_electrostatic_energy  + 1e-30);
}