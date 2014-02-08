#include "plasma/MacroState.h"

MacroState::MacroState(FILE *& InputDeck)
{
	MacroParameterizationFromFile initialization = MacroParameterizationFromFile(InputDeck);
	_plasma = initialization.get_plasma();
	_micro_state = std::unique_ptr<State>(new State(initialization));
	_simulation_time = _micro_state->get_simulation_time();
	_macro_iteration = std::make_shared<int>(0);

	// Using the Shay parameterization
	if (_plasma->get_number_of_populations() != 2)
	{
		throw std::runtime_error("There are " + std::to_string(_plasma->get_number_of_populations()) + " and not 2 populations as required to use the Shay parameterization.\n");
	}

	double vte = initialization.get_initial_thermal_vel(1); // Recover the electron thermal velocity
	_macro_dt = _plasma->get_dt() * _plasma->get_macro_to_micro_dt_ratio();

	_parameterization = std::unique_ptr<MacroParameterization>(new MacroParameterizationShay(initialization, vte));
	_parameterization->Initialize(*_micro_state);
}

std::ostream& operator<<( std::ostream& os, const MacroState& state)
{
	os << *state._micro_state; // Let's be lazy for now
	return os;
}


void MacroState::Step()
{
	// We suppose initially that the microstate has already been correctly loaded at the end of the previous step.
	_parameterization->Step(*_micro_state);
	/* By construction the Step() function should leave the PIC code correctly initialized for the next step */

	_micro_state->ComputeVelocityProfile();

	(*_macro_iteration)++;
	*_simulation_time = *_macro_iteration * _macro_dt;

}