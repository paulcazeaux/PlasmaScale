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

	double Te = initialization.get_initial_temperature(1); // Recover the electron temperature
	_number_of_microsteps = _plasma->get_number_of_microsteps();
	_macro_to_micro_dt_ratio = _plasma->get_macro_to_micro_dt_ratio();

	_parameterization = std::unique_ptr<MacroParameterization>(new MacroParameterizationShay(initialization, Te, _number_of_microsteps));
	_parameterization->RestrictAndPushback(_micro_state.get());
}

std::ostream& operator<<( std::ostream& os, const MacroState& state)
{
	os << *state._micro_state; // Let's be lazy for now
	return os;
}


void MacroState::Step()
{
	// We suppose initially that the microstate has already been loaded at the end of the previous step.

	// Step 1: Run the fine solver
	for (int i=0; i<_number_of_microsteps; i++)
	{
		_micro_state->Step();
		_parameterization->RestrictAndPushback(_micro_state.get());
	}
	// Step 2: Extrapolate using EPFI
	_parameterization->ExtrapolateAndLift(_macro_to_micro_dt_ratio);
	_micro_state->Load(*_parameterization);
	_micro_state->ComputeVelocityProfile();

	(*_macro_iteration)++;
	*_simulation_time = *_macro_iteration * _macro_to_micro_dt_ratio * _plasma->get_dt();

}