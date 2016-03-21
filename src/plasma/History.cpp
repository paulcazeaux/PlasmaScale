#include "History.h"

/* constuctor and destructor ============================================================ */
History::History(std::shared_ptr<const Plasma> plasma) :
		_plasma(plasma),
		_max_size(plasma->_max_size_history),
		_velocity_accumulation_interval(plasma->_velocity_accumulation_interval),
		_interval(1),
		_number_of_populations(plasma->_number_of_populations)
{
	_size = 0;
	_time_array.reserve(_max_size);
	_electrostatic_energy.reserve(_max_size);
	_kinetic_energy.reserve(_max_size);
	_total_energy.reserve(_max_size);

	_kinetic_energy_by_population = std::vector<std::vector<double>	>(_number_of_populations);
	_moment_by_population = std::vector<std::vector<double>	>(_number_of_populations);
	for (auto & it : _kinetic_energy_by_population)
		it.reserve(_max_size);
	for (auto & it : _moment_by_population)
		it.reserve(_max_size);
}

void History::Resize(int size)
{
	assert(size <= _max_size);

	_time_array.resize(size);
	_electrostatic_energy.resize(size);
	_kinetic_energy.resize(size);
	_total_energy.resize(size);

	for (auto & it : _moment_by_population)
		it.reserve(size);
	_size = size;
}

void History::Compute(const State& state)
{	
	assert(_size == _time_array.size());
	if (_size >= _max_size)	
		this->Comb();	/* Limit capacity reached : comb time histories */
	if (*(state._iteration) % _interval != 0) 
		return;				/* only accum every interval steps */

	_time_array.push_back(*(state._simulation_time));

	state.PushBackElectrostaticEnergy(_electrostatic_energy);
	state.PushBackKineticEnergy(_kinetic_energy, _kinetic_energy_by_population);
	state.PushBackMoment(_moment, _moment_by_population);

	_total_energy.push_back(_electrostatic_energy.back() + _kinetic_energy.back());
	_size++;
}

void History::Compute(const MacroState& state)
{
	this->Compute(*state._micro_state);
}

void History::Comb()
{
	for (int i=1, k=4; k < _size; i++, k+=4)
	{
		_time_array[i] = _time_array[k];
		_electrostatic_energy[i] = _electrostatic_energy[k];
		_kinetic_energy[i] = _kinetic_energy[k];
		_total_energy[i] = _total_energy[k];

		for (auto & it : _kinetic_energy_by_population)
			it[i] = it[k];		
		for (auto & it : _moment_by_population)
			it[i] = it[k];
	}
	_interval*=4;
	this->Resize(_size/4);
}

void 	History::SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics)
{
	diagnostics.emplace_back(new CurveDiagnostic("linlin", "Time", "Kinetic Energy(t)", 0, 670));
	for (int i=0; i < _number_of_populations; i++)
		diagnostics.back()->AddData(_time_array.data(), _kinetic_energy_by_population.at(i).data(), &_size, i);
	diagnostics.back()->AddData(_time_array.data(), _kinetic_energy.data(), &_size, 4);

	diagnostics.emplace_back(new CurveDiagnostic("linlog", "Time", "Field Energy(t)", 410, 670));
	diagnostics.back()->AddData(_time_array.data(), _electrostatic_energy.data(), &_size, 4);

	diagnostics.emplace_back(new CurveDiagnostic("linlin", "Time", "Total Energy(t)", 820, 670));
	diagnostics.back()->AddData(_time_array.data(), _total_energy.data(), &_size, 4);
}





