#include "History.h"

/* constuctor and destructor ============================================================ */
History::History(std::shared_ptr<const Plasma> plasma) :
		_plasma(plasma),
		_interval(1),
		_max_size(plasma->get_max_size()),
		_max_mode(plasma->get_max_mode()),
		_number_of_populations(plasma->get_number_of_populations()),
		_velocity_accumulation_interval(plasma->get_velocity_accumulation_interval())
{
	_size = 0;
	_time_array.reserve(_max_size);
	_electrostatic_energy.reserve(_max_size);
	_kinetic_energy.reserve(_max_size);
	_total_energy.reserve(_max_size);

	_electrostatic_energy_by_mode = std::vector<std::vector<double>	>(_max_mode);
	for (auto & it : _electrostatic_energy_by_mode)
		it.reserve(_max_size);

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
	for (auto & it : _electrostatic_energy_by_mode)
		it.resize(size);
	
	for (auto & it : _electrostatic_energy_by_mode)
		it.resize(size);
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

	(state._fields)->PushBackElectrostaticEnergy(_electrostatic_energy, _electrostatic_energy_by_mode);
	(state._populations)->PushBackKineticEnergy(_kinetic_energy, _kinetic_energy_by_population);
	(state._populations)->PushBackMoment(_moment, _moment_by_population);

	_total_energy.push_back(_electrostatic_energy.back() + _kinetic_energy.back());
	_size++;
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
		for (auto & it : _electrostatic_energy_by_mode)
			it[i] = it[k];
	}
	_interval*=4;
	this->Resize(_size/4);
}

void 	History::SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics)
{
	diagnostics.emplace_back(new CurveDiagnostic("linlin", "Time", "Kinetic Energy(t)", 800, 200));
	for (int i=0; i < _number_of_populations; i++)
		diagnostics.back()->AddData(_time_array.data(), _kinetic_energy_by_population.at(i).data(), &_size, i);
	diagnostics.back()->AddData(_time_array.data(), _kinetic_energy.data(), &_size, 4);

	diagnostics.emplace_back(new CurveDiagnostic("linlog", "Time", "Field Energy(t)", 700, 400));
	diagnostics.back()->AddData(_time_array.data(), _electrostatic_energy.data(), &_size, 4);

	char buffer[25];
	for (int m=0; m<_max_mode; m++)
	{
		sprintf(buffer, "Mode %d ESE", m+1);
		diagnostics.emplace_back(new CurveDiagnostic("linlog", "Time", buffer, 700, 400));
		diagnostics.back()->AddData(_time_array.data(), _electrostatic_energy_by_mode.at(m).data(), &_size, 4);
	}

	diagnostics.emplace_back(new CurveDiagnostic("linlin", "Time", "Total Energy(t)", 700, 600));
	diagnostics.back()->AddData(_time_array.data(), _total_energy.data(), &_size, 4);
}





