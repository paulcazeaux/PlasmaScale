#include "parameterization/MacroParameterizationFullPICtoHistogram.h"


MacroParameterizationFullPICtoHistogram::MacroParameterizationFullPICtoHistogram(MacroParameterizationFullPICtoHistogram &&parameterization) :
			MacroParameterization(std::move(parameterization)),
			_grid_size(parameterization._grid_size),
			_number_of_bins(parameterization._number_of_bins),
			_min_depth(std::move(parameterization._min_depth)),
			_max_depth(std::move(parameterization._max_depth)),
			_ion_vmax(std::move(parameterization._ion_vmax)),
			_histogram(std::move(parameterization._histogram))
		{}

MacroParameterizationFullPICtoHistogram& MacroParameterizationFullPICtoHistogram::operator=(MacroParameterizationFullPICtoHistogram &&parameterization)
{
	if (this == &parameterization)
		return *this;
	
	MacroParameterization::operator=(std::move(parameterization));
	_grid_size = parameterization._plasma->get_grid_size();
	_min_depth = std::move(parameterization._min_depth);
	_max_depth = std::move(parameterization._max_depth);
    _number_of_bins = std::pow(2, _max_depth);
	_ion_vmax = std::move(parameterization._ion_vmax);
	_histogram = std::move(parameterization._histogram);

	return *this;
}

MacroParameterizationFullPICtoHistogram::MacroParameterizationFullPICtoHistogram(MacroParameterization & parameterization, double ion_vmax) :
	MacroParameterization(std::move(parameterization)), _ion_vmax(ion_vmax)
{
	_grid_size = _plasma->get_macro_grid_size();
	_min_depth = _plasma->get_wavelet_cutoff();
	_max_depth = _plasma->get_wavelet_depth();
    _number_of_bins = std::pow(2, _max_depth);

	_histogram = std::vector<std::vector<double> >(_grid_size, std::vector<double>(_number_of_bins));
}

void MacroParameterizationFullPICtoHistogram::Initialize(State & state)
{
	this->ComputeVariables(state);
}

void MacroParameterizationFullPICtoHistogram::ComputeVariables(const State & state)
{
	int p = 1;
	std::vector<double>::iterator 	position 		= state.get_vector_of_position_arrays().at(p)->begin();
	std::vector<double>::iterator	velocity 		= state.get_vector_of_x_velocity_arrays().at(p)->begin();
	std::vector<double>::iterator 	weights 		= state.get_vector_of_weight_arrays().at(p)->begin();
	int size = state.get_vector_of_position_arrays().at(p)->size();
	double dv = (2.*_ion_vmax)/static_cast<double>(_number_of_bins);
	double dt = _plasma->get_dt();

    _grid_size = _plasma->get_grid_size();	
    double population_density = static_cast<double>(_grid_size)/static_cast<double>(size);

	_histogram.resize(_grid_size);
    for (auto & col : _histogram)
    {
        col.resize(_number_of_bins);
		std::fill(col.begin(), col.end(), 0.);
    }

	for (int i=0; i<size; i++)
	{
		double pos = position[i];
		double weight = weights[i];

		int xbin = _plasma->find_index_on_grid(pos);
		double right_weight =  _plasma->find_position_in_cell(pos) * weight;
		double left_weight = weight - right_weight;

		int vbin = static_cast<int>(velocity[i]/(dt*dv) + static_cast<double>(_number_of_bins/2));
		if (vbin < 0)
			vbin = 0;
		if (vbin > _number_of_bins-1)
			vbin = _number_of_bins-1;
		_histogram.at(xbin).at(vbin) += left_weight;
		if (xbin < _grid_size-1)
		{
			_histogram.at(xbin+1).at(vbin) += right_weight;
		}
		else
		{
			_histogram.at(0).at(vbin) += right_weight;
		}
	}
    
    for (auto & cell : _histogram)
    {
        for (auto & value : cell)
            value *= population_density;
    }

    while (_grid_size>_plasma->get_macro_grid_size())
    	this->Coarsen();
}

/*
void MacroParameterizationFullPICtoHistogram::ComputeVariables(const State & state)
{
	std::vector<double>::iterator 	position 		= state.get_vector_of_position_arrays().front()->begin();
	std::vector<double>::iterator	velocity 		= state.get_vector_of_x_velocity_arrays().front()->begin();
	std::vector<double>::iterator 	weights 		= state.get_vector_of_weight_arrays().front()->begin();
	int size = state.get_vector_of_position_arrays().front()->size();
	double dv = (2.*_ion_vmax)/static_cast<double>(_number_of_bins);
	double dt = _plasma->get_dt();

    _grid_size = _plasma->get_grid_size();	
    double population_density = static_cast<double>(_grid_size)/static_cast<double>(size);

	_histogram.resize(_grid_size);
    for (auto & col : _histogram)
    {
        col.resize(_number_of_bins);
		std::fill(col.begin(), col.end(), 0.);
    }

	for (int i=0; i<size; i++)
	{
		double pos = position[i];
		double weight = weights[i];

		int xbin = _plasma->find_index_on_grid(pos);
		double right_weight =  _plasma->find_position_in_cell(pos) * weight;
		double left_weight = weight - right_weight;

		int vbin = static_cast<int>(velocity[i]/(dt*dv) + static_cast<double>(_number_of_bins/2));
		if (vbin < 0)
			vbin = 0;
		if (vbin > _number_of_bins-1)
			vbin = _number_of_bins-1;
		_histogram.at(xbin).at(vbin) += left_weight;
		if (xbin < _grid_size-1)
		{
			_histogram.at(xbin+1).at(vbin) += right_weight;
		}
		else
		{
			_histogram.at(0).at(vbin) += right_weight;
		}
	}
    
    for (auto & cell : _histogram)
    {
        for (auto & value : cell)
            value *= population_density;
    }

    while (_grid_size>_plasma->get_macro_grid_size())
    	this->Coarsen();
}
*/

void MacroParameterizationFullPICtoHistogram::Coarsen()
{
	_grid_size /= 2;
	std::vector<double> dis0 = _histogram.at(0);
	for (int n=0; n<_grid_size-1; n++)
	{
		for (int i=0; i<_number_of_bins; i++)
			_histogram.at(n).at(i) = 0.25 * _histogram.at(2*n).at(i) + 0.5 * _histogram.at(2*n+1).at(i) + 0.25 * _histogram.at(2*n+2).at(i);
	}
	{
		int n=_grid_size-1;
		for (int i=0; i<_number_of_bins; i++)
			_histogram.at(n).at(i) = 0.25 * _histogram.at(2*n).at(i) + 0.5 * _histogram.at(2*n+1).at(i) + 0.25 * dis0.at(i);
	}

	_histogram.resize(_grid_size);
}


void MacroParameterizationFullPICtoHistogram::Step(State & state)
{	
	for (int i=0; i< _plasma->get_macro_to_micro_dt_ratio(); i++)
	{
		state.Step();
	}
	this->ComputeVariables(state);
}

void MacroParameterizationFullPICtoHistogram::SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics) {}


void MacroParameterizationFullPICtoHistogram::WriteData(std::fstream & fout)
{
	
	fout << "Histogram: " << _grid_size << " grid cells by " << _number_of_bins << " velocity bins" << std::endl;
    for (auto & hist : _histogram)
    {
        for (auto & value : hist)
        {
            fout << value << "\t";
        }
        fout << std::endl;
    }
}





