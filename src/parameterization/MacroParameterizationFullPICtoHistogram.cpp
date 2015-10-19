#include "parameterization/MacroParameterizationFullPICtoHistogram.h"


MacroParameterizationFullPICtoHistogram::MacroParameterizationFullPICtoHistogram(MacroParameterizationFullPICtoHistogram &&parameterization) :
			MacroParameterization(std::move(parameterization)),
			_grid_size(parameterization._grid_size),
			_number_of_bins(parameterization._number_of_bins),
			_min_depth(std::move(parameterization._min_depth)),
			_max_depth(std::move(parameterization._max_depth)),
			_ion_vmax(std::move(parameterization._ion_vmax)),
			_electron_thermal_vel(std::move(parameterization._electron_thermal_vel)),
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
	_electron_thermal_vel = std::move(parameterization._electron_thermal_vel);
	_histogram = std::move(parameterization._histogram);

	return *this;
}

MacroParameterizationFullPICtoHistogram::MacroParameterizationFullPICtoHistogram(MacroParameterization & parameterization, double electron_thermal_vel, double ion_vmax) :
	MacroParameterization(std::move(parameterization)), _ion_vmax(ion_vmax), _electron_thermal_vel(electron_thermal_vel)
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
	int p = 0;
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
		/*if (i%10 == 0)
			this->ResetElectrons(state);
			*/
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


void MacroParameterizationFullPICtoHistogram::ResetElectrons(State & state)
{
	int size = _grid_size;
	static std::vector<double> ion_density = std::vector<double>(size),
							   	electron_density = std::vector<double>(size), 
								ion_velocity = std::vector<double>(size),
								potential = std::vector<double>(size),
							   	update = std::vector<double>(size);
	static std::vector<double> residual = std::vector<double>(size),
								conj_dir = std::vector<double>(size),
								dir = std::vector<double>(size);
	double tol2 = 1e-20, iter_max = 10;

	/* Initialization */
	double dv = (2.*_ion_vmax)/static_cast<double>(_number_of_bins);
	for (int n=0; n<size; n++)
	{
		double mean = 0.;
		double dens = 0.;
		double v = -_ion_vmax + 0.5*dv;
		for (auto & val : _histogram.at(n))
		{
			dens += val;
			mean += val * v;
			v += dv;
		}
		ion_density.at(n) = dens;
		ion_velocity.at(n) = mean / dens;
		potential.at(n) = std::log(ion_density.at(n));
	}

	/* Newton's method loop */
	double scaling = - std::pow(_electron_thermal_vel/_plasma_pulsations.back(), 2.)/std::pow(_plasma->get_macro_dx(), 2.);
	for (int count=0; count<iter_max; count++)
	{
		/* Inner solve : CG algorithm */

		for (int i=0; i<size; i++)
		{
			electron_density.at(i) = std::exp(potential.at(i));
				/* Initial value for the CG algorithm : zero. */
			residual.at(i)	=  electron_density.at(i)  +  scaling * (
														  potential.at((i>0?i-1:size-1)) 
														+ potential.at((i<size-1?i+1:0)) 
														- 2*potential.at(i)
															     )
									- ion_density.at(i);
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
				dir.at(i) = (electron_density.at(i)-2*scaling)*conj_dir.at(i) 
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
	for (int i=0; i<size; i++)
	{
		electron_density.at(i) = std::exp(potential.at(i));
	}

	/* Refine velocity and density array */
	while (size < _plasma->get_grid_size())
	{
		electron_density.resize(2*size);
		ion_velocity.resize(2*size);
	
		for (int n=size-1; n>=0; n--)
		{
			electron_density.at(2*n+1) 	= electron_density.at(n);
			ion_velocity.at(2*n+1) 	= ion_velocity.at(n);
		}
		{
			electron_density.front() 	= 0.5*(electron_density.at(1) + electron_density.at(2*_grid_size-1));
			ion_velocity.front() 	= 0.5*(ion_velocity.at(1) + ion_velocity.at(2*_grid_size-1));
		}
		for (int n=1; n<size; n++)
		{
			electron_density.at(2*n) 	= 0.5*(electron_density.at(2*n-1) + electron_density.at(2*n+1));
			ion_velocity.at(2*n) 	= 0.5*(ion_velocity.at(2*n-1) + ion_velocity.at(2*n+1));
		}
		size *= 2;
	}

	/* Initialize quiet start arrays */
	std::vector<double> quiet_start_vel = std::vector<double>(_number_of_bins+1);
	std::vector<double> quiet_start_icdf = std::vector<double>(_number_of_bins+1);

	{
		double vmax = 5.0*_electron_thermal_vel;
		dv = 2.0*vmax/(static_cast<double>(_number_of_bins));
		quiet_start_vel.front() = -vmax;
		quiet_start_icdf.front() = 0.0;
		for (int i=1; i < _number_of_bins+1; i++)
		{
			double vv = ((static_cast<double>(i) - 0.5)*dv - vmax);
			quiet_start_vel.at(i) 	= static_cast<double>(i)*dv - vmax;
			quiet_start_icdf.at(i) = quiet_start_icdf.at(i-1) + std::exp(-0.5*vv*vv);
		}
		double norm = 1./quiet_start_icdf.back();
		for (int i=1; i<_number_of_bins+1; i++)
		{
			quiet_start_icdf.at(i) *= norm;
		}
	}

	/* Load the electrons */

	double L = _plasma->get_length();
	double dn = 1./static_cast<double>(size);

	double xs = 0.;			
	double fv = 0.5*dn;

	std::vector<double>::iterator position 		= state.get_vector_of_position_arrays().at(1)->begin();
	std::vector<double>::iterator velocity_x 	= state.get_vector_of_x_velocity_arrays().at(1)->begin();
	std::vector<double>::iterator velocity_y 	= state.get_vector_of_y_velocity_arrays().at(1)->begin();
	std::vector<double>::iterator weight 		= state.get_vector_of_weight_arrays().at(1)->begin();
	int population_size = *state.get_vector_of_sizes().at(1);

	for (int i=0; i<size; i++)
	{
		/* bit-reversed scrambling to reduce the correlation with the positions */
		double xsi = 0.5;
		xs -= 0.5;
		while (xs >= 0.0)
		{
			xsi *= 0.5;
			xs -= xsi;
		} 
		xs += 2.0*xsi;

		double pos = L*xs;
		int bin = _plasma->find_index_on_grid(pos);
		double cellpos = _plasma->find_position_in_cell(pos);

		auto it_vel = quiet_start_vel.begin();
		auto it_icdf = quiet_start_icdf.begin();

		while (fv >= *(it_icdf+1))
		{
			it_vel++;
			it_icdf++;
		}

		position[i] = pos;
		weight[i]  = Tools::EvaluateP1Function(electron_density, bin, cellpos);
		velocity_x[i] = Tools::EvaluateP1Function(ion_velocity, bin, cellpos) + (*it_vel + (*(it_vel+1) - *it_vel)*(fv - *it_icdf)/(*(it_icdf+1) - *it_icdf));
		fv += dn;
	}
	if (_cyclotronic_rotation_parameters.at(1) != 0.) 
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


