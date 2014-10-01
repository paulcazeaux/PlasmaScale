#include "representation/PUREHaarRepresentation.h"


PUREHaarRepresentation::PUREHaarRepresentation(std::shared_ptr<const Plasma> plasma, double vmax, int max_depth, int grid_size, int min_depth, int buffer) :
_plasma(plasma), _vmax(vmax), _max_depth(max_depth), _grid_size(grid_size), _min_depth(min_depth), _buffer(buffer)
{
	assert(_min_depth <= _max_depth);
	_number_of_bins = std::pow(2, max_depth);
	_dv = (2.*vmax)/static_cast<double>(_number_of_bins);

	_is_transformed = false;
	_histogram = std::vector<std::vector<double> >(_grid_size, std::vector<double>(_number_of_bins));
	_scaling_coefficients = std::vector<std::vector<double> >(_grid_size);
	_detail_coefficients = std::vector<std::vector<double> >(_grid_size);
	_mask = std::vector<std::vector<int> >(_grid_size, std::vector<int>(_number_of_bins, true));
}

void PUREHaarRepresentation::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights)
{
	double population_density = static_cast<double>(_grid_size)/static_cast<double>(size);
    double dt = _plasma->get_dt();
	this->Reset();

	for (int i=0; i<size; i++)
	{
		double pos = position[i];
		while (pos<0)
			pos += _plasma->get_length();
		double weight = weights[i];
		int xbin = _plasma->find_index_on_grid(pos);

		int vbin = static_cast<int>(velocity[i]/(dt*_dv) + static_cast<double>(_number_of_bins/2));
		if (vbin < 0)
			vbin = 0;
		if (vbin > _number_of_bins-1)
			vbin = _number_of_bins-1;

		_histogram.at(xbin).at(vbin) += weight;
	}
    
    for (auto & cell : _histogram)
    {
        for (auto & value : cell)
            value *= population_density;
    }
}

void PUREHaarRepresentation::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights,
								const double delay,
								const std::vector<double> & accfield)
{
	double population_density = static_cast<double>(_grid_size)/static_cast<double>(size);
    double dt = _plasma->get_dt();
	this->Reset();

	for (int i=0; i<size; i++)
	{
		double pos = position[i] + delay*velocity[i];
		while (pos<0)
			pos += _plasma->get_length();

		int xbin = _plasma->find_index_on_grid(pos);
		double cellpos = _plasma->find_position_in_cell(pos);
		
		double vel = velocity[i] + delay*Tools::EvaluateP1Function(accfield, xbin, cellpos);
		double weight = weights[i];

		int vbin = static_cast<int>(vel/(dt*_dv) + static_cast<double>(_number_of_bins/2));
		if (vbin < 0)
			vbin = 0;
		if (vbin > _number_of_bins-1)
			vbin = _number_of_bins-1;

		_histogram.at(xbin).at(vbin) += weight;
	}
    
    for (auto & cell : _histogram)
    {
        for (auto & value : cell)
            value *= population_density;
    }
}

void PUREHaarRepresentation::Load(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights)
{
	assert(size > 1);
	if (_is_transformed)
		this->Transform();

	static std::vector<double> vel = std::vector<double>(_number_of_bins+1);
	static std::vector<std::vector<double> > icdf = std::vector<std::vector<double> >(_grid_size);
	static std::vector<double> density = std::vector<double>(_grid_size);
	/* sanity checks */
	vel.resize(_number_of_bins+1);
	for (int i=0; i<=_number_of_bins; i++)
	{
		vel.at(i) = static_cast<double>(i)*_dv - _vmax;
	}
	icdf.resize(_grid_size);
	density.resize(_grid_size);

	for (int n=0; n<_grid_size; n++)
	{
		icdf.at(n).resize(_number_of_bins+1);
		icdf.at(n).front() = 0.;
		std::partial_sum(_histogram.at(n).begin(), _histogram.at(n).end(), icdf.at(n).begin()+1);

		density.at(n) = icdf.at(n).back();
		double norm = 1./density.at(n);
		for (int i=1; i<=_number_of_bins; i++)
		{
			icdf.at(n).at(i) *= norm;
		}
	}
	
	double mean_bin_size = static_cast<double>(size) / static_cast<double>(_grid_size);
	double dx = _plasma->get_dx();

	int bin_start_index = 0;
	for (int bin = 0; bin < _grid_size; bin++)
	{
		int bin_end_index = std::ceil(static_cast<double>(bin+1)*mean_bin_size-0.5);
		int bin_size = bin_end_index - bin_start_index;
		auto it_vel = vel.begin();
		auto it_icdf = icdf.at(bin).begin();
		double dn = 1./static_cast<double>(bin_size);
		
		double xs = 0.;
		for (int i=0; i<bin_size; i++)
		{
			double fv = (static_cast<double>(i) + 0.5)*dn;
			while (fv >= *(it_icdf+1)) 
			{
				it_vel++;
				it_icdf++;
			}
			/* bit-reversed scrambling to reduce the correlation with the positions */
			double xsi = 0.5;
			xs -= 0.5;
			while (xs >= 0.0)
			{
				xsi *= 0.5;
				xs -= xsi;
			} 
			xs += 2.0*xsi;
			// double xs = RandomTools::Generate_randomly_uniform(0., 1.);
			double cellpos = xs + 0.5/static_cast<double>(bin_size);

			position[bin_start_index+i] = (static_cast<double>(bin)+cellpos) * dx;
			velocity[bin_start_index+i] = (*it_vel + _dv*(fv - *it_icdf)/(*(it_icdf+1) - *it_icdf));
			weights[bin_start_index+i]   = density.at(bin);
		}
		bin_start_index = bin_end_index;
	}
    
    
}

void PUREHaarRepresentation::Coarsen()
{
	if (_is_transformed)
		this->Transform();
	_grid_size /= 2;
	std::vector<double> dis0 = _histogram.at(0);
	std::vector<int> mask0 = _mask.at(0);
	for (int n=0; n<_grid_size-1; n++)
	{
		for (int i=0; i<_number_of_bins; i++)
		{
			_histogram.at(n).at(i) = 0.25 * _histogram.at(2*n).at(i) + 0.5 * _histogram.at(2*n+1).at(i) + 0.25 * _histogram.at(2*n+2).at(i);
			_mask.at(n)[i] = _mask.at(2*n)[i] || _mask.at(2*n+1)[i] || _mask.at(2*n+2)[i];
		}
	}
	{
		int n=_grid_size-1;
		for (int i=0; i<_number_of_bins; i++)
		{
			_histogram.at(n).at(i) = 0.25 * _histogram.at(2*n).at(i) + 0.5 * _histogram.at(2*n+1).at(i) + 0.25 * dis0.at(i);
			_mask.at(n)[i] = _mask.at(2*n)[i] || _mask.at(2*n+1)[i] || mask0[i];
		}
	}

	_histogram.resize(_grid_size);
	_scaling_coefficients.resize(_grid_size);
	_detail_coefficients.resize(_grid_size);
	_mask.resize(_grid_size);
}


void PUREHaarRepresentation::Refine()
{
	if (_is_transformed)
		this->Transform();

	_histogram.resize(2*_grid_size);
	_mask.resize(2*_grid_size);
	for (int n=_grid_size; n<2*_grid_size; n++)
    {
		_histogram.at(n).resize(_number_of_bins);
        _mask.at(n).resize(_number_of_bins);
    }

	_scaling_coefficients.resize(2*_grid_size);
	_detail_coefficients.resize(2*_grid_size);

	for (int n=_grid_size-1; n>=0; n--)
	{
		for (int i=0; i<_number_of_bins; i++)
		{
			_histogram.at(2*n+1).at(i) 	= _histogram.at(n).at(i);
			_mask.at(2*n+1)[i]		= _mask.at(n)[i];
		}
	}
	{
		for (int i=0; i<_number_of_bins; i++)
		{
			_histogram.front().at(i) 	= 0.5*(_histogram.at(1).at(i) + _histogram.at(2*_grid_size-1).at(i));
			_mask.front()[i] 		= _mask.at(1)[i] || _mask.at(2*_grid_size-1)[i];
		}
	}
	for (int n=1; n<_grid_size; n++)
	{
		for (int i=0; i<_number_of_bins; i++)
		{
			_histogram.at(2*n).at(i) 	= 0.5*(_histogram.at(2*n-1).at(i) + _histogram.at(2*n+1).at(i));
			_mask.at(2*n)[i] 		= _mask.at(2*n-1)[i] || _mask.at(2*n+1)[i];
		}
	}
	_grid_size *= 2;
}

void PUREHaarRepresentation::Transform()
{
	if (_min_depth >= _max_depth)
	{
		/* Do nothing */
		std::cerr << "Warning, no transform when min_depth = max_depth for the Haar representation !" << std::endl;
		return;
	}
	if (_is_transformed)
	{
		/* Perform inverse transform */
		for (int n=0; n<_grid_size; n++)
		{
			HaarTools::InverseTransform(_histogram.at(n), _scaling_coefficients.at(n), _detail_coefficients.at(n), _max_depth, _min_depth);
		}
		_is_transformed = false;
	}
	else
	{
		/* Perform forward unnormalized transform */
		for (int n=0; n<_grid_size; n++)
		{
			HaarTools::ForwardTransform(_histogram.at(n), _scaling_coefficients.at(n), _detail_coefficients.at(n), _max_depth, _min_depth);
		}
		_is_transformed = true;
	}
}

void PUREHaarRepresentation::SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const double & thermal_velocity)
{
	// TODO
}

void PUREHaarRepresentation::SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const std::vector<double> & thermal_velocity)
{
	// TODO
}


void PUREHaarRepresentation::GetDensityVelocity(std::vector<double> & density, std::vector<double> & velocity) const
{
	density.resize(_grid_size);
	velocity.resize(_grid_size);

	assert(!_is_transformed);
	for (int n=0; n<_grid_size; n++)
	{
		double mean = 0.;
		double dens = 0.;
		double v = -_vmax + 0.5*_dv;
		for (int i=0; i<_number_of_bins; i++)
		{
			dens += _histogram.at(n).at(i);
			mean += _histogram.at(n).at(i) * v;
			v += _dv;
		}
		density.at(n) = dens;
		if (dens>0)
			velocity.at(n) = mean / dens;
		else
			velocity.at(n) = 0.;
	}
}

void PUREHaarRepresentation::GetDensityVelocityPressure(std::vector<double> & density, std::vector<double> & velocity, std::vector<double> & pressure) const
{
	density.resize(_grid_size);
	velocity.resize(_grid_size);
    pressure.resize(_grid_size);

	assert(!_is_transformed);
	for (int n=0; n<_grid_size; n++)
	{
		double mean = 0., dens = 0., sq = 0.;
		double v = -_vmax + 0.5*_dv;
		for (int i=0; i<_number_of_bins; i++)
		{
			dens += _histogram.at(n).at(i);
			mean += _histogram.at(n).at(i) * v;
			sq	 += _histogram.at(n).at(i) * v*v;
			v += _dv;
		}
		density.at(n) = dens;
		if (dens>0)
			velocity.at(n) = mean / dens;
		else
			velocity.at(n) = 0.;
		pressure.at(n) = sq - mean*velocity.at(n);
	}
}

void PUREHaarRepresentation::Denoise(int n_coef)
{
    		/* DWT */
    if (!_is_transformed)
		this->Transform();
	std::vector<double> help;
	for (int n=0; n<_grid_size; n++)
	{
		/* Find the correct threshold */
		help = _detail_coefficients.at(n);
		std::sort(help.begin(), help.end(), greater<double>());
		double thresh = help.at(n_coef-1);

		for (auto & coeff : _detail_coefficients.at(n))
		{
			if (abs(coeff) < thresh)
			{
				coeff = 0.;
			}
		}
	}
	this->Transform();
}

void PUREHaarRepresentation::Denoise(double thresh)
{
    		/* DWT */
    if (!_is_transformed)
		this->Transform();
	for (int n=0; n<_grid_size; n++)
	{
		for (auto & coeff : _detail_coefficients.at(n))
		{
			if (abs(coeff) < thresh)
			{
				coeff = 0.;
			}
		}
	}
	this->Transform();
}



void PUREHaarRepresentation::Cutoff(int depth)
{
    int cutoff = std::pow(2, depth);
    if (cutoff<_number_of_bins)
    {
    		/* DWT */
    	if (!_is_transformed)
			this->Transform();
        for (int n=0; n<_grid_size; n++)
        {
            std::fill(_detail_coefficients.at(n).begin()+cutoff, _detail_coefficients.at(n).end(), 0.);
        }

        /* iDWT */
		this->Transform();
    }
}

void PUREHaarRepresentation::PUREAdapt(const double intensity)
{
    if (intensity==0)
        return;
	if (_is_transformed)
		this->Transform();
	/* Empirical scaling : detect local empirical mean and variance */
	double alpha = intensity*Tools::ComputePoissonEmpiricalScaling(_histogram);

	/* Bandwise thresholding using PURE threshold optimization */

	for (int n=0; n<_grid_size; n++)
	{
		std::for_each(_histogram.at(n).begin(), _histogram.at(n).end(), [&](double& val) { val /= alpha;});
		HaarTools::ForwardTransform(_histogram.at(n), _scaling_coefficients.at(n), _detail_coefficients.at(n), _max_depth, _min_depth);
		HaarTools::PUREShrink(_scaling_coefficients.at(n), _detail_coefficients.at(n), _mask.at(n), _max_depth, _min_depth);
		HaarTools::InverseTransform(_histogram.at(n), _scaling_coefficients.at(n), _detail_coefficients.at(n), _max_depth, _min_depth);
		std::for_each(_histogram.at(n).begin(), _histogram.at(n).end(), [&](double& val) { val *= alpha;});
	}

	double count=0;
	for (auto& col : _mask)
	{
		for (int i = std::pow(2, _min_depth); i<_number_of_bins; i++)
		{
			if (col.at(i))
				count += 1;
		}
	}
	std::cout << "Compression ratio: " << 100*count / (_grid_size * (_number_of_bins - std::pow(2, _min_depth))) << "%"<< std::endl;
}

void PUREHaarRepresentation::CopyMask(const PUREHaarRepresentation& representation)
{
	if (this == &representation)
		return;

	assert(representation._grid_size == _grid_size);
	assert(representation._number_of_bins == _number_of_bins);
	for (int n=0; n<_grid_size; n++)
		std::copy(representation._mask.at(n).begin(), representation._mask.at(n).end(), _mask.at(n).begin());
    
}

void PUREHaarRepresentation::ResetMask()
{
	for (int n=0; n<_grid_size; n++)
	{
		std::fill(_mask.at(n).begin(), _mask.at(n).end(), true);
	}
}

void PUREHaarRepresentation::ApplyMask()
{
	for (int n=0; n<_grid_size; n++)
	{
		HaarTools::ForwardTransform(_histogram.at(n), _scaling_coefficients.at(n), _detail_coefficients.at(n), _max_depth, _min_depth);
		for (int i=0; i<_number_of_bins; i++)
		{
			if (!_mask.at(n).at(i))
				_detail_coefficients.at(n).at(i) = 0;
		}
		HaarTools::InverseTransform(_histogram.at(n), _scaling_coefficients.at(n), _detail_coefficients.at(n), _max_depth, _min_depth);
	}
}

void PUREHaarRepresentation::Reset()
{
	for (auto & hist : _histogram)
		std::fill(hist.begin(), hist.end(), 0.);
	_is_transformed = false;
}

void PUREHaarRepresentation::print(std::ostream& os) const
{
	os << "Histogram: " << _grid_size << " grid cells by " << _number_of_bins << " velocity bins" << std::endl;
    for (auto & hist : _histogram)
    {
        for (auto & value : hist)
        {
            os << value << "\t";
        }
        os << std::endl;
    }
}
