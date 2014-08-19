#include "representation/WaveletRepresentation.h"


WaveletRepresentation::WaveletRepresentation(std::shared_ptr<const Plasma> plasma, double vmax, int depth, int grid_size) :
_plasma(plasma), _vmax(vmax), _depth(depth), _grid_size(grid_size)
{
	_number_of_bins = std::pow(2, depth);
	_dv = (2.*vmax)/static_cast<double>(_number_of_bins);
	_filter = "db6";

	_is_transformed = false;
	_histogram = std::vector<std::vector<double> >(_grid_size, std::vector<double>(_number_of_bins));
	_coefficients = std::vector<std::vector<double> >(_grid_size);
	_flag = std::vector<std::vector<double> >(_grid_size);
	_length = std::vector<std::vector<int> >(_grid_size);
}


void WaveletRepresentation::Weigh(int size,
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

void WaveletRepresentation::Load(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights)
{
	assert(size > 1);
	if (_is_transformed)
		this->iDWT();

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
        //for (int i=0; i<_number_of_bins; i++)
        //    icdf.at(n).at(i+1) = icdf.at(n).at(i) + std::max(_histogram.at(n).at(i), 0.);
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

void WaveletRepresentation::Coarsen()
{
	if (_is_transformed)
		this->iDWT();
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
	_coefficients.resize(_grid_size);
	_flag.resize(_grid_size);
	_length.resize(_grid_size);
}


void WaveletRepresentation::Refine()
{
	if (_is_transformed)
		this->iDWT();

	_histogram.resize(2*_grid_size);
	for (int n=_grid_size; n<2*_grid_size; n++)
		_histogram.at(n).resize(_number_of_bins);

	_coefficients.resize(2*_grid_size);
	_flag.resize(2*_grid_size);
	_length.resize(2*_grid_size);

	for (int n=_grid_size-1; n>=0; n--)
	{
		for (int i=0; i<_number_of_bins; i++)
			_histogram.at(2*n+1).at(i) 	= _histogram.at(n).at(i);
	}
	{
		for (int i=0; i<_number_of_bins; i++)
			_histogram.front().at(i) 	= 0.5*(_histogram.at(1).at(i) + _histogram.at(2*_grid_size-1).at(i));
	}
	for (int n=1; n<_grid_size; n++)
	{
		for (int i=0; i<_number_of_bins; i++)
			_histogram.at(2*n).at(i) 	= 0.5*(_histogram.at(2*n-1).at(i) + _histogram.at(2*n+1).at(i));
	}
	_grid_size *= 2;
}

void WaveletRepresentation::SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const double & thermal_velocity)
{
	// TODO
}

void WaveletRepresentation::SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const std::vector<double> & thermal_velocity)
{
	// TODO
}


void WaveletRepresentation::GetDensityVelocity(std::vector<double> & density, std::vector<double> & velocity) const
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
		velocity.at(n) = mean / dens;
	}
}

void WaveletRepresentation::GetDensityVelocityPressure(std::vector<double> & density, std::vector<double> & velocity, std::vector<double> & pressure) const
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
		velocity.at(n) = mean / dens;
		pressure.at(n) = sq - mean*velocity.at(n);
	}
}


void WaveletRepresentation::DWT()
{
	for (int n=0; n<_grid_size; n++)
	{
		_coefficients.at(n).clear();
		_length.at(n).clear();
		_flag.at(n).clear();
		dwt(_histogram.at(n), _depth, _filter, _coefficients.at(n), _flag.at(n), _length.at(n));
	}
	_is_transformed = true;
}


void WaveletRepresentation::iDWT()
{
	for (int n=0; n<_grid_size; n++)
		idwt(_coefficients.at(n), _flag.at(n), _filter, _histogram.at(n), _length.at(n));
	_is_transformed = false;
}

void WaveletRepresentation::Denoise(int n_coef)
{
	if (_is_transformed)
	{
		std::vector<double> help;
		for (int n=0; n<_grid_size; n++)
		{
			/* Find the correct threshold */
			help = _coefficients.at(n);
			std::sort(help.begin(), help.end(), greater<double>());
			double thresh = help.at(n_coef-1);

			for (auto & coeff : _coefficients.at(n))
			{
				if (abs(coeff) < thresh)
				{
					coeff = 0.;
				}
			}
		}
	}
	else
	{
		this->DWT();
		std::vector<double> help;
		for (int n=0; n<_grid_size; n++)
		{
			/* Find the correct threshold */
			help = _coefficients.at(n);
			std::sort(help.begin(), help.end(), greater<double>());
			double thresh = help.at(n_coef-1);

			for (auto & coeff : _coefficients.at(n))
			{
				if (abs(coeff) < thresh)
				{
					coeff = 0.;
				}
			}
		}
		this->iDWT();
	}
}

void WaveletRepresentation::Denoise(double thresh)
{
	if (_is_transformed)
	{
		std::vector<double> help;
		for (int n=0; n<_grid_size; n++)
		{
			for (auto & coeff : _coefficients.at(n))
			{
				if (abs(coeff) < thresh)
				{
					coeff = 0.;
				}
			}
		}
	}
	else
	{
		this->DWT();
		std::vector<double> help;
		for (int n=0; n<_grid_size; n++)
		{
			for (auto & coeff : _coefficients.at(n))
			{
				if (abs(coeff) < thresh)
				{
					coeff = 0.;
				}
			}
		}
		this->iDWT();
	}
}



void WaveletRepresentation::Cutoff(int depth)
{
    int cutoff = std::pow(2, depth);
    if (cutoff<_number_of_bins)
    {
        if (_is_transformed)
        {
            std::vector<double> help;
            for (int n=0; n<_grid_size; n++)
            {
                std::fill(_coefficients.at(n).begin()+cutoff, _coefficients.at(n).end(), 0.);
            }
        }
        else
        {
        	/* DWT */
            for (int n=0; n<_grid_size; n++)
			{
				_coefficients.at(n).clear();
				_length.at(n).clear();
				_flag.at(n).clear();
				dwt(_histogram.at(n), _depth-depth, _filter, _coefficients.at(n), _flag.at(n), _length.at(n));
				assert(_length.at(n).at(1) == cutoff );
			}
            for (int n=0; n<_grid_size; n++)
            {
                std::fill(_coefficients.at(n).begin()+_length.at(n).at(1), _coefficients.at(n).end(), 0.);
            }

            /* iDWT */
			for (int n=0; n<_grid_size; n++)
				idwt(_coefficients.at(n), _flag.at(n), _filter, _histogram.at(n), _length.at(n));
        }
    }
}

void WaveletRepresentation::DiscardNegativeValues()
{
	for (auto & hist : _histogram)
	{
		//double corr = 0;
		for (auto & value : hist)
		{
			if (value < 0)
			{
				//corr -= value;
				value=0;
			}
		}
//		if (corr >0)
//		{
//			corr /= _number_of_bins;
//			for (auto & value : hist)
//				value += corr;
//		}
	}
}

void WaveletRepresentation::Reset()
{
	for (auto & hist : _histogram)
		std::fill(hist.begin(), hist.end(), 0.);
	_is_transformed = false;
}

void WaveletRepresentation::print(std::ostream& os) const
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
