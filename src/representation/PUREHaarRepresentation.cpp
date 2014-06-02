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

void PUREHaarRepresentation::Coarsen()
{
	if (_is_transformed)
		this->Transform();
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
	_scaling_coefficients.resize(_grid_size);
	_detail_coefficients.resize(_grid_size);
}


void PUREHaarRepresentation::Refine()
{
	if (_is_transformed)
		this->Transform();

	_histogram.resize(2*_grid_size);
	for (int n=_grid_size; n<2*_grid_size; n++)
		_histogram.at(n).resize(_number_of_bins);

	_scaling_coefficients.resize(2*_grid_size);
	_detail_coefficients.resize(2*_grid_size);

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
			int size = std::pow(2, _min_depth);
			std::copy(_scaling_coefficients.at(n).begin(), _scaling_coefficients.at(n).begin()+size, _histogram.at(n).begin());
			for (int j=_min_depth; j<_max_depth; j++)
			{
				int scale_begin = std::pow(2, j);
				for (int i=size-1; i>=0; i--)
				{
					_histogram.at(n).at(2*i+1) = 0.5 * (_histogram.at(n).at(i) - _detail_coefficients.at(n).at(scale_begin+i));
					_histogram.at(n).at(2*i)   = 0.5 * (_histogram.at(n).at(i) + _detail_coefficients.at(n).at(scale_begin+i));
				}
				size *= 2;
			}
		}
		_is_transformed = false;
	}
	else
	{
		/* Perform forward unnormalized transform */
		for (int n=0; n<_grid_size; n++)
		{
			int size = std::pow(2, _max_depth-1);
			_scaling_coefficients.resize(size);
			_detail_coefficients.resize(size);

			for (int i=size-1; i>size/2; i--)
			{
				_scaling_coefficients.at(n).at(i) = _histogram.at(n).at(2*i) + _histogram.at(n).at(2*i+1);
				_detail_coefficients.at(n).at(i)  = _histogram.at(n).at(2*i) - _histogram.at(n).at(2*i+1);
			}
			size /= 2;

			for (int j=_max_depth-1; j>_min_depth; j--)
			{
				for (int i=size-1; i>size/2; i--)
				{
					_scaling_coefficients.at(n).at(i) = _scaling_coefficients.at(n).at(2*i) + _scaling_coefficients.at(n).at(2*i+1);
					_detail_coefficients.at(n).at(i)  = _scaling_coefficients.at(n).at(2*i) - _scaling_coefficients.at(n).at(2*i+1);
				}
				size /= 2;
			}
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
		double v = _vmax + 0.5*_dv;
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
		velocity.at(n) = mean / dens;
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

void PUREHaarRepresentation::PUREAdapt()
{
	/* We rescale the data using the local density */
	static std::vector<double> scaling = std::vector<double>(_grid_size);
	for (int n=0; n<_grid_size; n++)
	{
		scaling.at(n) = std::accumulate(_histogram.at(n).begin(), _histogram.at(n).end(), 0.);
	}

	/* Empirical scaling : detect local empirical mean and variance */
	const int patch_size = 8;
	static std::vector<double> patch_mean = std::vector<double>(_grid_size*_number_of_bins/(patch_size*patch_size));
	static std::vector<double> patch_variance = std::vector<double>(_grid_size*_number_of_bins/(patch_size*patch_size));

	auto it_mean = patch_mean.begin();
	auto it_variance = patch_variance.begin();
	for (int n=0; n<_grid_size; n+=8)
	{
		for (int i=0; i<_number_of_bins; i+=8)
		{
			double mean=0., variance=0.;
			for (int np=n; np<n+patch_size; np++)
			{
				for (int ip=i; ip<i+patch_size; ip++)
				{
					double val=_histogram.at(np).at(ip)/scaling.at(n);
					mean += val;
					variance += val*val;
				}
			}
			mean /= static_cast<double>(patch_size*patch_size);
			variance = variance/static_cast<double>(patch_size*patch_size) - mean*mean;

			*it_mean++ 		= mean;
			*it_variance++ 	= variance;
		}
	}

	double alpha = Tools::EvaluateSlope(patch_mean, patch_variance);
	for (int n=0; n<_grid_size; n++)
	{
		scaling.at(n) *= alpha;
		std::for_each(_histogram.at(n).begin(), _histogram.at(n).end(), [&](double val) { val /= scaling.at(n);});
	}
		
	/* DWT */
	if (!_is_transformed)
		this->Transform();

	/* Bandwise thresholding using PURE-LET */
	static std::vector<double> theta0 = std::vector<double>(_grid_size/2), 
								theta1 = std::vector<double>(_grid_size/2),
								theta2 = std::vector<double>(_grid_size/2);
	Eigen::Vector3d A;
	Eigen::Vector3d C;
	Eigen::Matrix3d M;
	for (int j=_min_depth; j<_max_depth; j++)
	{
		for (int n=0; n<_grid_size; n++)
		{
			int scale_size = std::pow(2, j);

			/* Fill the predictors for the threshold */
			for (int i=0; i<scale_size; i++)
			{
				double detail = _detail_coefficients.at(n).at(i+scale_size);
				theta0.at(i) = detail;
				theta1.at(i) = (1 - std::exp(-detail*detail/(12*_scaling_coefficients.at(n).at(i+scale_size)))) * detail;
			}
			theta2.at(0) = _scaling_coefficients.at(n).at(2*scale_size-1) - _scaling_coefficients.at(n).at(1+scale_size);
			for (int i=1; i<scale_size-1; i++)
				theta2.at(i) = _scaling_coefficients.at(n).at(i-1+scale_size) - _scaling_coefficients.at(n).at(i+1+scale_size);
			theta2.at(scale_size-1) = _scaling_coefficients.at(n).at(2*scale_size-2) - _scaling_coefficients.at(n).at(scale_size);

			/* Compute the linear system used to determine the LET coefficients */
			M(0,0) = std::inner_product(theta0.begin(), theta0.end(), theta0.begin(), 0.);
			M(0,1) = std::inner_product(theta0.begin(), theta0.end(), theta1.begin(), 0.);
			M(1,0) = M(0,1);
			M(0,2) = std::inner_product(theta0.begin(), theta0.end(), theta2.begin(), 0.);
			M(2,0) = M(0,2);
			M(1,1) = std::inner_product(theta0.begin(), theta1.end(), theta1.begin(), 0.);
			M(1,2) = std::inner_product(theta0.begin(), theta1.end(), theta2.begin(), 0.);
			M(2,1) = M(1,2);
			M(2,2) = std::inner_product(theta0.begin(), theta2.end(), theta2.begin(), 0.);

			C(0)=0;
			C(1)=0; 
			C(2)=0;
			for (int i=0; i<scale_size; i++)
			{
				double detail = _detail_coefficients.at(n).at(i+scale_size);
				double scale = _scaling_coefficients.at(n).at(i+scale_size);

				C(0) += 0.5*detail*detail - scale;
				C(1) += (1 - std::exp(-(detail-1)*(detail-1)/(12*(scale-1))))*(detail-1)*(detail+scale)*0.5
				   	   +(1 - std::exp(-(detail+1)*(detail+1)/(12*(scale-1))))*(detail+1)*(detail-scale)*0.5;
				C(2) += detail*theta2.at(i);
			}

			/* Solve for the LET coefficients and assemble the denoised representation */
			A = M.ldlt().solve(C);
			std::cout << "depth: " << j << " ; " << "a_0=" << A(0) << " ; " << "a_1=" << A(1) << " ; " << "a_2=" << A(2) << std::endl;
			for (int i=0; i<scale_size; i++)
			{
				_detail_coefficients.at(n).at(i+scale_size) = A(0)*theta0.at(i) + A(1)*theta1.at(i) + A(2)*theta2.at(i);
			}
		}
	}

	this->Transform();

	/* Rescale the data */
	for (int n=0; n<_grid_size; n++)
	{
		std::for_each(_histogram.at(n).begin(), _histogram.at(n).end(), [&](double val) { val *= scaling.at(n);});
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
