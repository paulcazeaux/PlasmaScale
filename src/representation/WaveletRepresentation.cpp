#include "representation/WaveletRepresentation.h"


WaveletRepresentation::WaveletRepresentation(std::shared_ptr<const Plasma> plasma, double vmax, int depth, int grid_size) :
_plasma(plasma), _vmax(vmax), _depth(depth), _grid_size(grid_size)
{
	_number_of_bins = std::pow(2, depth);
	_dv = (2.*vmax)/static_cast<double>(_number_of_bins);
	/* R-Coiflet 6 
	_low_pass = std::vector<double> { -.000152916987, .000315249697, .001443474332,
									-.001358589300, -.007915890196, .006194347829,
									.025745731466, -.039961569717, -.049807716931,
    								.269094527854, .558133106629, .322997271647, 
    								-.040303265359, -.069655118535, .015323777973, 
    								.013570199856, -.002466300927, -.001196319329 };
	_high_pass =  std::vector<double> { -.001196319329, .002466300927,  .013570199856,
									-.015323777973, -.069655118535, .040303265359, 
    								.322997271647, -.558133106629,  .269094527854, 
									.049807716931, -.039961569717, -.025745731466, 
    								 .006194347829, .007915890196, -.001358589300,
    								 -.001443474332, .000315249697,  .000152916987};

	for (auto & v : _low_pass)
		v *= std::pow(2, 0.5);
	for (auto & v : _high_pass)
		v *= std::pow(2, 0.5);

	_filter_length = 18;
    */
    /* R-Coiflet 2 
    _low_pass = std::vector<double> {0.038580777747887, 
									-0.12696912539621, 
									-0.077161555495774, 
									0.60749164138568, 
									0.74568755893443, 
									0.22658426519707};
	_high_pass =  std::vector<double> {0.22658426519707,
									-0.74568755893443,
									0.60749164138568,
									0.077161555495774,
									-0.12696912539621,
									-0.038580777747887};
	_filter_length = 6;
	   */

	/* R-Coiflet 4 */
	_low_pass = std::vector<double> {0.0011945726958388,
									 -0.012845579755324,
									  0.024804330519353,
									  0.050023519962135,
									 -0.15535722285996,
									 -0.071638282295294,
									  0.57046500145033,
									  0.75033630585287,
									  0.28061165190244,
									 -0.0074103835186718,
									 -0.014611552521451,
									 -0.0013587990591632};
	_high_pass = std::vector<double> {-0.0013587990591632,
									  0.014611552521451,
									 -0.0074103835186718,
									 -0.28061165190244,
									  0.75033630585287,
									 -0.57046500145033,
									 -0.071638282295294,
									  0.15535722285996,
									  0.050023519962135,
									 -0.024804330519353,
									 -0.012845579755324,
									 -0.0011945726958388};
	_filter_length = 12;
	/*  */


	_is_transformed = false;
	_histogram = std::vector<std::vector<double> >(_grid_size, std::vector<double>(_number_of_bins));
	_coefficients = std::vector<std::vector<double> >(_grid_size);
	_length = std::vector<int>(_number_of_bins);
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

void WaveletRepresentation::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights,
								const double delay,
								const std::vector<double> & accfield)
{
	assert(accfield.size() == _grid_size);
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
		Tools::AssembleICDF(_histogram.at(n), icdf.at(n), density.at(n));
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
		double weight = density.at(bin);
		
		double xs = 0.;
		for (int i=0; i<bin_size; i++)
		{
			double fv = (static_cast<double>(i) + 0.5)*dn;

			while (weight*fv >= *(it_icdf+1)) 
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
			velocity[bin_start_index+i] = (*it_vel + _dv*(weight*fv - *it_icdf)/(*(it_icdf+1) - *it_icdf));
			weights[bin_start_index+i]   = weight;
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
}


void WaveletRepresentation::Refine()
{
	if (_is_transformed)
		this->iDWT();

	_histogram.resize(2*_grid_size);
	for (int n=_grid_size; n<2*_grid_size; n++)
		_histogram.at(n).resize(_number_of_bins);

	_coefficients.resize(2*_grid_size);

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
		if (dens>0)
			velocity.at(n) = mean / dens;
		else
			velocity.at(n) = 0.;
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
		if (dens>0)
			velocity.at(n) = mean / dens;
		else
			velocity.at(n) = 0.;
		pressure.at(n) = sq - mean*velocity.at(n);
	}
}

void WaveletRepresentation::GetVelocityMoment(double& moment)
{
	assert(!_is_transformed);
	moment = 0;
	for (int n=0; n<_grid_size; n++)
	{
		double v = -_vmax + 0.5*_dv;
		for (int i=0; i<_number_of_bins; i++)
		{
			moment += v*_histogram.at(n).at(i);
			v += _dv;
		}
	}
	moment /= _grid_size;
}


void WaveletRepresentation::DWT(const int J)
{
	static std::vector<double> tmp = std::vector<double>(_number_of_bins);
	_coefficients = _histogram;

	if (_length.size() < _number_of_bins)
		_length = std::vector<int>(_number_of_bins);
	int size = _number_of_bins;
	for (int j = 0 ; j<std::min(J, _depth); j++)
	{
		size /= 2;
		std::fill_n(_length.begin()+size, size, size);
	}
	std::fill_n(_length.begin(), size, size);

	for (int n=0; n<_grid_size; n++)
	{
		size = _number_of_bins;
		if (std::accumulate(_coefficients.at(n).begin(), _coefficients.at(n).end(), 0.) > 1e-10)
		{
			for (int j = 0 ; j<std::min(J, _depth); j++)
			{
				size = size/2;
				std::fill_n(tmp.begin(), 2*size, 0.);
				for (int i=0; i<size; i++)
				{
					for (int m= 0; m<_filter_length; m++)
					{
						int r = (2*i+m)&(2*size-1);
						tmp.at(i) += _low_pass.at(m)*_coefficients.at(n).at(r);
						tmp.at(size+i) += _high_pass.at(m)*_coefficients.at(n).at(r);
					}
				}
				std::copy(tmp.begin(), tmp.begin()+2*size, _coefficients.at(n).begin());
			}
		}
			
	}
	_is_transformed = true;
}


void WaveletRepresentation::iDWT()
{
	static std::vector<double> tmp = std::vector<double>(_number_of_bins);
	_histogram = _coefficients;

	for (int n=0; n<_grid_size; n++)
	{	
		int size = _length.at(0);
		while (size < _number_of_bins)
		{
			std::fill_n(tmp.begin(), 2*size, 0.);
			for (int i=0; i<size; i++)
			{
				for (int m=0; 2*m+1<_filter_length; m++)
				{
					int r = (i-m)&(size-1);
					tmp.at(2*i) += _low_pass.at(2*m)*_histogram.at(n).at(r) + _high_pass.at(2*m)*_histogram.at(n).at(size+r);
					tmp.at(2*i+1) += _low_pass.at(2*m+1)*_histogram.at(n).at(r) + _high_pass.at(2*m+1)*_histogram.at(n).at(size+r);
				}
			}
			std::copy(tmp.begin(), tmp.begin()+2*size, _histogram.at(n).begin());
			size *= 2;
		}
	}	
	std::fill_n(_length.begin(), _number_of_bins, _number_of_bins);
	_is_transformed = false;
}

void WaveletRepresentation::Denoise(const int n_coef)
{
	if (_is_transformed)
	{
		std::vector<double> help;
		for (int n=0; n<_grid_size; n++)
		{
			/* Find the correct threshold */
			help = _coefficients.at(n);
			std::sort(help.begin(), help.end(), std::greater<double>());
			double thresh = help.at(n_coef-1);

			for (auto & coeff : _coefficients.at(n))
			{
				if (std::abs(coeff) < thresh)
				{
					coeff = 0.;
				}
			}
		}
	}
	else
	{
		this->DWT(_depth);
		std::vector<double> help;
		for (int n=0; n<_grid_size; n++)
		{
			/* Find the correct threshold */
			help = _coefficients.at(n);
			std::sort(help.begin(), help.end(), std::greater<double>());
			double thresh = help.at(n_coef-1);

			for (auto & coeff : _coefficients.at(n))
			{
				if (std::abs(coeff) < thresh)
				{
					coeff = 0.;
				}
			}
		}
		this->iDWT();
	}
}

void WaveletRepresentation::Denoise(const double thresh, const int J)
{
    int L = std::pow(2, J);
	if (_is_transformed)
	{
		for (int n=0; n<_grid_size; n++)
		{
			for (int i = std::pow(2,J); i < _number_of_bins; i++)
			{
				if (std::abs(_coefficients.at(n).at(i)) < thresh)
				{
					_coefficients.at(n).at(i) = 0.;
				}
			}
		}
	}
	else
	{
		this->DWT(_depth - J);
		if (thresh == 0.)
			for (int n=0; n<_grid_size; n++)
                std::fill(_coefficients.at(n).begin()+L, _coefficients.at(n).end(), 0.);
        else
			for (int n=0; n<_grid_size; n++)
				{
					for (int i = L; i < _number_of_bins; i++)
					{
						if (std::abs(_coefficients.at(n).at(i)) < thresh)
						{
							_coefficients.at(n).at(i) = 0.;
						}
					}
				}
		this->iDWT();
	}
}



void WaveletRepresentation::Cutoff(int cutoff)
{	
    int L = std::pow(2, cutoff);
    if (L<_number_of_bins)
    {
        if (_is_transformed)
        {
            std::vector<double> help;
            for (int n=0; n<_grid_size; n++)
            {
                std::fill(_coefficients.at(n).begin()+L, _coefficients.at(n).end(), 0.);
            }
        }
        else
        {
        	/* DWT */
            this->DWT(_depth - cutoff);

            for (int n=0; n<_grid_size; n++)
                std::fill(_coefficients.at(n).begin()+L , _coefficients.at(n).end(), 0.);
            
            /* iDWT */
			this->iDWT();
        }
    }
}

void WaveletRepresentation::DiscardNegativeValues()
{
	for (auto & hist : _histogram)
	{
		for (auto & value : hist)
			if (value < 0) value=0;
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
	if (_is_transformed)
	{
		os << "Histogram: " << _grid_size/4 << " grid cells by " << _length.front() << " velocity bins" << std::endl;
	    for (int i=3*_grid_size/8; i<5*_grid_size/8; i++)
	    {
	        for (int j=0; j<_length.front(); j++)
	        {
	            os << std::scientific << _coefficients.at(i).at(j) << "\t";
	        }
	        os << std::endl;
	    }
	}
	else
	{
		os << "Histogram: " << _grid_size/4 << " grid cells by " << _number_of_bins << " velocity bins" << std::endl;
	    for (int i=3*_grid_size/8; i<5*_grid_size/8; i++)
	    {
	        for (int j=0; j<_length.front(); j++)
	        {
	            os << std::scientific << _histogram.at(i).at(j) << "\t";
	        }
	        os << std::endl;
	    }
	}
}
