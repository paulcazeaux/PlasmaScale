#include "representation/WaveletRepresentation.h"


WaveletRepresentation::WaveletRepresentation(std::shared_ptr<const Plasma> plasma, double vmin, double vmax, int depth, double reference_density) :
_plasma(plasma), _vmin(vmin), _vmax(vmax), _depth(depth), _grid_end(_plasma->_grid_end), _reference_density(reference_density)
{
	_number_of_bins = std::pow(2, depth);
	_dv = (vmax-vmin)/static_cast<double>(_number_of_bins);
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
	_histogram = std::vector<std::vector<double> >(_grid_end+1, std::vector<double>(_number_of_bins));
	_coefficients = std::vector<std::vector<double> >(_grid_end+1);
	_length = std::vector<int>(_number_of_bins);
}

void WaveletRepresentation::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights)
{
	this->Reset();

	for (int i=0; i<size; i++)
	{
		double pos = position[i];
		int xbin = _plasma->find_index_on_grid(pos);

		if (xbin>=0 && xbin<_grid_end)
		{
			double weight = weights[i];
			int vbin = static_cast<int>((velocity[i]/_plasma->_dt - _vmin)*_dv);
			if (vbin < 0)
				vbin = 0;
			if (vbin > _number_of_bins-1)
				vbin = _number_of_bins-1;

			_histogram.at(xbin).at(vbin) += weight;
		}
	}
    
    double n0 = 1./(_reference_density*_plasma->_dx);
    for (auto & cell : _histogram)
        for (auto & value : cell)
            value *= n0;
}

void WaveletRepresentation::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights,
								const double delay,
								const std::vector<double> & accfield)
{
	assert(accfield.size() == _grid_end+1);
	this->Reset();

	for (int i=0; i<size; i++)
	{
		double pos = position[i] + delay*velocity[i];
		if (pos>=0)
		{
			int xbin = _plasma->find_index_on_grid(pos);

			if (xbin>=0 && xbin<_grid_end)
			{
				double s = _plasma->find_position_in_cell(pos);
				double vel = velocity[i] + delay*Tools::EvaluateP1Function(accfield, xbin, s);

				double weight = weights[i];
				int vbin = static_cast<int>((vel/_plasma->_dt - _vmin)*_dv);
				if (vbin < 0)
					vbin = 0;
				if (vbin > _number_of_bins-1)
					vbin = _number_of_bins-1;

				_histogram.at(xbin).at(vbin) += weight;
			}
		}
		else
		{
			pos = - pos;

			int xbin = _plasma->find_index_on_grid(pos);

			if (xbin>=0 && xbin<_grid_end)
			{
				double s = _plasma->find_position_in_cell(pos);
				double vel = - velocity[i] + delay*Tools::EvaluateP1Function(accfield, xbin, s);

				double weight = weights[i];
				int vbin = static_cast<int>((vel/_plasma->_dt - _vmin)*_dv);
				if (vbin < 0)
					vbin = 0;
				if (vbin > _number_of_bins-1)
					vbin = _number_of_bins-1;

				_histogram.at(xbin).at(vbin) += weight;
			}
		}
	}

	double n0 = 1./(_reference_density*_plasma->_dx);
    for (auto & cell : _histogram)
        for (auto & value : cell)
            value *= n0;
}

void WaveletRepresentation::Load(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights)
{
	assert(size > 1);
	if (_is_transformed)
		this->iDWT();

	static std::vector<double> vel;
	static std::vector<double> icdf_pos;
	static std::vector<std::vector<double> > icdf_vel;
	static std::vector<double> density;
	vel.resize(_number_of_bins+1);
	icdf_pos.resize(_grid_end+1);
	icdf_vel.resize(_grid_end+1);

	/*********************************/
	/* Initialize quiet start arrays */	
	/*********************************/
		
	for (int i=0; i<=_number_of_bins; i++)
	{
		vel.at(i) = static_cast<double>(i)*_dv + _vmin;
	}
	for (int n=0; n<=_grid_end; n++)
	{
		double density;
		Tools::AssembleICDF(_histogram.at(n), icdf_vel.at(n), density);
		for (double & c: icdf_vel.at(n)) 
			c /= density;
	}

	icdf_pos.front() = 0;
	for (int i=1; i<=_grid_end; i++)
		icdf_pos.at(i) = icdf_pos.at(i-1)+density.at(i-1);

	for (auto & h : icdf_pos)
		h /= icdf_pos.back();

	/**************************/
	/* Quiet particle loading */	
	/**************************/

	double dn = 1./static_cast<double>(size);
	double fv = 0.5*dn;
	double xs = 0.;			

	for (int i=0; i<size; i++)
	{
		/************************************************************************/
		/* bit-reversed scrambling to reduce the correlation with the positions */
		/************************************************************************/
		
		double xsi = 0.5;
		xs -= 0.5;
		while (xs >= 0.0)
		{
			xsi *= 0.5;
			xs -= xsi;
		} 
		xs += 2.0*xsi;

		/**********************************/
    	/* Binary search for the position */
    	/**********************************/

		int index_down = 0;
		int index_up = _grid_end;
		while (index_up - index_down > 1)
		{
			int index_mid = (index_up + index_down)/2;
			if (xs < icdf_pos[index_mid])
				index_up = index_mid;
			else
				index_down = index_mid;
		}

		double xdown = icdf_pos[index_down];
		double xup = icdf_pos[index_up];
		double pos = _plasma->_dx*(index_down + (xs-xdown)/(xup-xdown));

		int bin = _plasma->find_index_on_grid(pos);

		/**********************************/
    	/* Binary search for the velocity */
    	/**********************************/

		index_down = 0;
		index_up = _number_of_bins;
		while (index_up - index_down > 1)
		{
			int index_mid = (index_up + index_down)/2;
			if (fv < icdf_vel.at(bin).at(index_mid))
				index_up = index_mid;
			else
				index_down = index_mid;
		}

		double fvdown = icdf_vel.at(bin).at(index_down);
		double fvup = icdf_vel.at(bin).at(index_up);

		position[i] = pos;
		velocity[i] = (vel.at(index_down) + _dv*(fv - fvdown)/(fvup - fvdown));
		weights[i]  = 1.;
	}
}

void WaveletRepresentation::Coarsen()
{
	_grid_end /= 2;

	if (!_is_transformed)
	{
		for (int i=0; i<_number_of_bins; i++)
		{
			_histogram.at(0).at(i) = 0.5*(_histogram.at(0).at(i) + _histogram.at(1).at(i));
			for (int n=1; n<_grid_end; n++)
					_histogram.at(n).at(i) = 0.25 * _histogram.at(2*n-1).at(i) + 0.5 * _histogram.at(2*n).at(i) + 0.25 * _histogram.at(2*n+1).at(i);

			_histogram.at(_grid_end).at(i) = 0.25 * _histogram.at(2*_grid_end-1).at(i) + 0.5 * _histogram.at(2*_grid_end).at(i);
		}
		_histogram.resize(_grid_end+1);
	}
	else
	{
		for (int i=0; i<_number_of_bins; i++)
		{
			_coefficients.at(0).at(i) = 0.5*(_coefficients.at(0).at(i) + _coefficients.at(1).at(i));

			for (int n=1; n<_grid_end; n++)
					_coefficients.at(n).at(i) = 0.25 * _coefficients.at(2*n-1).at(i) + 0.5 * _coefficients.at(2*n).at(i) + 0.25 * _coefficients.at(2*n+1).at(i);

			_coefficients.at(_grid_end).at(i) = 0.25 * _coefficients.at(2*_grid_end-1).at(i) + 0.5 * _coefficients.at(2*_grid_end).at(i);
		}
		_coefficients.resize(_grid_end+1);
	}
}


void WaveletRepresentation::Refine()
{
	if (!_is_transformed)
	{
		_histogram.resize(2*_grid_end+1);			
		for (int n=_grid_end; n>=0; n--)
				_histogram.at(2*n) = _histogram.at(n);

		for (int n=1; n<2*_grid_end; n+=2)
		{
			_histogram.at(n).resize(_number_of_bins);
			for (int i=0; i<_number_of_bins; i++)
				_histogram.at(n).at(i) 	= 0.5*(_histogram.at(n-1).at(i) + _histogram.at(n+1).at(i));
		}
	}
	else
	{
		_coefficients.resize(2*_grid_end+1);			
		for (int n=_grid_end; n>=0; n--)
				_coefficients.at(2*n+1) = _coefficients.at(n);

		for (int n=1; n<2*_grid_end; n+=2)
		{
			_coefficients.at(n).resize(_number_of_bins);
			for (int i=0; i<_number_of_bins; i++)
				_coefficients.at(n).at(i) 	= 0.5*(_coefficients.at(n-1).at(i) + _coefficients.at(n+1).at(i));
		}
	}
	
	_grid_end *= 2;
}

void WaveletRepresentation::GetDensityVelocity(std::vector<double> & density, std::vector<double> & velocity)
{
	density.resize(_grid_end+1);
	velocity.resize(_grid_end+1);

	if (_is_transformed) this->iDWT();

	for (int n=0; n<_grid_end; n++)
	{
		double mean = 0.;
		double dens = 0.;
		double v = _vmin + 0.5*_dv;
		for (double & val : _histogram.at(n))
		{
			dens += val;
			mean += val * v;
			v += _dv;
		}
		density.at(n) = dens;
		if (dens>0)
			velocity.at(n) = mean / dens;
		else
			velocity.at(n) = 0.;
	}
}

void WaveletRepresentation::GetDensityVelocityPressure(std::vector<double> & density, std::vector<double> & velocity, std::vector<double> & pressure)
{
	density.resize(_grid_end+1);
	velocity.resize(_grid_end+1);
    pressure.resize(_grid_end+1);

	if (_is_transformed) this->iDWT();

	for (int n=0; n<_grid_end; n++)
	{
		double dens = 0.;
		double mean = 0.;
		double sq = 0.;
		double v = _vmin + 0.5*_dv;
		for (double & val : _histogram.at(n))
		{
			dens += val;
			mean += val * v;
			sq	 += val * v*v;
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
	for (std::vector<double> col: _histogram)
	{
		double v = _vmin + 0.5*_dv;
		for (double val: col)
		{
			moment += v*val;
			v += _dv;
		}
	}
	moment *= _plasma->_dx;
}


void WaveletRepresentation::DWT(const int J)
{
	static std::vector<double> tmp = std::vector<double>(_number_of_bins);
	_coefficients = std::move(_histogram);

	if (_length.size() < _number_of_bins)
		_length = std::vector<int>(_number_of_bins);
	int size = _number_of_bins;
	for (int j = 0 ; j<std::min(J, _depth); j++)
	{
		size /= 2;
		std::fill_n(_length.begin()+size, size, size);
	}
	std::fill_n(_length.begin(), size, size);

	for (int n=0; n<=_grid_end; n++)
	{
		size = _number_of_bins;
		if (std::accumulate(_coefficients.at(n).begin(), _coefficients.at(n).end(), 0.) > 1e-10)
		{
			for (int j = 0 ; j<std::min(J, _depth); j++)
			{
				std::fill_n(tmp.begin(), size, 0.);
				size /=2;
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
	_histogram = std::move(_coefficients);

	for (int n=0; n<=_grid_end; n++)
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
			size *= 2;
			std::copy(tmp.begin(), tmp.begin()+size, _histogram.at(n).begin());
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
		for (int n=0; n<=_grid_end; n++)
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
		for (int n=0; n<=_grid_end; n++)
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
		for (int n=0; n<=_grid_end; n++)
		{
			for (int i = L; i < _number_of_bins; i++)
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
			for (int n=0; n<=_grid_end; n++)
                std::fill(_coefficients.at(n).begin()+L, _coefficients.at(n).end(), 0.);
        else
			for (int n=0; n<=_grid_end; n++)
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
            for (int n=0; n<=_grid_end; n++)
            {
                std::fill(_coefficients.at(n).begin()+L, _coefficients.at(n).end(), 0.);
            }
        }
        else
        {
        	/* DWT */
            this->DWT(_depth - cutoff);

            for (int n=0; n<=_grid_end; n++)
                std::fill(_coefficients.at(n).begin()+L , _coefficients.at(n).end(), 0.);
            
            /* iDWT */
			this->iDWT();
        }
    }
}

void WaveletRepresentation::DiscardNegativeValues()
{
	if (!_is_transformed)
		for (auto & col : _histogram)
			for (auto & value : col)
				if (value < 0) 
					value=0;
}

void WaveletRepresentation::Reset()
{
	if (_is_transformed)
	{
		_coefficients = std::move(_histogram);
		_is_transformed = false;
	}
	for (auto & col : _histogram)
		for (auto & val : col) 
			val = 0;
}

void WaveletRepresentation::print(std::ostream& os) const
{
	if (_is_transformed)
	{
		os << "Histogram: " << _grid_end/4 << " grid cells by " << _length.front() << " velocity bins" << std::endl;
	    for (int i=0; i<_grid_end/4; i++)
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
		os << "Histogram: " << _grid_end/4 << " grid cells by " << _number_of_bins << " velocity bins" << std::endl;
	    for (int i=0; i<_grid_end/4; i++)
	    {
	        for (int j=0; j<_length.front(); j++)
	        {
	            os << std::scientific << _histogram.at(i).at(j) << "\t";
	        }
	        os << std::endl;
	    }
	}
}
