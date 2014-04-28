#include "representation/WaveletRepresentation.h"


WaveletRepresentation(std::shared_ptr<const Plasma plasma, double vmax, double depth, int grid_size) :
_plasma(plasma), _vmax(vmax), _depth(depth), _grid_size(grid_size)
{
	_number_of_bins = std::pow(2, depth);
	_dv = (2.*vmax)/std::static_cast<double>(_number_of_bins);
	_filter = "db1";

	_is_transformed = false;
	_histogram = std::vector<std::vector<double> >(_grid_size);
	_coefficients = std::vector<std::vector<double> >(_grid_size);
	for (int n=0; n<_grid_size; n++)
	{
		_histogram.at(n).resize(_number_of_bins);
		_coefficients.at(n).resize(_number_of_bins);
	}
}


void WaveletRepresentation::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>:iterator  	velocity,
								std::vector<double>::iterator 	weight)
{
	double scaling = static_cast<double>(_grid_size)/(static_cast<double>(size)*_plasma->get_dt());
	this->Reset();

	for (int i=0; i<size; i++)
	{
		double pos = position[i];
		double weight = weight[i];
		int xbin = _plasma->find_index_on_grid(pos);

		double vel = velocity[i] * scaling;
		int vbin = static_cast<int>(vel/_dv) + _number_of_bins/2;
		if (vbin < 0)
			vbin = 0;
		if (vbin > _number_of_bins-1)
			vbin = _number_of_bins-1;

		_histogram.at(xbin).at(vbin) += weight;
	}
}

void WaveletRepresentation::Load(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>:iterator  	velocity,
								std::vector<double>::iterator 	weight)
{
	assert(size > 1);

	static std::vector<double> vel = std::vector<double>(_number_of_bins);
	static std::vector<double> icdf = std::vector<std::vector<double> >(_grid_size);
	static std::vector<double> density = std::vector<double>(_grid_size);
	/* sanity checks */
	vel.resize(_number_of_bins);			
	for (int i=0; i<=number_of_bins; i++)
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
	
	double mean_bin_size = static_cast<double>(population_size) / static_cast<double>(_grid_size) ;

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
			//double xs = RandomTools::Generate_randomly_uniform(-0.5, 0.5);
			double cellpos = xs + 0.5/static_cast<double>(bin_size);

			position[bin_start_index+i] = (static_cast<double>(bin)+cellpos) * dx;
			velocity[bin_start_index+i] = (*it_vel + _dv*(fv - *it_icdf)/(*(it_icdf+1) - *it_icdf));
			weight[bin_start_index+i]   = density;
		}
		bin_start_index = bin_end_index;
	}
}

void WaveletRepresentation::Coarsen()
{
	_grid_size /= 2;
	double dis0 = this->at(0);
	for (int n=0; n<_grid_size-1; n++)
	{
		for (int i=0; i<_number_of_bins; i++)
			this->at(n).at(i) = 0.25 * this->at(2*n).at(i) + 0.5 * this->at(2*n+1).at(i) + 0.25 * this->at(2*n+2).at(i);
	}
	{
		int n=size-1;
		for (int i=0; i<_number_of_bins; i++)
			this->at(n) = 0.25 * this->at(2*n).at(i) + 0.5 * this->at(2*n+1).at(i) + 0.25 * dis0.at(i);
	}
}


void WaveletRepresentation::Refine()
{
	for (int n=_grid_size-1; n>=0; n--)
	{
		for (int i=0; i<_number_of_bins; i++)
			this->at(2*n+1).at(i) 	= this->at(n).at(i);
	}
	{
		for (int i=0; i<_number_of_bins; i++)
			this->front() 			= 0.5*(this->at(1).at(i) + this->at(2*size-1).at(i));
	}
	for (int n=1; n<size; n++)
	{
		for (int i=0; i<_number_of_bins; i++)
			this->at(2*n).at(i) 	= 0.5*(this->at(2*n-1).at(i) + this->at(2*n+1).at(i));
	}
	_grid_size *= 2;
}

void WaveletRepresentation::DWT()
{
	for (int i=0; i<_number_of_bins; i++)
		std::dwt(_histogram.at(i), _depth, _filter, _coefficients, _flag, _length);
	_is_transformed = true;
}


void WaveletRepresentation::iDWT()
{
	for (int i=0; i<_number_of_bins; i++)
		std::idwt(_coefficients.at(i), _flag, _filter, _histogram, _length);
	_is_transformed = false;
}

void WaveletRepresentation::Denoise(int n_coef)
{
	if (_is_transformed)
	{
		double thresh;
		for (int i=0; i<_number_of_bins; i++)
		{
			std::findthresh(_coefficients.at(i), n_coef, thresh);
			for (auto & coeff : _coefficients.at(i))
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
		double thresh;
		for (int i=0; i<_number_of_bins; i++)
		{
			std::findthresh(_coefficients.at(i), n_coef, thresh);
			for (auto & coeff : _coefficients.at(i))
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

void WaveletRepresentation::Reset()
{
	std::fill(_histogram, 0.);
	_is_transformed = false;
}

void WaveletRepresentation::print(std::ostream& os)
{
}
