#include "representation/WaveletRepresentationP1.h"



WaveletRepresentationP1::WaveletRepresentationP1(std::shared_ptr<const Plasma> plasma, double vmax, int depth, int grid_size) :
WaveletRepresentation(plasma, vmax, depth, grid_size) {}

void WaveletRepresentationP1::Weigh(int size,
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
		double right_weight =  _plasma->find_position_in_cell(pos) * weight;
		double left_weight = weight - right_weight;

		int vbin = static_cast<int>(velocity[i]/(dt*_dv) + static_cast<double>(_number_of_bins/2));
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
}

void WaveletRepresentationP1::Load(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights)
{
	assert(size > 1);

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
		auto it_vel_left = vel.begin();
		auto it_icdf_left = icdf.at(bin).begin();
		auto it_vel_right = vel.begin();			
		auto it_icdf_right = icdf.at(bin+1 < _grid_size ? bin+1 : 0).begin();
		double dn = 1./static_cast<double>(bin_size);
		
		double xs = 0.;
		for (int i=0; i<bin_size; i++)
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
			//double xs = RandomTools::Generate_randomly_uniform(0., 1.);
			double cellpos = xs + 0.5/static_cast<double>(bin_size);

			double fv = (static_cast<double>(i) + 0.5)*dn;
			while (fv >= *(it_icdf_left+1)) 
			{
				it_vel_left++;
				it_icdf_left++;
			}
			double vel_left =  (*it_vel_left + _dv*(fv - *it_icdf_left)/(*(it_icdf_left+1) - *it_icdf_left));
			while (fv >= *(it_icdf_right+1)) 
			{
				it_vel_right++;
				it_icdf_right++;
			}
			double vel_right =  (*it_vel_right + _dv*(fv - *it_icdf_right)/(*(it_icdf_right+1) - *it_icdf_right));

			position[bin_start_index+i] = (static_cast<double>(bin)+cellpos) * dx;
			weights[bin_start_index+i]  = Tools::EvaluateP1Function(density, bin, cellpos);
			velocity[bin_start_index+i] = (1.-cellpos)*vel_left + cellpos*vel_right;

		}
		bin_start_index = bin_end_index;
	}
}
