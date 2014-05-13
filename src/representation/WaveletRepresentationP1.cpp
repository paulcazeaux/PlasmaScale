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
	}
	
	double L = _plasma->get_length();
	double dn = 1./static_cast<double>(size);

	double xs = 0.;			
	double fv = 0.5*dn;
	
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
		//double pos = RandomTools::Generate_randomly_uniform(0., L);
		double pos = L*xs;
		int bin = _plasma->find_index_on_grid(pos);
		double cellpos = _plasma->find_position_in_cell(pos);

		auto it_icdf_left = icdf.at(bin).begin();
		auto it_icdf_right = icdf.at(bin+1 < _grid_size ? bin+1 : 0).begin();
        double weight = Tools::EvaluateP1Function(density, bin, cellpos);

		int vindex=0;
		double fvdown = (1.-cellpos)*(*(it_icdf_left)) + cellpos*(*(it_icdf_right));
		double fvup = (1.-cellpos)*(*(++it_icdf_left)) + cellpos*(*(++it_icdf_right));
		while (weight*fv >= fvup)
		{
			vindex++;
			fvdown = fvup;
			fvup = (1.-cellpos)*(*(++it_icdf_left)) + cellpos*(*(++it_icdf_right));
		}

		position[i] = pos;
		weights[i]  = weight;
		velocity[i] = vel.at(vindex) + _dv*(weight*fv - fvdown)/(fvup - fvdown);
		fv += dn;
	}
}
