#include "representation/WaveletRepresentationP1.h"

WaveletRepresentationP1::WaveletRepresentationP1(std::shared_ptr<const Plasma> plasma, double vmax, int depth, int grid_size) :
WaveletRepresentation(plasma, vmax, depth, grid_size) {}

void WaveletRepresentationP1::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights)
{
	double population_density = static_cast<double>(_grid_size)/static_cast<double>(size);
    double idv = 1/(_dv*_plasma->get_dt());
    double vmin = static_cast<double>(_number_of_bins)/(2.*idv);

	this->Reset();

	for (int i=0; i<size; i++)
	{
		double pos = position[i];
		int xbin = _plasma->find_index_on_grid(pos);
		double cellpos = _plasma->find_position_in_cell(pos);

		int vbin = static_cast<int>((vmin+velocity[i])*idv)&(_number_of_bins-1); // Periodization in v

		_histogram.at(xbin).at(vbin) += (1.-cellpos)*weights[i];
		_histogram.at((xbin+1)&(_grid_size-1)).at(vbin) += cellpos*weights[i];

	}
    
    for (auto & cell : _histogram)
    {
        for (auto & value : cell)
            value *= population_density;
    }
}

void WaveletRepresentationP1::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights,
								const double delay,
								const std::vector<double> & accfield)
{
	double population_density = static_cast<double>(_grid_size)/static_cast<double>(size);
    double idv = 1/(_dv*_plasma->get_dt());
    double vmin = static_cast<double>(_number_of_bins)/(2.*idv);
	this->Reset();

	static std::vector<std::vector<double>::iterator> periodic_helper = std::vector<std::vector<double>::iterator>(_grid_size+1);
	for (int bin = 0; bin<_grid_size; bin++)
		periodic_helper.at(bin) = _histogram.at(bin).begin();
	periodic_helper.at(_grid_size) = _histogram.at(0).begin();

	for (int i=0; i<size; i++)
	{
		double vel = velocity[i];
		double pos = position[i] + delay*vel;

		int xbin = _plasma->find_index_on_grid(pos);
		double cellpos = _plasma->find_position_in_cell(pos);
		vel += (delay+.5)*Tools::EvaluateP1Function(accfield, xbin, cellpos);

		int vbin = static_cast<int>((vmin+vel)*idv)&(_number_of_bins-1); // Periodization in v

		periodic_helper.at(xbin)[vbin] += (1.-cellpos)*weights[i];
		periodic_helper.at(xbin+1)[vbin] += cellpos*weights[i];
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
		Tools::AssembleICDF(_histogram.at(n), icdf.at(n), density.at(n));
	}
	double L = _plasma->get_length();
	double dn = 1./static_cast<double>(size);

	double xs = 0.;			
	double fv = 0.5*dn;
	int index = 0;

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
	    double weight = Tools::EvaluateP1Function(density, bin, cellpos);

		auto it_icdf_left = icdf.at(bin).begin();
		auto it_icdf_right = icdf.at(bin+1 < _grid_size ? bin+1 : 0).begin();
	        
        	// Binary search
		int index_down = 0;
		int index_up = icdf.at(bin).size();
		while (index_up - index_down > 1)
		{
			int index_mid = (index_up + index_down)/2;
			if (weight*fv < (1.-cellpos)*it_icdf_left[index_mid] + cellpos*it_icdf_right[index_mid])
				index_up = index_mid;
			else
				index_down = index_mid;
		}

		double fvdown = (1.-cellpos)*it_icdf_left[index_down] + cellpos*it_icdf_right[index_down];
		double fvup = (1.-cellpos)*it_icdf_left[index_up] + cellpos*it_icdf_right[index_up];


		position[i] = pos;
		weights[i]  = weight;
		velocity[i] = vel.at(index_down) + _dv*(weight*fv - fvdown)/(fvup - fvdown);
		fv += dn;
	}
}
