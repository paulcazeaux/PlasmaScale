#include "representation/WaveletRepresentationP1.h"

WaveletRepresentationP1::WaveletRepresentationP1(std::shared_ptr<const Plasma> plasma, double vmin, double vmax, int depth, double reference_density) :
WaveletRepresentation(plasma, vmin, vmax, depth, reference_density) {}

void WaveletRepresentationP1::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights)
{
    double idv = 1/(_dv*_plasma->_dt);
    double vmin = _vmin*_plasma->_dt;
    double idx =  static_cast<double>(_grid_end)/_plasma->_length;

	this->Reset();

	for (int i=0; i<size; i++)
	{
		double pos = position[i] * idx;
		int xbin = static_cast<int>(pos);

		if (xbin>=0 && xbin<=_grid_end)
		{
			int vbin = static_cast<int>((velocity[i]-vmin)*idv);
			if (vbin<0) vbin = 0;
			if (vbin>=_number_of_bins) vbin = _number_of_bins-1;

			double s = pos - static_cast<double>(xbin);
			_histogram.at(xbin).at(vbin) += (1.-s)*weights[i];
			_histogram.at(xbin+1).at(vbin) += s*weights[i];
		}
	}
	for (double & c: _histogram.front())
		c *= 2;
	for (double & c: _histogram.back())
		c *= 2;
    
	double n0 = _grid_end/(_plasma->_length*_reference_density);
    for (auto & col : _histogram)
        for (auto & value : col)
            value *= n0;
}

void WaveletRepresentationP1::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights,
								const double delay,
								const std::vector<double> & accfield)
{
    double idv = 1/(_dv*_plasma->_dt);
    double vmin = _vmin*_plasma->_dt;
    double idx = 1./_plasma->_dx;
    int R = _plasma->_grid_end / _grid_end;
	this->Reset();

	for (int i=0; i<size; i++)
	{
		double pos = idx*(position[i] + delay*velocity[i]);
		double vel = velocity[i];
		if (pos < 0)
		{
			pos = -pos;
			vel = -vel;
		}

		int xbin = static_cast<int>(pos);
		if (xbin < _plasma->_grid_end)
		{
			vel += delay*Tools::EvaluateP1Function(accfield, xbin, pos - static_cast<double>(xbin));
			int vbin = static_cast<int>((vel-vmin)*idv);
			if (vbin<0) vbin = 0;
			if (vbin>=_number_of_bins) vbin = _number_of_bins-1;

			xbin /= R;
			double s = pos/R - xbin;
			_histogram.at(xbin)[vbin] += (1.-s)*weights[i];
			_histogram.at(xbin+1)[vbin] += s*weights[i];
		}
	}
	for (double & c: _histogram.front())
		c *= 2;
	for (double & c: _histogram.back())
		c *= 2;

	double n0 = _grid_end/(_plasma->_length*_reference_density);
    for (auto & cell : _histogram)
    {
        for (auto & value : cell)
            value *= n0;
    }
}

void WaveletRepresentationP1::Load(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights)
{
	assert(size > 1);


	static std::vector<double> vel;
	static std::vector<double> icdf_pos;
	static std::vector<std::vector<double> > icdf_vel;
	static std::vector<double> density;
	vel.resize(_number_of_bins+1);
	icdf_pos.resize(_grid_end+1);
	icdf_vel.resize(_grid_end+1);
	density.resize(_grid_end+1);

	/*********************************/
	/* Initialize quiet start arrays */	
	/*********************************/
	for (int i=0; i<=_number_of_bins; i++)
	{
		vel.at(i) = static_cast<double>(i)*_dv + _vmin;
	}


	for (int n=0; n<=_grid_end; n++)
	{
		Tools::AssembleICDF(_histogram.at(n), icdf_vel.at(n), density.at(n));
	}

	icdf_pos.front() = 0;
	for (int i=1; i<=_grid_end; i++)
		icdf_pos.at(i) = icdf_pos.at(i-1)+.5*(density.at(i-1)+density.at(i));
	for (auto & h : icdf_pos)
		h /= icdf_pos.back();

	/**************************/
	/* Quiet particle loading */	
	/**************************/

	double dn = 1./static_cast<double>(size);
	double fv = 0.5*dn;
	double xs = 0.;
	double x0 = 0.5*dn;

	for (int i=0; i<size; i++)
	{
		double pos = (x0+xs)*_plasma->_length;

		int bin = _plasma->find_index_on_grid(pos);
		double s = _plasma->find_position_in_cell(pos);
        double dens = Tools::EvaluateP1Function(density, bin, s);

		/* Deal with negligible contributions */
        if (dens < 1e-15)  
    	{
    		position[i] = pos; 
    		velocity[i] = 0; 
    		weights[i] = 0;
    	}
    	else
    	{
    		/**********************************/
	    	/* Binary search for the velocity */
	    	/**********************************/

			auto it_icdf_left = icdf_vel.at(bin).begin();
			auto it_icdf_right = icdf_vel.at(bin+1).begin();

			int index_down = 0;
			int index_up = _number_of_bins;
			while (index_up - index_down > 1)
			{
				int index_mid = (index_up + index_down)/2;
				if (dens*fv < (1.-s)*it_icdf_left[index_mid] + s*it_icdf_right[index_mid])
					index_up = index_mid;
				else
					index_down = index_mid;
			}

			double fvdown = (1.-s)*it_icdf_left[index_down] + s*it_icdf_right[index_down];
			double fvup = (1.-s)*it_icdf_left[index_up] + s*it_icdf_right[index_up];

			position[i] = pos;
			weights[i]  = dens;
			velocity[i] = vel.at(index_down) + _dv*(dens*fv - fvdown)/(fvup - fvdown);
    	}
        
		fv += dn;

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
	}

	double w0 = size / std::accumulate(weights, weights+size, 0.);
	for (int i=0; i<size; i++) weights[i] *= w0;
}
