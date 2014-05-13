#include "representation/MaxwellianRepresentationP1.h"


MaxwellianRepresentationP1::MaxwellianRepresentationP1(std::shared_ptr<const Plasma> plasma, int grid_size) :
MaxwellianRepresentation(plasma, grid_size) {}

void MaxwellianRepresentationP1::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights)
{
	double dt = _plasma->get_dt();
	double population_density = static_cast<double>(_grid_size)/static_cast<double>(size);
	this->Reset();

	static std::vector<double> velocitysq;
	velocitysq.resize(_grid_size);
	std::fill(velocitysq.begin(), velocitysq.end(), 0.);

	for (int i=0; i<size; i++)
	{		
		double pos = position[i];
		double weight = weights[i];
		double vel = velocity[i];

		int xbin = _plasma->find_index_on_grid(pos);
		double right_weight =  _plasma->find_position_in_cell(pos) * weight;
		double left_weight = weight - right_weight;

		_density.at(xbin) += left_weight;
		_velocity.at(xbin) += left_weight * vel;
		velocitysq.at(xbin) += left_weight * vel * vel;
		if (xbin < _grid_size-1)
		{
			_density.at(xbin+1) += right_weight;
			_velocity.at(xbin+1) += right_weight * vel;
			velocitysq.at(xbin+1) += right_weight * vel * vel;
		}
		else
		{
			_density.at(0) += right_weight;
			_velocity.at(0) += right_weight * vel;
			velocitysq.at(0) += right_weight * vel * vel;
		}
	}
	for (int n=0; n<_grid_size; n++)
	{
		_velocity.at(n) /= _density.at(n) * dt;
		_thermal_velocity.at(n) = std::sqrt( (velocitysq.at(n)/(_density.at(n)*dt*dt) - _velocity.at(n)*_velocity.at(n)));
		_density.at(n) *= population_density;
	}
}

void MaxwellianRepresentationP1::Load(int size,
				std::vector<double>::iterator 	position,
				std::vector<double>::iterator  	velocity, 
				std::vector<double>::iterator 	weights)
{
	assert(size > 1);

	double mean_bin_size = static_cast<double>(size) / static_cast<double>(_grid_size);
	double dx = _plasma->get_dx();

	int bin_start_index = 0;
	for (int bin = 0; bin < _grid_size; bin++)
	{
		int bin_end_index = std::ceil(static_cast<double>(bin+1)*mean_bin_size-0.5);
		int bin_size = bin_end_index - bin_start_index;

		auto it_vel = _quiet_start_vel.begin();
		auto it_icdf = _quiet_start_icdf.begin();

		double dn = 1./static_cast<double>(bin_size);
		double dv = *(it_vel+1) - *it_vel;

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
			while (fv >= *(it_icdf+1)) 
			{
				it_vel++;
				it_icdf++;
			}

			*(position+bin_start_index+i) = (static_cast<double>(bin)+cellpos) * dx;
			*(velocity+bin_start_index+i) = Tools::EvaluateP1Function(_velocity, bin, cellpos) + Tools::EvaluateP1Function(_thermal_velocity, bin, cellpos) * (*it_vel + dv*(fv - *it_icdf)/(*(it_icdf+1) - *it_icdf));
			*(weights+bin_start_index+i)  = Tools::EvaluateP1Function(_density, bin, cellpos);
		}
		bin_start_index = bin_end_index;
	}
}