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

void MaxwellianRepresentationP1::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights,
								const double delay,
								const std::vector<double> & accfield)
{
	double dt = _plasma->get_dt();
	double population_density = static_cast<double>(_grid_size)/static_cast<double>(size);
	this->Reset();

	static std::vector<double> velocitysq;
	velocitysq.resize(_grid_size);
	std::fill(velocitysq.begin(), velocitysq.end(), 0.);

	for (int i=0; i<size; i++)
	{		
		double pos = position[i] + delay*velocity[i];
		while (pos<0)
			pos += _plasma->get_length();

		int xbin = _plasma->find_index_on_grid(pos);
		double cellpos = _plasma->find_position_in_cell(pos);
		
		double vel = velocity[i] + delay*Tools::EvaluateP1Function(accfield, xbin, cellpos);
		double weight = weights[i];
		double right_weight =  cellpos * weight;
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

		auto it_vel = _quiet_start_vel.begin();
		auto it_icdf = _quiet_start_icdf.begin();

		while (fv >= *(it_icdf+1))
		{
			it_vel++;
			it_icdf++;
		}

		position[i] = pos;
		weights[i]  = Tools::EvaluateP1Function(_density, bin, cellpos);
		velocity[i] = Tools::EvaluateP1Function(_velocity, bin, cellpos) + Tools::EvaluateP1Function(_thermal_velocity, bin, cellpos) * (*it_vel + (*(it_vel+1) - *it_vel)*(fv - *it_icdf)/(*(it_icdf+1) - *it_icdf));
		fv += dn;
	}
}