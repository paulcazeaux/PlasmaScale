#include "representation/MaxwellianRepresentationP1.h"


MaxwellianRepresentationP1::MaxwellianRepresentationP1(std::shared_ptr<const Plasma> plasma, double reference_density) :
MaxwellianRepresentation(plasma, reference_density) {}

void MaxwellianRepresentationP1::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights)
{
	this->Reset();

	static std::vector<double> velocitysq;
	velocitysq.resize(_grid_end+1);
	std::fill(velocitysq.begin(), velocitysq.end(), 0.);

	for (int i=0; i<size; i++)
	{		
		double pos = position[i];
		double weight = weights[i];
		double vel = velocity[i];

		int xbin = _plasma->find_index_on_grid(pos);
		if (xbin >=0 && xbin < _grid_end)
		{
			double right_weight =  _plasma->find_position_in_cell(pos) * weight;
			double left_weight = weight - right_weight;

			_density.at(xbin) += left_weight;
			_velocity.at(xbin) += left_weight * vel;
			velocitysq.at(xbin) += left_weight * vel * vel;

			_density.at(xbin+1) += right_weight;
			_velocity.at(xbin+1) += right_weight * vel;
			velocitysq.at(xbin+1) += right_weight * vel * vel;
		}
	}
	
	double n0 = 1./(_reference_density*_plasma->_dx);
	double dt = _plasma->_dt;
	for (int n=0; n<=_grid_end; n++)
	{
		if (_density.at(n)>0)
		{
			_velocity.at(n) /= _density.at(n) * dt;
			_thermal_velocity.at(n) = std::sqrt( (velocitysq.at(n)/(_density.at(n)*dt*dt) - _velocity.at(n)*_velocity.at(n)));
			_density.at(n) *= n0;
		}
	}
	_density[0] *= 2;
	_density[_grid_end] *= 2;
}

void MaxwellianRepresentationP1::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights,
								const double delay,
								const std::vector<double> & accfield)
{
	double dt = _plasma->_dt;
	this->Reset();

	static std::vector<double> velocitysq;
	velocitysq.resize(_grid_end);
	std::fill(velocitysq.begin(), velocitysq.end(), 0.);

	for (int i=0; i<size; i++)
	{		
		double pos = position[i] + delay*velocity[i];
		if (pos>=0)
		{
			int xbin = _plasma->find_index_on_grid(pos);
		
			if (xbin < _grid_end)
			{
				double s = _plasma->find_position_in_cell(pos);
				
				double vel = velocity[i] + delay*Tools::EvaluateP1Function(accfield, xbin, s);
				double weight = weights[i];
				double right_weight =  s * weight;
				double left_weight = weight - right_weight;
	
				_density.at(xbin) += left_weight;
				_velocity.at(xbin) += left_weight * vel;
				velocitysq.at(xbin) += left_weight * vel * vel;
	
				_density.at(xbin+1) += right_weight;
				_velocity.at(xbin+1) += right_weight * vel;
				velocitysq.at(xbin+1) += right_weight * vel * vel;
			}
		}
		else
		{
			pos = -pos;
			int xbin = _plasma->find_index_on_grid(pos);
		
			if (xbin < _grid_end)
			{
				double s = _plasma->find_position_in_cell(pos);
				
				double vel = - velocity[i] + delay*Tools::EvaluateP1Function(accfield, xbin, s);
				double weight = weights[i];
				double right_weight =  s * weight;
				double left_weight = weight - right_weight;
	
				_density.at(xbin) += left_weight;
				_velocity.at(xbin) += left_weight * vel;
				velocitysq.at(xbin) += left_weight * vel * vel;
	
				_density.at(xbin+1) += right_weight;
				_velocity.at(xbin+1) += right_weight * vel;
				velocitysq.at(xbin+1) += right_weight * vel * vel;
			}
		}
	}

	double n0 = 1./(_reference_density*_plasma->_dx);
	for (int n=0; n<=_grid_end; n++)
	{
		if (_density.at(n)>0)
		{
			_velocity.at(n) /= _density.at(n) * dt;
			_thermal_velocity.at(n) = std::sqrt( (velocitysq.at(n)/(_density.at(n)*dt*dt) - _velocity.at(n)*_velocity.at(n)));
			_density.at(n) *= n0;
		}
	}
	_density[0] *= 2;
	_density[_grid_end] *= 2;
}

void MaxwellianRepresentationP1::Load(int size,
				std::vector<double>::iterator 	position,
				std::vector<double>::iterator  	velocity, 
				std::vector<double>::iterator 	weights)
{
	assert(size > 1);

	static std::vector<double> icdf;
	icdf.resize(_grid_end+1);
	icdf.front() = 0;
	for (int i=1; i<=_grid_end; i++)
		icdf.at(i) = icdf.at(i-1)+.5*(_density.at(i-1)+_density.at(i));
	for (auto & h : icdf)
		h /= icdf.back();

	double dv = _quiet_start_vel.at(1) - _quiet_start_vel.at(0);
	double dn = 1./static_cast<double>(size);	
	double fv = 0.5*dn;

	double xs = 0.;		
	int index = 0;

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

		/*****************/
    	/* Binary search */
    	/*****************/

		int index_down = 0;
		int index_up = _grid_end;
		while (index_up - index_down > 1)
		{
			int index_mid = (index_up + index_down)/2;
			if (xs < icdf[index_mid])
				index_up = index_mid;
			else
				index_down = index_mid;
		}

		double xdown = icdf[index_down];
		double xup = icdf[index_up];
		double pos = _plasma->_dx*(index_down + (xs-xdown)/(xup-xdown));

		int bin = _plasma->find_index_on_grid(pos);
		double s = _plasma->find_position_in_cell(pos);

		while (fv > _quiet_start_icdf.at(index+1)) ++index;

		position[i] = pos;
		weights[i]  = 1.;
		velocity[i] = Tools::EvaluateP1Function(_velocity, bin, s) + Tools::EvaluateP1Function(_thermal_velocity, bin, s) * (_quiet_start_vel.at(index) + dv*(fv - _quiet_start_icdf.at(index))/(_quiet_start_icdf.at(index+1) - _quiet_start_icdf.at(index)));
		fv += dn;
	}
}