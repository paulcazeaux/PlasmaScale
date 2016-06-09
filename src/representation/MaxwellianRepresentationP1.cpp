#include "representation/MaxwellianRepresentationP1.h"


MaxwellianRepresentationP1::MaxwellianRepresentationP1(std::shared_ptr<const Plasma> plasma, double reference_density) :
MaxwellianRepresentation(plasma, reference_density) {}

void MaxwellianRepresentationP1::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights)
{
    double idx =  static_cast<double>(_grid_end)/_plasma->_length;
	this->Reset();

	static std::vector<double> velocitysq;
	velocitysq.resize(_grid_end+1);
	std::fill(velocitysq.begin(), velocitysq.end(), 0.);

	for (int i=0; i<size; i++)
	{		
		double pos = position[i]*idx;
		double weight = weights[i];
		double vel = velocity[i];

		int xbin = static_cast<int>(pos);
		if (xbin >=0 && xbin < _grid_end)
		{
			double right_weight =  (pos - static_cast<double>(xbin)) * weight;
			double left_weight = weight - right_weight;

			_density.at(xbin) += left_weight;
			_velocity.at(xbin) += left_weight * vel;
			velocitysq.at(xbin) += left_weight * vel * vel;

			_density.at(xbin+1) += right_weight;
			_velocity.at(xbin+1) += right_weight * vel;
			velocitysq.at(xbin+1) += right_weight * vel * vel;
		}
	}
	
	double n0 = _grid_end/(_plasma->_length*_reference_density);
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
    double idx =  static_cast<double>(_grid_end)/_plasma->_length;
	this->Reset();

	static std::vector<double> velocitysq;
	velocitysq.resize(_grid_end);
	std::fill(velocitysq.begin(), velocitysq.end(), 0.);

	for (int i=0; i<size; i++)
	{		
		double pos = position[i] + delay*velocity[i];
		double vel = velocity[i];
		if (pos<0)
		{
			pos = -pos;
			vel = -vel;
		}

		int xbin = static_cast<int>(pos * idx);
		if (xbin < _grid_end)
		{
			double s = pos - static_cast<double>(xbin);
			
			vel += delay*Tools::EvaluateP1Function(accfield, _plasma->find_index_on_grid(pos), _plasma->find_position_in_cell(pos));
			double weight = weights[i];
			double right_weight =  (pos*idx - xbin) * weight;
			double left_weight = weight - right_weight;

			_density.at(xbin) += left_weight;
			_velocity.at(xbin) += left_weight * vel;
			velocitysq.at(xbin) += left_weight * vel * vel;

			_density.at(xbin+1) += right_weight;
			_velocity.at(xbin+1) += right_weight * vel;
			velocitysq.at(xbin+1) += right_weight * vel * vel;
		}
	}

	double n0 = _grid_end/(_plasma->_length*_reference_density);
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
	double x0 = .5/size;
	int index = 0;

	for (int i=0; i<size; i++)
	{
		double pos = (x0+xs)*_plasma->_length;

		int bin = _plasma->find_index_on_grid(pos);
		double s = _plasma->find_position_in_cell(pos);
        double dens = Tools::EvaluateP1Function(_density, bin, s);

		/* Guard against accidents */
        if (dens == 0.0)  s = 0.5;

		while (fv > _quiet_start_icdf.at(index+1)) ++index;

		position[i] = pos;
		weights[i]  = dens;
		velocity[i] = Tools::EvaluateP1Function(_velocity, bin, s) + Tools::EvaluateP1Function(_thermal_velocity, bin, s) * (_quiet_start_vel.at(index) + dv*(fv - _quiet_start_icdf.at(index))/(_quiet_start_icdf.at(index+1) - _quiet_start_icdf.at(index)));
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