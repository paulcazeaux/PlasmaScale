#include "representation/MaxwellianRepresentation.h"

std::vector<double> MaxwellianRepresentation::_quiet_start_vel;
std::vector<double> MaxwellianRepresentation::_quiet_start_icdf;

MaxwellianRepresentation::MaxwellianRepresentation(std::shared_ptr<const Plasma> plasma, double reference_density) :
_plasma(plasma), _grid_end(_plasma->_grid_end), _reference_density(reference_density)
{
	_density 		= std::vector<double>(_grid_end+1);
	_velocity 		= std::vector<double>(_grid_end+1);
	_thermal_velocity = std::vector<double>(_grid_end+1);
}

void MaxwellianRepresentation::InitializeQuietStartArrays(int number_of_bins)
{
	_quiet_start_vel.resize(number_of_bins);
	_quiet_start_icdf.resize(number_of_bins);

	double vmax = 5.0;
	double dv = 2.0*vmax/(static_cast<double>(number_of_bins-1));
	_quiet_start_vel.front() = -vmax;
	_quiet_start_icdf.front() = 0.0;
	for (int i=1; i < number_of_bins; i++)
	{
		double vv = ((static_cast<double>(i) - 0.5)*dv - vmax);
		_quiet_start_vel.at(i) 	= static_cast<double>(i)*dv - vmax;
		_quiet_start_icdf.at(i) = _quiet_start_icdf.at(i-1) + std::exp(-0.5*vv*vv);
	}
	double norm = 1./_quiet_start_icdf.back();
	for (auto & v: _quiet_start_icdf)
		v *= norm;
}

void MaxwellianRepresentation::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights)
{
	double dt = _plasma->_dt;
    double idx =  static_cast<double>(_grid_end)/_plasma->_length;
	double n0 = idx/_reference_density*_plasma->_dx;
	this->Reset();

	static std::vector<double> velocitysq;
	velocitysq.resize(_grid_end);
	std::fill(velocitysq.begin(), velocitysq.end(), 0.);

	for (int i=0; i<size; i++)
	{
		int xbin = static_cast<int>(position[i]*idx);
		double weight = weights[i];
		double vel = velocity[i];

		_density.at(xbin) += weight;
		_velocity.at(xbin) += weight * vel;
		velocitysq.at(xbin) += weight * vel * vel;
	}

	for (int n=0; n<=_grid_end; n++)
	{
		_velocity.at(n) /= _density.at(n) * dt;
		_thermal_velocity.at(n) = std::sqrt( (velocitysq.at(n)/(_density.at(n)*dt*dt) - _velocity.at(n)*_velocity.at(n)));
		_density.at(n) *= n0;
	}
}

void MaxwellianRepresentation::Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights,
								const double delay,
								const std::vector<double> & accfield)
{
	double dt = _plasma->_dt;
    double idx =  static_cast<double>(_grid_end)/_plasma->_length;
	double n0 = idx/_reference_density*_plasma->_dx;
	this->Reset();

	static std::vector<double> velocitysq;
	velocitysq.resize(_grid_end);
	std::fill(velocitysq.begin(), velocitysq.end(), 0.);

	for (int i=0; i<size; i++)
	{
		double pos 	= (position[i] + delay*velocity[i]) * idx;
		double vel 	= velocity[i];
		if (pos < 0)
		{
			pos = -pos;
			vel = -vel;	
		}
		
		int xbin = static_cast<int>(pos);
		if (xbin < _grid_end)
		{
			vel += delay*Tools::EvaluateP1Function(accfield, _plasma->find_index_on_grid(pos), _plasma->find_position_in_cell(pos));
			double weight 	= weights[i];

			_density.at(xbin) 	+= weight;
			_velocity.at(xbin) 	+= weight * vel;
			velocitysq.at(xbin) += weight * vel * vel;
		}
	}

	for (int n=0; n<=_grid_end; n++)
	{
		if (_density.at(n) > 0)
		{
			_velocity.at(n) /= _density.at(n) * dt;
			_thermal_velocity.at(n) = std::sqrt( (velocitysq.at(n)/(_density.at(n)*dt*dt) - _velocity.at(n)*_velocity.at(n)));
			_density.at(n) *= n0;	
		}
	}
}

void MaxwellianRepresentation::Load(int size,
				std::vector<double>::iterator 	position,
				std::vector<double>::iterator  	velocity, 
				std::vector<double>::iterator 	weights)
{
	assert(size > 1);
	
	static std::vector<double> icdf;
	icdf.resize(_grid_end+1);
	icdf.front() = 0;
	for (int i=1; i<=_grid_end; i++)
		icdf.at(i) = icdf.at(i-1)+_density.at(i-1);
	for (auto & h : icdf)
		h /= icdf.back();

	double dv = _quiet_start_vel.at(1) - _quiet_start_vel.at(0);
	double dn = 1./static_cast<double>(size);	
	double fv = 0.5*dn;

	double xs = 0.;
	double x0 = 0.5/size;		
	int index = 0;

	for (int i=0; i<size; i++)
	{
		double pos = (xs+x0)*_plasma->_length;

		int bin = _plasma->find_index_on_grid(pos);
		double s = _plasma->find_position_in_cell(pos);
        double dens = Tools::EvaluateP1Function(_density, bin, s);

		/* Guard against accidents */
        if (dens == 0.0)  s = 0.5;

		while (fv > _quiet_start_icdf.at(index+1)) ++index;

		position[i] = pos;
		velocity[i] = _velocity.at(bin) + _thermal_velocity.at(bin) * (_quiet_start_vel.at(index) + dv*(fv - _quiet_start_icdf.at(index))/(_quiet_start_icdf.at(index+1) - _quiet_start_icdf.at(index)));
		weights[i]  = dens;

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

void MaxwellianRepresentation::Coarsen()
{
	_grid_end /= 2;

	_density.at(0) = .5*(_density.at(0) + _density.at(1));
	_velocity.at(0) = .5*(_velocity.at(0) + _velocity.at(1));
	_thermal_velocity.at(0) = .5*(_thermal_velocity.at(0) + _thermal_velocity.at(1));
	for (int n=1; n<_grid_end; n++)
	{
		_density.at(n) = 0.25 * _density.at(2*n-1) + 0.5 * _density.at(2*n) + 0.25 * _density.at(2*n+1);
		_velocity.at(n) = 0.25 * _velocity.at(2*n-1) + 0.5 * _velocity.at(2*n) + 0.25 * _velocity.at(2*n+1);
		_thermal_velocity.at(n) = 0.25 * _thermal_velocity.at(2*n-1) + 0.5 * _thermal_velocity.at(2*n) + 0.25 * _thermal_velocity.at(2*n+1);
	}
	_density.at(_grid_end)  = 0.25 * _density.at(2*_grid_end-1) + 0.5 * _density.at(2*_grid_end);
	_velocity.at(_grid_end) = 0.25 * _velocity.at(2*_grid_end-1) + 0.5 * _velocity.at(2*_grid_end);
	_thermal_velocity.at(_grid_end) = 0.25 * _thermal_velocity.at(2*_grid_end-1) + 0.5 * _thermal_velocity.at(2*_grid_end);

	_density.resize(_grid_end+1);
	_velocity.resize(_grid_end+1);
	_thermal_velocity.resize(_grid_end+1);
}


void MaxwellianRepresentation::Refine()
{
	_density.resize(2*_grid_end+1);
	_velocity.resize(2*_grid_end+1);
	_thermal_velocity.resize(2*_grid_end+1);

	for (int n=_grid_end; n>=0; n--)
	{
		_density.at(2*n) 	= _density.at(n);
		_velocity.at(2*n) 	= _velocity.at(n);
		_thermal_velocity.at(2*n) 	= _thermal_velocity.at(n);
	}

	_grid_end *= 2;

	for (int n=1; n<_grid_end; n+=2)
	{
		_density.at(n) 			= 0.5*(_density.at(n-1) + _density.at(n+1));
		_velocity.at(n) 		= 0.5*(_velocity.at(n-1) + _velocity.at(n+1));
		_thermal_velocity.at(n) = 0.5*(_thermal_velocity.at(n-1) + _thermal_velocity.at(n+1));
	}
}

void MaxwellianRepresentation::SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const double & thermal_velocity)
{
	_grid_end = density.size()-1;
	_density = density;
	_velocity = velocity;
	_thermal_velocity.resize(_density.size());
	for (auto & vt: _thermal_velocity) vt = thermal_velocity;
}

void MaxwellianRepresentation::SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const std::vector<double> & thermal_velocity)
{
	_grid_end = density.size()-1;
	_density = density;
	_velocity = velocity;
	_thermal_velocity = thermal_velocity;
}

void MaxwellianRepresentation::GetDensityVelocity(std::vector<double> & density, std::vector<double> & velocity)
{
	density = _density;
	velocity = _velocity;
}

void MaxwellianRepresentation::GetDensityVelocityPressure(std::vector<double> & density, std::vector<double> & velocity, std::vector<double> & pressure)
{
	density = _density;
	velocity = _velocity;
	pressure.resize(_grid_end+1);
	for (int n=0; n<=_grid_end; n++)
		pressure.at(n) = _density.at(n) * std::pow(_thermal_velocity.at(n), 2.);
}

void MaxwellianRepresentation::Reset()
{
	for (auto & d: _density) d = 0;
	for (auto & v: _velocity) v = 0;
	for (auto & vt: _thermal_velocity) vt = 0;
}

void MaxwellianRepresentation::print(std::ostream& os) const

{
	os << "Density:" << std::endl;
	for (int i=0; i<_grid_end/4; i++)	os << _density[i] << "\t";
	os << std::endl << "Velocity:" << std::endl;
	for (int i=0; i<_grid_end/4; i++)	os << _velocity[i] << "\t";
	os << std::endl << "Pressure:" << std::endl;
	for (int i=0; i<_grid_end/4; i++)	os << _density.at(i) * std::pow(_thermal_velocity.at(i), 2.) << "\t";
	os << std::endl;
}

