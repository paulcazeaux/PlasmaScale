#include "representation/MaxwellianRepresentation.h"

std::vector<double> MaxwellianRepresentation::_quiet_start_vel;
std::vector<double> MaxwellianRepresentation::_quiet_start_icdf;

MaxwellianRepresentation::MaxwellianRepresentation(std::shared_ptr<const Plasma> plasma, int grid_size) :
_plasma(plasma), _grid_size(grid_size)
{
	_density 		= std::vector<double>(_grid_size);
	_velocity 		= std::vector<double>(_grid_size);
	_thermal_velocity = std::vector<double>(_grid_size);
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
	for (int i=1; i<number_of_bins; i++)
	{
		_quiet_start_icdf.at(i) *= norm;
	}
}

void MaxwellianRepresentation::Weigh(int size,
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
		int xbin = _plasma->find_index_on_grid(position[i]);
		double weight = weights[i];
		double vel = velocity[i];

		_density.at(xbin) += weight;
		_velocity.at(xbin) += weight * vel;
		velocitysq.at(xbin) += weight * vel * vel;
	}
	for (int n=0; n<_grid_size; n++)
	{
		_velocity.at(n) /= _density.at(n) * dt;
		_thermal_velocity.at(n) = std::sqrt( (velocitysq.at(n)/(_density.at(n)*dt*dt) - _velocity.at(n)*_velocity.at(n)));
		_density.at(n) *= population_density;
	}
}

void MaxwellianRepresentation::Weigh(int size,
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

		_density.at(xbin) += weight;
		_velocity.at(xbin) += weight * vel;
		velocitysq.at(xbin) += weight * vel * vel;
	}
	for (int n=0; n<_grid_size; n++)
	{
		_velocity.at(n) /= _density.at(n) * dt;
		_thermal_velocity.at(n) = std::sqrt( (velocitysq.at(n)/(_density.at(n)*dt*dt) - _velocity.at(n)*_velocity.at(n)));
		_density.at(n) *= population_density;
	}
}

void MaxwellianRepresentation::Load(int size,
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
			//double xs = RandomTools::Generate_randomly_uniform(0, 1.);
			double cellpos = xs + 0.5/static_cast<double>(bin_size);

			*(position+bin_start_index+i) = (static_cast<double>(bin)+cellpos) * dx;
			*(velocity+bin_start_index+i) = _velocity.at(bin) + _thermal_velocity.at(bin) * (*it_vel + dv*(fv - *it_icdf)/(*(it_icdf+1) - *it_icdf));
			*(weights+bin_start_index+i)  = _density.at(bin);
		}
		bin_start_index = bin_end_index;
	}
}

void MaxwellianRepresentation::Coarsen()
{
	_grid_size /= 2;
	double d0 = _density.at(0), v0 = _velocity.at(0), t0 = _thermal_velocity.at(0);
	for (int n=0; n<_grid_size-1; n++)
	{
		_density.at(n) = 0.25 * _density.at(2*n) + 0.5 * _density.at(2*n+1) + 0.25 * _density.at(2*n+2);
		_velocity.at(n) = 0.25 * _velocity.at(2*n) + 0.5 * _velocity.at(2*n+1) + 0.25 * _velocity.at(2*n+2);
		_thermal_velocity.at(n) = 0.25 * _thermal_velocity.at(2*n) + 0.5 * _thermal_velocity.at(2*n+1) + 0.25 * _thermal_velocity.at(2*n+2);
	}
	{
		int n=_grid_size-1;
		_density.at(n) = 0.25 * _density.at(2*n) + 0.5 * _density.at(2*n+1) + 0.25 * d0;
		_velocity.at(n) = 0.25 * _velocity.at(2*n) + 0.5 * _velocity.at(2*n+1) + 0.25 * v0;
		_thermal_velocity.at(n) = 0.25 * _thermal_velocity.at(2*n) + 0.5 * _thermal_velocity.at(2*n+1) + 0.25 * t0;
	}

	_density.resize(_grid_size);
	_velocity.resize(_grid_size);
	_thermal_velocity.resize(_grid_size);
}


void MaxwellianRepresentation::Refine()
{
	_density.resize(2*_grid_size);
	_velocity.resize(2*_grid_size);
	_thermal_velocity.resize(2*_grid_size);

	for (int n=_grid_size-1; n>=0; n--)
	{
		_density.at(2*n+1) 	= _density.at(n);
		_velocity.at(2*n+1) 	= _velocity.at(n);
		_thermal_velocity.at(2*n+1) 	= _thermal_velocity.at(n);
	}
	{
		_density.front() 	= 0.5*(_density.at(1) + _density.at(2*_grid_size-1));
		_velocity.front() 	= 0.5*(_velocity.at(1) + _velocity.at(2*_grid_size-1));
		_thermal_velocity.front() 	= 0.5*(_thermal_velocity.at(1) + _thermal_velocity.at(2*_grid_size-1));
	}
	for (int n=1; n<_grid_size; n++)
	{
		_density.at(2*n) 	= 0.5*(_density.at(2*n-1) + _density.at(2*n+1));
		_velocity.at(2*n) 	= 0.5*(_velocity.at(2*n-1) + _velocity.at(2*n+1));
		_thermal_velocity.at(2*n) 	= 0.5*(_thermal_velocity.at(2*n-1) + _thermal_velocity.at(2*n+1));
	}
	_grid_size *= 2;
}

void MaxwellianRepresentation::SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const double & thermal_velocity)
{
	_grid_size = density.size();
	_density = density;
	_velocity = velocity;
	_thermal_velocity.resize(_density.size());
	std::fill(_thermal_velocity.begin(), _thermal_velocity.end(), thermal_velocity);
}

void MaxwellianRepresentation::SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const std::vector<double> & thermal_velocity)
{
	_grid_size = density.size();
	_density = density;
	_velocity = velocity;
	_thermal_velocity = thermal_velocity;
}

void MaxwellianRepresentation::GetDensityVelocity(std::vector<double> & density, std::vector<double> & velocity) const
{
	density = _density;
	velocity = _velocity;
}

void MaxwellianRepresentation::GetDensityVelocityPressure(std::vector<double> & density, std::vector<double> & velocity, std::vector<double> & pressure) const
{
	density = _density;
	velocity = _velocity;
	pressure.resize(_grid_size);
	for (int n=0; n<_grid_size; n++)
		pressure.at(n) = _density.at(n) * std::pow(_thermal_velocity.at(n), 2.);
}

void MaxwellianRepresentation::Reset()
{
	std::fill(_density.begin(), _density.end(), 0.);
	std::fill(_velocity.begin(), _velocity.end(), 0.);
	std::fill(_thermal_velocity.begin(), _thermal_velocity.end(), 0.);
}

void MaxwellianRepresentation::print(std::ostream& os) const

{
	os << "Density:" << std::endl;
	for (const double & density : _density)
		os << density << "\t";
	os << std::endl << "Velocity:" << std::endl;
	for (const double & velocity : _velocity)
		os << velocity << "\t";
	os << std::endl << "Pressure:" << std::endl;
	for (int n=0; n<_grid_size; n++)
		os << _density.at(n) * std::pow(_thermal_velocity.at(n), 2.) << "\t";
	os << std::endl;
}

