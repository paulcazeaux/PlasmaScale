#include "representation/MaxwellianRepresentation.h"


MaxwellianRepresentation(std::shared_ptr<const Plasma plasma) :
_plasma(plasma)
{
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

void MaxwellianRepresentation::Weigh(double weight, double velocity)
{
	_density += weight;
	_velocity += weight * velocity;
	_velocitysq += weight * std::pow(velocity, 2.0);
}

void MaxwellianRepresentation::NormalizeAndFinalize(double population_density, double dt)
{
	_velocity /= _density * dt;
	_velocitysq /= dt*dt;
	_density *= population_density;
	_pressure = population_density * (_velocitysq - _velocity*_velocity);
}

void MaxwellianRepresentation::ComputeThermalVelocity()
{
	_thermal_velocity = std::sqrt(_pressure / _density);
}

void MaxwellianRepresentation::LoadBin(int bin, int bin_size,
				std::vector<double>::iterator 	position,
				std::vector<double>:iterator  	velocity, 
				std::vector<double>::iterator 	weight)
{
	assert(bin_size > 1);

	static std::vector<double> vel = std::vector<double>(_number_of_bins);
	static std::vector<double> icdf = std::vector<double>(_number_of_bins+1);
	/* sanity check */

	vel.resize(_number_of_bins);		
	for (int i=0; i<=number_of_bins; i++)
	{
		vel.at(i) = static_cast<double>(i)*_dv - _vmax;
	}

	icdf.resize(_number_of_bins+1);
	icdf.front() = 0.;
	for (int i=1; i<=number_of_bins; i++)
	{
		icdf.at(i) = icdf.at(i-1) + _histogram.at(i-1);
	}
	double density = icdf.back();
	double norm = 1./density;
	for (int i=1; i<=_number_of_bins; i++)
	{
		icdf.at(i) *= norm;
	}

	auto it_vel = vel.begin();
	auto it_icdf = icdf.begin();
	double dn = 1./static_cast<double>(bin_size);
	double dx = _plasma->get_dx();

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
		//double xs = RandomTools::Generate_randomly_uniform(-0.5, 0.5);
		double cellpos = xs + 0.5/static_cast<double>(bin_size) - 0.5;

		*(position+i) = (static_cast<double>(bin)+cellpos) * dx;
		*(velocity+i) = _velocity + _thermal_velocity * (*it_vel + _dv*(fv - *it_icdf)/(*(it_icdf+1) - *it_icdf));
		*(weight+i)   = density;
	}
	bin_start_index = bin_end_index;
}


void MaxwellianRepresentation::Reset()
{
	_density = 0.;
	_velocity = 0.;
	_velocitysq = 0.;
	_thermal_velocity = 0.;
}

void MaxwellianRepresentation::print(std::ofstream& os)
{
}
