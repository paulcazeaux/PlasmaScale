#include "representation/MaxwellianRepresentationP1.h"


MaxwellianRepresentationP1(std::shared_ptr<std::vector<WaveletRepresentation>> representations, int bin) :
_representations(representations), _bin(bin)
{
	_size = _representations.size();
}

void MaxwellianRepresentationP1::InitializeQuietStartArrays(int number_of_bins)
{
	MaxwellianRepresentation::InitializeQuietStartArrays(number_of_bins);
}

void MaxwellianRepresentationP1::Weigh(double weight, double velocity)
{
	_representations.at(_bin)->Weigh(weight, velocity);
}

void MaxwellianRepresentation::Weigh(double weight, double velocity, double cell_pos)
{
	double right_weight =  cell_pos * weight;
	double left_weight = weight - right_weight;

	_representations->at(bin).Weigh(left_weight, velocity);
	if (bin < _size-1)
	{
		_representations->at(bin+1).Weigh(right_weight, velocity);
	}
	else
	{
		_representations->at(0).Weigh(right_weight, velocity);
	}
}

void MaxwellianRepresentationP1::NormalizeAndFinalize(double population_density, double dt)
{
	_representations.at(_bin)->NormalizeAndFinalize(density, dt);
}

void MaxwellianRepresentationP1::ComputeThermalVelocity()
{
	_representations.at(_bin)->ComputeThermalVelocity();
}

void MaxwellianRepresentationP1::LoadBin(int bin, int bin_size,
				std::vector<double>::iterator 	position,
				std::vector<double>:iterator  	velocity, 
				std::vector<double>::iterator 	weight)
{
	assert(bin = _bin);
	assert(bin_size > 1);

	static std::vector<double> vel = std::vector<double>(_number_of_bins);
	static std::vector<double> icdf = std::vector<double>(_number_of_bins+1);
	/* sanity check */

	vel.resize(_number_of_bins);
	vel.front() = -vmax;				
	for (int i=1; i<=number_of_bins; i++)
	{
		vel.at(i) = static_cast<double>(i)*_dv - _vmax;
	}

	icdf.resize(_number_of_bins+1);
	for (int i=0; i<=number_of_bins; i++)
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
		double cellpos = xs + 0.5/static_cast<double>(bin_size);

		*(position+i) = (static_cast<double>(bin)+cellpos) * dx;
		*(velocity+i) = _velocity + _thermal_velocity * (*it_vel + _dv*(fv - *it_icdf)/(*(it_icdf+1) - *it_icdf));
		*(weight+i)   = density;
	}
	bin_start_index = bin_end_index;
}


void MaxwellianRepresentationP1::Reset()
{
	_representations.at(_bin)->Reset();
}

void MaxwellianRepresentationP1::print(std::ostream& os)
{
	for (auto & representation : *_representations)
		os = representation.print(os);
}