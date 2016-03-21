/* HEADER MaxwellianRepresentation ===============================================================
 * author: Paul Cazeaux
 * date: 2014-04-25
 *
 * description:
 *  -> store the wavelet representation for the velocity distribution of a set of particles
 *  -> implement some useful operations: DWT and its inverse, denoising, particle weighting
 *  -> Implement some linear operations: addition, multiplication by scalar.
 *  
 * ============================================================================= */

#ifndef DEF_MAXWELLIANREPRESENTATION
#define DEF_MAXWELLIANREPRESENTATION

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cassert>
#include "plasma/Plasma.h"
#include "representation/Representation.h"

class MaxwellianRepresentation : public Representation
{
	protected:
		/* class members ======================================================================== */
		/* pointer on object */
		std::shared_ptr<const Plasma>		_plasma;
		/* parameter */
		int 								_grid_end;
		double 								_reference_density;
		/* values */
		std::vector<double> 				_density;
		std::vector<double> 				_velocity;
		std::vector<double> 				_thermal_velocity;

		/* Arrays describing a normalized Maxwellian distribution with variance 1 */
		static std::vector<double> 			_quiet_start_vel;
		static std::vector<double> 			_quiet_start_icdf;

	public:
		/* constructor  ========================================================================= */
		MaxwellianRepresentation() {}
		MaxwellianRepresentation(std::shared_ptr<const Plasma> plasma, double reference_density);

		/* static method */
		static void InitializeQuietStartArrays(int number_of_bins);

		/* methods */
		virtual void Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weight);
		virtual void Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weight,
								const double delay,
								const std::vector<double> & accfield);

		virtual void Load(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity, 
								std::vector<double>::iterator 	weight);
		virtual void Coarsen();
		virtual void Refine();

		virtual void SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const double & thermal_velocity);
		virtual void SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const std::vector<double> & thermal_velocity);
		virtual void GetDensityVelocity(std::vector<double> & density, std::vector<double> & velocity);
		virtual void GetDensityVelocityPressure(std::vector<double> & density, std::vector<double> & velocity, std::vector<double> & pressure);

		virtual void Reset();
		virtual void print(std::ostream& os) const;
		virtual int  get_grid_end() const 	{ return _grid_end;	}
		virtual void set_grid_end(const int new_grid_end)
		{
			_density.resize(new_grid_end+1);
			_velocity.resize(new_grid_end+1);
			_thermal_velocity.resize(new_grid_end+1);
			_grid_end = new_grid_end;
		}

		/* operator overload */
		MaxwellianRepresentation& operator=(const MaxwellianRepresentation& rhs)
		{
			_plasma = rhs._plasma;
			_grid_end = rhs._grid_end;
			_reference_density = rhs._reference_density;

			_density = std::vector<double>(rhs._density.begin(), rhs._density.begin()+_grid_end+1);
			_density = std::vector<double>(rhs._velocity.begin(), rhs._velocity.begin()+_grid_end+1);
			_density = std::vector<double>(rhs._thermal_velocity.begin(), rhs._thermal_velocity.begin()+_grid_end+1);

			return *this;
		}

		MaxwellianRepresentation& operator+=(const MaxwellianRepresentation& rhs)
		{
			assert((_grid_end==rhs._grid_end)&&(_reference_density=rhs._reference_density));
			for (int n=0; n<=_grid_end; n++)
			{
				_density.at(n) += rhs._density.at(n);
				_velocity.at(n) += rhs._velocity.at(n);

				/* We add the pressures, not the thermal velocities */
				double p;
				p = _density.at(n) * std::pow(_thermal_velocity.at(n), 2.) 
					+ rhs._density.at(n) * std::pow(rhs._thermal_velocity.at(n), 2.);
				_thermal_velocity.at(n) = std::sqrt(p / _density.at(n));
			}
			return *this;
		}

		MaxwellianRepresentation& operator-=(const MaxwellianRepresentation& rhs)
		{
			assert((_grid_end==rhs._grid_end)&&(_reference_density=rhs._reference_density));
			for (int n=0; n<=_grid_end; n++)
			{
				_density.at(n) -= rhs._density.at(n);
				_velocity.at(n) -= rhs._velocity.at(n);

				/* We add the pressures, not the thermal velocities */
				double p;
				p = _density.at(n) * std::pow(_thermal_velocity.at(n), 2.) 
					- rhs._density.at(n) * std::pow(rhs._thermal_velocity.at(n), 2.);
				_thermal_velocity.at(n) = std::sqrt(p / _density.at(n));
			}
			return *this;
		}

		MaxwellianRepresentation& operator*=(const double lambda)
		{
			for (int n=0; n<=_grid_end; n++)
			{
				_density.at(n) *= lambda;
				_velocity.at(n) *= lambda;
				_thermal_velocity.at(n) *= lambda;
			}

			return *this;
		}

};


inline MaxwellianRepresentation operator+(MaxwellianRepresentation lhs, const MaxwellianRepresentation& rhs)
{
	lhs += rhs;
	return lhs;
}

inline MaxwellianRepresentation operator-(MaxwellianRepresentation lhs, const MaxwellianRepresentation& rhs)
{
	lhs -= rhs;
	return lhs;
}

inline MaxwellianRepresentation operator*(MaxwellianRepresentation lhs, const double lambda)
{
	lhs *= lambda;
	return lhs;
}

inline MaxwellianRepresentation operator*(const double lambda, MaxwellianRepresentation rhs)
{
	rhs *= lambda;
	return rhs;
}

#endif