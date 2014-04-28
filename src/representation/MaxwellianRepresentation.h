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
#include "plasma/Plasma.h"

class MaxwellianRepresentation : public Representation
{
	private:
		/* class members ======================================================================== */
		/* pointer on object */
		std::shared_ptr<const Plasma>		_plasma;

		/* values */
		double 								_density;
		double 								_velocity;
		double 								_velocitysq;
		double 								_pressure;
		double 								_thermal_velocity;

		/* Arrays describing a normalized Maxwellian distribution with variance 1 */
		static std::vector<double> 			_quiet_start_vel;
		static std::vector<double> 			_quiet_start_icdf;

	public:
		/* constructor  ========================================================================= */
		MaxwellianRepresentation() {};
		MaxwellianRepresentation(std::shared_ptr<const Plasma> plasma, double vmax, int depth);

		/* static method */
		static void InitializeQuietStartArrays(int number_of_bins);

		/* methods */
		virtual void Weigh(double weight, double velocity);
		void FinalizeWeigh();
		void ComputeThermalVelocity();

		virtual void LoadBin(int bin, int bin_size,
								std::vector<double>::iterator 	position,
								std::vector<double>:iterator  	velocity,
								std::vector<double>::iterator 	weight);

		virtual void Reset();
		virtual void print(std::ostream& os) const;
		virtual int get_grid_size() const 	{ return _grid_size;	}

		/* operator overload */
		MaxwellianRepresentation& operator+=(const MaxwellianRepresentation& rhs)
		{
			_density += rhs._density;
			_velocity += rhs._velocity;
			_velocitysq += rhs._velocitysq;
			_pressure += rhs._pressure;
			_thermal_velocity += rhs._thermal_velocity;

			return *this;
		}

		MaxwellianRepresentation& operator*=(const double lambda)
		{
			_density *= lambda;
			_velocity *= lambda;
			_velocitysq  *= lambda;
			_pressure *= lambda;
			_thermal_velocity *= lambda;

			return *this;
		}

};

#endif