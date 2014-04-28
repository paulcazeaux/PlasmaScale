/* HEADER MaxwellianRepresentationP1 ===============================================================
 * author: Paul Cazeaux
 * date: 2014-04-25
 *
 * description:
 *  -> store the wavelet representation for the velocity distribution of a set of particles
 *  -> implement some useful operations: DWT and its inverse, denoising, particle weighting
 *  -> Implement some linear operations: addition, multiplication by scalar.
 *  
 * ============================================================================= */

#ifndef DEF_MAXWELLIANREPRESENTATIONP1
#define DEF_MAXWELLIANREPRESENTATIONP1

#include <iostream>
#include <string>
#include <vector>
#include "plasma/Plasma.h"

class MaxwellianRepresentationP1 : public Representation
{
	private:
		/* class members ======================================================================== */
		/* pointer on object */
		std::shared_ptr<std::vector<MaxwellianRepresentation> >	_representations;
		/* Relevant position in the array */
		int 													_bin;
		
		int 													_size;

	public:
		/* constructor  ========================================================================= */
		MaxwellianRepresentationP1(std::shared_ptr<std::vector<MaxwellianRepresentation> > representations, int bin);

		/* static method for practical construction */
		static void MakeMaxwellianRepresentationP1Vector(std::shared_ptr<std::vector<MaxwellianRepresentation> > representations, std::vector<MaxwellianRepresentationP1> & result)
		{
			result.resize(0);
			for (int bin=0; bin<representations.size(); bin++)
				result.emplace_back(std::make_shared<MaxwellianRepresentationP1>(representations, bin));
		}

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
		MaxwellianRepresentationP1& operator+=(const MaxwellianRepresentationP1& rhs)
		{
			assert(_bin == rhs._bin);
			_representations->at(_bin) += rhs._representations->at(_bin);
			return *this;
		}

		MaxwellianRepresentationP1& operator*=(const double lambda)
		{
			_representations->at(_bin) *= lambda;
			return *this;
		}

};

#endif