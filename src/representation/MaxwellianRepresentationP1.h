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
#include <algorithm>
#include "tools/Tools.h"
#include "plasma/Plasma.h"
#include "representation/MaxwellianRepresentation.h"

class MaxwellianRepresentationP1 : public MaxwellianRepresentation
{
	public:
		/* constructor  ========================================================================= */
		MaxwellianRepresentationP1() {}
		MaxwellianRepresentationP1(std::shared_ptr<const Plasma> plasma, int grid_size);

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

};


inline MaxwellianRepresentationP1 operator+(MaxwellianRepresentationP1 lhs, const MaxwellianRepresentationP1& rhs)
{
	lhs += rhs;
	return lhs;
}

inline MaxwellianRepresentationP1 operator-(MaxwellianRepresentationP1 lhs, const MaxwellianRepresentationP1& rhs)
{
	lhs -= rhs;
	return lhs;
}

inline MaxwellianRepresentationP1 operator*(MaxwellianRepresentationP1 lhs, const double lambda)
{
	lhs *= lambda;
	return lhs;
}

inline MaxwellianRepresentationP1 operator*(const double lambda, MaxwellianRepresentationP1 rhs)
{
	rhs *= lambda;
	return rhs;
}

#endif