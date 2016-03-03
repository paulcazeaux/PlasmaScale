/* HEADER PUREHaarRepresentationP1 ===============================================================
 * author: Paul Cazeaux
 * date: 2014-04-25
 *
 * description:
 *  -> store the wavelet representation for the velocity distribution of a set of particles
 *  -> implement some useful operations: DWT and its inverse, denoising, particle weighting
 *  -> Implement some linear operations: addition, multiplication by scalar.
 *  
 * ============================================================================= */

#ifndef DEF_PUREHaarREPRESENTATIONP1
#define DEF_PUREHaarREPRESENTATIONP1

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "tools/Tools.h"
#include "plasma/Plasma.h"
#include "representation/PUREHaarRepresentation.h"

class PUREHaarRepresentationP1 : public PUREHaarRepresentation
{
	public:
		/* constructor  ========================================================================= */
		PUREHaarRepresentationP1() {}
		PUREHaarRepresentationP1(std::shared_ptr<const Plasma> plasma, double vmax, int max_depth, int grid_size, int min_depth = 1, int buffer = 1);

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

inline PUREHaarRepresentationP1 operator+(PUREHaarRepresentationP1 lhs, const PUREHaarRepresentationP1& rhs)
{
	lhs += rhs;
	return lhs;
}

inline PUREHaarRepresentationP1 operator-(PUREHaarRepresentationP1 lhs, const PUREHaarRepresentationP1& rhs)
{
	lhs -= rhs;
	return lhs;
}

inline PUREHaarRepresentationP1 operator*(PUREHaarRepresentationP1 lhs, const double lambda)
{
	lhs *= lambda;
	return lhs;
}

inline PUREHaarRepresentationP1 operator*(const double lambda, PUREHaarRepresentationP1 rhs)
{
	rhs *= lambda;
	return rhs;
}

#endif