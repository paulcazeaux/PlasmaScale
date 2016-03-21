/* HEADER WaveletRepresentationP1 ===============================================================
 * author: Paul Cazeaux
 * date: 2014-04-25
 *
 * description:
 *  -> store the wavelet representation for the velocity distribution of a set of particles
 *  -> implement some useful operations: DWT and its inverse, denoising, particle weighting
 *  -> Implement some linear operations: addition, multiplication by scalar.
 *  
 * ============================================================================= */

#ifndef DEF_WAVELETREPRESENTATIONP1
#define DEF_WAVELETREPRESENTATIONP1

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "tools/Tools.h"
#include "plasma/Plasma.h"
#include "representation/WaveletRepresentation.h"

class WaveletRepresentationP1 : public WaveletRepresentation
{
	public:
		/* constructor  ========================================================================= */
		WaveletRepresentationP1() {};
		WaveletRepresentationP1(std::shared_ptr<const Plasma> plasma, double vmin, double vmax, int depth, double reference_density);

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

inline WaveletRepresentationP1 operator+(WaveletRepresentationP1 lhs, const WaveletRepresentationP1& rhs)
{
	lhs += rhs;
	return lhs;
}

inline WaveletRepresentationP1 operator-(WaveletRepresentationP1 lhs, const WaveletRepresentationP1& rhs)
{
	lhs -= rhs;
	return lhs;
}

inline WaveletRepresentationP1 operator*(WaveletRepresentationP1 lhs, const double lambda)
{
	lhs *= lambda;
	return lhs;
}

inline WaveletRepresentationP1 operator*(const double lambda, WaveletRepresentationP1 rhs)
{
	rhs *= lambda;
	return rhs;
}

#endif