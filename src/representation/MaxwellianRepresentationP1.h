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
#include "plasma/Plasma.h"
#include "representation/MaxwellianRepresentation.h"

class MaxwellianRepresentationP1 : public MaxwellianRepresentation
{
	public:
		/* constructor  ========================================================================= */
		MaxwellianRepresentationP1(std::shared_ptr<const Plasma> plasma, int grid_size);

		/* methods */
		virtual void Weigh(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weight);

		virtual void Load(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity, 
								std::vector<double>::iterator 	weight);

};

#endif