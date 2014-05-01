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
		WaveletRepresentationP1(std::shared_ptr<const Plasma> plasma, double vmax, int depth, int grid_size);

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