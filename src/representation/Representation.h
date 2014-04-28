/* HEADER Representation ===============================================================
 * author: Paul Cazeaux
 * date: 2014-04-25
 *
 * description:
 *  -> Template for the various representations possible for each population
 *  -> Provides the necessary virtual functions for the parameterization algorithm to work
 *  
 * ============================================================================= */

#ifndef DEF_REPRESENTATION
#define DEF_REPRESENTATION

#include <iostream>
#include <string>
#include <vector>
#include "plasma/Plasma.h"

class Representation
{
	public:
		/* friends */
		friend std::ostream & operator<<( std::ostream& os, const WaveletRepresentation& representation)
		{
			representation.print(os);
			return os;
		}

		/* methods */

		virtual void Weigh(int size,
				std::vector<double>::iterator 	position,
				std::vector<double>:iterator  	velocity,
				std::vector<double>::iterator 	weight);
		virtual void Reset();

		virtual void Load(int size,
				std::vector<double>::iterator 	position,
				std::vector<double>:iterator  	velocity,
				std::vector<double>::iterator 	weight);
		virtual void Coarsen();
		virtual void Refine();

		virtual void print(std::ostream& os) const;
		virtual int get_grid_size() const;

		/* operator overload */
		virtual WaveletRepresentation& operator+=(const WaveletRepresentation& rhs);

		virtual WaveletRepresentation& operator*=(const double lambda);
};

inline Representation operator+(Representation lhs, const Representation& rhs)
{
	lhs += rhs;
	return lhs;
}

inline Representation operator*(Representation lhs, const double lambda)
{
	lhs *= lambda;
	return lhs;
}

inline Representation operator*(const double lambda, Representation rhs)
{
	rhs *= lambda;
	return rhs;
}



#endif