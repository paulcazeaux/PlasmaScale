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
#include "tools/Tools.h"
#include "tools/RandomTools.h"

class Representation
{
	public:
		/* friends */
		friend std::ostream & operator<<( std::ostream& os, const Representation& representation)
		{
			representation.print(os);
			return os;
		}

		/* methods */
		virtual void Weigh(int size,
				std::vector<double>::iterator 	position,
				std::vector<double>::iterator  	velocity,
				std::vector<double>::iterator 	weight) = 0;
		virtual void Weigh(int size,
				std::vector<double>::iterator 	position,
				std::vector<double>::iterator  	velocity,
				std::vector<double>::iterator 	weight,
				const double delay,
				const std::vector<double> & accfield) = 0;
		virtual void Reset() = 0;

		virtual void Load(int size,
				std::vector<double>::iterator 	position,
				std::vector<double>::iterator  	velocity,
				std::vector<double>::iterator 	weight) = 0;
		virtual void Coarsen() = 0;
		virtual void Refine() = 0;

		virtual void SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const double & thermal_velocity) = 0;
		virtual void SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const std::vector<double> & thermal_velocity) = 0;
		virtual void GetDensityVelocity(std::vector<double> & density, std::vector<double> & velocity) = 0;
		virtual void GetDensityVelocityPressure(std::vector<double> & density, std::vector<double> & velocity, std::vector<double> & pressure) = 0;

		virtual void print(std::ostream& os) const = 0;
		virtual int get_grid_end() const = 0;
		virtual void set_grid_end(const int new_grid_end) = 0;
};
#endif