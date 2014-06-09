/* HEADER Tools ====================================================================================
 * author: Paul Cazeaux
 * date: 2013-11-28
 * 
 * description:
 *  -> P1 interpolation
 *  -> Slope fitting by mean-squares
 *  -> ICDF computation with possibly negative histogram values
 *  -> Computation of an empirical scaling for Poisson distribution fitting
 * 
 * ============================================================================================== */

#ifndef DEF_PLASMASCALE_TOOLS
#define DEF_PLASMASCALE_TOOLS

/* includes ===================================================================================== */
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <cassert>
#include <utility>

class Tools
{

	public:
		static void DisplayTitle();

		template<typename T>
		static T EvaluateP1Function(const std::vector<T> & values, const int bin, const double cellpos)
		{
			assert((cellpos>=0) && (cellpos<=1));
			if (bin+1 < values.size())
				return (1-cellpos)*values.at(bin) + cellpos*values.at(bin+1);
			else
				return (1-cellpos)*values.at(bin) + cellpos*values.front();
		}

		template<typename T>
		static T EvaluateSlope(const std::vector<T> & values )
		{
			int size = values.size();
			assert(size>0);
			T sum = values.at(0), first_moment = 0 * sum;
			for (int n=1; n<size; n++)
			{
				sum				+= values.at(n);
				first_moment 	+= static_cast<double>(n)*values.at(n);
			}
			return 6./static_cast<double>(size*(size+1)) * (2./static_cast<double>(size-1)*first_moment - sum);
		}

		template<typename T>
		static T EvaluateSlope(const std::vector<T> & values, const int size)
		{
			assert(size>0);
			T sum = values.at(0), first_moment = 0 * sum;
			for (int n=1; n<size; n++)
			{
				sum				+= values.at(n);
				first_moment 	+= static_cast<double>(n)*values.at(n);
			}
			return 6./static_cast<double>(size*(size+1)) * (2./static_cast<double>(size-1)*first_moment - sum);
		}

		template<typename T>
		static T EvaluateSlope(const std::vector<T> & points, const std::vector<T> & values)
		{
			assert(points.size()==values.size());
			int size = values.size();
			assert(size>0);
			T covar = points.at(0)*values.at(0), var = points.at(0)*points.at(0);
			for (int n=1; n<size; n++)
			{
				covar	+= points.at(n)*values.at(n);
				var		+= points.at(n)*points.at(n);
			}
			return covar/var;
		}

		static void AssembleICDF(const std::vector<double>& histogram, std::vector<double>& icdf, double& density);
		static double ComputePoissonEmpiricalScaling(const std::vector<std::vector<double> >& histogram, const int& patch_size = 8);
};


 #endif