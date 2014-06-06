/* HEADER Tools ====================================================================================
 * author: Paul Cazeaux
 * date: 2013-11-28
 * 
 * description:
 *  -> FFT routines
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
		/* constructor and destructor =========================================================== */
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

		static void AssembleICDF(const std::vector<double>& histogram, std::vector<double>& icdf, double& density)
		{
			int nbins = histogram.size();
			icdf.resize(nbins+1);
			icdf.front() = 0.;
			std::partial_sum(histogram.begin(), histogram.end(), icdf.begin()+1);
			density = icdf.back();

			/* Ensure that the ICDF is a monotonic function */
			double oldval=0.;
			for (int i=1; i<=nbins; i++)
			{
				if (icdf.at(i) < oldval)
				{
					int ind1 = i-1, ind2 = i;
					double high = oldval, low = icdf.at(i);
					if (high>density) high=density;

					/* The vector icdf is monotonic up to index i-1 */
                    if (ind2<nbins)
                        while (icdf.at(++ind2) < high)
                        {
                            if (icdf.at(ind2) < low)
                                low = icdf.at(ind2);
                        }
					if (low<0.) low=0.;
                    if (ind1>0)
                        while (icdf.at(--ind1) > low) {}
					/* double val1 = icdf.at(ind1), dn = (icdf.at(ind2)-val1)/static_cast<double>(ind2-ind1);
					for (int j=1; j<=ind2-ind1; j++)
					{
						icdf.at(ind1+j) = val1 + static_cast<double>(j)*dn;
					}*/
                    for (int j=ind1+1; j<i; j++)
                        icdf.at(j) = low;
				}
				oldval = icdf.at(i);
			}

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

		static void PUREShrink(std::vector<double>& scaling, std::vector<double>& detail, std::vector<bool>& pattern, const int& maxdepth, const int& mindepth = 0.);

};


 #endif