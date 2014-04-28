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
#include <vector>
#include <cassert>

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


		template<typename T>
		static T EvaluateSlope(const std::vector<T> & values )
		{
			int size = values.size();
			T sum = 0, first_moment = 0;
			for (int n=0; n<size; n++)
			{
				sum				+= values.at(n);
				first_moment 	+= static_cast<double>(n)*values.at(n);
			}
			return 6./static_cast<double>(size*(size+1)) * (2./static_cast<double>(size-1)*first_moment - sum);
		}
};


 #endif