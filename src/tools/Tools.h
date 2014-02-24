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

		static double EvaluateP1Function(const std::vector<double> & values, const int number, const double cellpos);
		static double EvaluateSlope(const std::vector<double> & values );
};


 #endif