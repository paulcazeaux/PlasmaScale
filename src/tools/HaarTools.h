/* HEADER HaarTools ====================================================================================
 * author: Paul Cazeaux
 * date: 2014-06-07
 * 
 * description:
 *  -> Haar Transforms
 *  -> Denoising with PUREShrink
 * 
 * ============================================================================================== */

#ifndef DEF_PLASMASCALE_HAARTOOLS
#define DEF_PLASMASCALE_HAARTOOLS

/* includes ===================================================================================== */
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <cassert>
#include <utility>

class HaarTools
{

	public:
		static void ForwardTransform(const std::vector<double>& hist, std::vector<double>& scaling, std::vector<double>& detail, const int& maxdepth, const int& mindepth = 0.);
		static void InverseTransform(std::vector<double>& hist, std::vector<double>& scaling, const std::vector<double>& detail, const int& maxdepth, const int& mindepth = 0.);
		static void PUREShrink(std::vector<double>& scaling, std::vector<double>& detail, std::vector<bool>& pattern, const int& maxdepth, const int& mindepth = 0.);
};


 #endif