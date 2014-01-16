/* HEADER RandomTools ====================================================================================
 * author: Paul Cazeaux
 * date: 2013-11-28
 * 
 * description:
 *  -> FFT routines
 * 
 * ============================================================================================== */

#ifndef DEF_PLASMASCALE_RANDOMTOOLS
#define DEF_PLASMASCALE_RANDOMTOOLS

/* includes ===================================================================================== */
#include <iostream>
#include <random>

class RandomTools
{
	public:
		/* Uniform random generator */
		template<typename T>
		static T Generate_randomly_uniform(const T start, const T end);

		/* Normal distributed random generator */
		static double Generate_randomly_normal(const double mean, const double stddev);

};


 #endif