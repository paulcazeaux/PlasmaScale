
#include "Tools.h"

void Tools::DisplayTitle()
{
	std::cout << std::endl << std::endl << "PlasmaScale - Based on ES1 - Electrostatic 1 Dimensional Code" << std::endl << std::endl;
	//printf("version 4.1\n");
	//printf("(c) Copyright 1987-92 Regents of the University of California\n");
	//printf("Plasma Theory and Simulation Group\n");
	//printf("University of California - Berkeley\n");
}//

void Tools::AssembleICDF(std::vector<double>& histogram, std::vector<double>& icdf, double& density)
{
	// std::vector<double> tmp = icdf; // Debug mode
	// try{
	int nbins = histogram.size();
	static std::vector<double> icdf_reverse;
	icdf.resize(nbins+1);
	icdf.front() = 0.;
	std::partial_sum(histogram.begin(), histogram.end(), icdf.begin()+1);
	density = icdf.back();

	if (density <= 0)
	{
		/* Give up reconstruction */
		std::fill(icdf.begin(), icdf.end(), 0.);
		density = 0;
		return;
	}

	/* First, we take care of negative values and overshots */
	int i0 = 0, i1 = nbins, imax = 0;
	for (int i=1; i<nbins; i++)
	{
		if (icdf.at(i) < 0) i0 = i;
		if (icdf.at(i) > density && i1 > i) i1 = i;
		if (icdf.at(i) > .5*density && icdf.at(i-1) < .5*density) imax = i;
	}

	if (i1 < i0) /* Accidents happen */
	{
		double df = (density - 0.)/(i0 - i1);
		for (int j=0; j<= i1; j++) icdf.at(j) = 0.;
		for (int j=i1+1; j<i0; j++) icdf.at(j) = (j-i1)*df;
		for (int j=i0; j<=nbins; j++) icdf.at(j) = density;
		i0 = 0; i1 = nbins;
	}

	if (i0 > 0)
	{
		int ind1 = i0 - 1;
		int ind2 = i0 + 1;
		double S = 0;
		for (int j = 1; j < ind2; j++) 
			S += icdf.at(j);
		double df = 2.*S/static_cast<double>((ind2-ind1)*(ind2-ind1-1));
		while (df<0 || (ind2-1-ind1)*df > icdf.at(ind2))
		{
			if (ind1 > 1  && icdf.at(ind2) <= density)
				ind1--;
			else if (ind2 < imax)
			{
				S += icdf.at(ind2);
				ind2++;
				ind1 = i0-1;
			}
			else // There is something wrong with the distribution. Give up momentum conservation
       		{
				ind1 = i0-1;
       			ind2 = i0+1;
	   			df = std::min(density, icdf.at(ind2))/(ind2-ind1);
       			break;
       		}

			df = 2.*S/static_cast<double>((ind2-ind1)*(ind2-ind1-1));
		}
		for (int j=0; j<ind1+1; j++)
			icdf.at(j) = 0;
		for (int j = ind1+1; j<ind2; j++)
			icdf.at(j) = (j-ind1)*df;
	}

	if (i1 < nbins)
	{
		int ind1 = i1 - 1;
		int ind2 = i1 + 1;
		double S = 0;
		for (int j = ind1+1; j < nbins; j++) 
			S += icdf.at(j);
		double df = 2.*((nbins-ind1-1)*density-S)/static_cast<double>((ind2-ind1)*(ind2-ind1-1));
		while (df<0 || (ind2-1-ind1)*df > density - icdf.at(ind1))
		{
			if (ind2 < nbins+1 && icdf.at(ind1) >= 0)
				ind2++;
			else if (ind1 > imax)
			{
				S += icdf.at(ind1);
				ind1--;
				ind2 = i1+1;
			}
			else
			{
				ind1 = i1-1;
				ind2 = i1+1;
				df =  (std::max(icdf.at(ind1), 0.) - density)/(ind1-ind2);
				break;
			}
			df = 2.*((nbins-ind1-1)*density-S)/static_cast<double>((ind2-ind1)*(ind2-ind1-1));
		}
		for (int j = ind1+1; j<ind2; j++)
			icdf.at(j) = density + (j-ind2)*df;
		for (int j=ind2; j<nbins; j++)
			icdf.at(j) = density;
	}

	double lowval=0.;
	int i = 1;
	/* Ensure that the ICDF is a monotonic function */
	while (i < imax)
	{
		if (lowval > icdf.at(i))
		{
			/* The vector icdf is monotonic up to index i-1 */
			int ind1 = i-1, ind2 = i+1;
			double minval = icdf.at(i);

            while (lowval > icdf.at(ind2))
            {
            	ind2++;
                if (minval > icdf.at(ind2)) 
               	{
               		minval = icdf.at(ind2);
               		i = ind2;
               	}
            }
            while (minval < icdf.at(ind1)) ind1--;
            double S = 0;
            for (int j = ind1+1; j<i; j++)
            	S += icdf.at(j);
            double df = 2.*(S - (i - ind1 - 1)*icdf.at(ind1))/static_cast<double>((i-ind1)*(i-ind1-1));
            while (df<0 || df*(i-1-ind1) > icdf.at(i) - icdf.at(ind1))
            {
           		if (i<imax)
           		{
	            	S += icdf.at(i);
            		i++;
           			df = 2.*(S - (i - ind1 - 1)*icdf.at(ind1))/static_cast<double>((i-ind1)*(i-ind1-1));
           		}
           		else // There is something wrong with the distribution. Give up momentum conservation
           		{
           			i = ind2;
           			df = (icdf.at(i) - icdf.at(ind1))/(i-ind1);
           			break;
           		}

            }
            for (int j=ind1+1; j<i; j++)
               icdf.at(j) = icdf.at(ind1) + df*(j-ind1);
		}
		lowval = icdf.at(i);
		i++;
	}
	double highval = density;
	i=nbins-1;
	while (i>imax)
	{
		if (highval < icdf.at(i))
		{        
			int ind1 = i+1, ind2 = i-1;
			double maxval = icdf.at(i); 
			while (highval < icdf.at(ind2))
			{
				ind2--;
				if (maxval < icdf.at(ind2))
				{
					maxval = icdf.at(ind2);
					i = ind2;
				}
			}
			while (maxval > icdf.at(ind1)) ind1++;  

			double S = 0;
			for (int j=i+1; j<ind1-1; j++)
				S += icdf.at(j);
			double df = 2/((ind1-i)*(ind1-i-1))*((ind1 - i -1) * icdf.at(ind1) - S);
			while (df < 0 || df*(ind1-i-1) > icdf.at(ind1) - icdf.at(i))
			{
				if (i>imax)
				{
					S += icdf.at(i);
					i--;
					df = 2/((ind1-i)*(ind1-i-1))*((ind1 - i -1) * icdf.at(ind1) - S);
				}
				else // There is something wrong with the distribution. Give up momentum conservation
				{
           			i = ind2;
           			df = (icdf.at(ind1) - icdf.at(i))/(ind1-i);
           			break;
				}
			}
            for (int j=ind1-1; j>i; j--)
               icdf.at(j) = icdf.at(ind1) + df*(j-ind1);
		}

		highval = icdf.at(i);
		i--;
	}
	for (int i=0; i<nbins; i++)
		histogram.at(i) = icdf.at(i+1) - icdf.at(i);
// }
// catch(...){ 
// std::cout << "Failure in the ICDF reconstruction\n"; 
// icdf= tmp;


// int nbins = histogram.size();
// static std::vector<double> icdf_reverse;
// icdf.resize(nbins+1);
// icdf.front() = 0.;
// std::partial_sum(histogram.begin(), histogram.end(), icdf.begin()+1);
// density = icdf.back();
// std::cout << "1:\n";
// for (double & s: icdf) std::cout << s << "\t"; std::cout << "\n\n"; 

// if (density <= 0)
// {
// 	/* Give up reconstruction */
// 	std::fill(icdf.begin(), icdf.end(), 0.);
// 	density = 0;
// 	return;
// }

// /* First, we take care of negative values and overshots */
// int i0 = 0, i1 = nbins, imax = 0;
// for (int i=1; i<nbins; i++)
// {
// 	if (icdf.at(i) < 0) i0 = i;
// 	if (icdf.at(i) > density && i1 > i) i1 = i;
// 	if (icdf.at(i) > .5*density && icdf.at(i-1) < .5*density) imax = i;
// }

// if (i1 < i0) /* Accidents happen */
// {
// 	double df = (density - 0.)/(i0 - i1);
// 	for (int j=0; j<= i1; j++) icdf.at(j) = 0.;
// 	for (int j=i1+1; j<i0; j++) icdf.at(j) = (j-i1)*df;
// 	for (int j=i0; j<=nbins; j++) icdf.at(j) = density;
// 	i0 = 0; i1 = nbins;
// }

// if (i0 > 0)
// {
// 	int ind1 = i0 - 1;
// 	int ind2 = i0 + 1;
// 	double S = 0;
// 	for (int j = 1; j < ind2; j++) 
// 		S += icdf.at(j);
// 	double df = 2.*S/static_cast<double>((ind2-ind1)*(ind2-ind1-1));
// 	while (df<0 || (ind2-1-ind1)*df > icdf.at(ind2))
// 	{
// 		if (ind1 > 1  && icdf.at(ind2) <= density)
// 			ind1--;
// 		else if (ind2 < imax)
// 		{
// 			S += icdf.at(ind2);
// 			ind2++;
// 			ind1 = i0-1;
// 		}
// 		else // There is something wrong with the distribution. Give up momentum conservation
//    		{
// 			ind1 = i0-1;
//    			ind2 = i0+1;
//    			df = std::min(density, icdf.at(ind2))/(ind2-ind1);
//    			break;
//    		}

// 		df = 2.*S/static_cast<double>((ind2-ind1)*(ind2-ind1-1));
// 	}
// 	for (int j=0; j<=ind1; j++)
// 		icdf.at(j) = 0;
// 	for (int j = ind1+1; j<ind2; j++)
// 		icdf.at(j) = (j-ind1)*df;
// 	std::cout << "2:\t ind1:\t" << ind1 << "\ti0:\t" << i0 << "\tind2:\t" << ind2 << "\n";
// }
// for (double & s: icdf) std::cout << s << "\t"; std::cout << "\n\n"; 

// if (i1 < nbins)
// {
// 	int ind1 = i1 - 1;
// 	int ind2 = i1 + 1;
// 	double S = 0;
// 	for (int j = ind1+1; j < nbins; j++) 
// 		S += icdf.at(j);
// 	double df = 2.*((nbins-ind1-1)*density-S)/static_cast<double>((ind2-ind1)*(ind2-ind1-1));
// 	while (df<0 || (ind2-1-ind1)*df > density - icdf.at(ind1))
// 	{
// 		if (ind2 < nbins+1 && icdf.at(ind1) >= 0)
// 			ind2++;
// 		else if (ind1 > imax)
// 		{
// 			S += icdf.at(ind1);
// 			ind1--;
// 			ind2 = i1+1;
// 		}
// 		else
// 		{
// 			ind1 = i1-1;
// 			ind2 = i1+1;
// 			df =  (std::max(icdf.at(ind1), 0.) - density)/(ind1-ind2);
// 			break;
// 		}
// 		df = 2.*((nbins-ind1-1)*density-S)/static_cast<double>((ind2-ind1)*(ind2-ind1-1));
// 	}
// 	for (int j = ind1+1; j<ind2; j++)
// 		icdf.at(j) = density + (j-ind2)*df;
// 	for (int j=ind2; j<nbins; j++)
// 		icdf.at(j) = density;
// std::cout << "3:\t ind1:\t" << ind1 << "\ti1:\t" << i1 << "\tind2:\t" << ind2 << "\n";
// }
// for (double & s: icdf) std::cout << s << "\t"; std::cout << "\n\n"; 

// double lowval=0.;
// int i = 1;
// /* Ensure that the ICDF is a monotonic function */
// while (i < imax)
// {
// 	if (lowval > icdf.at(i))
// 	{
// 		/* The vector icdf is monotonic up to index i-1 */
// 		int ind1 = i-1, ind2 = i+1;
// 		double minval = icdf.at(i);

//         while (lowval > icdf.at(ind2))
//         {
//         	ind2++;
//             if (minval > icdf.at(ind2)) 
//            	{
//            		minval = icdf.at(ind2);
//            		i = ind2;
//            	}
//         }
//         while (minval < icdf.at(ind1)) ind1--;
//         double S = 0;
//         for (int j = ind1+1; j<i; j++)
//         	S += icdf.at(j);
//         double df = 2.*(S - (i - ind1 - 1)*icdf.at(ind1))/static_cast<double>((i-ind1)*(i-ind1-1));
//         while (df<0 || df*(i-1-ind1) > icdf.at(i) - icdf.at(ind1))
//         {
//        		if (i<imax)
//        		{
//             	S += icdf.at(i);
//         		i++;
//        			df = 2.*(S - (i - ind1 - 1)*icdf.at(ind1))/static_cast<double>((i-ind1)*(i-ind1-1));
//        		}
//        		else // There is something wrong with the distribution. Give up momentum conservation
//        		{
//        			i = ind2;
//        			df = (icdf.at(i) - icdf.at(ind1))/(i-ind1);
//        			break;
//        		}

//         }
//         for (int j=ind1+1; j<i; j++)
//            icdf.at(j) = icdf.at(ind1) + df*(j-ind1);
// 	}
// 	lowval = icdf.at(i);
// 	i++;
// }
// std::cout << "4:\n";
// for (double & s: icdf) std::cout << s << "\t"; std::cout << "\n\n"; 
// double highval = density;
// i=nbins-1;
// while (i>imax)
// {
// 	if (highval < icdf.at(i))
// 	{        
// 		int ind1 = i+1, ind2 = i-1;
// 		double maxval = icdf.at(i); 
// 		while (highval < icdf.at(ind2))
// 		{
// 			ind2--;
// 			if (maxval < icdf.at(ind2))
// 			{
// 				maxval = icdf.at(ind2);
// 				i = ind2;
// 			}
// 		}
// 		while (maxval > icdf.at(ind1)) ind1++;  

// 		double S = 0;
// 		for (int j=i+1; j<ind1-1; j++)
// 			S += icdf.at(j);
// 		double df = 2/((ind1-i)*(ind1-i-1))*((ind1 - i -1) * icdf.at(ind1) - S);
// 		while (df < 0 || df*(ind1-i-1) > icdf.at(ind1) - icdf.at(i))
// 		{
// 			if (i>imax)
// 			{
// 				S += icdf.at(i);
// 				i--;
// 				df = 2/((ind1-i)*(ind1-i-1))*((ind1 - i -1) * icdf.at(ind1) - S);
// 			}
// 			else // There is something wrong with the distribution. Give up momentum conservation
// 			{
//        			i = ind2;
//        			df = (icdf.at(ind1) - icdf.at(i))/(ind1-i);
//        			break;
// 			}
// 		}
//         for (int j=ind1-1; j>i; j--)
//            icdf.at(j) = icdf.at(ind1) + df*(j-ind1);
// 	}

// 	highval = icdf.at(i);
// 	i--;
// }
// std::cout << "5:\n";
// for (double & s: icdf) std::cout << s << "\t"; std::cout << "\n\n"; 
// for (int i=0; i<nbins; i++)
// 	histogram.at(i) = icdf.at(i+1) - icdf.at(i);

// 	exit(1);
// }
}