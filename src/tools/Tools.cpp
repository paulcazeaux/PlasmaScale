
#include "Tools.h"

void Tools::DisplayTitle()
{
	std::cout << std::endl << std::endl << "PlasmaScale - Based on ES1 - Electrostatic 1 Dimensional Code" << std::endl << std::endl;
	//printf("version 4.1\n");
	//printf("(c) Copyright 1987-92 Regents of the University of California\n");
	//printf("Plasma Theory and Simulation Group\n");
	//printf("University of California - Berkeley\n");
}//

void Tools::AssembleICDF(const std::vector<double>& histogram, std::vector<double>& icdf, double& density)
{
	int nbins = histogram.size();
	static std::vector<double> icdf_reverse;
	icdf.resize(nbins+1);
	icdf.front() = 0.;
	std::partial_sum(histogram.begin(), histogram.end(), icdf.begin()+1);
	density = icdf.back();

	/* First, we take care of negative values and overshots */
	int i0 = 0, i1 = nbins, imax = 0;
	for (int i=1; i<nbins; i++)
	{
		if (icdf.at(i) < 0) i0 = i;
		if (icdf.at(i) > density && i1 > i) i1 = i;
		if (icdf.at(i) > .5*density && icdf.at(i-1) < .5*density) imax = i;
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
			if (ind1 > 1)
				ind1--;
			else
			{
				S += icdf.at(ind2);
				ind2++;
				ind1 = i0-1;
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
			if (ind2 < nbins+1)
				ind2++;
			else
			{
				S += icdf.at(ind1);
				ind1--;
				ind2 = i1+1;
			}
			df = 2.*((nbins-ind1-1)*density-S)/static_cast<double>((ind2-ind1)*(ind2-ind1-1));
		}
		for (int j = ind1+1; j<ind2; j++)
			icdf.at(j) = density + (j-ind2)*df;
		for (int j=ind2; j<nbins; j++)
			icdf.at(j) = density;
	}

	/* Ensure that the ICDF is a monotonic function */
	double lowval=0.;
	int i = 1;
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
            	S += icdf.at(i);
            	i++;
           		df = 2.*(S - (i - ind1 - 1)*icdf.at(ind1))/static_cast<double>((i-ind1)*(i-ind1-1));

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
				S += icdf.at(i);
				i--;
				df = 2/((ind1-i)*(ind1-i-1))*((ind1 - i -1) * icdf.at(ind1) - S);
			}
		}

		highval = icdf.at(i);
		i--;
	}
}


double Tools::ComputePoissonEmpiricalScaling(const std::vector<std::vector<double> >& histogram, const int& patch_size)
{
	/* Empirical scaling : detect local empirical mean and variance */
	int grid_size = histogram.size();
	int number_of_bins = histogram.front().size();
	static std::vector<double> patch_mean = std::vector<double>(grid_size*number_of_bins/(patch_size*patch_size));
	static std::vector<double> patch_variance = std::vector<double>(grid_size*number_of_bins/(patch_size*patch_size));
    
	auto it_mean = patch_mean.begin();
	auto it_variance = patch_variance.begin();
	for (int n=0; n<grid_size; n+=8)
	{
		for (int i=0; i<number_of_bins; i+=8)
		{
			double mean=0., variance=0.;
			for (int np=n; np<n+patch_size; np++)
			{
				for (int ip=i; ip<i+patch_size; ip++)
				{
					double val=histogram.at(np).at(ip);
					mean += val;
					variance += val*val;
				}
			}
			mean /= static_cast<double>(patch_size*patch_size);
			variance = variance/static_cast<double>(patch_size*patch_size) - mean*mean;
            
			*it_mean++ 		= mean;
			*it_variance++ 	= variance;
		}
	}
    
	return Tools::EvaluateSlope(patch_mean, patch_variance);
}