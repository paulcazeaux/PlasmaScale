
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
	for (auto & val : icdf)
	{
		if (val < 0) val = 0;
		if (val > density) val = density;
	}

	/* Ensure symmetry of the operation */
	icdf_reverse.resize(nbins+1);
	std::copy(icdf.rbegin(), icdf.rend(), icdf_reverse.begin());


	/* Ensure that the ICDF is a monotonic function */
	double oldval=0.;
	for (int i=1; i<=nbins; i++)
	{
		if (icdf.at(i) < oldval)
		{
			int ind1 = i-1, ind2 = i;
			double high = oldval, low = icdf.at(i);

			/* The vector icdf is monotonic up to index i-1 */
            while (ind2<nbins && icdf.at(++ind2) < high)
            {
                if (icdf.at(ind2) < low) low = icdf.at(ind2);
            }
            while (ind1>0 && icdf.at(--ind1) > low) {}
			double val1 = icdf.at(ind1);
			double dn = (icdf.at(ind2)-val1)/static_cast<double>(ind2-ind1);
			for (int j=1; j<ind2-ind1; j++)
			{
				icdf.at(ind1+j) = val1 + static_cast<double>(j)*dn;
			}
            //for (int j=ind1+1; j<i; j++)
            //    icdf.at(j) = low;
		}
		oldval = icdf.at(i);
	}

	/* Same operation in reverse for the reverse ICDF */

				/* Ensure that the ICDF is a monotonic function */
	oldval=density;
	for (int i=1; i<=nbins; i++)
	{
		if (icdf_reverse.at(i)>density) 
		{
			icdf_reverse.at(i)=density;
		}
		if (icdf_reverse.at(i) > oldval)
		{
			int ind1 = i-1, ind2 = i;
			double low = oldval, high = icdf_reverse.at(i);
			/* The vector icdf_reverse is monotonic up to index i-1 */
            while (ind2<nbins && icdf_reverse.at(++ind2) > low)
            {
                if (icdf_reverse.at(ind2) > high)
                    high = icdf_reverse.at(ind2);
            }
            while (ind1>0 && icdf_reverse.at(--ind1) < high) {}
			double val1 = icdf_reverse.at(ind1);
			double dn = (icdf_reverse.at(ind2)-val1)/static_cast<double>(ind2-ind1);
			for (int j=1; j<ind2-ind1; j++)
			{
				icdf_reverse.at(ind1+j) = val1 + static_cast<double>(j)*dn;
			}
            // for (int j=ind1+1; j<i; j++)
            //     icdf_reverse.at(j) = high;
		}
		oldval = icdf_reverse.at(i);
	}

	/* Finish by averaging the two results */
	auto it_icdf_reverse = icdf_reverse.rbegin();
	for (auto & val : icdf)
	{
		val = .5 * (val + *it_icdf_reverse++);
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