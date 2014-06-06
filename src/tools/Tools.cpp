
#include "Tools.h"

void Tools::DisplayTitle()
{
	std::cout << std::endl << std::endl << "PlasmaScale - Based on ES1 - Electrostatic 1 Dimensional Code" << std::endl << std::endl;
	//printf("version 4.1\n");
	//printf("(c) Copyright 1987-92 Regents of the University of California\n");
	//printf("Plasma Theory and Simulation Group\n");
	//printf("University of California - Berkeley\n");
}//

void Tools::PUREShrink(std::vector<double>& scaling, std::vector<double>& detail, std::vector<bool>& pattern, const int& maxdepth, const int& mindepth)
{
	assert(scaling.size() == std::pow(2, maxdepth));
	static std::vector<double> breaks;
	static std::vector<std::pair<int,int> > index;
	static std::vector<double> theta;
	theta.resize(std::pow(2, maxdepth-1));
	breaks.reserve(std::pow(2, maxdepth-1));
	index.reserve(std::pow(2, maxdepth-1));


	for (int j=mindepth; j<maxdepth; j++)
	{
		int scale_size = std::pow(2, j);
		breaks.resize(0);
		index.resize(0);

        int b=0;
		for (int i=scale_size; i<2*scale_size; i++)
		{
			if (scaling.at(i) != 0)
			{
				breaks.push_back( std::abs(detail.at(i)) / std::sqrt(std::abs(scaling.at(i))) );
				index.emplace_back(i,b++);
			}
		}
		std::sort(index.begin(), index.end(), [&] (const std::pair<int,int>& a, const std::pair<int,int>& b) { return breaks.at(a.second) < breaks.at(b.second); });

		/* Find the value of the threshold minizing an approximation of the PURE estimate */

		double c1=0., c2=0., PURE=0., tmp;
        double a = breaks.at(0);
		for (int i=1; i<index.size(); i++)
		{
			int prev = index.at(i-1).first;
            
			c1 += std::abs(scaling.at(prev));
			c2 += 2*(std::abs(scaling.at(prev)) - std::pow(detail.at(prev), 2.));
			tmp = c1 * std::pow(breaks.at(index.at(i).second), 2.) + c2;
			if (tmp<PURE)
			{
				a = breaks.at(index.at(i).second);
				PURE = tmp;
			}
		}
		c2 += 2*(std::abs(scaling.at(index.back().first)) - std::pow(detail.at(index.back().first), 2.));
		if (c2 < PURE)
        {
			a=0;
            PURE = c2;
        }

		/* Update the detail coefficients and the sparsity pattern */

		for (int i=scale_size; i<2*scale_size; i++)
		{
			double d = std::abs(detail.at(i)) - a*std::sqrt(std::abs(scaling.at(i)));
			if (d>0)
			{
				detail.at(i) = d * (detail.at(i) > 0 ? 1 : -1);
				pattern.at(i) = true;
			}
			else
			{
				detail.at(i) = 0.;
				pattern.at(i) = false;
			}
		}
	}

}