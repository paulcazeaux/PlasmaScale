
#include "Tools.h"

void Tools::DisplayTitle()
{
	printf("\n\nES1 - Electrostatic 1 Dimensional Code\n");
	printf("version 4.1\n");
	printf("(c) Copyright 1987-92 Regents of the University of California\n");
	printf("Plasma Theory and Simulation Group\n");
	printf("University of California - Berkeley\n");
}

double Tools::EvaluateSlope(const std::vector<double> & values)
{
	int size = values.size();
	double sum = 0, first_moment = 0;
	for (int n=0; n<size; n++)
	{
		sum			+= values.at(n);
		first_moment 	+= static_cast<double>(n)*values.at(n);
	}

	return 6./static_cast<double>(size*(size+1)) * (2./static_cast<double>(size-1)*first_moment - sum);
}