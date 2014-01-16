#include "RandomTools.h"

template<>
int RandomTools::Generate_randomly_uniform<int>(const int start, const int end)
{
	static std::random_device rd;
	static std::mt19937_64 gen( rd() );
	std::uniform_int_distribution<int> dis(start, end);
	return dis(gen);
}

template<>
double RandomTools::Generate_randomly_uniform<double>(const double start, const double end)
/* randomly genrate a double between start and end */
{
	static std::random_device rd;
	static std::mt19937_64 gen( rd() );
	std::uniform_real_distribution<double> dis(start, end);
	return dis(gen);
}

double RandomTools::Generate_randomly_normal(const double mean, const double stddev)
{
	static std::random_device rd;
	static std::mt19937_64 gen( rd() );
	std::normal_distribution<double> dis(mean, stddev);
	return dis(gen);
}