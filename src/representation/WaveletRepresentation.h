/* HEADER WaveletRepresentation ===============================================================
 * author: Paul Cazeaux
 * date: 2014-04-25
 *
 * description:
 *  -> store the wavelet representation for the velocity distribution of a set of particles
 *  -> implement some useful operations: DWT and its inverse, denoising, particle weighting
 *  -> Implement some linear operations: addition, multiplication by scalar.
 *  
 * ============================================================================= */

#ifndef DEF_WAVELETREPRESENTATION
#define DEF_WAVELETREPRESENTATION

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cassert>
#include "plasma/Plasma.h"
#include "representation/Representation.h"

#include <chrono>

/* Declarations */
 
class WaveletRepresentation : public Representation
{
	protected:
		/* class members ======================================================================== */
		/* pointer on object */
		std::shared_ptr<const Plasma>		_plasma;

		/* parameters */
		double 								_vmin;
		double 								_vmax;
		int 								_depth;
		int 								_number_of_bins;
		int 								_grid_end;
		double 								_reference_density;
		double 								_dv;
		std::vector<double> 				_low_pass;
		std::vector<double> 				_high_pass;
		int 								_filter_length;

		bool 								_is_transformed;

		/* grid */
		std::vector<std::vector<double> >	_histogram;
		std::vector<std::vector<double> > 	_coefficients;
		std::vector<int>					_length;

	public:
		/* constructor  ========================================================================= */
		WaveletRepresentation() {};
		WaveletRepresentation(std::shared_ptr<const Plasma> plasma, double vmin, double vmax, int depth, double reference_density);

		/* methods */
		virtual void Weigh(int size,
				std::vector<double>::iterator 	position,
				std::vector<double>::iterator  	velocity,
				std::vector<double>::iterator 	weight);
		virtual void Weigh(int size,
				std::vector<double>::iterator 	position,
				std::vector<double>::iterator  	velocity,
				std::vector<double>::iterator 	weight,
				const double delay,
				const std::vector<double> & accfield);
		virtual void Load(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights);
		virtual void Coarsen();
		virtual void Refine();

		virtual void SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const double & thermal_velocity) { /* TODO */ };
		virtual void SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const std::vector<double> & thermal_velocity) { /* TODO */ };
		virtual void GetDensityVelocity(std::vector<double> & density, std::vector<double> & velocity);
		virtual void GetDensityVelocityPressure(std::vector<double> & density, std::vector<double> & velocity, std::vector<double> & pressure);
		virtual void GetVelocityMoment(double& moment);

		void DWT(const int J);
		void iDWT();
		void Denoise(const int n_coef);
		void Denoise(const double thresh, const int J);
		void Cutoff(const int depth);
		void DiscardNegativeValues();

		virtual void Reset();
		virtual void print(std::ostream& os) const;
		virtual int get_grid_end() const 	{ return _grid_end;	}
		virtual void set_grid_end(const int new_grid_end)
		{
			if (!_is_transformed)
			{
				_histogram.resize(new_grid_end+1);
				for (int n=_grid_end+1; n<=new_grid_end; n++)
					_histogram.at(n).resize(_number_of_bins);
			}
			else
			{
				_coefficients.resize(new_grid_end+1);
				for (int n=_grid_end+1; n<=new_grid_end; n++)
					_coefficients.at(n).resize(_number_of_bins);
			}
			_grid_end = new_grid_end;
		}

		/* operator overload */
		WaveletRepresentation& operator=(const WaveletRepresentation & rhs)
		{
			_plasma = rhs._plasma;
			_vmin  = rhs._vmin;
			_vmax = rhs._vmax;
			_depth = rhs._depth;
			_number_of_bins = rhs._number_of_bins;
			_grid_end = rhs._grid_end;
			_reference_density = rhs._reference_density;
			_dv = rhs._dv;
			_low_pass = rhs._low_pass;
			_high_pass = rhs._high_pass;
			_filter_length = rhs._filter_length;

			_is_transformed = rhs._is_transformed;

			_histogram.resize(_grid_end+1);
			_coefficients.resize(_grid_end+1);

            if (_is_transformed)
            {
                std::copy(rhs._coefficients.begin(), rhs._coefficients.end(), _coefficients.begin());
                std::copy(rhs._length.begin(), rhs._length.end(), _length.begin());
            }
            else
            {
                std::copy(rhs._histogram.begin(), rhs._histogram.end(), _histogram.begin());
            }
			return *this;
		}

		WaveletRepresentation& operator+=(const WaveletRepresentation& rhs)
		{
			assert((_grid_end==rhs._grid_end)&&(_number_of_bins==rhs._number_of_bins)&&(_vmax==rhs._vmax)&&(_vmin==rhs._vmin)&&(_reference_density=rhs._reference_density));

			if (_is_transformed)
			{
				for (int n=0; n<=_grid_end; n++)
				{
					for (int i=0; i<_number_of_bins; i++)
					{
						_coefficients.at(n).at(i) += rhs._coefficients.at(n).at(i);
					}
				}
			}
			else
			{
				for (int n=0; n<=_grid_end; n++)
				{
					for (int i=0; i<_number_of_bins; i++)
					{
						_histogram.at(n).at(i) += rhs._histogram.at(n).at(i);
					}
				}
			}
			return *this;
		}

		WaveletRepresentation& operator-=(const WaveletRepresentation& rhs)
		{
			assert((_grid_end==rhs._grid_end)&&(_number_of_bins==rhs._number_of_bins)&&(_vmax==rhs._vmax)&&(_vmin==rhs._vmin)&&(_reference_density=rhs._reference_density));

			if (_is_transformed)
			{
				for (int n=0; n<=_grid_end; n++)
				{
					for (int i=0; i<_number_of_bins; i++)
					{
						_coefficients.at(n).at(i) -= rhs._coefficients.at(n).at(i);
					}
				}
			}
			else
			{
				for (int n=0; n<=_grid_end; n++)
				{
					for (int i=0; i<_number_of_bins; i++)
					{
						_histogram.at(n).at(i) -= rhs._histogram.at(n).at(i);
					}
				}
			}
			return *this;
		}

		WaveletRepresentation& operator*=(const double lambda)
		{
			if (_is_transformed)
			{
				for (int n=0; n<=_grid_end; n++)
				{
					for (int i=0; i<_number_of_bins; i++)
					{
						_coefficients.at(n).at(i) *= lambda;
					}
				}
			}
			else
			{
				for (int n=0; n<=_grid_end; n++)
				{
					for (int i=0; i<_number_of_bins; i++)
					{
						_histogram.at(n).at(i) *= lambda;
					}
				}
			}
			return *this;
		}
};


inline WaveletRepresentation operator+(WaveletRepresentation lhs, const WaveletRepresentation& rhs)
{
	lhs += rhs;
	return lhs;
}

inline WaveletRepresentation operator-(WaveletRepresentation lhs, const WaveletRepresentation& rhs)
{
	lhs -= rhs;
	return lhs;
}

inline WaveletRepresentation operator*(WaveletRepresentation lhs, const double lambda)
{
	lhs *= lambda;
	return lhs;
}

inline WaveletRepresentation operator*(const double lambda, WaveletRepresentation rhs)
{
	rhs *= lambda;
	return rhs;
}

#endif