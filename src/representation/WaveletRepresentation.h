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
#include "wavelets/wavelet2s.h"
#include "plasma/Plasma.h"
#include "representation/Representation.h"

class WaveletRepresentation : public Representation
{
	protected:
		/* class members ======================================================================== */
		/* pointer on object */
		std::shared_ptr<const Plasma>		_plasma;

		/* parameters */
		double 								_vmax;
		int 								_depth;
		int 								_number_of_bins;
		int 								_grid_size;
		double 								_dv;
		string 								_filter;

		bool 								_is_transformed;

		/* grid */
		std::vector<std::vector<double> >	_histogram;
		std::vector<std::vector<double> > 	_coefficients;
		std::vector<std::vector<double> > 	_flag;
		std::vector<std::vector<int> >		_length;

	public:
		/* constructor  ========================================================================= */
		WaveletRepresentation() {};
		WaveletRepresentation(std::shared_ptr<const Plasma> plasma, double vmax, int depth, int grid_size);

		/* methods */
		virtual void Weigh(int size,
				std::vector<double>::iterator 	position,
				std::vector<double>::iterator  	velocity,
				std::vector<double>::iterator 	weight,
				const double delay = 0.);
		virtual void Load(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights);
		virtual void Coarsen();
		virtual void Refine();

		virtual void SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const double & thermal_velocity);
		virtual void SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const std::vector<double> & thermal_velocity);
		virtual void GetDensityVelocity(std::vector<double> & density, std::vector<double> & velocity) const;
		virtual void GetDensityVelocityPressure(std::vector<double> & density, std::vector<double> & velocity, std::vector<double> & pressure) const;

		void DWT();
		void iDWT();
		void Denoise(int n_coef);
		void Denoise(double thresh);
		void Cutoff(int depth);
		void DiscardNegativeValues();

		virtual void Reset();
		virtual void print(std::ostream& os) const;
		virtual int get_grid_size() const 	{ return _grid_size;	}
		virtual void set_grid_size(const int new_grid_size)
		{
			_histogram.resize(new_grid_size);
			for (int n=_grid_size; n<new_grid_size; n++)
				_histogram.at(n).resize(_number_of_bins);

			_coefficients.resize(new_grid_size);
			_flag.resize(new_grid_size);
			_length.resize(new_grid_size);
			_grid_size = new_grid_size;
		}

		/* operator overload */
		WaveletRepresentation& operator=(const WaveletRepresentation & rhs)
		{
			_plasma = rhs._plasma;
			_vmax = rhs._vmax;
			_depth = rhs._depth;
			_number_of_bins = rhs._number_of_bins;
			_grid_size = rhs._grid_size;
			_dv = rhs._dv;
			_filter = rhs._filter;

			_is_transformed = rhs._is_transformed;

			_histogram.resize(_grid_size);
			_coefficients.resize(_grid_size);
			_flag.resize(_grid_size);
			_length.resize(_grid_size);

            if (_is_transformed)
            {
                std::copy(rhs._coefficients.begin(), rhs._coefficients.end(), _coefficients.begin());
                std::copy(rhs._flag.begin(), rhs._flag.end(), _flag.begin());
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
			assert((_grid_size==rhs._grid_size)&&(_number_of_bins==rhs._number_of_bins)&&(_vmax==rhs._vmax));

			if (_is_transformed)
			{
				for (int n=0; n<_grid_size; n++)
				{
					for (int i=0; i<_number_of_bins; i++)
					{
						_coefficients.at(n).at(i) += rhs._coefficients.at(n).at(i);
					}
				}
			}
			else
			{
				for (int n=0; n<_grid_size; n++)
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
			assert((_grid_size==rhs._grid_size)&&(_number_of_bins==rhs._number_of_bins)&&(_vmax==rhs._vmax));

			if (_is_transformed)
			{
				for (int n=0; n<_grid_size; n++)
				{
					for (int i=0; i<_number_of_bins; i++)
					{
						_coefficients.at(n).at(i) -= rhs._coefficients.at(n).at(i);
					}
				}
			}
			else
			{
				for (int n=0; n<_grid_size; n++)
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
				for (int n=0; n<_grid_size; n++)
				{
					for (int i=0; i<_number_of_bins; i++)
					{
						_coefficients.at(n).at(i) *= lambda;
					}
				}
			}
			else
			{
				for (int n=0; n<_grid_size; n++)
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