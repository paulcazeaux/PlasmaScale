/* HEADER PUREHaarRepresentation ===============================================================
 * author: Paul Cazeaux
 * date: 2014-04-25
 *
 * description:
 *  -> store the wavelet representation for the velocity distribution of a set of particles
 *  -> implement some useful operations: DWT and its inverse, denoising, particle weighting
 *  -> Implement some linear operations: addition, multiplication by scalar.
 *  
 * ============================================================================= */

#ifndef DEF_PUREHAARREPRESENTATION
#define DEF_PUREHAARREPRESENTATION

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <eigen3/Eigen/Dense>
#include "wavelets/wavelet2s.h"
#include "plasma/Plasma.h"
#include "tools/Tools.h"
#include "tools/HaarTools.h"
#include "representation/Representation.h"

class PUREHaarRepresentation : public Representation
{
	protected:
		/* class members ======================================================================== */
		/* pointer on object */
		std::shared_ptr<const Plasma>		_plasma;

		/* parameters */
		double 								_vmax;
		int 								_max_depth;
		int 								_number_of_bins;
		int 								_grid_size;
		double 								_dv;

		bool 								_is_transformed;

		/* Parameters for adapativity */
		int 								_min_depth;
		int 								_buffer;

		/* grid */
		std::vector<std::vector<double> >	_histogram;
		std::vector<std::vector<double> > 	_scaling_coefficients;
		std::vector<std::vector<double> >   _detail_coefficients;

		std::vector<std::vector<int> >  	_mask;

	public:
		/* constructor  ========================================================================= */
		PUREHaarRepresentation() {};
		PUREHaarRepresentation(std::shared_ptr<const Plasma> plasma, double vmax, int max_depth, int grid_size, int min_depth = 1, int buffer = 1);

		/* methods */
		virtual void Weigh(int size,
				std::vector<double>::iterator 	position,
				std::vector<double>::iterator  	velocity,
				std::vector<double>::iterator 	weight);
		virtual void Load(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>::iterator  	velocity,
								std::vector<double>::iterator 	weights);
		virtual void Coarsen();
		virtual void Refine();
		void 		 Transform();

		virtual void SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const double & thermal_velocity);
		virtual void SetAdiabaticValues(const std::vector<double> & density, const std::vector<double> & velocity, const std::vector<double> & thermal_velocity);
		virtual void GetDensityVelocity(std::vector<double> & density, std::vector<double> & velocity) const;
		virtual void GetDensityVelocityPressure(std::vector<double> & density, std::vector<double> & velocity, std::vector<double> & pressure) const;

		void Denoise(int n_coef);
		void Denoise(double thresh);
		void Cutoff(int depth);
		void PUREAdapt(const double intensity = 1.);

		virtual void Reset();
		virtual void print(std::ostream& os) const;
		virtual int get_grid_size() const 	{ return _grid_size;	}
		virtual void set_grid_size(const int new_grid_size)
		{
			_histogram.resize(new_grid_size);
			_mask.resize(new_grid_size);
			for (int n=_grid_size; n<new_grid_size; n++)
			{
				_histogram.at(n).resize(_number_of_bins);
				_mask.at(n).resize(_number_of_bins);
			}

			_scaling_coefficients.resize(new_grid_size);
			_detail_coefficients.resize(new_grid_size);
			_grid_size = new_grid_size;
		}

		/* operator overload */
		PUREHaarRepresentation& operator=(const PUREHaarRepresentation & rhs)
		{
			_plasma = rhs._plasma;
			_vmax = rhs._vmax;
			_max_depth = rhs._max_depth;
			_number_of_bins = rhs._number_of_bins;
			_grid_size = rhs._grid_size;
			_dv = rhs._dv;

			_is_transformed = rhs._is_transformed;

			_histogram.resize(_grid_size);
			_scaling_coefficients.resize(_grid_size);
			_detail_coefficients.resize(_grid_size);
			_mask.resize(_grid_size);

            if (_is_transformed)
            {
                std::copy(rhs._scaling_coefficients.begin(), rhs._scaling_coefficients.end(), _scaling_coefficients.begin());
                std::copy(rhs._detail_coefficients.begin(), rhs._detail_coefficients.end(), _detail_coefficients.begin());
            }
            else
            {
                std::copy(rhs._histogram.begin(), rhs._histogram.end(), _histogram.begin());
            }                
            std::copy(rhs._mask.begin(), rhs._mask.end(), _mask.begin());

			return *this;
		}

		PUREHaarRepresentation& operator+=(const PUREHaarRepresentation& rhs)
		{
			assert((_grid_size==rhs._grid_size)&&(_number_of_bins==rhs._number_of_bins)&&(_vmax==rhs._vmax));

			if (_is_transformed)
			{
				for (int n=0; n<_grid_size; n++)
				{
					for (int i=0; i<_number_of_bins; i++)
					{
						_scaling_coefficients.at(n).at(i) += rhs._scaling_coefficients.at(n).at(i);
						_detail_coefficients.at(n).at(i) += rhs._detail_coefficients.at(n).at(i);
						_mask.at(n)[i] = _mask.at(n)[i] || rhs._mask.at(n)[i];
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
						_mask.at(n)[i] = _mask.at(n)[i] || rhs._mask.at(n)[i];
					}
				}
			}
			return *this;
		}

		PUREHaarRepresentation& operator-=(const PUREHaarRepresentation& rhs)
		{
			assert((_grid_size==rhs._grid_size)&&(_number_of_bins==rhs._number_of_bins)&&(_vmax==rhs._vmax));

			if (_is_transformed)
			{
				for (int n=0; n<_grid_size; n++)
				{
					for (int i=0; i<_number_of_bins; i++)
					{
						_scaling_coefficients.at(n).at(i) -= rhs._scaling_coefficients.at(n).at(i);
						_detail_coefficients.at(n).at(i) -= rhs._detail_coefficients.at(n).at(i);
						_mask.at(n)[i] = _mask.at(n)[i] || rhs._mask.at(n)[i];

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
						_mask.at(n)[i] = _mask.at(n)[i] || rhs._mask.at(n)[i];

					}
				}
			}
			return *this;
		}

		PUREHaarRepresentation& operator*=(const double lambda)
		{
			if (_is_transformed)
			{
				for (int n=0; n<_grid_size; n++)
				{
					for (int i=0; i<_number_of_bins; i++)
					{
						_scaling_coefficients.at(n).at(i) *= lambda;
						_detail_coefficients.at(n).at(i) *= lambda;
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


inline PUREHaarRepresentation operator+(PUREHaarRepresentation lhs, const PUREHaarRepresentation& rhs)
{
	lhs += rhs;
	return lhs;
}

inline PUREHaarRepresentation operator-(PUREHaarRepresentation lhs, const PUREHaarRepresentation& rhs)
{
	lhs -= rhs;
	return lhs;
}

inline PUREHaarRepresentation operator*(PUREHaarRepresentation lhs, const double lambda)
{
	lhs *= lambda;
	return lhs;
}

inline PUREHaarRepresentation operator*(const double lambda, PUREHaarRepresentation rhs)
{
	rhs *= lambda;
	return rhs;
}

#endif