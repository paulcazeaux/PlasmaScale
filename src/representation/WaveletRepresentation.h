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
#include "plasma/Plasma.h"
#include "wavelets/wavelet2s.h"

class WaveletRepresentation : public Representation
{
	private:
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
		std::vector<std::vector<double> >				_histogram;
		std::vector<std::vector<double> > 				_coefficients;
		std::vector<std::vector<double> > 				_flag;
		std::vector<std::vector<int> >					_length;

	public:
		/* constructor  ========================================================================= */
		WaveletRepresentation() {};
		WaveletRepresentation(std::shared_ptr<const Plasma> plasma, double vmax, int depth, int grid_size);

		/* methods */
		virtual void Weigh(int size,
				std::vector<double>::iterator 	position,
				std::vector<double>:iterator  	velocity,
				std::vector<double>::iterator 	weight);
		virtual void Load(int size,
								std::vector<double>::iterator 	position,
								std::vector<double>:iterator  	velocity,
								std::vector<double>::iterator 	weight);
		virtual void Coarsen();
		virtual void Refine();
		void DWT();
		void iDWT();
		void Denoise();

		virtual void Reset();
		virtual void print(std::ostream& os) const;
		virtual int get_grid_size() const 	{ return _grid_size;	}

		/* operator overload */
		WaveletRepresentation& operator+=(const WaveletRepresentation& rhs)
		{
			assert((_grid_size==rhs._grid_size)&&(_number_of_bins==rhs._number_of_bins)&&(_vmax==rhs._vmax));

			if _is_transformed
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

		WaveletRepresentation& operator*=(const double lambda)
		{
			if _is_transformed
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

#endif