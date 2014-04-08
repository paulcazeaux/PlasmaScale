/* HEADER Field ===============================================================
 * author: Paul Cazeaux
 * date: 2013-11-28
 *
 * description:
 *  -> store the values of a field (potential, energy, density...) in the plasma
 *  -> Particle weighting
 *  
 * ============================================================================= */

#ifndef DEF_PLASMASCALE_FIELD
#define DEF_PLASMASCALE_FIELD

#include <iostream>
#include <string>
#include <vector>
#include <fftw3.h>
#include "plasma/Plasma.h"

class Field
{
	friend class PlasmaFields;

	private:
		/* class members ======================================================= */
		/* pointer on object */
		std::shared_ptr<const Plasma>		_plasma;

		/* grid */
		const int 							_size;

		/* values */
		double*								_values;

		fftw_plan 							_fft_forward;
		fftw_plan							_fft_inverse;
	public:
		/* constuctor and destructor =========================================== */
		Field(std::shared_ptr<const Plasma> plasma);
		~Field();

		/* getters */

		double * get_ptr()			const   {return _values;	}
		int 	 get_size()			const 	{return _size; 		}

		/* methods ==============================================================*/

		void Reset();
		void WeighParticle(const double position, const double charge);

		/* FFTs =================================================================*/

		void FFT()
		{
			fftw_execute(_fft_forward);
		}

		void iFFT()
		{
			fftw_execute(_fft_inverse);
		}


};


 #endif