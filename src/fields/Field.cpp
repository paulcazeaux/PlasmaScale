#include "Field.h"

Field::Field(std::shared_ptr<const Plasma> plasma) : 
	_plasma(plasma), _size(plasma->get_grid_size())
/* constructor from Plasma */
{
	_values = fftw_alloc_real(_size + 1);

	_fft_forward = fftw_plan_r2r_1d(_size, _values, _values, 
						FFTW_R2HC, FFTW_MEASURE);
	_fft_inverse = fftw_plan_r2r_1d(_size, _values, _values, 
						FFTW_HC2R, FFTW_MEASURE);
}


/* Destructor */
Field::~Field()
{
	fftw_destroy_plan(_fft_forward);
	fftw_destroy_plan(_fft_inverse);
	fftw_free(_values);
}



/* Reset operator */
void Field::Reset()
{
	for (int i=0; i < _size; i++)
	{
		_values[i] = 0;
	}
}

/* Add a charge from a particle on the charge density field */
void Field::WeighParticle(const double position, const double charge)
{
		static double dx 	= _plasma->get_dx();
		int bin 			= _plasma->find_index_on_grid(position);
		double cellpos		= _plasma->find_position_in_cell(position);

		_values[bin] += (1. - cellpos) * charge / dx;
		if (bin < _size - 1)
			_values[bin+1] += cellpos * charge / dx;
		else
			_values[0] += cellpos * charge / dx;
}
