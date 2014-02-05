#include "PlasmaFields.h"



/* Constructor =============================================================== */

PlasmaFields::PlasmaFields(std::shared_ptr<const Plasma> plasma,
					std::shared_ptr<double> simulation_time,
					std::shared_ptr<int> iteration)
					: 	_plasma(plasma),
					  	_simulation_time(simulation_time),
					  	_iteration(iteration),

						_charge(Field(plasma)),
						_electrical_field(Field(plasma)),
						_potential(Field(plasma)),

						_length		(plasma->get_length()),
						_epsilon	(plasma->get_epsilon()),
						_E0			(plasma->get_E0()),
						_w0			(plasma->get_w0()),
						_grid_size  (plasma->get_grid_size()),
						_dx			(plasma->get_dx()),
				_filter_parameter_1	(plasma->get_filter_parameter_1()),
				_filter_parameter_2	(plasma->get_filter_parameter_2()),
						_max_mode	(plasma->get_max_mode())
{

	_electrostatic_energy_by_mode = std::vector<double>(_max_mode);
	_electrostatic_energy_total  = 0.;
	/* Initialize the filter and the filtered Poisson kernel */
	std::vector<double> filtered_kernel 	= std::vector<double>(_grid_size/2 + 1		);
	std::vector<double> filter 				= std::vector<double>(_grid_size/2 + 1, 1.	);
	double kdx2, kperp2 = std::pow(AKPERP * plasma->get_inverse_of_particle_radius(), 2), temp, temp1;

	for (int k = 1; k < _grid_size/2 + 1; k++) 
	{
		kdx2 = k * (M_PI/_grid_size);
		if ((_filter_parameter_1 != 0.0) || (_filter_parameter_2 != 0.0)) 
		{
			if (k < _grid_size)
			{
				temp = std::pow( std::tan(kdx2), 4);
			}
			else
			{
				temp = 1e20;
			}
			temp1 = std::pow(std::sin(kdx2), 2);

			filter.at(k) = std::exp(_filter_parameter_1*temp1 - _filter_parameter_2*temp);
		}

		filtered_kernel.at(k) = filter.at(k) * filter.at(k) 
								* _epsilon / ( kperp2 + 4.0*std::pow( std::sin(kdx2)/ _dx, 2) ) ;
	}
	_filtered_kernel 	= std::unique_ptr<const std::vector<double> > (new std::vector<double>( filtered_kernel	));
	_filter 			= std::unique_ptr<const std::vector<double> > (new std::vector<double>( filter 			));
}

void PlasmaFields::ComputeAndFilter()
{
	int half_grid_size 			= _grid_size / 2;
	double * charge_array 		= _charge.get_ptr();
	double * potential_array	= _potential.get_ptr();
	double * e_field_array 		= _electrical_field.get_ptr();
	
	/* Rescale charge */
	for (int j=0; j < _grid_size; j++) 
	{
		charge_array[j] *=  _dx;
	}

	/* Fourier transform into frequency space */
	_charge.FFT();

	charge_array[0] = 0.; // neutralizing background
	potential_array[0] = 0.;

	/* Compute the filtered Fourier transform of the potential */
	double eses = 0;
	for (int k = 1; k < half_grid_size; k++) 
	{
		int sym_k 				= _grid_size - k;
		potential_array[k] 		= _filtered_kernel->at(k)*charge_array[k    ]; 
		potential_array[sym_k]	= _filtered_kernel->at(k)*charge_array[sym_k]; 
		eses 	   		   		+= charge_array[k    ]*potential_array[k    ];
		eses			   		+= charge_array[sym_k]*potential_array[sym_k];
	}
	potential_array[half_grid_size] = _filtered_kernel->at(half_grid_size) * charge_array[half_grid_size];
	eses += 0.5 * charge_array[half_grid_size] * potential_array[half_grid_size];

	/* Compute the electrostatic energy, total and by mode */

	_electrostatic_energy_total = eses/_length;

	for (int k = 1; k <= _max_mode; k++) 
	{
		int sym_k = _grid_size - k;
		_electrostatic_energy_by_mode[k-1] = (charge_array[k]*potential_array[k] + charge_array[sym_k]*potential_array[sym_k])/_length;
		if (k == sym_k)
		{
			_electrostatic_energy_by_mode[k-1] *= 0.25;
		}
	}

	/* Filter the charge density */
	for (int k=1; k < half_grid_size; k++)
	{
		int sym_k 		 = _grid_size - k;
		charge_array[k]		  	*= _filter->at(k); 
		potential_array[sym_k] 	*= _filter->at(k); 
	}
	charge_array[half_grid_size] *= _filter->at(half_grid_size);  

	/* Normalization */
	for (int k=0; k < _grid_size ;k++)
	{
		charge_array[k]    /= _length;
		potential_array[k] /= _length;
	}

	/* Back to normal space */
	_charge.iFFT();
	_potential.iFFT();

	/* Compute the electrical field */

	double external_e_field = _E0*cos(_w0*(*_simulation_time));
	for (int j=1; j<_grid_size-1; j++)
	{
		e_field_array[j] 			= (potential_array[j-1]-potential_array[j+1])/(2.*_dx) + external_e_field;
	}
	e_field_array[0] 				= (potential_array[_grid_size - 1]-potential_array[1])/(2.*_dx) + external_e_field;
	e_field_array[_grid_size -1] 	= (potential_array[_grid_size - 2]-potential_array[0])/(2.*_dx) + external_e_field;
}

void PlasmaFields::WeighParticle(double position, double charge)
{
	_charge.WeighParticle(position, charge);
}

void PlasmaFields::PushBackElectrostaticEnergy(std::vector<double>& electrostatic_energy, std::vector<std::vector<double>	>& electrostatic_energy_by_mode) const
{
	electrostatic_energy.push_back(_electrostatic_energy_total  + 1e-30);
	for (int m = 0; m<_max_mode; m++)
	{
		electrostatic_energy_by_mode.at(m).push_back(_electrostatic_energy_by_mode[m] + 1e-30);
	}
}