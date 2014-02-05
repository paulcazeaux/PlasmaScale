/* HEADER  Plasma ====================================================================================
 * author: Paul Cazeaux
 * date: 2013-11-29
 * 
 * description:
 *  -> store the main parameters of the plasma (timestep, length...)
 *	-> act as a library
 * 
 * ============================================================================================== */

#ifndef DEF_PLASMASCALE_PLASMA
#define DEF_PLASMASCALE_PLASMA

/* includes ===================================================================================== */
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <memory>
#include <algorithm>



class Plasma
{
	private:
		/* class members ======================================================================== */

		/* Timestep and spatial grid step */
		double		 								_dt;
		double 										_dx;

		/* Parameters for extrapolation by EPFI */
		int 										_number_of_microsteps;
		int 										_macro_to_micro_dt_ratio;

		/* size */
		const double 								_length;
		double										_inverse_of_particle_radius;

		/* grid */
		std::unique_ptr<int> 						_grid_size;
		std::unique_ptr<int> 						_macro_grid_size;
		std::unique_ptr<std::vector<double> >		_x_grid;
		std::unique_ptr<std::vector<double> >		_k_grid;
		int 										_highest_mode;

		/* Plasma parameters */
		const double								_epsilon;
		const double								_rho0;
		double										_E0;
		double										_w0;

		/* Simulation parameters */
		const int 									_number_of_populations;

		const double								_filter_parameter_1;
		const double								_filter_parameter_2;

		/* Diagnostics parameters */
		const int									_max_size_history;
		const int 									_max_mode;
		int 										_velocity_accumulation_interval;
		bool										_profiling_active;

		/* singleton ============================================================================ */
		Plasma(const Plasma&);
		Plasma& operator=(const Plasma&);
		Plasma& operator=(Plasma);

	public:
		/* constuctor and destructor ============================================================ */
		Plasma(double 	length 						= 6.283185307,
				double	dt 							= 0.2,
				int number_of_microsteps 			= 50,
				int macro_to_micro_dt_ratio 		= 200,
				double epsilon 						= 1.0,
				double la							= 0.,
				double rho0							= 1.,
				double E0 							= 0,
				double w0 							= 0,
				int number_of_populations 			= 1,
				int grid_size 						= 128,
				int macro_grid_size 				= 32,
				int velocity_accumulation_interval 	= 1,
				int max_mode 						= 1,
				double filter_parameter_1			= 0.,
				double filter_parameter_2 			= 0.,
				int max_size_history 				= 4096 );
		Plasma(FILE * Inputdeck);
		~Plasma() {}

		/* operator << ========================================================================== */

		friend std::ostream& operator<<( std::ostream& os, const Plasma& plasma);

		/* getters ============================================================================== */
		/* Timestep and spatial grid step */
		double			get_dt()				const {	return _dt;			}
		double 			get_dx()				const {	return _dx;			}

		/* Parameters for extrapolation by EPFI */
		int 	get_number_of_microsteps()		const { return _number_of_microsteps;		}
		int 	get_macro_to_micro_dt_ratio()	const { return _macro_to_micro_dt_ratio;	}								

		/* Length */
		const double 	get_length()			const {	return _length;		}

		/* grid */
		const int 		get_grid_size()			const { return *_grid_size;	}
		const int 		get_macro_grid_size()	const { return *_macro_grid_size;	}
		int *			get_grid_size_ptr()		const { return _grid_size.get(); }
		double *		get_x_grid_ptr()		const {	return _x_grid->data();	}
		double *		get_k_grid_ptr()		const {	return _k_grid->data();	}

		/* Plasma parameters */
		const double 	get_inverse_of_particle_radius()	const {	return _inverse_of_particle_radius;	}
		const double	get_epsilon()						const {	return _epsilon;	}
		const double	get_rho0()							const {	return _rho0;		}
		double			get_E0()							const {	return _E0;			}
		double			get_w0()							const {	return _w0;			}

		/* Simulation parameters */
		const int		get_number_of_populations()			const { return _number_of_populations;}

		const double	get_filter_parameter_1()			const {	return _filter_parameter_1;			}
		const double	get_filter_parameter_2()			const {	return _filter_parameter_2;			}

		/* Diagnostics parameters */
		const double	get_max_size()						const { return _max_size_history;	}
		const int 		get_max_mode()						const { return _max_mode;			}
		const int	get_velocity_accumulation_interval()	const {	return _velocity_accumulation_interval;	}
		const bool		is_velocity_profiling_active()		const { return _profiling_active;	}
		
		/* setters ============================================================================== */
		/* Timestep */
		void			set_dt(	double dt )		{	_dt = dt;	}

		/* Plasma parameters */

		void			set_E0( double E0 )		{	_E0 = E0;	}
		void			set_w0( double w0 )		{	_w0 = w0;	}	

		/* Multiscale parameters */
		void 	set_macro_grid_size( int macro_grid_size)	{ *_macro_grid_size = macro_grid_size;	}
		
		/* methods ============================================================================== */
		int 			find_index_on_grid( double x)		const {return static_cast<int>(std::floor(x/_dx)) %  *_grid_size;	}
		double 			find_position_in_cell(double x) 	const {double xg = x/_dx; return xg - std::floor(xg);	}
};


 #endif