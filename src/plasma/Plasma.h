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
#include <eigen3/Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat;

class Plasma
{
	private:
		/* Poisson operator data */
		SpMat										_M;
		Eigen::SimplicialCholesky<SpMat>			_solver;

		std::unique_ptr<std::vector<double>	>		_x_grid;
		std::unique_ptr<std::vector<double>	>		_macro_x_grid;
		std::unique_ptr<int>						_grid_end_ptr;
		std::unique_ptr<int> 						_macro_grid_end_ptr;


	public:
		/* class members ======================================================================== */

		/* Timestep and spatial grid step */
		const double		 						_dt;
		const double 								_dx;
		const double								_macro_dx;

		/* Parameters for extrapolation by EFPI */
		const int 									_number_of_microsteps;
		const int 									_macro_to_micro_dt_ratio;

		const bool 									_record_microsteps;
		const bool 									_use_full_PIC;

		/* size */
		const double 								_length;

		/* grid */
		const int 				 					_grid_end;
		const int 			 						_macro_grid_end;
		const int									_depth;
		const int									_cutoff;
		const int									_intensity;

		/* Plasma parameters */
		const double								_epsilon;
		const double								_E0;
		const double								_w0;



		/* Simulation parameters */
		const int 									_number_of_populations;

		/* Diagnostics parameters */
		const int									_max_size_history;
		int 										_velocity_accumulation_interval;
		bool										_profiling_active;

		/* singleton ============================================================================ */
		Plasma(const Plasma&);
		Plasma& operator=(const Plasma&);
		Plasma& operator=(Plasma);

		/* constuctor and destructor ============================================================ */
		Plasma(double 	length 						= 6.283185307,
				double	dt 							= 0.2,
				int number_of_microsteps 			= 50,
				int macro_to_micro_dt_ratio 		= 200,
				double epsilon 						= 1.0,
				double E0 							= 0,
				double w0 							= 0,
				int number_of_populations 			= 1,
				int grid_end 						= 128,
				int macro_grid_end 					= 32,
				int velocity_accumulation_interval 	= 1,
				int depth 							= 5,
				int cutoff 							= 5,
				double intensity 					= 1.,
				int max_size_history 				= 4096,
				bool use_full_PIC					= true,
				bool record_microsteps 				= false);
		Plasma(FILE * Inputdeck);
		~Plasma() {}

		/* operator << ========================================================================== */

		friend std::ostream& operator<<( std::ostream& os, const Plasma& plasma);

		/* getters ============================================================================== */


		/* grid pointers for the diagnostics */
		int *			get_grid_end_ptr()		const { return _grid_end_ptr.get(); }
		int *		get_macro_grid_end_ptr()	const { return _macro_grid_end_ptr.get(); }
		double *		get_x_grid_ptr()		const {	return _x_grid->data();	}
		double *		get_macro_x_grid_ptr()	const {	return _macro_x_grid->data(); }	
		
		/* methods ============================================================================== */
		int 			find_index_on_grid(double x)		const {return static_cast<int>(std::floor(x/_dx));	}
		double 			find_position_in_cell(double x) 	const {double xg = x/_dx; return xg - std::floor(xg);	}
		void 			Poisson_Solver(std::vector<double> & charge, std::vector<double> & potential) const;
};


 #endif