
/* HEADER MacroParameterizationFromFile ===============================================================
 * author: Paul Cazeaux
 * date: 2014-01-27
 *
 * description:
 *  -> store the parameters of a collection of populations of particles modeling the plasma
 *	-> Particle loading & pushing
 *  
 * ============================================================================= */

#ifndef DEF_PLASMASCALE_MACROPARAMETERIZATIONFROMFILE
#define DEF_PLASMASCALE_MACROPARAMETERIZATIONFROMFILE

#include "particles/MacroParameterization.h"

/* Declarations */

class MacroParameterizationFromFile : public MacroParameterization
{
	private:
		/* class members ======================================================================== */
		/* pointer to objects */

			// Initialization
		std::unique_ptr<std::vector<int> >		_group_sizes;
		std::unique_ptr<std::vector<double> >	_mean_velocities;
		std::unique_ptr<std::vector<double> >	_quiet_start_exponents;
		std::unique_ptr<std::vector<double> >	_quiet_mean_temperatures;
		std::unique_ptr<std::vector<double> >	_random_mean_temperatures;

			// Perturbation
		std::unique_ptr<std::vector<int> >		_modes;
		std::unique_ptr<std::vector<double> >	_density_amplitudes;
		std::unique_ptr<std::vector<double> >	_density_phases;
		std::unique_ptr<std::vector<double> >	_velocity_amplitudes;
		std::unique_ptr<std::vector<double> >	_velocity_phases;

			// Diagnostics
		std::unique_ptr<std::vector<int> >		_bin_numbers;
		std::unique_ptr<std::vector<double> >	_upper_velocities;
		std::unique_ptr<std::vector<double> >	_lower_velocities;

	public:

		MacroParameterizationFromFile(FILE *& InputDeck, int macro_grid_size);
		
		virtual void 	Load(State * state);

		double 			get_initial_temperature(int population_index);

		virtual bool 	HaveVelocityDiagnostics()	{return true;	}
		virtual double 	GetBinStart(int index);
		virtual	double	GetBinWidth(int index);
		virtual int		GetNumberOfBins(int index);
};

#endif