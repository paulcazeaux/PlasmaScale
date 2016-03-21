
/* HEADER MacroParameterization ===============================================================
 * author: Paul Cazeaux
 * date: 2014-01-27
 *
 * description:
 *  -> store the parameters of a collection of populations of particles modeling the plasma
 *	-> Particle loading & pushing
 *  
 * ============================================================================= */

#ifndef DEF_PLASMASCALE_MACROPARAMETERIZATION
#define DEF_PLASMASCALE_MACROPARAMETERIZATION

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "plasma/Plasma.h"
#include "tools/Diagnostic.h"

/* Forward declarations */
class State;

/* Declarations */

class MacroParameterization
{
	protected:
		/* class members ======================================================================== */
		/* pointer to objects */
		std::shared_ptr<const Plasma>	_plasma;

		int								_number_of_populations;

		/* Properties */
		std::vector<int>				_population_sizes;

		std::vector<double>				_reference_densities;
		std::vector<double>				_unit_charges;
		std::vector<double>				_unit_masses;
		std::vector<double>				_plasma_pulsations;
		std::vector<double> 			_thermal_velocities;

		/* Velocity bounds */
		std::vector<int>			_bin_numbers;
		std::vector<double>			_upper_velocities;
		std::vector<double>			_lower_velocities;

	public:
		/* constuctor and destructor ============================================================ */

		MacroParameterization() {}
		virtual ~MacroParameterization() {}

		/* move constuctor and assignment ======================================================= */

		MacroParameterization(MacroParameterization &&parameterization);
		MacroParameterization& operator=(MacroParameterization &&parameterization);

		/* getters ============================================================================== */

		std::shared_ptr<const Plasma> get_plasma() 				const	{ return _plasma;	}
		int		get_number_of_populations() 					const	{ return _number_of_populations;	}

		int		get_population_size(int index)					const 	{ return _population_sizes.at(index);	}
		double 	get_unit_charge(int index)						const 	{ return _unit_charges.at(index);		}
		double 	get_unit_mass(int index)						const 	{ return _unit_masses.at(index);		}

		/* virtual methods ====================================================================== */

		virtual void Load(State & state) const = 0;

		virtual void SetAccField(State & state) {}

		virtual void Initialize(State &) {}

		double 	get_initial_thermal_vel(int population_index)	const {	 return _thermal_velocities.at(population_index); }

		/* If implemented, the Step() function should leave the PIC micro state in a correctly initialized state */
		virtual void Step(State &) {}

        virtual bool 	HaveVelocityDiagnostics()       const;
		virtual double 	GetBinStart(int index)			const;
		virtual double 	GetBinEnd(int index)			const;
		virtual	double	GetBinWidth(int index)			const;
		virtual int		GetNumberOfBins(int index)		const;

		virtual void SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics) {}
		virtual void WriteData(std::fstream & fout) {}
		virtual void RecordMicroSteps(std::fstream & fout) {}
};

#endif