
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

		std::vector<double>				_unit_charges;
		std::vector<double>				_unit_masses;
		std::vector<double>				_plasma_pulsations;
		std::vector<double>				_cyclotronic_rotation_parameters;

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
		double 	get_cyclotronic_rotation_parameter(int index)	const 	{ return _cyclotronic_rotation_parameters.at(index);	}

		/* virtual methods ====================================================================== */

		virtual void Load(State & state) const 
		{
			std::cout << "A Load function should be implemented for your choice of parameterization !" << std::endl;
		}
		virtual void SetAccField(State & state) {}

		virtual void Initialize(State &) {}

		/* If implemented, the Step() function should leave the PIC micro state in a correctly initialized state */
		virtual void Step(State &) {}

		virtual bool HaveVelocityDiagnostics()					const 	{return false;	}
		virtual double 	GetBinStart(int)						const 	{return 0.;		}
		virtual double 	GetBinEnd(int)							const 	{return 0.;		}
		virtual	double	GetBinWidth(int)						const 	{return 0.;		}
		virtual	int		GetNumberOfBins(int)					const 	{return 0;		}

		virtual void SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics) {}
		virtual void WriteData(std::fstream & fout) {}
		virtual void RecordMicroSteps(std::fstream & fout) {}
};

#endif