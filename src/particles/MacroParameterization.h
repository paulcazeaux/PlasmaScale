
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

#include <iostream>
#include <string>
#include <vector>

#include "plasma/Plasma.h"
#include "plasma.State.h"

/* Declarations */

class MacroParameterization
{
	protected:
		/* class members ======================================================================== */
		/* pointer to objects */
		std::shared_ptr<const Plasma>			_plasma;

		int										_number_of_populations;

		/* Properties */

		std::unique_ptr<std::vector<int> >		_population_sizes;

		std::unique_ptr<std::vector<double> >	_unit_charges;
		std::unique_ptr<std::vector<double> >	_unit_masses;
		std::unique_ptr<std::vector<double> >	_cyclotronic_rotation_parameters;

	public:
		/* constuctor and destructor ============================================================ */

		virtual MacroParameterization() {}
		virtual ~MacroParameterization() {}

		/* move constuctor and assignment ======================================================= */

		MacroParameterization(MacroParameterization &&parameterization);
		MacroParameterization& operator=(MacroParameterization &&parameterization);

		/* getters ============================================================================== */

		std::shared_ptr<const Plasma> get_plasma() 				const	{ return _plasma;	}
		int		get_number_of_populations() 					const	{ return _number_of_populations;	}

		int		get_population_size(int index)					const 	{ return _population_sizes->at(index);	}
		double 	get_unit_charge(int index)						const 	{ return _unit_charges->at(index);		}
		double 	get_unit_mass(int index)						const 	{ return _unit_masses->at(index);		}
		double 	get_cyclotronic_rotation_parameter(int index)	const 	{ return _cyclotronic_rotation_parameters->at(index);	}

		/* virtual methods ====================================================================== */

		virtual void Load(State * state);

		virtual void RestrictAndPushback(State * state_) {}
		virtual void ExtrapolateAndLoad() {}

		virtual bool Have_Velocity_Diagnostics()	{return false;	}
		virtual double 	GetBinStart()				{return 0.;		}
		virtual	double	GetBinWidth()				{return 0.;		}
		virtual	int		GetNumberOfBins()			{return 0;		}
};

#endif