/* HEADER  State ====================================================================================
 * author: Paul Cazeaux
 * date: 2013-11-29
 * 
 * description:
 *  -> store the state of the simulation
 * 
 * ============================================================================================== */

#ifndef DEF_PLASMASCALE_STATE
#define DEF_PLASMASCALE_STATE

/* includes ===================================================================================== */
#include <iostream>
#include <vector>
#include <memory>

#include "plasma/Plasma.h"
#include "fields/PlasmaFields.h"
#include "particles/PopulationOfParticles.h"
#include "parameterization/MacroParameterization.h"
#include "tools/CurveDiagnostic.h"
#include "tools/ScatterDiagnostic.h"

/* Forward declarations */
class MacroParameterization;

class State
{
	friend class History;

	private:
		/* class members ======================================================================== */

		/* pointer on object */
		std::shared_ptr<const Plasma>							_plasma;
		std::unique_ptr<PlasmaFields>							_fields;
		std::vector<std::unique_ptr<PopulationOfParticles> >	_populations;

		/* members */
		const int												_number_of_populations;
		std::shared_ptr<double>									_simulation_time;
		std::shared_ptr<int>									_iteration;

		/* singleton ============================================================================ */
		State(const State&);
		State& operator=(const State&);
		State& operator=(State);

	public:
		/* constuctor and destructor ============================================================ */
		State(	std::shared_ptr<const Plasma> plasma );
		State(	const MacroParameterization & parameterization);
		~State() {}

		/* getter */
		std::shared_ptr<double> get_simulation_time()	const 	{ return _simulation_time;	}
		double 	get_kinetic_energy(int population_index) 	const { return _populations.at(population_index)->get_kinetic_energy(); }
		int get_number_of_particles(int population_index) const { return _populations.at(population_index)->get_size();	};
		double get_total_weight(int population_index) const { return _populations.at(population_index)->get_total_weight();	};

		/* getters for the diagnostics ========================================================== */
		std::vector<std::vector<double> * >	get_vector_of_position_arrays() const;
		std::vector<std::vector<double> * >	get_vector_of_velocity_arrays() const;
		std::vector<std::vector<double> * > get_vector_of_weight_arrays() const;
		std::vector<int *>					get_vector_of_sizes() const;

		std::vector<std::vector<double> * > get_vector_of_bin_arrays() const;
		std::vector<std::vector<double> * > get_vector_of_velocity_profiles() const;
		std::vector<int *>					get_vector_of_number_of_bins() const;

		/* operator ============================================================================= */
		friend std::ostream& operator<<( std::ostream& os, const State& state);

		/* method =============================================================================== */

		void set_new_number_of_particles(const int population_index, const int new_population_size);
		void Load(const MacroParameterization & parameterization);
		void Prepare(const bool toggle_half_step = true);
		void Prepare(const int population_index, const bool toggle_half_step = true);

		void Reset()
		{
			for (auto & population : _populations)
			{
				population->Reset();
			}
		}

		void Weigh()
		{
			_fields->ResetCharge();
			for (auto & population : _populations)
			{
				population->Weigh(*_fields);
			}
			
			/* Account for BCs */
			_fields->_charge.front() 	*= 2.;
			_fields->_charge.back()		*= 2.;
		}

		void Accelerate(double factor = 1.0)
		{
			for (auto & population : _populations)
			{
				population->Accelerate(*_fields, factor);
			}
		}

		void Move()
		{
			for (auto & population : _populations)
			{
				population->Move();
			}
			this->Weigh();
		}

		void Step();
		void GetEField(std::vector<double> & Efield);
		void SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics);

		void ComputeVelocityProfile();
		void ComputeVelocityProfile(std::vector<std::vector<double> >& profile_by_population);

		void PushBackElectrostaticEnergy(std::vector<double>& electrostatic_energy) const;
		void PushBackKineticEnergy(std::vector<double> &kinetic_energy, std::vector<std::vector<double> >& kinetic_energy_by_population) const;
		void PushBackMoment(std::vector<double> &moment, std::vector<std::vector<double> >& moment_by_population) const;
};


 #endif