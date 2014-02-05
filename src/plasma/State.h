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
#include "tools/CurveDiagnostic.h"
#include "tools/ScatterDiagnostic.h"

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

		std::shared_ptr<double> get_simulation_time()	const { return _simulation_time;	}

		/* getters for the diagnostics ========================================================== */
		std::vector<std::vector<double> * >	get_vector_of_position_arrays();
		std::vector<std::vector<double> * >	get_vector_of_x_velocity_arrays();
		std::vector<std::vector<double> * >	get_vector_of_y_velocity_arrays();
		std::vector<std::vector<double> * > get_vector_of_weight_arrays();
		std::vector<bool> 					get_vector_of_magnetizations();
		std::vector<int *>					get_vector_of_sizes();

		std::vector<std::vector<double> * > get_vector_of_bin_arrays();
		std::vector<std::vector<double> * > get_vector_of_velocity_profiles();
		std::vector<int *>					get_vector_of_number_of_bins();

		/* operator ============================================================================= */
		friend std::ostream& operator<<( std::ostream& os, const State& state);

		/* method =============================================================================== */
		void Load(const MacroParameterization & parameterization);

		void Reset()
		{
			for (auto & population : _populations)
			{
				population->Reset();
			}
		}

		void Prepare()
		{
			_fields.ComputeAndFilter();
			for (auto & population : _populations)
			{
				population->Prepare(*_fields);
			}

		}

		void Weigh()
		{
			_fields->ResetCharge();
			for (auto & population : _populations)
			{
				population->Weigh(*_fields);
			}
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
		void SetupDiagnostics(std::vector<std::unique_ptr<Diagnostic> > &diagnostics);

		void ComputeVelocityProfile();

		void PushBackKineticEnergy(std::vector<double> &kinetic_energy, std::vector<std::vector<double> >	&kinetic_energy_by_population);
		void PushBackMoment(std::vector<double> &moment, std::vector<std::vector<double> >	&moment_by_population);
};


 #endif