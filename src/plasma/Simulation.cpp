#include "Simulation.h"


/* Initialization ============================================= */
void Simulation::Setup(int argc, char ** argv)
{
	FILE * InputDeck;

	if (!argc>1) 
	{
		InputDeck = std::fopen("es1data","r");
	}
	else 
	{
		InputDeck = std::fopen(argv[2],"r");
	}
	
	if (!InputDeck)
	{
		std::printf("\nCan't find input file %s\n",argv[1]);
		std::printf("\nCorrect syntax is: ES1 -i file.inp\n");
		std::exit(-1);
	}

	_state  = std::unique_ptr<MacroState> (new MacroState(InputDeck));
	std::fclose(InputDeck);

	_history = std::unique_ptr<History> (new History(_state->get_plasma()));

	_state->SetupDiagnostics(_diagnostics);
	_history->SetupDiagnostics(_diagnostics);

	_is_initialized = true;
}

/* methods ============================================================= */

void Simulation::Step()
{
	_state->Step();
	_history->Compute(*_state);
}

void Simulation::InitWindows()
{
	assert(this->is_initialized());

	/*********************************************/
	/* Set up each window structure              */
	/*********************************************/
	for (auto & it : _diagnostics)
		it->InitWindow();

	XGStart();
}