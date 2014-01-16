#include "Simulation.h"


/* Initialization ============================================= */
void Simulation::Setup(int argc, char ** argv)
{
	/* Parameters for the plasma */
	double length, dt, epsilon, la, e0, w0;
	int number_of_populations, grid_size, nt, iw, ec, velocity_accumulation_interval, max_mode;	// nt, iw, ec are unused
	double filter_parameter_1, filter_parameter_2;

	char a_char[80];
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
	/* read lines until we get to numbers */
	
	while (std::fscanf(InputDeck,"%d %lg %lg %d %d %lg %d", &number_of_populations, &length, &dt, &nt, &max_mode, &la, &velocity_accumulation_interval) <7)
	{
		std::fscanf(InputDeck, "%s", a_char);
	}
	/* note: la is l/a */

	
	while (std::fscanf(InputDeck," %d %d %d %lg %lg %lg %lg %lg", &grid_size, &iw, &ec, &epsilon, &filter_parameter_1, &filter_parameter_2, &e0, &w0) < 8)
	{
		std::fscanf(InputDeck, "%s", a_char);
	}

	if(velocity_accumulation_interval<0) 
	{ 
		std::printf("\nError:  accum can't be negative! \n"); exit(1);
	}

	printf(" nsp = %2d     l = %8.5f \n", number_of_populations, length);
	printf(" dt = %4.5f    nt = %4d \n",dt,nt);
	printf(" ng = %5d   iw = %2d   ec = %2d  accum = %4d\n", grid_size, iw, ec, velocity_accumulation_interval);
	printf(" epsi = %4.2f  a1 = %4.2f  a2 = %4.2f \n",epsilon, filter_parameter_1, filter_parameter_2);

	_plasma = std::make_shared<const Plasma>(length, dt, 
				epsilon, la, 0.0, e0, w0,
				number_of_populations, grid_size, velocity_accumulation_interval, max_mode,
				filter_parameter_1, filter_parameter_2);

	_state  = std::unique_ptr<State> (new State(_plasma, InputDeck));
	_history = std::unique_ptr<History> (new History(_plasma));

	std::fclose(InputDeck);

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