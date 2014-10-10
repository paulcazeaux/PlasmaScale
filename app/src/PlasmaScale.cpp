/* 
**	ES1 - ELECTROSTATIC 1 DIMENSIONAL PLASMA SIMULATION
**
*/

#include <memory>
#include "plasma/Simulation.h"
#include "xgrafix.h"
#include "tools/Tools.h"

Simulation PlasmaDevice;

int main(int argc, char **argv)
{
	Tools::DisplayTitle();

	/*********************************************************/
	/*	 Read input file, and initialize params and vars. 	 */
	/*********************************************************/
	PlasmaDevice.Setup(argc,argv);
	if (argv[argc-2][0] == '-' && argv[argc-2][1] == 'n' && argv[argc-2][2] == 'm')
		XGInit(argc-2, argv, PlasmaDevice.get_simulation_time_ptr());
	else
		XGInit(argc, argv, PlasmaDevice.get_simulation_time_ptr());
	PlasmaDevice.InitWindows();

	XGStart();
	return 0;
}

extern "C"
void XGMainLoop()
{
	PlasmaDevice.Step();
}

extern "C"
void Dump(char *filename) 
{	
	std::fstream fout;
	fout.open(filename, std::fstream::app);
	PlasmaDevice.WriteData(fout);
	fout.close();
}

extern "C"
void Quit() {}