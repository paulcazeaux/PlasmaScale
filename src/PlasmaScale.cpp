/* 
**	ES1 - ELECTROSTATIC 1 DIMENSIONAL PLASMA SIMULATION
**
*/

#include "PlasmaScale.h"

Simulation PlasmaDevice;

int main(int argc, char **argv)
{
	Tools::DisplayTitle();

	/*********************************************************/
	/*	 Read input file, and initialize params and vars. 	 */
	/*********************************************************/
	PlasmaDevice.Setup(argc,argv);

	XGInit(argc, argv, PlasmaDevice.get_simulation_time_ptr());
	PlasmaDevice.InitWindows();

	XGStart();
	return 0;
}