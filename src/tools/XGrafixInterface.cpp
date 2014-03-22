#include "XGrafixInterface.h"

void XGMainLoop()
{
	PlasmaDevice.Step();
}

void Dump(char *filename) 
{	
	std::fstream fout;
	fout.open(filename, std::fstream::app);

	PlasmaDevice.WriteData(fout);

	fout.close();
}

void Quit() {}