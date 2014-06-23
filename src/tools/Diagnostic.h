/* HEADER Diagnostic ===============================================================
 * author: Paul Cazeaux
 * date: 2013-11-28
 *
 * description:
 *  -> Generic diagnostic interface with the XGrafix library
 *  -> User-friendly interface
 *  
 * ============================================================================= */


#ifndef DEF_PLASMASCALE_DIAGNOSTIC
#define DEF_PLASMASCALE_DIAGNOSTIC

#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <vector>

#include "xgrafix.h"

class Diagnostic
{
	private:
		std::string _type; 		// one of "linlin", "linlog", "loglin", "loglog"
		std::string _x_label;
		std::string _y_label;
		std::string _state;		// one of "open" or "closed"
		int    _ul_x, _ul_y;	// position of upper left corner
		double _x_scaling;	// scaling factor for x array
		double _y_scaling;	// scaling factor for y array
		bool   _auto_rescale_x, _auto_rescale_y;
		double _x_min, _x_max;
		double _y_min, _y_max;

 	protected:
		std::vector<std::pair<double *, double*> >	_datasets;
		std::vector<int *>							_sizes;
		std::vector<int>							_colors;

	public:
		Diagnostic(std::string type, std::string x_label, std::string y_label, int ul_x, int ul_y, double x_scaling, double y_scaling, 
				bool auto_rescale_x, bool auto_rescale_y, double x_min, double x_max, double y_min, double y_max, std::string state = std::string("closed"));
		~Diagnostic() {}
		void AddData(std::vector<double> * x_data, std::vector<double> * y_data, int * length, int color);
		void AddData(double * x_data, double * y_data, int * length, int color);
		void InitWindow();
		void MenuItem();
		virtual void Window() {}
};

 #endif