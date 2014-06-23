/* HEADER ScatterDiagnostic ===============================================================
 * author: Paul Cazeaux
 * date: 2013-11-28
 *
 * description:
 *  -> Scatterplot diagnostic interface with the XGrafix library
 *  -> User-friendly interface
 *  
 * ============================================================================= */


#ifndef DEF_PLASMASCALE_SCATTERDIAGNOSTIC
#define DEF_PLASMASCALE_SCATTERDIAGNOSTIC

#include <cassert>
#include "Diagnostic.h"


 class ScatterDiagnostic : public Diagnostic
 {
 	public:
 		ScatterDiagnostic(const char * type, const char * x_label, const char * y_label, int ul_x, int ul_y, 
 				double x_scaling = 1.0, double y_scaling = 1.0, 
				bool auto_rescale_x = true, bool auto_rescale_y = true, 
				double x_min = 0.0, double x_max = 0.0, double y_min = 0.0, double y_max = 0.0,
				std::string state = std::string("closed"));
 		~ScatterDiagnostic() {}

 		void Window();
 };


#endif