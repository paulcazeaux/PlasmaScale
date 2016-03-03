/* HEADER CurveDiagnostic ===============================================================
 * author: Paul Cazeaux
 * date: 2013-11-28
 *
 * description:
 *  -> Curve plot diagnostic interface with the XGrafix library
 *  -> User-friendly interface
 *  
 * ============================================================================= */


#ifndef DEF_PLASMASCALE_CURVEDIAGNOSTIC
#define DEF_PLASMASCALE_CURVEDIAGNOSTIC

#include <cassert>
#include "Diagnostic.h"


 class CurveDiagnostic : public Diagnostic
 {
 	public:
 		CurveDiagnostic(const char * type, const char * x_label, const char * y_label, int ul_x, int ul_y, 
 				double x_scaling = 1.0, double y_scaling = 1.0, 
				bool auto_rescale_x = true, bool auto_rescale_y = true, 
				double x_min = 0.0, double x_max = 0.0, double y_min = 0.0, double y_max = 0.0,
				std::string state = std::string("closed"));
 		~CurveDiagnostic() {}

 		void Window();
 };


#endif