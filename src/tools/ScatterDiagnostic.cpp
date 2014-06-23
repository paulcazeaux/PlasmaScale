#include "ScatterDiagnostic.h"

ScatterDiagnostic::ScatterDiagnostic(const char * type, const char * x_label, const char * y_label, int ul_x, int ul_y,
				double x_scaling /* = 1.0 */, double y_scaling /* = 1.0 */, 
				bool auto_rescale_x /* = true */, bool auto_rescale_y /* = true */, 
				double x_min /* = 0.0 */, double x_max /* = 0.0 */, double y_min /* = 0.0 */, double y_max /* = 0.0 */,
				std::string state /* = std::string("closed") */) :
				Diagnostic(type, x_label, y_label, ul_x%881, (ul_y%676)+125, x_scaling, y_scaling, 
				auto_rescale_y, auto_rescale_y, x_min, x_max, y_min, y_max, 
				state) {}

void ScatterDiagnostic::Window()
{
	for (int i=0; i<_datasets.size(); i++)
	{
		XGScat2D(_datasets[i].first, _datasets[i].second, _sizes[i], _colors[i]); 
	}
}
