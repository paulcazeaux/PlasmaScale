#include "Diagnostic.h"

Diagnostic::Diagnostic(std::string type, std::string x_label, std::string y_label, int ul_x, int ul_y, double x_scaling, double y_scaling, 
				bool auto_rescale_x, bool auto_rescale_y, double x_min, double x_max, double y_min, double y_max):
				_type(type), _x_label(x_label), _y_label(y_label),
				_ul_x(ul_x), _ul_y(ul_y),
				_x_scaling(x_scaling), _y_scaling(y_scaling),
				_auto_rescale_x(auto_rescale_x), _auto_rescale_y(auto_rescale_y),
				_x_min(x_min), _x_max(x_max), _y_min(y_min), _y_max(y_max)
			{
				_state = std::string("closed");
			}

void Diagnostic::AddData(std::vector<double> * x_data, std::vector<double> * y_data, int * length, int color)
{
	assert(x_data->size() == y_data->size());
	assert(x_data->size() == *length);

	_datasets.push_back(std::pair<double *, double*>(x_data->data(), y_data->data()));
	_sizes.push_back(length);
	_colors.push_back(color);
}

void Diagnostic::AddData(double * x_data, double * y_data, int * length, int color)
{
	// Potentially unsafe
	_datasets.push_back(std::pair<double *, double*>(x_data, y_data));
	_sizes.push_back(length);
	_colors.push_back(color);
}

void Diagnostic::InitWindow()
{
	this->MenuItem();
	this->Window();
}
void Diagnostic::MenuItem()
{
	XGSet2D(_type.c_str(), _x_label.c_str(), _y_label.c_str(), _state.c_str(), _ul_x, _ul_y,
		 _x_scaling, _y_scaling, _auto_rescale_x, _auto_rescale_y,
		 _x_min, _x_max, _y_min, _y_max);
}