#define XLFunction_cpp
#include "XLFunction.h"

XLFunction::XLFunction (const char* n_name, int n_type, const char* n_helptext, const char* n_result_helptext, int n_result_type) :
	name (n_name), type (n_type), arguments (), helptext (n_helptext), result_helptext (n_result_helptext), result_type (n_result_type)
{
}

// EOF %M%
	