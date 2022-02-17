// DExpression.cpp: implementation of the DExpression class.
//
//////////////////////////////////////////////////////////////////////

#include "kexpression.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

KDoubleExpression::UnaryFunctionMap KDoubleExpression::UnaryMathFunctions = KDoubleExpression::InitUnary();
KDoubleExpression::BinaryFunctionMap KDoubleExpression::BinaryMathFunctions = KDoubleExpression::InitBinary();


double KDoubleExpression::max (double a, double b) {return (a>b) ? a : b;}
double KDoubleExpression::min (double a, double b) {return (a<b) ? a : b;}

KDoubleExpression::UnaryFunctionMap KDoubleExpression::InitUnary() 
{
	UnaryFunctionMap m;
	m["ACOS"] = acos; 
	m["ASIN"] = asin; 
	m["ATAN"] = atan; 
	m["COS"] = cos; 
	m["COSH"] = cosh; 
	m["EXP"] = exp; 
	m["FABS"] = fabs; 
	m["LOG"] = log; 
	m["LOG10"] = log10; 
	m["SIN"] = sin; 
	m["SINH"] = sinh; 
	m["TAN"] = tan; 
	m["TANH"] = tanh; 
	m["SQRT"] = sqrt; 
	m["CEIL"] = ceil; 
	m["FLOOR"] = floor; 
	return m;
}

KDoubleExpression::BinaryFunctionMap KDoubleExpression::InitBinary() 
{
	BinaryFunctionMap m;
	m["ATAN2"] = atan2;	
	m["FMOD"] = fmod;	
	m["MAX"] = max;
	m["MIN"] = min;
	m["POW"] = pow;	
	return m;
}
