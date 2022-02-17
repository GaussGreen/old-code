
#ifndef __INTNORMAL_H__
#define __INTNORMAL_H__

#include "INTSTL.h"
#include "INTUtilities.h"


//----------------------------------------------------------------------------------------------

// Normal density function
 double Ndf(double x);

// Error function
 double erf(double x);

// Complementary error function
 double erfc(double x);

// Normal cumulative density function, Implementation 1
 double Ncdf1(double x);

// Normal cumulative density function, Implementation 2
 double Ncdf2(double x);

// Normal inverse cumulative density function, Implementation 1
 double Nicdf1(double u);

// Normal inverse cumulative density function, Implementation 2
 double Nicdf2(double u);

//----------------------------------------------------------------------------------------------


#endif





