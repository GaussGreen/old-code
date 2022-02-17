/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file long_double.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_LONG_DOUBLE_H
#define _GP_CF_LONG_DOUBLE_H


#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>




CC_BEGIN_NAMESPACE(ARM)

class long_double
{
	public:
	double logarithm;
	int sign;
	long_double(const long_double& d):logarithm(d.logarithm),sign(d.sign){}
	long_double(const double& x)
	{
		logarithm=log(fabs(x));
		sign=(int)(x/fabs(x));
	}
	long_double(const long double& x)
	{
		logarithm=log(fabs(x));
		sign=(int)(x/fabs(x));
	}
	long_double(const int& x)
	{
		logarithm=log(fabs((double)x));
		sign=(int)(x/fabs((double)x));
	}
	long_double(const double& x, const int& m)
	{
		logarithm=x;
		sign=m;
	}
	long_double()
	{
		logarithm=0;
		sign=1;
	}
	double todouble() const
	{
		return sign*exp(logarithm);
	}
	
};


long_double operator*(const long_double& a,const long_double& b);
long_double operator/(const long_double& a,const long_double& b);
long_double logLD(const long_double& a);
long_double expLD(const long_double& a);
long_double expLD( const double& a);
long_double powLD(const long_double& a,const long_double& b);
long_double powLD(const long_double& a,const  double& b);
long_double powLD(const double& a,const  double& b);

bool operator==(const long_double& a ,const long_double& b);
bool operator!=(const long_double& a ,const long_double& b);
bool operator<=(const long_double& a ,const long_double& b);
bool operator>=(const long_double& a ,const long_double& b);
bool operator>(const long_double& a ,const long_double& b);
bool operator<(const long_double& a ,const long_double& b);


bool operator==(const double& a ,const long_double& b);
bool operator!=(const double& a ,const long_double& b);
bool operator<=(const double& a ,const long_double& b);
bool operator>=(const double& a ,const long_double& b);
bool operator>(const double& a ,const long_double& b);
bool operator<(const double& a ,const long_double& b);

bool operator==(const long_double& a ,const double& b);
bool operator!=(const long_double& a ,const double& b);
bool operator<=(const long_double& a ,const double& b);
bool operator>=(const long_double& a ,const double& b);
bool operator>(const long_double& a ,const double& b);
bool operator<(const long_double& a ,const double& b);

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/



