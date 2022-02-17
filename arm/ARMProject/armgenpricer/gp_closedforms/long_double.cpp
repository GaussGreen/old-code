/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file long_double.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>

#include "gpclosedforms/long_double.h"

#include "expt.h"   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)


long_double operator*(const long_double& a,const long_double& b)
{
	return long_double(a.logarithm+b.logarithm,a.sign*b.sign);
}

long_double operator/(const long_double& a,const long_double& b)
{
	return long_double(a.logarithm-b.logarithm,a.sign*b.sign);
}

long_double powLD(const long_double& a,const long_double& b)
{
	if(a.sign!=1)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"pow: negative base" );
	}
	return long_double(a.logarithm*b.todouble(),1);
}

long_double powLD(const long_double& a,const double& b)
{
	if(a.sign!=1)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"pow: negative base" );
	}
	return long_double(a.logarithm*b,1);
}

long_double powLD(const double& a,const double& b)
{
	if(a<0)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"pow: negative base" );
	}
	double val=log(a)*b;
	return long_double(val,1);
}


long_double logLD(const long_double& a)
{
	if(a.sign!=1)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"pow: negative base" );
	}
	int signlog=-1;
	if (a.logarithm >=1) signlog=1;
	return long_double(log(fabs(a.logarithm)),signlog);
}

long_double expLD(const long_double& a)
{
	return long_double(a.todouble(),1);
}
long_double expLD(const double& a)
{
	return long_double(a,1);
}


bool operator==(const long_double& a ,const long_double& b)
{
	if((a.logarithm==b.logarithm) && (a.sign==b.sign)) return true;
		else return false;
}


bool operator!=(const long_double& a ,const long_double& b)
{
	if((a.logarithm==b.logarithm) && (a.sign==b.sign)) return false;
		else return true;
}


bool operator<=(const long_double& a ,const long_double& b)
{
	if((a.sign==1) && (b.sign!=1))
	{
		return false;
	}
	if((a.sign!=1) && (b.sign==1))
	{
		return true;
	}
	if((a.sign==1) && (b.sign==1))
	{
		return (a.logarithm<=b.logarithm);
	}
	if((a.sign!=1) && (b.sign!=1))
	{
		return (a.logarithm>=b.logarithm);
	}
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"operator<=(long_double a ,long_double b): should not reach here" );
	return false;
}

bool operator>=(const long_double& a ,const long_double& b)
{
	if((a.sign==1) && (b.sign!=1))
	{
		return true;
	}
	if((a.sign!=1) && (b.sign==1))
	{
		return false;
	}
	if((a.sign==1) && (b.sign==1))
	{
		return (a.logarithm>=b.logarithm);
	}
	if((a.sign!=1) && (b.sign!=1))
	{
		return (a.logarithm<=b.logarithm);
	}
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"operator>=(long_double a ,long_double b): should not reach here" );
	return false;
}

bool operator>(const long_double& a ,const long_double& b)
{
	if((a.sign==1) && (b.sign!=1))
	{
		return true;
	}
	if((a.sign!=1) && (b.sign==1))
	{
		return false;
	}
	if((a.sign==1) && (b.sign==1))
	{
		return (a.logarithm>b.logarithm);
	}
	if((a.sign!=1) && (b.sign!=1))
	{
		return (a.logarithm<b.logarithm);
	}
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"operator>=(long_double a ,long_double b): should not reach here" );
	return false;
}

bool operator<(const long_double& a ,const long_double& b)
{
	if((a.sign==1) && (b.sign!=1))
	{
		return false;
	}
	if((a.sign!=1) && (b.sign==1))
	{
		return true;
	}
	if((a.sign==1) && (b.sign==1))
	{
		return (a.logarithm<b.logarithm);
	}
	if((a.sign!=1) && (b.sign!=1))
	{
		return (a.logarithm>b.logarithm);
	}
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"operator<=(long_double a ,long_double b): should not reach here" );
	return false;
}



bool operator==(const double& a ,const long_double& b)
{
	long_double al(a);
	return al==b;
}
bool operator!=(const double& a ,const long_double& b)
{
	long_double al(a);
	return al!=b;
}

bool operator<=(const double& a ,const long_double& b)
{
	long_double al(a);
	return al<=b;
}

bool operator>=(const double& a ,const long_double& b)
{
	long_double al(a);
	return al>=b;
}

bool operator>(const double& a ,const long_double& b)
{
	long_double al(a);
	return al>b;
}

bool operator<(const double& a ,const long_double& b)
{
	long_double al(a);
	return al<b;
}

bool operator==(const long_double& a ,const double& b)
{
	long_double bl(b);
	return a==bl;
}
bool operator!=(const long_double& a ,const double& b)
{
	long_double bl(b);
	return a!=bl;
}


bool operator<=(const long_double& a ,const double& b)
{
	long_double bl(b);
	return a<=bl;
}

bool operator>=(const long_double& a ,const double& b)
{
	long_double bl(b);
	return a>=bl;
}

bool operator>(const long_double& a ,const double& b)
{
	long_double bl(b);
	return a>bl;
}
bool operator<(const long_double& a ,const double& b)
{
	long_double bl(b);
	return a<bl;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

