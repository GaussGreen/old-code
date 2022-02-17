//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : LessEqualGreaterEps.cpp
//
//   Description : Comparisons with epsilon (un)tolerance as the objects
//
//   Date        : March 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/LessEqualGreaterEps.hpp"
#include "edginc/ModelException.hpp"

DRLIB_BEGIN_NAMESPACE

///////////////////////////////////////////////////////////////////////////////
strictLessEps::strictLessEps(double gEpsilon)
	:	EPSILON(gEpsilon)
{
	if (EPSILON < 0.0)
	{
		throw ModelException("EPSILON < 0", "strictLessEps");
	}
}

bool strictLessEps::operator()(double g1, double g2) const
{
	return ((g1 + EPSILON) < g2);
}
///////////////////////////////////////////////////////////////////////////////
looseLessEps::looseLessEps(double gEpsilon)
	:	EPSILON(gEpsilon)
{
	if (EPSILON < 0.0)
	{
		throw ModelException("EPSILON < 0", "looseLessEps");
	}
}

bool looseLessEps::operator()(double g1, double g2) const
{
	return ( g1 <= (g2 + EPSILON) );
}

///////////////////////////////////////////////////////////////////////////////
equalEps::equalEps(double gEpsilon)
	:	EPSILON(gEpsilon)
{
	if (EPSILON < 0.0)
	{
		throw ModelException("EPSILON < 0", "equalEps");
	}
}

bool equalEps::operator()(double g1, double g2) const
{
	return (((g1 + EPSILON) > g2) && (g1 < (g2+ EPSILON)));
};

///////////////////////////////////////////////////////////////////////////////
looseGreaterEps::looseGreaterEps(double gEpsilon)
	:	EPSILON(gEpsilon)
{
	if (EPSILON < 0.0)
	{
		throw ModelException("EPSILON < 0", "looseGreaterEps");
	}
}

bool looseGreaterEps::operator()(double g1, double g2) const
{
	return ( (g1 + EPSILON) >= g2 );
}

///////////////////////////////////////////////////////////////////////////////
strictGreaterEps::strictGreaterEps(double gEpsilon)
	:	EPSILON(gEpsilon)
{
	if (EPSILON < 0.0)
	{
		throw ModelException("EPSILON < 0", "strictGreaterEps");
	}
}

bool strictGreaterEps::operator()(double g1, double g2) const
{
	return (g1 > (g2 + EPSILON));
};

///////////////////////////////////////////////////////////////////////////////

DRLIB_END_NAMESPACE
