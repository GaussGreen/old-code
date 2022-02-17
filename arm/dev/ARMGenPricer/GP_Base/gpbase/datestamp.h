/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file datestamp.h
 *
 *  \brief date stamp to avoid old xlls
 *	\author  E Benhamou
 *	\version 1.0
 *	\date FEBRUARY 2004
 */

#ifndef _INGPBASE_DATESTAMP_H
#define _INGPBASE_DATESTAMP_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
#include "env.h"
#include "port.h"

#include <iostream>
CC_USING_NS (std,ostream)

#include <vector>
CC_USING_NS (std,vector)

CC_BEGIN_NAMESPACE( ARM )

/// \struct ARM_TheMostSimpleDate
/// very light structure to get the month,day,year of a date..
struct ARM_TheMostSimpleDate
{
	int itsMonth;
	int itsDay;
	int itsYear;
};

/// function to print an ARM_TheMostSimpleDate
ostream& operator<<(ostream& os, const ARM_TheMostSimpleDate& date );

/// compute the month difference between two ARM_TheMostSimpleDates
double Compute_MonthDiff( const ARM_TheMostSimpleDate& lhs, const ARM_TheMostSimpleDate& rhs );

///  computes the current date
ARM_TheMostSimpleDate ComputeCurrentDate();

/// vector of System months name
vector<string> ComputeSystemMonthsName();

/// Computes the compile date
ARM_TheMostSimpleDate ComputeCompileDate();
 
/// global constant for current date
extern const ARM_TheMostSimpleDate ARM_CurrentDate;

/// global constant for compile date
extern const ARM_TheMostSimpleDate ARM_CompileDate;



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
