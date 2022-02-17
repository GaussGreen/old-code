/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file datestamp.cpp
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date February 2004
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

#include "gpbase/datestamp.h"
#include "gpbase/ostringstream.h"

#include <ctime>
#include <cstdlib>

#include <sstream>
CC_USING_NS (std,istringstream)

#include <string>
CC_USING_NS (std,string)

#include <algorithm>
CC_USING_NS (std,find)

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Routine: operator<<
///	Returns: ostream& 
///	Action : Enables to print ARM_TheMostSimpleDate
////////////////////////////////////////////////////
ostream& operator<<(ostream& os, const ARM_TheMostSimpleDate& date )
{
	vector<string> monthsName( ComputeSystemMonthsName() );
	os << monthsName[date.itsMonth] << "/" << date.itsDay << "/" << date.itsYear;
	return os;
}


////////////////////////////////////////////////////
///	Routine: Compute_MonthDiff
///	Returns: double
///	Action : difference of month between two ARM_TheMostSimpleDates
////////////////////////////////////////////////////

double Compute_MonthDiff( const ARM_TheMostSimpleDate& lhs, const ARM_TheMostSimpleDate& rhs )
{
	return lhs.itsMonth-rhs.itsMonth + (lhs.itsYear-rhs.itsYear)*12.0+ (lhs.itsDay-rhs.itsDay)/30.0;
}



////////////////////////////////////////////////////
///	Routine: ComputeCurrentDate
///	Returns: ARM_TheMostSimpleDate 
///	Action : Computes the current date at init time
////////////////////////////////////////////////////

ARM_TheMostSimpleDate ComputeCurrentDate()
{
	time_t aclock;
	time( &aclock );							/// Get time in seconds
	struct tm* runTime = localtime( &aclock );	///	Convert time to struct

	ARM_TheMostSimpleDate date = { runTime->tm_mon, runTime->tm_mday, runTime->tm_year +1900};

	/// although runTime is a pointor, you should not delete it!
	/// and there is no memory leak!
	return date;
}


////////////////////////////////////////////////////
///	Routine: Init_ComputeSystemMonthsName
///	Returns: vector<string> 
///	Action : Computes the vector of all the system months
////////////////////////////////////////////////////
vector<string> ComputeSystemMonthsName()
{
	vector<string> result(12);
	struct tm time = { 0,0,0,0,0,0};
	for(size_t i=0;i<12;++i)
	{
		char output[20];
		time.tm_mon = i;
		strftime(output,20,"%b", &time);
		result[i] = output;
	}
	return result;
}




////////////////////////////////////////////////////
///	Routine: ComputeCurrentDate
///	Returns: ARM_TheMostSimpleDate
///	Action : Computes the compile date at init time
////////////////////////////////////////////////////

ARM_TheMostSimpleDate ComputeCompileDate()
{
	/// get compile Time Stamp
	CC_Ostringstream os;
	os << __DATE__ << " " __TIME__;
	istringstream compileDateInput(os.str());
	int day;
	int year;
	int month = 0;
	string monthString;
	compileDateInput >> monthString >> day >> year;
	vector<string> monthsName( ComputeSystemMonthsName() );
	std::vector<string>::iterator found = find( monthsName.begin(), monthsName.end(), monthString );

	if( found != monthsName.end())
		month = found -monthsName.begin();

	ARM_TheMostSimpleDate date = { month, day, year };
	return date;
}


///////////////////////////////////////////////////
///	Variable : ARM_CurrentDate
///	Meaning  : it is the current date (run time date)
////////////////////////////////////////////////////
extern const ARM_TheMostSimpleDate ARM_CurrentDate = ComputeCurrentDate();


///////////////////////////////////////////////////
///	Variable : ARM_CompileDate
///	Meaning  : it is the compile date (compile time date)
////////////////////////////////////////////////////
extern const ARM_TheMostSimpleDate ARM_CompileDate = ComputeCompileDate();



CC_END_NAMESPACE()

///-----------------------------------------------------------------------------
/*---- End of file ----*/
