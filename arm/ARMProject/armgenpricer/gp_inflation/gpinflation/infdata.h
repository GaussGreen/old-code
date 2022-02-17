/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: infdata.h,v $
 * Revision 1.11  2003/11/20 15:21:20  ebenhamou
 * added publish lag
 *
 * Revision 1.10  2003/09/17 18:12:07  ebenhamou
 * move ARM_Constants to armglob.cpp
 *
 * Revision 1.9  2003/09/05 07:20:32  ebenhamou
 * namespace constant
 *
 * Revision 1.8  2003/08/26 11:48:31  ebenhamou
 * added copyright
 *
 * Revision 1.7  2003/08/25 07:20:43  ebenhamou
 * avoid crash when index not defined
 *
 * Revision 1.6  2003/08/05 08:28:26  ebenhamou
 * change Getits into Get
 * same for Set
 *
 * Revision 1.5  2003/07/29 08:51:40  ebenhamou
 * change following code review of MAB
 *
 * Revision 1.3  2003/07/18 08:59:29  ebenhamou
 * more explicit terminology
 *
 * Revision 1.1  2003/07/16 07:01:45  ebenhamou
 * Initial revision
 *
 *
 */


	
#ifndef _INGPINFLATION_INFDATA_H
#define _INGPINFLATION_INFDATA_H

#include <string>
#include "gpbase/port.h"


CC_BEGIN_NAMESPACE( ARM )


/// for namespace reason
/// use std::string in the header

/// Inflation definition
/// this should move to armglob
/// but for the time being 
/// plugged here to avoid
/// another painful recompilation


/// for each index
/// provide a ticker corresponding to Bloomberg, Reuters
/// ccy and calendar monthly interpolation default method
/// daily interpolation default method
/// and the Day Count Fraction between the fixing and the date 
/// used to compute linear interpolation

typedef struct 
{
    char*	ticker;
    char*	ccy;
	char*	calendar;
	long	monthInterp;
	long	dailyInterp;
	char*	DCFLag;
	char*	ResetLag;
	long	DCFMonthly;
	long	DCFDaily;
	long	ExtrapolType;
	char*   PublishLag;
} Inf_Table_Row;



/// the various accessors to the data
/// the data are encapsulated
/// only accessors are public
class InfData{
private:
	static Inf_Table_Row InfTable[];
public:
	static int GetIndexPos( const char* dataName, const char* itsIndexName );
	static char* GetCurrency( const char* itsIndexName );
	static char* GetCalendar( const char* itsIndexName );
	static long GetMonthInterpolation( const char* itsIndexName );
	static long GetDailyInterpolation( const char* itsIndexName );
	static std::string GetDCFLag( const char* itsIndexName );
	static std::string GetResetLag( const char* itsIndexName );
	static long GetDCFMonthly( const char* itsIndexName );
	static long GetDCFDaily( const char* itsIndexName );
	static long GetExtrapolType( const char* itsIndexName );
	static std::string GetPublishLag(const char* itsIndexName );
};


typedef struct
{
	long	methodFlag;
	char*	methodName;
} interpMappingRow;


/// to encapsulate data for interpolation mapping
class InfInterp{
private:
	static interpMappingRow MappingInterpTable[100];
public:
	static char* GetMappingName( long methodFlag );
};

CC_END_NAMESPACE()

#endif
