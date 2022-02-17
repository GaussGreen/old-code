/*
 * $Log: infdata.cpp,v $
 * Revision 1.11  2003/11/20 15:21:37  ebenhamou
 * added publish lag
 *
 * Revision 1.10  2003/09/29 08:17:55  ebenhamou
 * using macro CC_USING_NS
 *
 * Revision 1.9  2003/09/22 13:52:20  ebenhamou
 * using __ARM_NO_NAMESPACE instead of std
 *
 * Revision 1.8  2003/09/02 17:47:47  ebenhamou
 * change default
 *
 * Revision 1.7  2003/08/25 07:20:31  ebenhamou
 * avoid crash when index not defined
 *
 * Revision 1.6  2003/08/05 08:28:44  ebenhamou
 * change Getits into Get same for Set
 *
 * Revision 1.5  2003/07/29 08:52:09  ebenhamou
 * change following code review of MAB
 *
 * Revision 1.3  2003/07/18 08:59:13  ebenhamou
 * more explicit terminology
 *
 * Revision 1.1  2003/07/16 07:01:52  ebenhamou
 * Initial revision
 *
 *
 *
 */


#include "gpinflation/infdata.h"

#include <glob/expt.h>	// neceassary for exception throwing
CC_USING_NS( std, pair )

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class      : InfInterp
///	Description: table to Get the name for an interpolation method
////////////////////////////////////////////////////
interpMappingRow InfInterp::MappingInterpTable[100] =
{
	/// methodFlag			methodName
	{	K_CPILINEAR,		"CPILINEAR"			}, // LINEAR sur les CPI
	{	K_CPISTEPWISE,		"CPISTEPWISE"		}, // StepWise sur les CPI
	{	K_ZCLINEAR,			"ZCLINEAR"			}, // Linear sur les ZC
	{	K_ZCCTFWD,			"ZCCTFWD"			}, // Utilise Ct Fwd ZC
	{	K_CPISTEPWISESTART,	"CPISTEPWISESTART"	}, // StepWiseStart sur les CPI
	{	K_CPISTEPWISEEND,	"CPISTEPWISEEND"	}, // StepWiseEnd sur les CPI
	{	K_CPISTEPWISEMIDDLE,"CPISTEPWISEMIDDLE" }, // StepWiseMiddle sur les CPI


	// extrapolation method
	{	K_LASTTWO,			"LASTTWO"			}, // Use the last two points
	{	K_MIDDLE,			"MIDDLE"			}, // Use second second last and last
	{	K_FIRSTTWO,			"FIRSTTWO"			}, // Use second second last and second last
};


////////////////////////////////////////////////////
///	Class  : InfInterp
///	Routine: GetMappingName
///	Returns: char*
///	Action : general function to find the position of an index
////////////////////////////////////////////////////
char* InfInterp::GetMappingName( long methodFlag )
{
	int nbofRows = sizeof( MappingInterpTable ) / sizeof( MappingInterpTable[0] );

    for(int i=0; i<nbofRows; i++)
		if( MappingInterpTable[i].methodFlag == methodFlag )
			return MappingInterpTable[i].methodName;
	
	char buffer[20];
	sprintf( "%s", buffer, methodFlag );
    // if not found
    throw 	Exception(__LINE__, __FILE__, ERR_OBJECT_UNK,
		strcat("Unknown indexName ... unable to find method for ", buffer ) );
}



////////////////////////////////////////////////////
///	Class      : InfInterp
///	Description: table for all default values for inflation mkt data
/// see the header for meaning of the fields
////////////////////////////////////////////////////
Inf_Table_Row InfData::InfTable[] = 
{
	/// CPI index	Ccy		Cal		MonthInterp		DailyInterp	 DCFLag	ResetLag DCFMonthly		DCFDaily		ExtrapolType	PublishLag
	  { "CPALEMU",	"EUR", 	"EUR",	K_CPILINEAR,	K_CPILINEAR, "3m", "0m",	 KACTUAL_REAL,	KACTUAL_REAL,	K_LASTTWO,		"2m" },
	  { "EMU",		"EUR", 	"EUR",	K_CPILINEAR,	K_CPILINEAR, "3m", "0m",	 KACTUAL_REAL,	KACTUAL_REAL,	K_LASTTWO,		"2m" },
	  { "EMUT",		"EUR", 	"EUR",	K_CPILINEAR,	K_CPILINEAR, "3m", "0m",	 KACTUAL_REAL,	KACTUAL_REAL,	K_LASTTWO,		"2m" },	
	  { "IFRF",		"EUR",	"EUR",	K_CPILINEAR,	K_CPILINEAR, "3m", "0m",	 KACTUAL_REAL,	KACTUAL_REAL,	K_LASTTWO,		"2m" },	

	  // DEFAULT HAS to be the last row
	  { "DEFAULT",	"EUR",	"EUR",	K_CPILINEAR,	K_CPILINEAR, "3m", "0m",	 K30_360,		KACTUAL_ACTUAL,	K_LASTTWO,		"2m" },
};	

////////////////////////////////////////////////////
///	Class  : InfData
///	Routine: GetIndexPos
///	Returns: int
///	Action : general function to find the position of an index
////////////////////////////////////////////////////
int InfData::GetIndexPos( const char* dataName,
	const char* itsIndexName )
{
	int nbofRows = sizeof( InfTable ) / sizeof( InfTable[0] );
	int i;

    for(i=0; i<nbofRows; i++)
		if( strcmp( InfTable[i].ticker, itsIndexName ) == 0 )
			return i;
	
    /// if not found
	/// if flag _MUST_FIND_INDEX, throw an exception
	/// otherwise just take default
	#ifdef _MUST_FIND_INDEX
		char msg[100];
		sprintf( msg, "Unknown indexName %s .. unable to find its %s", dataName );
		throw 	Exception(__LINE__, __FILE__, ERR_OBJECT_UNK, msg );
	#endif
	
	/// DEFAULT is the last row
	return nbofRows-1;
}


////////////////////////////////////////////////////
///	Class  : InfData
///	Routine: GetCurrency
///	Returns: char* 
///	Action : Get the currency of the index
////////////////////////////////////////////////////
char* InfData::GetCurrency( 
	const char* itsIndexName )
{
	int pos = GetIndexPos( "currency", itsIndexName );
	return InfTable[pos].ccy;
}


////////////////////////////////////////////////////
///	Class  : InfData
///	Routine: GetCalendar
///	Returns: char* 
///	Action : Get the calendar of the index
////////////////////////////////////////////////////

char* InfData::GetCalendar(
	const char* itsIndexName )
{
	int pos = GetIndexPos( "calendar", itsIndexName );
	return InfTable[pos].calendar;
}


////////////////////////////////////////////////////
///	Class  : InfData
///	Routine: GetMonthInterpolation
///	Returns: long
///	Action : Get Default monthly interpolation
////////////////////////////////////////////////////

long InfData::GetMonthInterpolation(
	const char* itsIndexName )
{
	int pos = GetIndexPos( "month Interpolation", itsIndexName );
	return InfTable[pos].monthInterp;
}

////////////////////////////////////////////////////
///	Class  : InfData
///	Routine: GetDailyInterpolation
///	Returns: long
///	Action : Get Default daily interpolation
////////////////////////////////////////////////////

long InfData::GetDailyInterpolation(
	const char* itsIndexName )
{
	int pos = GetIndexPos( "dailyInterp", itsIndexName );
	return InfTable[pos].dailyInterp;
}

////////////////////////////////////////////////////
///	Class  : InfData
///	Routine: GetDCFLag
///	Returns: string
///	Action : Get Default DCF time gap
////////////////////////////////////////////////////

string InfData::GetDCFLag(
	const char* itsIndexName )
{
	int pos = GetIndexPos( "DCF Lag", itsIndexName );
	return InfTable[pos].DCFLag;
}


////////////////////////////////////////////////////
///	Class  : InfData
///	Routine: GetResetLag
///	Returns: string
///	Action : Get Default DCF time gap as a string
////////////////////////////////////////////////////

string InfData::GetResetLag(
	const char* itsIndexName )
{
	int pos = GetIndexPos( "Reset Lag", itsIndexName );
	return InfTable[pos].ResetLag;
}


////////////////////////////////////////////////////
///	Class  : InfData
///	Routine: GetDCFMonthly
///	Returns: long
///	Action : Get Default DCF Monthly
////////////////////////////////////////////////////
long InfData::GetDCFMonthly(
	const char* itsIndexName )
{
	int pos = GetIndexPos( "DCFMonth", itsIndexName );
	return InfTable[pos].DCFMonthly;
}


////////////////////////////////////////////////////
///	Class  : InfData
///	Routine: GetDCFDaily
///	Returns: long
///	Action : Get Default DCF Daily
////////////////////////////////////////////////////

long InfData::GetDCFDaily(
	const char* itsIndexName )
{
	int pos = GetIndexPos( "DCFDaily", itsIndexName );
	return InfTable[pos].DCFDaily;
}


////////////////////////////////////////////////////
///	Class  : InfData
///	Routine: GetExtrapolType
///	Returns: long
///	Action : Get Default ExtrapolType
////////////////////////////////////////////////////

long InfData::GetExtrapolType(
	const char* itsIndexName )
{
	int pos = GetIndexPos( "ExtrapolType", itsIndexName );
	return InfTable[pos].ExtrapolType;
}


////////////////////////////////////////////////////
///	Class  : InfData
///	Routine: GetPublishLag
///	Returns: string
///	Action : Get Default oublishing lag
////////////////////////////////////////////////////

string InfData::GetPublishLag(
	const char* itsIndexName )
{
	int pos = GetIndexPos( "PublishLag", itsIndexName );
	return InfTable[pos].PublishLag;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------*/
/*---- End Of File ----*/

