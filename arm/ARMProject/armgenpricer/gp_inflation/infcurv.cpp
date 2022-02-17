/*
 * $Log: infcurv.cpp,v $
 * Revision 1.29  2004/05/11 10:46:01  jpriaudel
 * function added by mingzhi
 *
 * Revision 1.28  2003/11/25 10:12:41  rguillemot
 * Inflation Season Manager
 *
 * Revision 1.27  2003/11/19 19:59:39  ebenhamou
 * use CC_DISTANCE
 *
 * Revision 1.26  2003/11/06 13:14:06  ebenhamou
 * change to have day display in date in view
 *
 * Revision 1.25  2003/10/27 07:17:16  ebenhamou
 *  change for .Net compatibility
 *
 * Revision 1.24  2003/10/22 10:02:16  ebenhamou
 * better check on pointor
 *
 * Revision 1.23  2003/09/29 09:26:11  ebenhamou
 * char date should be at least 20 to avoid crash
 *
 * Revision 1.21  2003/09/26 17:05:19  ebenhamou
 * version with namespace handled by macro in port.h
 *
 * Revision 1.20  2003/09/23 17:32:57  ebenhamou
 * accessor to MonthyInterpType
 *
 * Revision 1.19  2003/09/22 13:53:02  ebenhamou
 * more strict using directive
 *
 * Revision 1.18  2003/09/11 10:42:22  ebenhamou
 * GetCPIIndexValue
 *
 * Revision 1.17  2003/09/05 07:19:28  ebenhamou
 * added more constant
 *
 * Revision 1.16  2003/08/27 07:50:57  ebenhamou
 * using namespace
 *
 * Revision 1.15  2003/08/22 15:37:37  ebenhamou
 * added exception
 *
 * Revision 1.14  2003/08/21 08:44:05  ebenhamou
 * inherit from ARM_ZERO_CURVE
 *
 * Revision 1.13  2003/08/18 11:06:54  ebenhamou
 * remove namespace as not supported by unix yet
 *
 * Revision 1.12  2003/08/15 07:08:15  ebenhamou
 * use CountYearsWithoutException for extrapolation purpose
 *
 * Revision 1.11  2003/08/12 08:11:27  ebenhamou
 * cutoff date is infidx ref date
 *
 * Revision 1.10  2003/08/05 08:32:51  ebenhamou
 *  keep inf curve derived from arm_object ... use ARM_Currency
 *
 * Revision 1.7  2003/07/18 08:59:40  ebenhamou
 * more explicit terminology
 *
 * Revision 1.5  2003/07/16 07:00:56  ebenhamou
 * version with monthly and daily interp
 *
 * Revision 1.4  2003/06/30 17:14:12  ebenhamou
 * dos2unix
 *
 * Revision 1.3  2003/06/30 16:11:27  ebenhamou
 * discounting version
 *
 * Revision 1.2  2003/06/23 12:39:24  ebenhamou
 * dos2unix
 *
 * Revision 1.1  2003/06/20 17:15:57  ebenhamou
 * Initial revision
 *
 *
 */

#include <glob/firsttoinc.h>
#include "gpinflation/infcurv.h"
#include "gpinflation/infdata.h"
#include "gpinflation/resetmanager.h"
#include "gpinflation/seasonmanager.h"
#include "gpinflation/infidx.h"

/// standard libs
#include <utility>
#include <stdarg.h>
#include <algorithm>
#include <functional>

///	flag for interp
#include <ccy/currency.h>


CC_USING_NS( std, string )
CC_USING_NS( std, vector )
//CC_USING_NS_T( std, vector, string )

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class   : ARM_InfCurv
///	Constant: minCPIVal
/// to make the difference between a CPI value
/// we assume that any CPI Value is always above MINCPIVal
////////////////////////////////////////////////////
const double ARM_InfCurv::minCPIVal	=	90.0;


////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: ConvertString2Double
///	Returns: double
///	Action : Conversion routine for dates
///				returns the julian date for a string 
///				corresponding to either 
////////////////////////////////////////////////////

double ARM_InfCurv::ConvertString2Double( const string& s, const char* calendar ) const
{
	const string maturity( "yYMmWwDd" );
	const double JULIANDATEMIN = 2000000.0;
	const double JULIANDATEADD = 2415019.0;

	/// two cases
	/// either it is already a Julian date or or an excel date
	/// and we just take it
	if( s.find_first_of( maturity ) == string::npos )
	{
		double datedble=  atof(s.c_str() );
		if( datedble < 0.0 )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid Date");

		if( datedble < JULIANDATEMIN )
			datedble += JULIANDATEADD;
		return datedble;
	}

	/// or we have a maturity
	/// and we add to itsCPIIndexDate which is stored as a double
	else
	{	
		/// since the function AddPeriod is not const
		/// we are forced to use a tempDate

		/// by definition, we add to the CPIrefDate
		ARM_Date tempDate = itsCPIIndexDate;
		tempDate.AddPeriod( s, calendar );
		return tempDate.GetJulian();
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: ConvertZCRate2CPI
///	Returns: double
///	Action : Conversion of a Zero coupon Rate to corresponding year term date
////////////////////////////////////////////////////

double ARM_InfCurv::ConvertZCRate2CPI( double ZCRate, double expiry ) const
{
	return( pow(1.0+ZCRate / CC_NS( ARM_Constants, rateBase ), expiry ) * itsCPIIndexValue );
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: GetCPIIndexValue
///	Returns: double
///	Action : get the cpi index value (accessor)
////////////////////////////////////////////////////

double ARM_InfCurv::GetCPIIndexValue() const
{
	return itsCPIIndexValue;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: ConvertCPI2ZCRate
///	Returns: double
///	Action : Conversion of a CPI Rate to  a ZC Rate
////////////////////////////////////////////////////

double ARM_InfCurv::ConvertCPI2ZCRate( double CPIVal, double expiry) const
{
	if( expiry < 0 )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid Expiry");
	else 
	{
		if( expiry == 0 )
			return -111111;
		else
			return( (pow(CPIVal/itsCPIIndexValue, 1.0 / expiry ) - 1.0 ) * CC_NS( ARM_Constants, rateBase ) );
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: ARM_InfCurv
///	Returns: 
///	Action : Constructor
///			 Tailor made Constructor
///			 in the sense that it transforms inputs
///			 find expiry terms from string and convert
///			 ZC rate to CPI Value
////////////////////////////////////////////////////

ARM_InfCurv::ARM_InfCurv(
	const ARM_Date& asOfDate,	
	const string& indexName,
	double CPIIndexValue,
	const ARM_Date& CPIIndexDate,
	const vector<string>& MktTerms,
	const vector<double>& MktValues,

	long MonthlyInterpType,
	long DailyInterpType,
	long DCFMonthly,
	long DCFDaily,
	long ExtrapolType,

	ARM_ResetManager* ResetManager,
	ARM_SeasonalityManager* SeasonalityManager)
: 
	itsAsOfDate(asOfDate),
	itsIndexName( indexName ), 
	itsCPIIndexValue( CPIIndexValue ), 
	itsCPIIndexDate( CPIIndexDate ),
	itsMonthlyInterpType( MonthlyInterpType == -1? 
		InfData::GetMonthInterpolation( itsIndexName.c_str() ) : MonthlyInterpType ), 
	itsDailyInterpType( DailyInterpType == -1 ?
		InfData::GetDailyInterpolation( itsIndexName.c_str() ) : DailyInterpType  ),
	itsDCFMonthly( DCFMonthly == -1?
		InfData::GetDCFMonthly( itsIndexName.c_str() ): DCFMonthly  ),
	itsDCFDaily( DCFDaily== -1?
		InfData::GetDCFDaily( itsIndexName.c_str() ): DCFDaily  ),
	itsExtrapolType( ExtrapolType == -1?
		InfData::GetExtrapolType( itsIndexName.c_str() ): ExtrapolType ),

	itsDCFLag( InfData::GetDCFLag( itsIndexName.c_str() ) ),
	itsResetLag( InfData::GetResetLag( itsIndexName.c_str() ) ),
	itsCalendar( InfData::GetCalendar( itsIndexName.c_str() ) ),
	itsMktTerms( MktTerms ),
	itsMktValues( MktValues ),
	itsResetManager( NULL ),
	itsSeasonalityManager( NULL ),
	itsCurrency( NULL ),
	itsMonthlyCPIValue()
{
	SetName(ARM_INFCURV);
	BuildCurve( itsMktTerms, itsMktValues );

	/// to avoid memory leak due to exception 
	/// does all the new Vector after building the curve
	itsResetManager = ResetManager == NULL ? NULL: (ARM_ResetManager*) ResetManager->Clone();
	itsSeasonalityManager = SeasonalityManager == NULL ? NULL: (ARM_SeasonalityManager*) SeasonalityManager->Clone();

	/// SetCurrency
	itsCurrency = new ARM_Currency( InfData::GetCurrency( itsIndexName.c_str() ) ); 
}



////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: ARM_InfCurv
///	Returns: 
///	Action : Default constructor
////////////////////////////////////////////////////
ARM_InfCurv::ARM_InfCurv()
:
	itsResetManager( NULL ),
	itsSeasonalityManager( NULL ),
	itsCurrency( NULL )
{
	SetName(ARM_INFCURV);
}




////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: BuildCurve
///	Returns: 
///	Action : Set to build the curve from data
///			 function to rebuild the curve
////////////////////////////////////////////////////

void ARM_InfCurv::BuildCurve(
	const vector<string>& MktTerms,
	const vector<double>& MktValues )
{
	/// convert string to Julian Date
	/// and compute the ZCRate and CPI
	vector<double>::const_iterator iterCPIInput	= MktValues.begin();
	vector<string>::const_iterator iterExpInput	= MktTerms.begin();

	/// initialise the vector
	itsExpiryTermsVec	= vector<double>( MktTerms.size() );
	itsCPIValuesVec		= vector<double>( MktValues.size() );
	itsZCValuesVec		= vector<double>( MktValues.size() );

	vector<double>::iterator iterExpiry	= itsExpiryTermsVec.begin();
	vector<double>::iterator iterCPI	= itsCPIValuesVec.begin();
	vector<double>::iterator iterZC		= itsZCValuesVec.begin();
	
	/// check wheter it is a CPI rate or a ZC Rate
	/// based on the fact that the input is lower than MinCPIVal
	/// and fill the appropriate vector
	for( iterCPIInput; iterCPIInput< MktValues.end(); ++iterCPIInput, ++iterExpInput, 
			++iterExpiry, ++iterCPI, ++iterZC  )
	{
		*iterExpiry = ConvertString2Double( *iterExpInput, itsCalendar.c_str() );
		double expiry = CountYearsWithoutException( itsDCFMonthly, itsCPIIndexDate.GetJulian(), *iterExpiry );

		if( *iterCPIInput < minCPIVal )
		{
			*iterCPI = ConvertZCRate2CPI( *iterCPIInput, expiry );
			*iterZC  = *iterCPIInput;
		}
		else
		{
			*iterCPI = *iterCPIInput;
			*iterZC  = ConvertCPI2ZCRate( *iterCPIInput,expiry );
		}
	}
	
	/// test that expiries are sorted
	vector<double> sortedVec = itsExpiryTermsVec;

	/// we are forced to use the std prefix to avoid ambiguity!!
	/// do not remove it
	std::sort( sortedVec.begin(), sortedVec.end(), less<double>() );

	if( sortedVec != itsExpiryTermsVec )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Expiries are not sorted");

	itsCutOffDate = ARM_Date( *itsExpiryTermsVec.begin() );

	SetBucketEndPeriod(*(iterExpiry-1));  /// ajoute par Mingzhi
}



////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : code to factorise the copy of all members data
////////////////////////////////////////////////////

void ARM_InfCurv::CopyNoCleanUp( const ARM_InfCurv& srcInfCurv )
{
	itsExpiryTermsVec	= srcInfCurv.itsExpiryTermsVec;
	itsCPIValuesVec		= srcInfCurv.itsCPIValuesVec;
	itsZCValuesVec		= srcInfCurv.itsZCValuesVec;

	itsIndexName		= srcInfCurv.itsIndexName;
	itsAsOfDate			= srcInfCurv.itsAsOfDate;
	itsCutOffDate		= srcInfCurv.itsCutOffDate;
	itsCPIIndexValue	= srcInfCurv.itsCPIIndexValue;
	itsCPIIndexDate		= srcInfCurv.itsCPIIndexDate;

	itsMonthlyInterpType= srcInfCurv.itsMonthlyInterpType;
	itsExtrapolType		= srcInfCurv.itsExtrapolType;
	itsDailyInterpType	= srcInfCurv.itsDailyInterpType;
	itsMonthlyCPIValue	= srcInfCurv.itsMonthlyCPIValue;
	itsDCFMonthly		= srcInfCurv.itsDCFMonthly;
	itsDCFDaily			= srcInfCurv.itsDCFDaily;

	itsDCFLag			= srcInfCurv.itsDCFLag;
	itsResetLag			= srcInfCurv.itsResetLag;
	itsCalendar			= srcInfCurv.itsCalendar;
	itsMktTerms			= srcInfCurv.itsMktTerms;	
	itsMktValues		= srcInfCurv.itsMktValues;

	itsBucketStartPeriod= srcInfCurv.itsBucketStartPeriod;
	itsBucketEndPeriod  = srcInfCurv.itsBucketEndPeriod;

	itsCurrency			= srcInfCurv.itsCurrency ?
		static_cast<ARM_Currency*>( srcInfCurv.itsCurrency->Clone() )
	:	NULL;

	/// copy of the reset manager
	itsResetManager	= srcInfCurv.itsResetManager ?
		static_cast<ARM_ResetManager*>( srcInfCurv.itsResetManager->Clone() )
	:	NULL;

	itsSeasonalityManager	= srcInfCurv.itsSeasonalityManager ?
		static_cast<ARM_SeasonalityManager*>( srcInfCurv.itsSeasonalityManager->Clone() )
	:	NULL;


}


////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: CleanUp
///	Returns: 
///	Action : code to factorise the copy of all members data
////////////////////////////////////////////////////

void ARM_InfCurv::CleanUp()
{
	delete itsCurrency;
	itsCurrency = NULL;
	delete itsResetManager;
	itsResetManager = NULL;
	delete itsSeasonalityManager;
	itsResetManager = NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////

ARM_InfCurv::~ARM_InfCurv()
{
	CleanUp();
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: Copy constructor
///	Returns: 
///	Action : Copy constructor
////////////////////////////////////////////////////

ARM_InfCurv::ARM_InfCurv (const ARM_InfCurv& srcInfCurv) 
:	ARM_Object(srcInfCurv)
{
	CopyNoCleanUp( srcInfCurv );
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: Operator = 
///	Returns: 
///	Action : Operator = 
////////////////////////////////////////////////////

ARM_InfCurv & ARM_InfCurv::operator = (const ARM_InfCurv& srcInfCurv)
{
	if( this != &srcInfCurv )
	{
		(*this).ARM_Object::operator = (srcInfCurv);
		// assign to all data members
		CleanUp();
		CopyNoCleanUp( srcInfCurv );
	}
	return *this;
}



////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: Clone
///	Returns: 
///	Action : Clone
////////////////////////////////////////////////////
ARM_Object* ARM_InfCurv::Clone(void)
{
	return new ARM_InfCurv(*this);
}



////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: CPIOneDateInterpolateAndStore
///	Returns: 
///	Action : CPIOneDateInterpolateAndStore
///		for the monthly points
///		interpolate and store the point
///		with the following possibilities
///		there is no check that the date is a monthly
///		date for performance reason
///		there is no interpFlag as the interp method 
///		is stored in the curve
///		
///		- K_CPILINEAR:	linear on the CPI
///		- K_ZCLINEAR:	linear on the Zero Coupon rate
///		- K_ZCCTFWD:		assuming a constant forward ZC rate
///		
///		Be aware that this function returns only CPI
///		hence ZC Rates have to be converted!
////////////////////////////////////////////////////

double ARM_InfCurv::CPIOneDateInterpolateAndStore( const ARM_Date& Date,
	long MonthlyInterpType )
{
	double expiry = Date.GetJulian();

	/// take the default?
	if( MonthlyInterpType == -1 )
		MonthlyInterpType = itsMonthlyInterpType;

	/// already computed?
	doubleDoubleMap::iterator iter = itsMonthlyCPIValue.find( expiry );
	if( iter != itsMonthlyCPIValue.end() )
		return  (*iter).second;

	/// Otherwise we need to compute the value
	// find the iterator after Date
	vector<double>::const_iterator iterExpiryBegin = itsExpiryTermsVec.begin();
	vector<double>::const_iterator iterExpiry = std::find_if(
		iterExpiryBegin, 
		( vector<double>::const_iterator ) itsExpiryTermsVec.end(), 
		CC_NS( std, bind2nd )( CC_NS( std, greater)<double>(), expiry) );

	/// find the same position in the CPI vector or ZC according to the 
	/// Monthly interpolation method
	vector<double>::const_iterator iterCPI;


	/// linear interpolation on the CPI value, so take CPI values
	switch( MonthlyInterpType )
	{
		/// use the CPI values
		case K_CPILINEAR: case K_ZCCTFWD:
			iterCPI = itsCPIValuesVec.begin() + CC_NS( std, CC_DISTANCE)( iterExpiryBegin, iterExpiry );
			break;
		/// use ZC Rate values
		case K_ZCLINEAR:  
			iterCPI = itsZCValuesVec.begin() + CC_NS( std, CC_DISTANCE)( iterExpiryBegin, iterExpiry );
			break;

		default:
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"Interpolation method not supported");
	}

	double ValueAfter;	/// next point X value
	double ExpiryAfter; /// next point Y value
	double ValueBefore;	/// previous point X value
	double ExpiryBefore;/// previous point Y Value

	/// special cases
	/// after the last point
	if( iterCPI == itsCPIValuesVec.end() )
	{
		/// extrapolation case
		/// either the last two points
		/// or the last one and the last last last one
		/// or the last last one and the last last last last one
		switch( itsExtrapolType )
		{
			case K_LASTTWO:
			{
				ValueAfter	= *(--iterCPI);
				ExpiryAfter = *(--iterExpiry);
				ValueBefore = *(--iterCPI);
				ExpiryBefore= *(--iterExpiry);
				break;
			}
			case K_MIDDLE:
			{
				ValueAfter	= *(--iterCPI);
				ExpiryAfter = *(--iterExpiry);
				--iterCPI;
				--iterExpiry;
				ValueBefore = *(--iterCPI);
				ExpiryBefore= *(--iterExpiry);
				break;
			}
			case K_FIRSTTWO:
			{
				--iterCPI;
				--iterExpiry;
				ValueAfter	= *(--iterCPI);
				ExpiryAfter = *(--iterExpiry);
				ValueBefore = *(--iterCPI);
				ExpiryBefore= *(--iterExpiry);
				break;
			}
			default:
		      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                         "Extrapolation method not supported");
		}
	}
	// standard case
	else
	{
		ValueAfter	= *iterCPI;
		ExpiryAfter = *iterExpiry;
		ValueBefore	= *(--iterCPI);
		ExpiryBefore= *(--iterExpiry);
	}

	double result;

	switch( MonthlyInterpType )
	{
		/// linear interpolation of the CPI
		case K_CPILINEAR: 
		{	
			/// the weight is given by the nb of days between ExpiryBefore and expiry
			/// over the nb of days between ExpiryBefore and ExpiryAfter
			double weight	= CountYearsWithoutException( itsDCFMonthly, expiry, ExpiryAfter ) 
				/ CountYearsWithoutException( itsDCFMonthly, ExpiryBefore, ExpiryAfter );
			result = weight*ValueBefore+(1.0- weight)*ValueAfter;
			break;
		}

		/// linear on ZC rate
		case K_ZCLINEAR:
		{	
			double weight	= CountYearsWithoutException( itsDCFMonthly, expiry, ExpiryAfter )
				/ CountYearsWithoutException( itsDCFMonthly, ExpiryBefore, ExpiryAfter );
			double ZCResult = weight*ValueBefore+(1.0- weight)*ValueAfter;

			double T = CountYearsWithoutException( itsDCFMonthly, itsCPIIndexDate.GetJulian(), expiry );
			result = ConvertZCRate2CPI( ZCResult, T );
			break;
		}

		/// in this particular case, we compute the constant forward
		/// and assume that it remains the same
		/// (1+X2)^T2 / (1+X1)^T1 = (1+CtFwd)^(T2-T1);
		/// returns (1+X1)^T1 * (1+CtFwd)^(T-T1)

		case K_ZCCTFWD:
		{
			double T1		= CountYearsWithoutException( itsDCFMonthly, itsCPIIndexDate.GetJulian(), ExpiryBefore );
			double T2		= CountYearsWithoutException( itsDCFMonthly, itsCPIIndexDate.GetJulian(), ExpiryAfter );
			double T		= CountYearsWithoutException( itsDCFMonthly, itsCPIIndexDate.GetJulian(), expiry );

			double CtFwd	= pow( ValueAfter/ ValueBefore, 1.0/(T2-T1) ) - 1.0;
			result = ValueBefore * pow( 1.0+CtFwd, T-T1);
			break;
		}

		/// other cases are error
		default:
	      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                         "Interpolation method not supported");
	}

	// If there is a seasonality manager we add a season adjustment
	if (itsSeasonalityManager)
	{
		itsSeasonalityManager->SeasonSpreadCorrect( itsAsOfDate, result, expiry, ExpiryBefore, ExpiryAfter);
	}
	
	/// store the new value
	itsMonthlyCPIValue.insert( pair<const double,double>( expiry, result) );
	return result;
}




////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: GetMonthBoundingDates
///	Returns: pair<ARM_Date,ARM_Date> 
///	Action : Get the value of the CPI for the beginning of the month
////////////////////////////////////////////////////

pair<ARM_Date,ARM_Date> ARM_InfCurv::GetMonthBoundingDates( const ARM_Date& CPIDate ) const
{
	pair<double,double> result;

	// get the two 1st months between the CPIDate
	ARM_Date CurrentMonth = CPIDate;
	int nbDays	= CurrentMonth.GetDay()-1 ;
	CurrentMonth.AddDays( -nbDays );

	/// get the next month
	ARM_Date NextMonth = CurrentMonth;
	NextMonth.AddMonths( 1 );

	return pair<ARM_Date,ARM_Date>( CurrentMonth, NextMonth );
};



////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: CPIInterpolate
///	Returns: double 
///	Action : 
///  The curve defines a way to interpolate the CPI each month
///  based on the method given in the curve
///  then we do an interpolation between 2 monthly points with
///  the following possibilities
/// 
///  -StepwiseStart		returns value at the start of the period
///  -StepwiseMid			returns mid of start and end value 
/// 						of the period
///  -Linear				returns linear interpolation based on the other
/// 						expiry of the CPI value
/// 
///  in the case of linear interpolation, because the 
///  weight can be computed on different dates as 
///  the one of the reset we take as an argument 
///  - the CPIDate corresponding to the reset
///  - the DCFDate corresponding to the computation of the weight
////////////////////////////////////////////////////

double ARM_InfCurv::CPIInterpolate(
		const ARM_Date& CPIDate,
		const ARM_Date& DCFDate,
		long DailyInterpType,
		double weight )
{
	/// use the default dailyInterp?
	if( DailyInterpType == -1 )
		DailyInterpType = itsDailyInterpType;

	/// historical reset?
	if( CPIDate < itsCutOffDate )
	{
		if( !itsResetManager )
		{
			char strDate[20];
			char strDate2[20];
			ARM_Date tempDate = CPIDate;
			tempDate.JulianToStrDateDay( strDate );
			tempDate = itsCutOffDate;
			tempDate.JulianToStrDateDay( strDate2 );
			char msg[80];
			sprintf( msg, "date %s before the first date %s with no resetmanager ... You should probably use a reset manager", strDate, strDate2 );
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
		}
		else
		{
			pair<ARM_Date,ARM_Date> BoundingMonthDate = GetMonthBoundingDates( CPIDate );
			ARM_Date MonthStart		= BoundingMonthDate.first;
			ARM_Date MonthNext		= BoundingMonthDate.second;
			
			const string MKT = "INF";
			const double monthStartCPI = itsResetManager->GetReset( MonthStart.GetJulian(), itsIndexName, MKT );

			
			// first checks that we have the correct date
			if( CPIDate == MonthStart )
				return monthStartCPI;

			/// second if it is CPILINEAR allow to compute the linear interpolation
			if( DailyInterpType == K_CPILINEAR )
			{
				const double monthNextCPI = itsResetManager->GetReset( MonthNext.GetJulian(), itsIndexName, MKT );
				
				pair<ARM_Date,ARM_Date> BoundingDCFMonthDate = GetMonthBoundingDates( DCFDate );
				
				double linearWeight = CountYearsWithoutException( itsDCFDaily, DCFDate, BoundingDCFMonthDate .second )
					/ CountYearsWithoutException( itsDCFDaily, BoundingDCFMonthDate .first, BoundingDCFMonthDate .second );

				return	linearWeight*monthStartCPI+(1.0-linearWeight)*monthNextCPI;
			}
			else
			{
				/// other cases are exception
				char msg[255];
				char dateChar[20];
				ARM_Date tempDate = CPIDate;
				tempDate.JulianToStrDateDay( dateChar );
				sprintf(msg, "could not find the date %s for the index %s in the reset data", 
					dateChar, itsIndexName.c_str() );
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					msg );
			}
		}
	}

	pair<ARM_Date,ARM_Date> BoundingMonthDate = GetMonthBoundingDates( CPIDate );
	ARM_Date MonthStart		= BoundingMonthDate.first;
	ARM_Date MonthNext		= BoundingMonthDate.second;

	/// we want to be sure to have the exact point in the curve
	if( CPIDate == MonthStart )
		return CPIOneDateInterpolateAndStore( MonthStart );


	double CPIValue = 0.0;

	switch( DailyInterpType )
	{
		case K_CPISTEPWISESTART:
			CPIValue = CPIOneDateInterpolateAndStore( MonthStart );
			break;
		case K_CPISTEPWISEEND:
			CPIValue = CPIOneDateInterpolateAndStore( MonthNext );
			break;
		case K_CPISTEPWISEMIDDLE:
			CPIValue = 0.5 * CPIOneDateInterpolateAndStore( MonthStart )
				   +0.5 * CPIOneDateInterpolateAndStore( MonthNext );
			break;
		case K_CPISTEPWISE:
			CPIValue = (1.0-weight)* CPIOneDateInterpolateAndStore( MonthStart )
				   + weight * CPIOneDateInterpolateAndStore( MonthNext );
			break;
		case K_CPILINEAR:
		{
			pair<ARM_Date,ARM_Date> BoundingMonthDate = GetMonthBoundingDates( DCFDate );
			
			double linearWeight = CountYearsWithoutException( itsDCFDaily, DCFDate, BoundingMonthDate.second )
				/ CountYearsWithoutException( itsDCFDaily, BoundingMonthDate.first, BoundingMonthDate.second );

			CPIValue = linearWeight * CPIOneDateInterpolateAndStore( MonthStart )
				   +(1.0-linearWeight) * CPIOneDateInterpolateAndStore( MonthNext );
			break;
		}
		case K_ZCLINEAR:
		{
			/// Computation of the weight
			pair<ARM_Date,ARM_Date> BoundingMonthDate = GetMonthBoundingDates( DCFDate );

			double linearWeight = CountYearsWithoutException( itsDCFDaily, DCFDate, BoundingMonthDate.second )
				/ CountYearsWithoutException( itsDCFDaily, BoundingMonthDate.first, BoundingMonthDate.second );

			/// computation of the ZC Rate
			double StartExpiry = CountYearsWithoutException( itsDCFMonthly, itsCPIIndexDate, MonthStart );
			double ZCStart = ConvertCPI2ZCRate( CPIOneDateInterpolateAndStore( MonthStart ),
				StartExpiry );

			double NextExpiry = CountYearsWithoutException( itsDCFMonthly, itsCPIIndexDate, MonthNext );
			double ZCNext = ConvertCPI2ZCRate( CPIOneDateInterpolateAndStore( MonthNext ),
				NextExpiry );

			double expiry =  CountYearsWithoutException( itsDCFMonthly, itsCPIIndexDate, CPIDate );
			CPIValue = ConvertZCRate2CPI( linearWeight*ZCStart+(1.0-linearWeight)*ZCNext,
				expiry );
			break;
		}

		case K_ZCCTFWD:
		{
			double ValueBefore = CPIOneDateInterpolateAndStore( MonthStart );
			double T1 = CountYearsWithoutException( itsDCFMonthly, itsCPIIndexDate, MonthStart );

			double ValueAfter = CPIOneDateInterpolateAndStore( MonthNext );
			double T2 = CountYearsWithoutException( itsDCFMonthly, itsCPIIndexDate, MonthNext );

			double CtFwd = pow( ValueAfter/ ValueBefore, 1.0/(T2-T1))-1.0;

			double T = CountYearsWithoutException( itsDCFMonthly, itsCPIIndexDate, CPIDate );
			CPIValue = ValueBefore * pow( 1.0+CtFwd, T-T1 );
			break;
		}
		default:
	      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                         "Interpolation method not supported");
	}

	return CPIValue;
}



////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: CPIInterpolate
///	Returns: double 
///	Action :  Version for the interface
///			with dates and lag as string
////////////////////////////////////////////////////

double ARM_InfCurv::CPIInterpolate(	const ARM_Date& date, 
	const string& DCFLag,
	long DailyInterpType,
	const string& ResetLag, 
	double weight )
{
	// use default DCFLag ?
	string DCFLagNew = DCFLag == "-1" ? itsDCFLag : DCFLag;
	/// use default ResetLag ?
	string ResetLagNew = ResetLag == "-1" ? itsResetLag : ResetLag;

	ARM_Date CPIDate = date;
	ARM_Date DCFDate = date;
	CPIDate.AddPeriod( ResetLagNew, itsCalendar.c_str() );
	DCFDate.AddPeriod( DCFLagNew, itsCalendar.c_str() );
	return CPIInterpolate( CPIDate, DCFDate, DailyInterpType, weight );
}



////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: ZCRateInterpolate
///	Returns: double 
///	Action :  Version on ZCRate
///				with dates and lag as string
///				the choice is to compute everything in terms of 
///				CPI and translate as ZC Rate
////////////////////////////////////////////////////

double ARM_InfCurv::ZCRateInterpolate( const ARM_Date& date, 
	const string& DCFLag,
	long DailyInterpType,
	const string& ResetLag, 
	double weight )
{
	double CPI = CPIInterpolate( date, DCFLag, DailyInterpType, ResetLag, weight );

	/// use default ResetLag ?
	string ResetLagNew = ResetLag == "-1" ? itsResetLag : ResetLag;

	ARM_Date CPIDate = date;
	CPIDate.AddPeriod( ResetLagNew, itsCalendar.c_str() );

	pair<ARM_Date,ARM_Date> BoundingMonthDate = GetMonthBoundingDates( CPIDate );
	ARM_Date MonthStart	= BoundingMonthDate.first;

	double expiry = CountYearsWithoutException( itsDCFMonthly, itsCPIIndexDate, MonthStart )
		+ CountYearsWithoutException( itsDCFDaily, MonthStart, date );

	return ConvertCPI2ZCRate( CPI, expiry );
}



////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: View
///	Returns: void
///	Action : View the details of the curve
////////////////////////////////////////////////////

void ARM_InfCurv::View(char* id, FILE* ficOut)
{
	/// to be consistent with other interface
	/// need to use fprintf
    FILE* fOut;
    char fOutName[40];
    char strDate[20];
	
	/// do we have already a file opened?
    if ( ficOut == NULL )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;
	
	/// printing of the asofDate 
    fprintf(fOut, "\n\n\t =====> Inflation Curve \n\n");
	ARM_Date tempDate = itsAsOfDate;
    tempDate.JulianToStrDateDay(strDate);

	/// AsOfDate
	fprintf(fOut, "\n\n\t AsOfDate  : %s \n\n", strDate);
	/// Index
	fprintf(fOut, "\t Index     : %s\n", itsIndexName.c_str() );
	
	/// index Date and Ref
	tempDate = itsCPIIndexDate;
	tempDate.JulianToStrDateDay( strDate );
	fprintf(fOut, "\n\t Index Date: %s \t Index Ref : %.3lf\n", strDate, itsCPIIndexValue);

	/// currency and calendar
	fprintf(fOut, "\t Currency  : %s \t\t Calendar  : %s\n", itsCurrency->GetCcyName(), itsCalendar.c_str() );
	
	/// Interpolation spec
	fprintf(fOut, "\t Monthly   : %s \t\t Daily     : %s\n",
		InfInterp::GetMappingName( itsMonthlyInterpType ),
		InfInterp::GetMappingName( itsDailyInterpType ) );

	// Day count fraction
	fprintf(fOut, "\t DCFMonhtly: %s \t DCFDaily  : %s\n",
		DCGetName( itsDCFMonthly ),
		DCGetName( itsDCFDaily) );

	/// Lag CPI and DCF
	fprintf(fOut, "\t DCFLag    : %s \t\t ResetLag  : %s\n", itsDCFLag.c_str(), itsResetLag.c_str() );

	/// extrapolation
	fprintf(fOut, "\t Extrapol  : %s\n",
		InfInterp::GetMappingName( itsExtrapolType ) );

	/// points of the curve
	/// CPI date, Nb of days CPI rate and ZC Rate
	int sz = itsExpiryTermsVec.size();
	fprintf(fOut, "\n\n Date \t\t\t YearFrac \t CPI \t\t\t ZCRate \n\n");
	int  i;

	for( i = 0; i < sz; ++i )
	{
		tempDate = ARM_Date( itsExpiryTermsVec[i] );
		double yearFrac = CountYearsWithoutException( itsDCFMonthly, itsCPIIndexDate, tempDate );
		tempDate.JulianToStrDateDay(strDate);
		fprintf(fOut, " %s \t %.4lf \t %.3lf \t\t %.3lf \n",
			strDate, yearFrac, itsCPIValuesVec[i], itsZCValuesVec[i] );
	}			


	if( itsResetManager )
	{
		fprintf( fOut, "\n\n\n\n" );
		itsResetManager->View( id,fOut);
	}

	if( itsSeasonalityManager )
	{
		fprintf( fOut, "\n\n\n\n" );
		itsSeasonalityManager->View( id,fOut);
	}
	
	if ( ficOut == NULL )
		fclose(fOut);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCurv
///	Routine: GenerateShiftCurve
///	Returns: ARM_InfCurv*
///	Action : Function to Generate a shiftedCurve
////////////////////////////////////////////////////

ARM_InfCurv* ARM_InfCurv::GenerateShiftCurve(
	const vector<string>& mktTerms, 
    const vector<double>& epsilon )
{
	ARM_InfCurv* newCrv = NULL;

	vector<string> newMktTerms = itsMktTerms;
	vector<double> newMktValues = itsMktValues;

	int szinp = epsilon.size();

	int k=0;
	int l=0;

	double tempCPIIndexValue = GetCPIIndexValue();
	string sCPIIndexDate;
	char buffer[20];
	sprintf( buffer, "%f", GetCPIIndexDate().GetJulian() );
	sCPIIndexDate = buffer;

	for (k=0; k<newMktTerms.size(); k++)
	{
		for (l=0; l<szinp; l++)
		{
			if (mktTerms[l] == newMktTerms[k])
			{
				newMktValues[k] += epsilon[l];

				// Si c'est le CPI de référence, on le bumpe aussi
				if (mktTerms[l] == sCPIIndexDate)
					tempCPIIndexValue += epsilon[l];
				
				l=szinp;
			}
		}
	}

	newCrv = new ARM_InfCurv(GetAsOf(),
							 GetInfIdxName(),
							 tempCPIIndexValue,
							 GetCPIIndexDate(),
							 newMktTerms,
							 newMktValues,
							 GetMonthlyInterpType(),
							 GetDailyInterpType(),
							 GetDCFMonthly(),
							 GetDCFDaily(),
							 GetExtrapolType(),
							 GetResetManager(),
							 GetSeasonalityManager());

	return newCrv;
}

double ARM_InfCurv::DiscountFunction(double yearTerm )
{
	ARM_Date tmpDate = itsCPIIndexDate;
	tmpDate.AddYears( yearTerm );
	return itsCPIIndexValue / CPIInterpolate( tmpDate, tmpDate );
}

/// returns derivatives of discount price
double ARM_InfCurv::D1DiscountFunction(double yearTerm)
{
	const double epsilon = 0.00001;
	return ( DiscountFunction(yearTerm+epsilon)-DiscountFunction(yearTerm) )/ epsilon;
}

/// accessor to the monthly interpolation type
long ARM_InfCurv::GetMonthlyInterpType() const
{
	return itsMonthlyInterpType;
}

/// accessor to the daily interpolation type
long ARM_InfCurv::GetDailyInterpType() const
{
	return itsDailyInterpType;
}

long ARM_InfCurv::GetDCFMonthly() const
{
	return itsDCFMonthly;
}

long ARM_InfCurv::GetDCFDaily() const
{
	return itsDCFDaily;
}

long ARM_InfCurv::GetExtrapolType() const
{
	return itsExtrapolType;
}

/// set accessor for the reset manager
void ARM_InfCurv::SetResetManager(ARM_ResetManager* ResetManager)
{
	delete itsResetManager;
	itsResetManager = static_cast<ARM_ResetManager*>(ResetManager->Clone());
}


/// set accessor for the Seasonality manager
void ARM_InfCurv::SetSeasonalityManager(ARM_SeasonalityManager* SeasonalityManager)
{
	delete itsSeasonalityManager;
	itsSeasonalityManager = static_cast<ARM_SeasonalityManager*>(SeasonalityManager->Clone());
}


/// create the corresponding inflation index
ARM_InfIdx* ARM_InfCurv::GetInfIdx() const
{
	return new ARM_InfIdx( itsIndexName );
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------*/
/*---- End Of File ----*/

