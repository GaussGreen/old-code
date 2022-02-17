/*! \file ARM_local_infcurv.cpp
 *
 *  \brief file for the local interface of inflation objects
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date June 2003
 */

#include "firstToBeIncluded.h"
#include "ARM_local_gp_inflation.h"
#include "ARM_local_class.h"
#include "ARM_local_persistent.h"
#include "ARM_result.h"
#include "ARM_local_glob.h"
#include "ARM_local_gp_genericaddin.h"
#include "ARM_local_wrapper.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <ARM\local_xlarm\XL_local_xlarm_common.h>
#include "ARM_local_etoolkit.h"
#include "ARM_local_parsexml.h"

#include "glob\dates.h"
#include "ccy\currency.h"
#include "crv\volcube.h"
#include "crv\correlmanager.h"
#include "util\fromto.h"

#include <gpinfra\argconvdefault.h>

#include "GP_Base\gpbase\typedef.h"
#include "GP_Base\gpbase\\argconvdefault.h"
#include "GP_Base\gpbase\datestrip.h"
#include "GP_Base\gpbase\gplinalgtypedef.h"
#include "GP_Base\gpbase\gplinalgconvert.h"
#include "GP_Base\gpbase\gpvector.h"
#include "GP_Base\gpbase\stringmanip.h"
#include "GP_Help\gphelp\crmcookies.h"
#include "GP_Base\gpbase\cloneutilityfunc.h"
#include "GP_Base\gpbase\countedptr.h"				//ARM_CountedPtr

#include "GP_Inflation\gpinflation\infcurv.h"
#include "GP_Inflation\gpinflation\infdata.h"
#include "GP_Inflation\gpinflation\resetmanager.h"
#include "GP_Inflation\gpinflation\sparsevolcube.h"
#include "GP_Inflation\gpinflation\seasonmanager.h"
#include "GP_Inflation\gpinflation\infidx.h"
#include "GP_Inflation\gpinflation\fixzc.h"
#include "GP_Inflation\gpinflation\infleg.h"
#include "GP_Inflation\gpinflation\infcorridorleg.h"
#include "GP_Inflation\gpinflation\infcapfloor.h"
#include "GP_Inflation\gpinflation\infcurvmodel.h"
#include "GP_Inflation\gpinflation\infbsmodel.h"
#include "GP_Inflation\gpinflation\infbssmiledmodel.h"
#include "GP_Inflation\gpinflation\implicitcorrel.h"
#include "GP_Inflation\gpinflation\livretACurv.h"
#include "GP_Inflation\gpinflation\infswopfactory.h"
#include "GP_Inflation\gpinflation\infmultibsmodel.h"
#include "GP_Inflation\gpinflation\infcapfloorrielyield.h"
#include "GP_Inflation\gpinflation\infhybriddigital.h"
#include "GP_Inflation\gpinflation\argconvdefault.h"
#include "GP_Inflation\gpinflation\infhybmodel.h"

#include "GP_Inflation\gpinflation\hybridinfleg.h"
#include "GP_Inflation\gpinflation\hybridinfpricer.h"

#include "GP_Inflation\gpinflation\infpayoffvisitor.h"
#include "GP_Inflation\gpinflation\infoptionspreadvisitor.h"
#include "GP_Inflation\gpinflation\infdoubledigitalvisitor.h"
#include "GP_Inflation\gpinflation\infhybridpayoffvisitor.h"
#include "GP_Inflation\gpinflation\infcorridorvisitor.h"

#include "GP_Inflation\gpinflation\infmodelvisitor.h"
#include "GP_Inflation\gpinflation\infhkvisitor.h"
#include "GP_Inflation\gpinflation\infbilogvisitor.h"

#include <gpinfra\mktdatamanagerrep.h>
#include <gpinfra\mktdatamanager.h>

#include "GP_Models\gpmodels\EqModVolSto.h"

using namespace ARM ;


static ARM_CPI_REF	RefIFRF[] = {
                       { 37347, 37438, 105.4},
                       { 37377, 37456, 105.6},
                       { 37408, 37487, 105.5},
                       { 37438, 37530, 105.5},
                       { 37469, 37557, 105.8},
                       { 37500, 37578, 106.0},
                       { 37530, 37606, 106.2},
                       { 37561, 37634, 106.2},
                       { 37591, 37657, 106.3},
                       { 37622, 37685, 106.3},
                       { 37653, 37712, 107.1},
                       { 37681, 37753, 107.5},
                       { 37712, 37796, 107.4}
                    };


static ARM_CPI_REF RefEMU[] = {
                       { 37347, 37438, 111.0},
                       { 37377, 37456, 111.2},
                       { 37408, 37487, 111.1},
                       { 37438, 37530, 110.9},
                       { 37469, 37547, 111.0},
                       { 37500, 37578, 111.4},
                       { 37530, 37606, 111.6},
                       { 37561, 37630, 111.5},
                       { 37591, 37657, 111.7},
                       { 37622, 37699, 111.5},
                       { 37653, 37733, 112.0},
                       { 37681, 37753, 112.5},
                       { 37712, 37796, 112.7}
                    };

static ARM_CPI_REF RefEMUT[] = {
                       { 37347, 37438, 111.4},
                       { 37377, 37456, 111.5},
                       { 37408, 37487, 111.5},
                       { 37438, 37526, 111.3},
                       { 37469, 37547, 111.4},
                       { 37500, 37578, 111.7},
                       { 37530, 37606, 112.0},
                       { 37561, 37630, 111.9},
                       { 37591, 37657, 112.1},
                       { 37622, 37693, 112.1},
                       { 37653, 37733, 112.5},
                       { 37681, 37753, 113.1},
                       { 37712, 37796, 113.2}
                    };

/// function to convert date into gap if it is not really a gap but rather a 
///	date! uses a calendar
int GetGapOrJulianDate( int gapInput )
{
	/// corresponding to 1st Jan 1995
	const int DateMin = 34700;
	if(gapInput<DateMin)
		return gapInput;
	else
	{
		/// else it is a date
		char  charDate[11];
		Local_XLDATE2ARMDATE( gapInput, charDate );
		ARM_Date date( charDate );	
		return (int) date.GetJulian();
	}
}


/*!
 * Function to create an inflation curve
 */
extern long ARMLOCAL_InfCurv_Create(
			double					asOfDate,
			const CCString&			indexName,
			double					CPIIndexValue,
			double					CPIIndexDate,
			const VECTOR<CCString>&	maturities,
			const VECTOR<double>&	values,
			long					MonthlyInterpType,
			long					DailyInterpType,
			long					DCFMonthly,
			long					DCFDaily,
			long					ExtrapolType,
			long					resetManagerId,
			long					seasonManagerId,
			ARM_result&				result, 
			long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) 
		||	!sameVectorSize( maturities, "Maturities", values, "Index Values", result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_InfCurv* infCurv = NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Inf Curve" );

		/// date conversions
		char  charDate[11];
		Local_XLDATE2ARMDATE( asOfDate, charDate );
		ARM_Date startDate( charDate );

		Local_XLDATE2ARMDATE( CPIIndexDate, charDate );
		ARM_Date CPIIndexDate( charDate );
			
		/// to convert to a vector of double
		/// we are forced to use a temp variable
		vector<double> pdvalues( values.size());
		int i;

		for(i = 0; i < values.size(); i++)
			pdvalues[i]	= values[i];

		/// we also need to do a conversion from 
		/// VECTOR<CCString> to vector<string>
		vector<string> maturitiesV( maturities.size() );
		
		/// convert dates to Julian type
		string maturity( "yYMmWwDd" );
		string s;
		char buffer[20];

		for( i = 0; i < maturities.size(); i++)
		{
			s = maturities[i];
			if( s.find_first_of( maturity ) == string::npos )
			{
				double julianDate = atof( s.c_str() );
				sprintf( buffer, "%f", XLDateToJulian( julianDate ) );
				s = buffer;
			}
			maturitiesV[i] = s;
		}

		ARM_ResetManager* resetManager = NULL;
		if( !GetObjectFromIdwNull( &resetManager, resetManagerId, ARM_RESETMANAGER ) )
		{
			result.setMsg ("ARM_ERR: reset manager is not of a good type");
			return ARM_KO;
		};

		ARM_SeasonalityManager* seasonManager = NULL;

		if( !GetObjectFromIdwNull( &seasonManager, seasonManagerId, ARM_SEASONMANAGER ) )
		{
			result.setMsg ("ARM_ERR: season manager is not of a good type");
			return ARM_KO;
		};


		infCurv	= new ARM_InfCurv( 
			startDate,
			CCSTringToSTLString( indexName ),
			CPIIndexValue,
			CPIIndexDate,
			maturitiesV,
			pdvalues,
			MonthlyInterpType,
			DailyInterpType,
			DCFMonthly,
			DCFDaily,
			ExtrapolType,
			resetManager,
			seasonManager);

		/// assign object
		if( !assignObject( infCurv, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		delete infCurv;

		x.DebugPrint();
		
		ARM_RESULT();
	}
}



/*!
 * Function to interpolate the curve
 * it can either do a CPI interpolation
 * or a ZC Rate depending on the 
 * type of what
 */
extern long ARMLOCAL_InfCurv_Interp(
	long			idCurve,
	double			CPIDateDble, 
	ARM_result&		result,
	const CCString&	tmpDCFLag,
	long			DailyInterpType,
	const CCString&	tmpResetLag,
	double			weight,
	instrument		what )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");

	try
	{
		/// date conversions
		char  charDate[11];
		Local_XLDATE2ARMDATE( CPIDateDble, charDate );
		ARM_Date CPIDate( charDate );

		ARM_InfCurv* crv = (ARM_InfCurv *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		/// An inflation curve should be of the base type ZERO_CURVE
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(crv, ARM_INFCURV ) == 0)
		{
			result.setMsg ("ARM_ERR: Inf Curve is not of a good type");
			return ARM_KO;
		}

		double dResult;
		string DCFLag	= CCSTringToSTLString( tmpDCFLag );
		string ResetLag	= CCSTringToSTLString( tmpResetLag );

		switch( what )
		{
			case CPI:
				dResult = crv->CPIInterpolate( CPIDate, DCFLag, 
					DailyInterpType, ResetLag, weight );
				break;
			case ZCRate:
				dResult = crv->ZCRateInterpolate( CPIDate, DCFLag, 
					DailyInterpType, ResetLag, weight );
				break;
			default:
	           throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
		        "Invalid instrument to interpolate");
		}

		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


/*!
 * Function to create an inflation leg
 * using date strip objects
 */
long ARMLOCAL_InfIdx_Create(
			const CCString& indexName,
			const CCString& resetLag,
			const CCString& DCFLag,
			long ccyId,
			long infCurveId,
			ARM_result&	result,
			long objId 	)
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_CountedPtr<ARM_Currency> ccy ;
	ARM_CountedPtr<ARM_InfCurv> infCurve ;
	ARM_InfIdx* infIdx	= NULL;
	bool createdNewCcy	= false;


	try
	{
		ARM_Currency* ccyPtr		= ccy.operator->();
		ARM_InfCurv* infCurvePtr	= infCurve.operator->();
		if( !GetObjectFromId( &ccyPtr, ccyId, ARM_CURRENCY) )
		{
			result.setMsg ("ARM_ERR: Currency is not of a good type");
			return ARM_KO;
		}

		if( !GetObjectFromId( &infCurvePtr, infCurveId, ARM_INFCURV  ) )
		{
			result.setMsg ("ARM_ERR: forward inf curve is not of a good type");
			return ARM_KO;
		}

		/// new inflation index
		infIdx	= new ARM_InfIdx( ccy,infCurve,
			CCSTringToSTLString( indexName ),
			CCSTringToSTLString( resetLag ),
			CCSTringToSTLString( DCFLag ) );

		/// assign object
		if( !assignObject( infIdx, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		delete infIdx;
		
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


/*!
 * Function to create an inflation swap leg
 */
extern long ARMLOCAL_InfLegwDateStrip_Create(
	const CCString& indexName,
	long rcvOrPay,
	int interpType,
	double multiple,
    double constant,
	long finalNotionalType,
	long numStripDateId,
	long denomStripDateId,
	ARM_result&	result,
	long objId 	)
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_InfLeg* infLeg			= NULL;

	try
	{
		/// numStripDate is read from the global table
		/// hence should not be deleted
		ARM_DateStrip* numStripDate		= NULL;
		ARM_DateStrip* denomStripDate	= NULL;

		/// if the id is given by ARM_NULL_OBJECT
		/// then returns a null pointor
		if( !GetObjectFromIdwNull( &numStripDate, numStripDateId, ARM_DATESTRIP ) )
		{
			result.setMsg ("ARM_ERR: Numerator Date Strip is not of a good type");
			return ARM_KO;
		};

		/// if the id is given by ARM_NULL_OBJECT
		/// then returns a null pointor
		if( !GetObjectFromIdwNull( &denomStripDate, denomStripDateId, ARM_DATESTRIP ) )
		{
			result.setMsg ("ARM_ERR: Denominator Date Strip is not of a good type" );
			return ARM_KO;
		};
		
		/// new inflation index
		infLeg	= new ARM_InfLeg( 
			CCSTringToSTLString( indexName ),
			rcvOrPay,
			interpType,
			multiple,
			constant,
			finalNotionalType,
			numStripDate,
			denomStripDate );

		/// assign object
		if( !assignObject( infLeg, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		delete infLeg;

		x.DebugPrint();
		
		ARM_RESULT();
	}
}


/* 
 * function to create a year to year inflation leg
 * 
 */
extern long ARMLOCAL_InfLegAllInputs_Create(
	double startDate, 
	double endDate,
	const CCString& indexName,
	int swapType,
	int rcvOrPay,
	int interpType,
	double multiple,
	double CoMultiple,	
	double constant,
	int resetFreq,
	int dayCount,
	const CCString& resetCalendar,
	int fwdRule,
	int intRule,
	int stubRule,
	int resetNumGap,
	int resetDenomGap,
	int payFreq,
	int payGap,
	const CCString& payCalendar,
	int adjFirstDate,
	int finalNotionalType,
	double firstReset,
	ARM_result&	result,
	long objId 	)
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_InfLeg* infLeg = NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Generic Inflation Leg" );

		/// dateConversion
		char  startDateChar[11];
		char  endDateChar[11];

		Local_XLDATE2ARMDATE( startDate, startDateChar );
		Local_XLDATE2ARMDATE( endDate, endDateChar );

		char* resetCal	= resetCalendar.c_str();
		char* payCal	= payCalendar.c_str();

		resetNumGap		= GetGapOrJulianDate( resetNumGap );
		resetDenomGap	= GetGapOrJulianDate( resetDenomGap );
		payGap			= ARM_InfLeg::GetGapFromGapOrDate( GetGapOrJulianDate( payGap ), payCal,	(ARM_Date) startDateChar);

		/// new inflation index
		infLeg	= new ARM_InfLeg( 
			(ARM_Date) startDateChar,
			(ARM_Date) endDateChar,
			CCSTringToSTLString( indexName ),
			swapType,
			rcvOrPay,
			interpType,
			multiple,
			CoMultiple,
			constant,
			resetFreq,
			dayCount,
			resetCal,
			fwdRule,
			intRule,
			stubRule,
			resetNumGap,
			resetDenomGap,
			payFreq,
			payGap,
			payCal,
			adjFirstDate,
			finalNotionalType,
			firstReset );

		delete resetCal;
		delete payCal;

		/// assign object
		if( !assignObject( infLeg, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		delete infLeg;
		
		x.DebugPrint();
		
		ARM_RESULT();
	}

}

//////
long ARM_HybridInfIrLeg_CreateFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	string tmp;

	ARM_GenericParams*	genericParams	= GetGenericParams();
	ARM_HybridInfIrLeg*	hybridInfIrLeg	= NULL;


	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "HybridInfIrLeg Create" );

		char asOfDate	[20];
		char startDate	[20];
		char endDate	[20];

		Local_XLDATE2ARMDATE(	genericParams->GetParamValue("AsOfDate").GetDouble(),	asOfDate);
		Local_XLDATE2ARMDATE(	genericParams->GetParamValue("StartDate").GetDouble(), startDate);
		Local_XLDATE2ARMDATE(	genericParams->GetParamValue("EndDate").GetDouble(),	endDate);

		string	resetCal		= genericParams->GetParamValue("ResetCal").GetString();
		stringToUpper(resetCal);

		string	payCal			= genericParams->GetParamValue("PayCal").GetString();
		stringToUpper(payCal);

		string	mainIndex		= genericParams->GetParamValue("MainIndex").GetString();
		stringToUpper(mainIndex);

		string	subIndex		= genericParams->GetParamValue("SubIndex").GetString();
		stringToUpper(subIndex);

		string	supIndex		= genericParams->GetParamValue("SupIndex").GetString();
		stringToUpper(supIndex);
		
		tmp	= genericParams->GetParamValue("ResetFreq").GetString();
		stringToUpper(tmp);
		int	resetFreq	= ARM_ArgConv_LgNameFrequency.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("PayFreq").GetString();
		stringToUpper(tmp);
		int	payFreq	= ARM_ArgConv_LgNameFrequency.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("ResetTiming").GetString();
		stringToUpper(tmp);
		int	resetTiming	= ARM_ArgConv_Timing.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("PayTiming").GetString();
		stringToUpper(tmp);
		int	payTiming	= ARM_ArgConv_Timing.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("IntRule").GetString();
		stringToUpper(tmp);
		int	irIntRule	= ARM_ArgConv_IntRule.GetNumber(tmp);
	
		tmp	= genericParams->GetParamValue("StubRule").GetString();
		stringToUpper(tmp);
		int	stubRule	= ARM_ArgConv_StubRules.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("FwdRule").GetString();
		stringToUpper(tmp);
		int	fwdRule	= ARM_ArgConv_FwdRule.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("AdjFirstRule").GetString();
		stringToUpper(tmp);
		int	adjFirstRule	= ARM_ArgConv_IntRule.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("dayCount").GetString();
		stringToUpper(tmp);
		int	dayCount	= ARM_ArgConv_LgNameDayCount.GetNumber(tmp);

		stringToUpper(tmp);
		int	cpnDayCount	= ARM_ArgConv_LgNameDayCount.GetNumber(tmp);

		string	currency		=	genericParams->GetParamValue("Currency").GetString().c_str();
		int		infInterType	=	ARM_ArgConv_InterpolInfType.GetNumber(genericParams->GetParamValue("InfInterType").GetString().c_str());

		int		resetGap		=	genericParams->GetParamValue("ResetGap").GetDouble();
		int		resetNumGap		=	genericParams->GetParamValue("ResetNumGap").GetDouble();
		int		resetDemGap		=	genericParams->GetParamValue("ResetDemGap").GetDouble();
		int		payGap			=	genericParams->GetParamValue("PayGap").GetDouble();

		int		finNotioType	=	K_NX_NONE;
		int		firstReset		=	GETDEFAULTVALUE;
		int		decompFreq		=	K_COMP_PROP;
		int		compMeth		=	K_COMP_PROP;
		int		decompPriceFlag	=	1;
		int		maturity		=	-1;
		string	refDate			=	"NULL";

		// Compilation of all inputs
		
		map<string,ARM_Date>	mDate;
		map<string,string>		mString;
		map<string,ARM_Curve*>	mCurve;
		map<string,int>			mInt;
		map<string,double>		mDouble;
		
		mDate.insert(pair<string, ARM_Date>			(	"asOfDate"		,	asOfDate		)	);
		mDate.insert(pair<string, ARM_Date>			(	"startDate"		,	startDate		)	);
		mDate.insert(pair<string, ARM_Date>			(	"endDate"		,	endDate			)	);

		mString.insert(pair<string, string>			(	"mainIndex"		,	mainIndex		)	);
		mString.insert(pair<string, string>			(	"subIndex"		,	subIndex		)	);
		mString.insert(pair<string, string>			(	"supIndex"		,	supIndex		)	);
		mString.insert(pair<string, string>			(	"currency"		,	currency		)	);
		mString.insert(pair<string, string>			(	"resetCal"		,	resetCal		)	);
		mString.insert(pair<string, string>			(	"payCal"		,	payCal			)	);
		mString.insert(pair<string, string>			(	"refDate"		,	refDate			)	);

		mInt.insert(pair<string, int>				(	"resetFreq"		,	resetFreq		)	);
		mInt.insert(pair<string, int>				(	"resetGap"		,	resetGap		)	);
		mInt.insert(pair<string, int>				(	"resetTiming"	,	resetTiming		)	);
		mInt.insert(pair<string, int>				(	"resetNumGap"	,	resetNumGap		)	);
		mInt.insert(pair<string, int>				(	"resetDemGap"	,	resetDemGap		)	);
		mInt.insert(pair<string, int>				(	"payFreq"		,	payFreq			)	);
		mInt.insert(pair<string, int>				(	"payGap"		,	payGap			)	);
		mInt.insert(pair<string, int>				(	"payTiming"		,	payTiming		)	);
		mInt.insert(pair<string, int>				(	"dayCount"	,		dayCount		)	);
		mInt.insert(pair<string, int>				(	"infInterType"	,	infInterType	)	);		
		mInt.insert(pair<string, int>				(	"fwdRule"		,	fwdRule			)	);
		mInt.insert(pair<string, int>				(	"intRule"		,	irIntRule		)	);
		mInt.insert(pair<string, int>				(	"stubRule"		,	stubRule		)	);
		mInt.insert(pair<string, int>				(	"adjFirstRule"	,	adjFirstRule	)	);
		mInt.insert(pair<string, int>				(	"finNotioType"	,	finNotioType	)	);
		mInt.insert(pair<string, int>				(	"firstReset"	,	firstReset		)	);
		mInt.insert(pair<string, int>				(	"decompFreq"	,	decompFreq		)	);
		mInt.insert(pair<string, int>				(	"compMeth"		,	compMeth		)	);
		mInt.insert(pair<string, int>				(	"decompPriceFlag",	decompPriceFlag	)	);
		mInt.insert(pair<string, int>				(	"maturity"		,	maturity		)	);

		hybridInfIrLeg = 	new ARM_HybridInfIrLeg(	mDate, mString, mInt);
		if ( !assignObject( hybridInfIrLeg, result, objId ) )
			return ARM_KO; 
		else
			return ARM_OK; 
	}

	catch(Exception& x)
	{
		delete	hybridInfIrLeg;

		x.DebugPrint();
		ARM_RESULT();
	}
}



long ARM_HybridInfIrPayOff_CreateFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	string tmp;

	ARM_GenericParams*		genericParams		= GetGenericParams();
	ARM_InfHybridPayOff*	hybridInfIrPayOff	= NULL;


	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "HybridInfIrPayOff Create" );

		ARM::ARM_MAP_Curve cpnCurve;
		ARM::ARM_MAP_Curve optCurve;


//==> MainCpn
		string		mainCpnName	= genericParams->GetParamValue("MainCpnName").GetString();
		stringToUpper(mainCpnName);	
		ARM_Curve*	mainCpnCoef	= NULL;
		mainCpnCoef = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("MainCpnCoef").GetObjectId()));
		if (!mainCpnCoef)	{
			result.setMsg ("ARM_ERR: nominal should be a curve");
			return ARM_KO;
		}
		else{
			pair<string,ARM_GP_CurvePtr> p( mainCpnName,ARM_GP_CurvePtr(mainCpnCoef) );
			cpnCurve.insert(p);
		}

//==> SubCpn
		string		subCpnName	= genericParams->GetParamValue("SubCpnName").GetString();
		stringToUpper(subCpnName);	
		ARM_Curve*	subCpnCoef	= NULL;
		subCpnCoef = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("SubCpnCoef").GetObjectId()));
		if (!subCpnCoef)	{
			result.setMsg ("ARM_ERR: nominal should be a curve");
			return ARM_KO;
		}
		else{
			pair<string,ARM_GP_CurvePtr> p( subCpnName,ARM_GP_CurvePtr(subCpnCoef) );
			cpnCurve.insert(p);
		}

//==> SupCpn
		string		supCpnName	= genericParams->GetParamValue("SupCpnName").GetString();
		stringToUpper(supCpnName);	
		ARM_Curve*	supCpnCoef	= NULL;
		supCpnCoef = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("SupCpnCoef").GetObjectId()));
		if (!supCpnCoef)	{
			result.setMsg ("ARM_ERR: nominal should be a curve");
			return ARM_KO;
		}
		else{
			pair<string,ARM_GP_CurvePtr> p( supCpnName,ARM_GP_CurvePtr(supCpnCoef) );
			cpnCurve.insert(p);
		}

//==> CstCpn
		string		cstCpnName	= "NO";	
		ARM_Curve*	cstCpnCoef	= NULL;
		cstCpnCoef = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("CstCpnCoef").GetObjectId()));
		if (!cstCpnCoef)	{
			result.setMsg ("ARM_ERR: nominal should be a curve");
			return ARM_KO;
		}
		else{
			pair<string,ARM_GP_CurvePtr> p( cstCpnName,ARM_GP_CurvePtr(cstCpnCoef) );
			cpnCurve.insert(p);
		}

//==> MainOpt
		string		mainOptName	= genericParams->GetParamValue("MainOptName").GetString();
		stringToUpper(mainOptName);	
		ARM_Curve*	mainOptCoef	= NULL;
		mainOptCoef = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("MainOptCoef").GetObjectId()));
		if (!mainOptCoef)	{
			result.setMsg ("ARM_ERR: nominal should be a curve");
			return ARM_KO;
		}
		else{
			pair<string,ARM_GP_CurvePtr> p( mainOptName,ARM_GP_CurvePtr(mainOptCoef) );
			optCurve.insert(p);
		}

//==> SubOpt
		string		subOptName	= genericParams->GetParamValue("SubOptName").GetString();
		stringToUpper(subOptName);	
		ARM_Curve*	subOptCoef	= NULL;
		subOptCoef = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("SubOptCoef").GetObjectId()));
		if (!subOptCoef)	{
			result.setMsg ("ARM_ERR: nominal should be a curve");
			return ARM_KO;
		}
		else{
			pair<string,ARM_GP_CurvePtr> p( subOptName,ARM_GP_CurvePtr(subOptCoef) );
			optCurve.insert(p);
		}

//==> SupOpt
		string		supOptName	= genericParams->GetParamValue("SupOptName").GetString();
		stringToUpper(supOptName);	
		ARM_Curve*	supOptCoef	= NULL;
		supOptCoef = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("SupOptCoef").GetObjectId()));
		if (!supOptCoef)	{
			result.setMsg ("ARM_ERR: nominal should be a curve");
			return ARM_KO;
		}
		else{
			pair<string,ARM_GP_CurvePtr> p( supOptName,ARM_GP_CurvePtr(supOptCoef) );
			optCurve.insert(p);
		}

//==> CstOpt
		string		cstOptName	= "NO";
		ARM_Curve*	cstOptCoef	= NULL;
		cstOptCoef = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("CstOptCoef").GetObjectId()));
		if (!cstOptCoef)	{
			result.setMsg ("ARM_ERR: nominal should be a curve");
			return ARM_KO;
		}
		else{
			pair<string,ARM_GP_CurvePtr> p( cstOptName,ARM_GP_CurvePtr(cstOptCoef) );
			optCurve.insert(p);
		}

		hybridInfIrPayOff = new ARM_InfHybridPayOff(	cpnCurve,  optCurve);
		if ( !assignObject( hybridInfIrPayOff, result, objId ) )
			return ARM_KO; 
		else
			return ARM_OK; 
	}

	catch(Exception& x)
	{
		delete	hybridInfIrPayOff;

		x.DebugPrint();
		ARM_RESULT();
	}
}



/////
long ARM_Corridor_InfIr_CreateFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	string tmp;

	ARM_GenericParams*	genericParams	= GetGenericParams();
	ARM_InfCorridorLeg*	infCorridorLeg	= NULL;

	ARM_Curve*	notional	= NULL;
	ARM_Curve*	irLeverage	= NULL;
	ARM_Curve*	infLeverage	= NULL;
	ARM_Curve*	constant	= NULL;
	ARM_Curve*	multipleUp	= NULL;
	ARM_Curve*	multipleDown= NULL;
	ARM_Curve*	rangeUp		= NULL;
	ARM_Curve*	rangeDown	= NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Corridor Inf/Ir Leg Create" );

		char asOfDate	[20];
		char startDate	[20];
		char endDate	[20];

		Local_XLDATE2ARMDATE(	genericParams->GetParamValue("AsOfDate").GetDouble(),	asOfDate);
		Local_XLDATE2ARMDATE(	genericParams->GetParamValue("StartDate").GetDouble(), startDate);
		Local_XLDATE2ARMDATE(	genericParams->GetParamValue("EndDate").GetDouble(),	endDate);

		string	resetCal		= genericParams->GetParamValue("ResetCal").GetString();
		stringToUpper(resetCal);

		string	payCal			= genericParams->GetParamValue("PayCal").GetString();
		stringToUpper(payCal);

		string	infIndex		= genericParams->GetParamValue("InfIndex").GetString();
		stringToUpper(infIndex);
		
		tmp	= genericParams->GetParamValue("ResetFreq").GetString();
		stringToUpper(tmp);
		int	resetFreq	= ARM_ArgConv_LgNameFrequency.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("PayFreq").GetString();
		stringToUpper(tmp);
		int	payFreq	= ARM_ArgConv_LgNameFrequency.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("ResetTiming").GetString();
		stringToUpper(tmp);
		int	resetTiming	= ARM_ArgConv_Timing.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("PayTiming").GetString();
		stringToUpper(tmp);
		int	payTiming	= ARM_ArgConv_Timing.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("IrIntRule").GetString();
		stringToUpper(tmp);
		int	irIntRule	= ARM_ArgConv_IntRule.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("InfIntRule").GetString();
		stringToUpper(tmp);
		int	infIntRule	= ARM_ArgConv_IntRule.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("StubRule").GetString();
		stringToUpper(tmp);
		int	stubRule	= ARM_ArgConv_StubRules.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("FwdRule").GetString();
		stringToUpper(tmp);
		int	fwdRule	= ARM_ArgConv_FwdRule.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("AdjFirstRule").GetString();
		stringToUpper(tmp);
		int	adjFirstRule	= ARM_ArgConv_IntRule.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("IrDayCount").GetString();
		stringToUpper(tmp);
		int	irDayCount	= ARM_ArgConv_LgNameDayCount.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("InfDayCount").GetString();
		stringToUpper(tmp);
		int	infDayCount	= ARM_ArgConv_LgNameDayCount.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("CpnDayCount").GetString();
		stringToUpper(tmp);
		int	cpnDayCount	= ARM_ArgConv_LgNameDayCount.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("RecOrPay").GetString();
		stringToUpper(tmp);
		int	recOrPay	= ARM_ArgConv_RcvOrPay.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("IsIrCritera").GetString();
		stringToUpper(tmp);
		int	irCritera	= ARM_ArgConv_YesNo.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("IsModulable").GetString();
		stringToUpper(tmp);
		int	isModulable	= ARM_ArgConv_YesNo.GetNumber(tmp);

		tmp	= genericParams->GetParamValue("Model").GetString();
		stringToUpper(tmp);
		int	model	= ARM_ArgConv_Copula.GetNumber(tmp);

		notional = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("Notional").GetObjectId()));
		if (!notional)
		{
			result.setMsg ("ARM_ERR: nominal should be a curve");
			return ARM_KO;
		}

		irLeverage = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("IrLeverage").GetObjectId()));
		if (!irLeverage)
		{
			result.setMsg ("ARM_ERR: IrLeverage factor should be a curve");
			return ARM_KO;
		}

		infLeverage = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("InfLeverage").GetObjectId()));
		if (!infLeverage)
		{
			result.setMsg ("ARM_ERR: InfLeverage factor should be a curve");
			return ARM_KO;
		}

		constant = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("Constant").GetObjectId()));
		if (!constant)
		{
			result.setMsg ("ARM_ERR: Constant factor should be a curve");
			return ARM_KO;
		}
		multipleUp = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("MultipleUp").GetObjectId()));
		if (!multipleUp)
		{
			result.setMsg ("ARM_ERR: Multiple Up factor should be a curve");
			return ARM_KO;
		}
		multipleDown = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("MultipleDown").GetObjectId()));
		if (!multipleDown)
		{
			result.setMsg ("ARM_ERR: Multiple down factor should be a curve");
			return ARM_KO;
		}
		rangeUp = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("RangeUp").GetObjectId()));
		if (!rangeUp)
		{
			result.setMsg ("ARM_ERR: Range up factor should be a curve");
			return ARM_KO;
		}

		rangeDown = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("RangeDown").GetObjectId()));
		if (!rangeDown)
		{
			result.setMsg ("ARM_ERR: Range down factor should be a curve");
			return ARM_KO;
		}

		int		irIndex			= ARM_ArgConv_IndexType.GetNumber(genericParams->GetParamValue("IrIndex").GetString().c_str());

		string	currency		=	genericParams->GetParamValue("Currency").GetString().c_str();
		int		infInterType	=	ARM_ArgConv_InterpolInfType.GetNumber(genericParams->GetParamValue("InfInterType").GetString().c_str());

		int		resetGap		=	genericParams->GetParamValue("ResetGap").GetDouble();
		int		resetNumGap		=	genericParams->GetParamValue("ResetNumGap").GetDouble();
		int		resetDemGap		=	genericParams->GetParamValue("ResetDemGap").GetDouble();
		int		payGap			=	genericParams->GetParamValue("PayGap").GetDouble();

		double	epsilon			=	genericParams->GetParamValue("Epsilon").GetDouble();
		double	nbGaussLeg		=	genericParams->GetParamValue("NbGaussLeg").GetDouble();
		double	ptGaussLeg		=	genericParams->GetParamValue("PtGaussLeg").GetDouble();

		int		finNotioType	=	K_NX_NONE;
		int		firstReset		=	GETDEFAULTVALUE;
		int		decompFreq		=	K_COMP_PROP;
		int		compMeth		=	K_COMP_PROP;
		int		decompPriceFlag	=	1;
		int		maturity		=	-1;
		string	refDate			=	"NULL";

		// Compilation of all inputs
		
		map<string,ARM_Date>	mDate;
		map<string,string>		mString;
		map<string,ARM_Curve*>	mCurve;
		map<string,int>			mInt;
		map<string,double>		mDouble;
		
		mDate.insert(pair<string, ARM_Date>			(	"asOfDate"		,	asOfDate		)	);
		mDate.insert(pair<string, ARM_Date>			(	"startDate"		,	startDate		)	);
		mDate.insert(pair<string, ARM_Date>			(	"endDate"		,	endDate			)	);

		mString.insert(pair<string, string>			(	"infIndex"		,	infIndex		)	);
		mString.insert(pair<string, string>			(	"currency"		,	currency		)	);
		mString.insert(pair<string, string>			(	"resetCal"		,	resetCal		)	);
		mString.insert(pair<string, string>			(	"payCal"		,	payCal			)	);
		mString.insert(pair<string, string>			(	"refDate"		,	refDate			)	);

		mCurve.insert(pair<string, ARM_Curve*>		(	"notional"		,	notional		)	);
		mCurve.insert(pair<string, ARM_Curve*>		(	"irLeverage"	,	irLeverage		)	);
		mCurve.insert(pair<string, ARM_Curve*>		(	"infLeverage"	,	infLeverage		)	);
		mCurve.insert(pair<string, ARM_Curve*>		(	"constant"		,	constant		)	);
		mCurve.insert(pair<string, ARM_Curve*>		(	"multipleUp"	,	multipleUp		)	);
		mCurve.insert(pair<string, ARM_Curve*>		(	"multipleDown"	,	multipleDown	)	);
		mCurve.insert(pair<string, ARM_Curve*>		(	"rangeUp"		,	rangeUp			)	);
		mCurve.insert(pair<string, ARM_Curve*>		(	"rangeDown"		,	rangeDown		)	);

		mInt.insert(pair<string, int>				(	"infInterType"	,	infInterType	)	);
		mInt.insert(pair<string, int>				(	"irIndex"		,	irIndex			)	);
		mInt.insert(pair<string, int>				(	"resetFreq"		,	resetFreq		)	);
		mInt.insert(pair<string, int>				(	"resetGap"		,	resetGap		)	);
		mInt.insert(pair<string, int>				(	"resetTiming"	,	resetTiming		)	);
		mInt.insert(pair<string, int>				(	"resetNumGap"	,	resetNumGap		)	);
		mInt.insert(pair<string, int>				(	"resetDemGap"	,	resetDemGap		)	);
		mInt.insert(pair<string, int>				(	"payFreq"		,	payFreq			)	);
		mInt.insert(pair<string, int>				(	"payGap"		,	payGap			)	);
		mInt.insert(pair<string, int>				(	"payTiming"		,	payTiming		)	);
		mInt.insert(pair<string, int>				(	"cpnDayCount"	,	cpnDayCount		)	);	
		mInt.insert(pair<string, int>				(	"irDayCount"	,	irDayCount		)	);
		mInt.insert(pair<string, int>				(	"infDayCount"	,	infDayCount		)	);
		mInt.insert(pair<string, int>				(	"fwdRule"		,	fwdRule			)	);
		mInt.insert(pair<string, int>				(	"irIntRule"		,	irIntRule		)	);
		mInt.insert(pair<string, int>				(	"infIntRule"	,	infIntRule		)	);
		mInt.insert(pair<string, int>				(	"stubRule"		,	stubRule		)	);
		mInt.insert(pair<string, int>				(	"adjFirstRule"	,	adjFirstRule	)	);
		mInt.insert(pair<string, int>				(	"recOrPay"		,	recOrPay		)	);
		mInt.insert(pair<string, int>				(	"irCritera"		,	irCritera		)	);
		mInt.insert(pair<string, int>				(	"isModulable"	,	isModulable		)	);
		mInt.insert(pair<string, int>				(	"finNotioType"	,	finNotioType	)	);
		mInt.insert(pair<string, int>				(	"firstReset"	,	firstReset		)	);
		mInt.insert(pair<string, int>				(	"decompFreq"	,	decompFreq		)	);
		mInt.insert(pair<string, int>				(	"compMeth"		,	compMeth		)	);
		mInt.insert(pair<string, int>				(	"decompPriceFlag",	decompPriceFlag	)	);
		mInt.insert(pair<string, int>				(	"maturity"		,	maturity		)	);

		mDouble.insert(pair<string, double>			(	"epsilon"		,	epsilon			)	);
		mDouble.insert(pair<string, double>			(	"nbGaussLeg"	,	nbGaussLeg		)	);
		mDouble.insert(pair<string, double>			(	"ptGaussLeg"	,	ptGaussLeg		)	);

		infCorridorLeg = new ARM_InfCorridorLegWithModel<ARM_GaussCopula>(	mDate, mString, mCurve, mInt, mDouble);

		if ( !assignObject( infCorridorLeg, result, objId ) )
			return ARM_KO; 
		else
			return ARM_OK; 
	}

	catch(Exception& x)
	{
		delete	infCorridorLeg;

		x.DebugPrint();
		ARM_RESULT();
	}
}


long ARMLOCAL_infYCmod(	long idZeroCurve, 
						long idFwdInfCurve, 
						ARM_result& result, 
						long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_InfCurvModel* myModel = NULL;

	try
	{
		ARM_ZeroCurve* zeroCurveModel = NULL;
		if( !GetObjectFromId( &zeroCurveModel, idZeroCurve, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: forecast fwd curve is not of a good type");
			return ARM_KO;
		};

		ARM_InfCurv* fwdInfModel = NULL;

		if( !GetObjectFromId( &fwdInfModel, idFwdInfCurve, ARM_INFCURV ) )
		{
			result.setMsg ("ARM_ERR: forward inf curve is not of a good type");
			return ARM_KO;
		};

		myModel = new ARM_InfCurvModel( zeroCurveModel, fwdInfModel );

		/// assign object
		if( !assignObject( myModel, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
		
	catch(Exception& x)
	{
		delete myModel;
		
		x.DebugPrint();
		
		ARM_RESULT();
	}
}



/*! 
 * function to create a fix ZC leg
 * 
 */
extern long ARMLOCAL_fixZC(
	double startDate, 
	double endDate,
	double fixRate,
	int rcvOrPay,
	int dayCount,
	int intRule,
	int stubRule,
	int payGap,
	const CCString& payCalendar,
	long ccyId,
	ARM_result& result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_FixZC* fixZCLeg = NULL;
	ARM_Currency* ccy	= NULL;

	try
	{
		/// dateConversion
		char  startDateChar[11];
		char  endDateChar[11];

		Local_XLDATE2ARMDATE( startDate, startDateChar );
		Local_XLDATE2ARMDATE( endDate, endDateChar );

		if( ccyId == ARM_NULL_OBJECT )
			ccy = ARM_DEFAULT_CURRENCY;
		else
			if( !GetObjectFromId( &ccy, ccyId, ARM_CURRENCY ) )
			{
				result.setMsg ("ARM_ERR: ccy is not of a good type");
				return ARM_KO;
			};

		char* payCal = payCalendar.c_str();
		
		if( strcmp( payCal, GETDEFAULTVALUESTR ) == 0 )
		{
			delete payCal;
			payCal = ARM_DEFAULT_CURRENCY->GetCcyName();
		}

		payGap	= ARM_InfLeg::GetGapFromGapOrDate( GetGapOrJulianDate(payGap), payCal,	(ARM_Date) startDateChar);

		/// new inflation index
		fixZCLeg = new ARM_FixZC( 
			(ARM_Date) startDateChar,
			(ARM_Date) endDateChar,
			fixRate,
			rcvOrPay,
			dayCount,
			intRule,
			stubRule,
			payGap,
			payCal,
			ccy );

		/// assign object
		if( !assignObject( fixZCLeg, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		delete fixZCLeg;
		
		x.DebugPrint();
		
		ARM_RESULT();
	}

}


/*! 
 * function to create an inflation reset manager
 */
extern long ARMLOCAL_GetData(
	const VECTOR<CCString>& vecCCString,
	int nbrows,
	int nbcolumns,
	ARM_result& result,
	long objId )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");

	ARM_ResetManager* resetManager = NULL;

	try
	{
		vector<string> vecString( vecCCString.size() );

		int i;

		for( i=0; i<vecCCString.size(); i++ )
			vecString[i] = vecCCString[i];

		int decimal, sign;

		// conversion de date Excel en date Julienne
		for( i=2; i<nbrows; i++ )
		{
			if (strcmp(vecString[i*nbcolumns].c_str(),"") != 0)
			{
				char sDate[11];
				Local_XLDATE2ARMDATE(atof(vecCCString[i*nbcolumns]),sDate);
				vecString[i*nbcolumns] = ecvt(((ARM_Date)sDate).GetJulian(),7, &decimal, &sign);
			}
		}

		resetManager = new ARM_ResetManager( vecString, nbrows, nbcolumns );

		if( !assignObject( resetManager, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}

	catch(Exception& x)
	{
		delete resetManager;

		x.DebugPrint();

		ARM_RESULT();
	}
}



/*! 
 * function to create a sparse vol cube
 * a sparse vol cube is an input to create a vol cube
 *
 * this function creates a sparse vol cube and fill it
 *
 */
extern long ARMLOCAL_SparseVolCube_CreateNFill(
	double asOfDate,
	double lastKnownDate,
	const CCString& indexName,
	long dim1Type,
	const VECTOR<CCString>& dim1Value,
	long dim2Type,
	double dim2Value,
	const VECTOR<double>& strikes,
	const VECTOR<double>& vols,
	long strikeType,
	long volType,
	ARM_result& result,
	long objId )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	std::vector<double>& vDim1	= NULL;
	std::vector<double>& vStrikes	= NULL;
	ARM_Matrix* vVols		= NULL;

	int dim1_size			= dim1Value.size ();
	int strikes_size		= strikes.size ();
	int vols_size			= vols.size();

	if(vols_size != dim1_size * strikes_size )
	{
		result.setMsg ("ARM_ERR: check your volatility matrix dimension");
		return ARM_KO;
	}

	ARM_SparseVolCube* sparseVolCube= NULL;

	CCString msg ("");

	try
	{
		char startDate[11];
		Local_XLDATE2ARMDATE(asOfDate,startDate);
		ARM_Date asOfARMDate( startDate );

		char lastKnownDateChar[11];
		Local_XLDATE2ARMDATE(lastKnownDate,lastKnownDateChar);
		ARM_Date lastKnownDateARMDate( lastKnownDateChar );

		
		ARM_Vector* vtmpDim1	= CreateARMMatuVectorFromStrVECTOR(dim1Value);
		vDim1					= To_pARM_GP_Vector( vtmpDim1 );
		delete vtmpDim1;
		vStrikes= CreateARMGPVectorFromVECTOR(strikes);
		vVols	= CreateARMMatrixFromVECTOR( vols, dim1_size, strikes_size );

		/// to avoid memory leak!
		char* indexNameChar = indexName.c_str();
		string indexNameStr( indexNameChar );
		delete indexNameChar;

		sparseVolCube	= new ARM_SparseVolCube( 
			asOfARMDate,
			lastKnownDateARMDate,
			indexNameStr,
			dim1Type,
			vDim1,
			dim2Type,
			dim2Value,
			vStrikes,
			vVols,
			strikeType,
			volType );

		delete vStrikes;
		vStrikes = NULL;
		delete vVols;
		vVols = NULL;

		if( !assignObject( sparseVolCube, result, objId ) )
		{
			/// the assignObjet explicitly delete the sparseVolCube
			/// hence not required to delete it here
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}

	catch(Exception& x)
	{
		delete sparseVolCube;
		delete vStrikes;
		delete vVols;

		x.DebugPrint();

		ARM_RESULT();
	}
}



/*! 
 * this functions fills an already existing sparse vol cube
 */
extern long ARMLOCAL_SparseVolCube_Fill(
	long dim1Type,
	const VECTOR<CCString>& dim1,
	long dim2Type,
	double dim2Value,
	const VECTOR<double>& strikes,
	const VECTOR<double>& vols,
	long sparseVolCubeId,
	ARM_result& result,
	long objId )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	std::vector<double>& vDim1		= NULL;
	std::vector<double>& vStrikes	= NULL;
	ARM_Matrix* vVols				= NULL;

	int dim1_size					= dim1.size ();
	int strikes_size				= strikes.size ();
	int vols_size					= vols.size();
	ARM_SparseVolCube* sparseVolCube = NULL;

	if(vols_size != dim1_size * strikes_size )
	{
		result.setMsg ("ARM_ERR: check your volatility matrix dimension");
		return ARM_KO;
	}


	CCString msg ("");

	try
	{
		ARM_Vector* vtmpDim1	= CreateARMMatuVectorFromStrVECTOR(dim1);
		vDim1					= To_pARM_GP_Vector(vtmpDim1);
		delete vtmpDim1;

		vStrikes= CreateARMGPVectorFromVECTOR(strikes);
		vVols	= CreateARMMatrixFromVECTOR( vols, dim1_size, strikes_size );

		ARM_SparseVolCube* prevSparseVolCube = NULL;
		if( !GetObjectFromId( &prevSparseVolCube, sparseVolCubeId, ARM_SPARSE_VOL_CUBE ) )
		{
			result.setMsg ("ARM_ERR: sparse vol cube is not of a good type");
			return ARM_KO;
		};

		sparseVolCube = (ARM_SparseVolCube*) prevSparseVolCube->Clone();

		sparseVolCube->AddVolCurves( dim1Type, vDim1, dim2Type, dim2Value, vStrikes, vVols );
	
		delete vStrikes;
		vStrikes = NULL;
		delete vVols;
		vVols = NULL;

		if( !assignObject( sparseVolCube, result, objId ) )
		{
			/// the assignObjet explicitly delete the sparseVolCube
			/// hence not required to delete it here
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}

	}

	catch(Exception& x)
	{
		delete sparseVolCube;
		delete vStrikes;
		delete vVols;

		x.DebugPrint();

		ARM_RESULT();
	}
}






/*! 
 * this functions fills an already existing sparse vol cube
 */
extern long ARMLOCAL_VolCubeFromSparseVolCube(
	long sparseVolCubeId,
	ARM_result& result,
	long objId )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");

	ARM_VolCube* volCube = NULL;

	try
	{
		ARM_SparseVolCube* sparseVolCube	= NULL;
		if( !GetObjectFromId( &sparseVolCube, sparseVolCubeId, ARM_SPARSE_VOL_CUBE ) )
		{
			result.setMsg ("ARM_ERR: sparse vol cube is not of a good type");
			return ARM_KO;
		};

		volCube = sparseVolCube->ConvertToVolCube();
	
		if( !assignObject( volCube, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}

	}

	catch(Exception& x)
	{
		delete volCube;

		x.DebugPrint();

		ARM_RESULT();
	}
}




/* !
 * function to create an inflation BS model
 */
long ARMLOCAL_infBSMod(	double asOfDateDble,
						long idZeroCurve, 
						long idFwdInfCurve, 
						long idVolCurve,
						long idCorrelManager,
						long idIRBSModel,
						long infSwoptCurveId,
						long IRSwoptCurveId,
						ARM_result& result, 
						long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_InfBSModel* myModel = NULL;

	try
	{
		char  charDate[11];
		Local_XLDATE2ARMDATE( asOfDateDble, charDate );
		ARM_Date asOfDate( charDate );

		ARM_ZeroCurve* zeroCurveModel = NULL;
		if( !GetObjectFromId( &zeroCurveModel, idZeroCurve, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: forecast fwd curve is not of a good type");
			return ARM_KO;
		};

		ARM_InfCurv* fwdInfModel = NULL;
		if( !GetObjectFromId( &fwdInfModel, idFwdInfCurve, ARM_INFCURV ) )
		{
			result.setMsg ("ARM_ERR: forward inf curve is not of a good type");
			return ARM_KO;
		};

		ARM_VolCurve* volCurve = NULL;
		if( !GetObjectFromId( &volCurve, idVolCurve, ARM_VOL_CURVE ) )
		{
			result.setMsg ("ARM_ERR: volatility curve is not of a good type");
			return ARM_KO;
		};

		ARM_CorrelManager* correlManager = NULL;
		if( !GetObjectFromIdwNull( &correlManager, idCorrelManager, ARM_CORRELMANAGER ) )
		{
			result.setMsg ("ARM_ERR: correl manager is not of a good type");
			return ARM_KO;
		};


		ARM_BSModel* IRBSModel = NULL;
		if( !GetObjectFromIdwNull( &IRBSModel, idIRBSModel, ARM_BSMODEL ) &&
			!GetObjectFromIdwNull( &IRBSModel, idIRBSModel, ARM_BSSMILEDMODEL ) 
			)
		{
			result.setMsg ("ARM_ERR: IR BS Model is not of a good type");
			return ARM_KO;
		};

		ARM_VolCurve* infSwoptCurve = NULL;
		if( !GetObjectFromIdwNull( &infSwoptCurve, infSwoptCurveId, ARM_VOL_CURVE ) )
		{
			result.setMsg ("ARM_ERR: volatility curve is not of a good type");
			return ARM_KO;
		};

		ARM_VolCurve* IRSwoptCurve = NULL;
		if( !GetObjectFromIdwNull( &IRSwoptCurve , IRSwoptCurveId, ARM_VOL_CURVE ) )
		{
			result.setMsg ("ARM_ERR: volatility curve is not of a good type");
			return ARM_KO;
		};

		ARM_VolCube* IRSwoptCube = NULL;
		if( GetObjectFromId(&IRSwoptCube, IRSwoptCurveId, ARM_VOL_CUBE)==true )
		{
			myModel = new ARM_InfBSModel( asOfDate, zeroCurveModel, fwdInfModel, volCurve, 
			correlManager, IRBSModel, infSwoptCurve, IRSwoptCube );
		}	
			else 
		{
			myModel = new ARM_InfBSModel( asOfDate, zeroCurveModel, fwdInfModel, volCurve, 
			correlManager, IRBSModel, infSwoptCurve, IRSwoptCurve );
		};

		

		/// assign object
		if( !assignObject( myModel, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
		
	catch(Exception& x)
	{
		delete myModel;
		
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


long ARMLOCAL_infBSSmiledModel(		double		C_AsOfDate,
									long		C_DiscountCurvId,
									long		C_InfFwdCurvId,
									long		C_VolSigmaId,
									long		C_VolNuId,
									long		C_VolRhoId,
									long		C_VolBetaId,
									long		C_VolAtmIrId,
									long		C_CorrelId,
									long		C_CorrelAdjId,
									ARM_result& result, 
									long		objId ){
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_InfBSSmiledModel* myModel = NULL;

	try
	{
		char  charDate[11];
		Local_XLDATE2ARMDATE( C_AsOfDate, charDate );
		ARM_Date asOfDate( charDate );

		ARM_ZeroCurve* zeroCurve = NULL;
		if( !GetObjectFromId( &zeroCurve, C_DiscountCurvId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: forecast fwd curve is not of a good type");
			return ARM_KO;
		};

		ARM_InfCurv* fwdInfModel = NULL;
		if( !GetObjectFromId( &fwdInfModel, C_InfFwdCurvId, ARM_INFCURV ) )
		{
			result.setMsg ("ARM_ERR: forward inf curve is not of a good type");
			return ARM_KO;
		};

		ARM_VolCurve* volSigma = NULL;
		if( !GetObjectFromId( &volSigma, C_VolSigmaId, ARM_VOL_CURVE ) )
		{
			result.setMsg ("ARM_ERR Sigma: volatility curve is not of a good type");
			return ARM_KO;
		};

		ARM_VolCurve* volRho = NULL;
		if( !GetObjectFromId( &volRho, C_VolRhoId, ARM_VOL_CURVE ) )
		{
			result.setMsg ("ARM_ERR Rho: volatility curve is not of a good type");
			return ARM_KO;
		};

		ARM_VolCurve* volNu = NULL;
		if( !GetObjectFromId( &volNu, C_VolNuId, ARM_VOL_CURVE ) )
		{
			result.setMsg ("ARM_ERR Nu: volatility curve is not of a good type");
			return ARM_KO;
		};

		ARM_VolCurve* volBeta = NULL;
		if( !GetObjectFromId( &volBeta, C_VolBetaId, ARM_VOL_CURVE ) )
		{
			result.setMsg ("ARM_ERR Beta: volatility curve is not of a good type");
			return ARM_KO;
		};

		ARM_VolCurve* volAtmIr = NULL;
		if( !GetObjectFromId( &volAtmIr, C_VolAtmIrId, ARM_VOL_CURVE ) )
		{
			result.setMsg ("ARM_ERR Atm ir: volatility curve is not of a good type");
			return ARM_KO;
		};

		ARM_VolCurve* correl = NULL;
		if( !GetObjectFromId( &correl, C_CorrelId, ARM_VOL_CURVE ) )
		{
			result.setMsg ("ARM_ERR Correl: volatility curve is not of a good type");
			return ARM_KO;
		};

		ARM_VolCurve* correlAdj = NULL;
		GetObjectFromId( &correlAdj, C_CorrelAdjId, ARM_VOL_CURVE ) ;

		myModel = new ARM_InfBSSmiledModel( asOfDate, zeroCurve, fwdInfModel, volAtmIr, correl, correlAdj, volSigma, volNu, volRho, volBeta);
	
		/// assign object
		if( !assignObject( myModel, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
		
	catch(Exception& x)
	{
		delete myModel;
		
		x.DebugPrint();
		
		ARM_RESULT();
	}
}



/* !
 * function to create a year to year inflation leg
 * 
 */
extern long ARMLOCAL_InfCapFloor_Create(
	double startDate,
	double endDate,
	const CCString& indexName,
	int capOrFloor,
	double strike,
	double leverage,
	double spread,
	int swapType,
	int rcvOrPay,
	int interpType,
	int resetFreq,
	int dayCount,
	const CCString& resetCalendar,
	int fwdRule,
	int intRule,
	int stubRule,
	int resetGap,
	int payFreq,
	int payGap,
	const CCString& payCalendar,
	int adjFirstDate,
	double firstReset,
	long ccyId,
	ARM_result&	result,
	long objId 	)
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_InfCapFloor* infCapFloor = NULL;
	bool createdNewCcy	= false;
	ARM_Currency* ccy	= NULL;

	try
	{
		/// gets the currency 
		/// get the defaulted one if not provided
		if ( ccyId == ARM_NULL_OBJECT )
		{
			createdNewCcy = true;
			ccy = new ARM_Currency( InfData::GetCurrency( indexName.c_str() ) );
		}
		else
		{
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}

		/// dateConversion
		char  startDateChar[11];
		char  endDateChar[11];

		Local_XLDATE2ARMDATE( startDate, startDateChar );
		Local_XLDATE2ARMDATE( endDate, endDateChar );

		char* resetCal	= resetCalendar.c_str();
		char* payCal	= payCalendar.c_str();

		resetGap		= ARM_InfLeg::GetGapFromGapOrDate( GetGapOrJulianDate(resetGap), resetCal, (ARM_Date) startDateChar);
		payGap			= ARM_InfLeg::GetGapFromGapOrDate( GetGapOrJulianDate(payGap), payCal,	(ARM_Date) startDateChar);

		/// new inflation index
		infCapFloor	= new ARM_InfCapFloor( 
			(ARM_Date) startDateChar,
			(ARM_Date) endDateChar,
			CCSTringToSTLString( indexName ),
			capOrFloor,
			strike,
			leverage,
			spread,
			swapType,
			rcvOrPay,
			interpType,
			resetFreq,
			dayCount,
			resetCal,
			fwdRule,
			intRule,
			stubRule,
			resetGap,
			payFreq,
			payGap,
			payCal,
			adjFirstDate,
			firstReset );

		delete resetCal;
		delete payCal;

		/// assign object
		if( !assignObject( infCapFloor, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		if( createdNewCcy )
			delete ccy;

		delete infCapFloor;
		
		x.DebugPrint();
		
		ARM_RESULT();
	}

}


/*! 
 * function to create a correlation matrix
 */
extern long ARMLOCAL_CorrelMat_Create(
	double	asOfDate,
	const VECTOR<CCString>& X,
	const VECTOR<CCString>& Y,
	const VECTOR<double>&	Z,
	ARM_result& result,
	long objId )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");

	ARM_CorrelMatrix* correlMatrix = NULL;
	ARM_Vector* vX;
	ARM_Vector* vY;
	ARM_Matrix*	vZ;

	try
	{
		char startDate[11];
		Local_XLDATE2ARMDATE( asOfDate,startDate);
		ARM_Date asOfARMDate( startDate );

		vX		= CreateARMMatuVectorFromStrVECTOR(X);
		vY		= CreateARMMatuVectorFromStrVECTOR(Y);
		vZ		= CreateARMMatrixFromVECTOR( Z, vX->GetSize(), vY->GetSize() );

		correlMatrix = new ARM_CorrelMatrix( asOfARMDate, vX, vY, vZ );

		if( !assignObject( correlMatrix, result, objId ) )
		{
			delete vX;
			delete vY;
			delete vZ;
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}

	catch(Exception& x)
	{
		delete vX;
		delete vY;
		delete vZ;
		delete correlMatrix;

		x.DebugPrint();

		ARM_RESULT();
	}
}




/*! 
 * function to create a correlation manager
 */
extern long ARMLOCAL_CorrelManager_Create(
	const CCString& mktTag, 
	const CCString& intraMktTag,
	double	asOfDate,
	const VECTOR<CCString>& X,
	const VECTOR<CCString>& Y,
	const VECTOR<double>& Z,
	ARM_result& result,
	long objId )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");

	ARM_CorrelManager* correlManager = NULL;
	ARM_Vector* vX;
	ARM_Vector* vY;
	ARM_Matrix*	vZ;

	try
	{
		char startDate[11];
		Local_XLDATE2ARMDATE(asOfDate,startDate);
		ARM_Date asOfARMDate( startDate );

		vX		= CreateARMMatuVectorFromStrVECTOR(X);
		vY		= CreateARMMatuVectorFromStrVECTOR(Y);
		vZ		= CreateARMMatrixFromVECTOR( Z, vX->GetSize(), vY->GetSize() );

		correlManager = new ARM_CorrelManager( CCSTringToSTLString( mktTag ), CCSTringToSTLString( intraMktTag ), asOfARMDate, vX, vY, vZ );

		if( !assignObject( correlManager, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}

	catch(Exception& x)
	{
		delete correlManager;

		x.DebugPrint();

		ARM_RESULT();
	}
}






/*! 
 * function to create a correlation manager
 */
extern long ARMLOCAL_CorrelManager_CreateFromMat(
	const CCString& mktTag, 
	const CCString& intraMktTag,
	long CorrelMatId,
	ARM_result& result,
	long objId )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");
	ARM_CorrelManager* correlManager = NULL;

	try
	{
		ARM_CorrelMatrix* correlMat = NULL;
		if( !GetObjectFromId( &correlMat, CorrelMatId, ARM_CORRELMATRIX ) )
		{
			result.setMsg ("ARM_ERR: correl matrix is not of a good type");
			return ARM_KO;
		};

		correlManager = new ARM_CorrelManager( CCSTringToSTLString( mktTag ), CCSTringToSTLString( intraMktTag ), correlMat );

		if( !assignObject( correlManager, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}

	catch(Exception& x)
	{
		delete correlManager;

		x.DebugPrint();

		ARM_RESULT();
	}
}



long ARMLOCAL_CreateGenCorrelManager(VECTOR<CCString>& mktTags,
								     VECTOR<CCString>& intraMktTags,
									 vector<long>& correlVolIds, 
									 ARM_result& result, 
									 long objId)
{
	long CorrelManagerId;

	int mkt_size = mktTags.size ();
	int intraMkt_size = intraMktTags.size ();
	int correlId_size = correlVolIds.size();

    ARM_VolLInterpol* CorrelCurve = NULL;
    ARM_VolFlat* CorrelFlat = NULL;

	ARM_CorrelManager* CorrelManager = NULL;
	ARM_CorrelManager* CorrelManagerOrg = NULL;

	//verify the size
	if(! (mkt_size==intraMkt_size && mkt_size==correlId_size) || mkt_size == 0)
	{
		result.setMsg ("ARM_ERR: check your matrix dimension");
		return ARM_KO;		
	}


	CCString msg ("");

	try
	{

		char vMktTags[200][20];
		char vIntraMktTags[200][20];


		if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
		{
			result.setMsg ("ARM_ERR: Pb with accessing objects");
			return ARM_KO;
		}

		////////// Create CorrelManager ////////////
		sprintf(vMktTags[0], (const char *) mktTags[0]);
		sprintf(vIntraMktTags[0], (const char*) intraMktTags[0]);


		CorrelCurve = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(correlVolIds[0]);
		
		if ( CorrelCurve->GetName() != ARM_VOL_LIN_INTERPOL ) // Flat Curve
		{
			CorrelFlat = (ARM_VolFlat*) LOCAL_PERSISTENT_OBJECTS->GetObject(correlVolIds[0]);

			ARM_CorrelMatrix* mat = new ARM_CorrelMatrix(CorrelFlat);

			CorrelManager= new ARM_CorrelManager(vMktTags[0], vIntraMktTags[0], mat);

			delete mat;
		}
		else  //Curve interpol
		{
			ARM_CorrelMatrix* mat = new ARM_CorrelMatrix(CorrelCurve);

			CorrelManager = new ARM_CorrelManager(vMktTags[0], vIntraMktTags[0], mat);
					
			delete mat;
		}

		//////////Fill CorrelManager//////////////
		for(int i = 1; i < mkt_size; i++)
		{
			sprintf(vMktTags[i], (const char*)mktTags[i]);
			sprintf(vIntraMktTags[i], (const char*)intraMktTags[i]);

			CorrelCurve = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(correlVolIds[i]);

			if ( CorrelCurve->GetName() != ARM_VOL_LIN_INTERPOL )
			{
				CorrelFlat = (ARM_VolFlat*) LOCAL_PERSISTENT_OBJECTS->GetObject(correlVolIds[i]);
			
                if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(CorrelCurve, ARM_VOL_CURVE) == 0 )
				{
					result.setMsg("ARM_ERR: CorrelCurve  is not of the right type");

					return(ARM_KO);
				}

				ARM_CorrelMatrix* mat = new ARM_CorrelMatrix(CorrelFlat);

				CorrelManager->Fill(vMktTags[i], vIntraMktTags[i], mat);

				delete mat;
			}
			else
			{
				ARM_CorrelMatrix* mat = new ARM_CorrelMatrix(CorrelCurve);

				CorrelManager->Fill(vMktTags[i], vIntraMktTags[i], mat);

				delete mat;
			}
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			CorrelManagerId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) CorrelManager);

			if ( CorrelManagerId == RET_KO )
			{
				if (CorrelManager)
					delete CorrelManager;
				CorrelManager = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");				
				
                return(ARM_KO);
			}

			result.setLong(CorrelManagerId);

			return ARM_OK;			
		}
		else
		{
			CorrelManagerOrg = (ARM_CorrelManager *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(CorrelManagerOrg, ARM_CORRELMANAGER) == 1)
			{
				if (CorrelManagerOrg)
				{
					delete CorrelManagerOrg;
					CorrelManagerOrg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*) CorrelManager, objId);

				return(ARM_OK);
			}
			else
			{
				if (CorrelManager)
					delete CorrelManager;
				CorrelManager = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}
		
	catch(Exception& x)
	{
/*
		if (CorrelCurve)
  		   delete CorrelCurve;
		CorrelCurve = NULL;

		if (CorrelFlat)
			delete CorrelFlat;
		CorrelFlat = NULL;
*/

		if (CorrelManager)
			delete CorrelManager;
		CorrelManager = NULL;
/*
		if (CorrelManagerOrg)
			delete CorrelManagerOrg;
		CorrelManagerOrg = NULL;
*/

		ARM_RESULT();
	}

}


/*! 
 * function to fill a correlation manager
 */
extern long ARMLOCAL_CorrelManager_Fill(
	const CCString& mktTag,
	const CCString& intraMktTag,
	const VECTOR<CCString>& X,
	const VECTOR<CCString>& Y,
	const VECTOR<double>& Z,
	long correlManagerId,
	ARM_result& result,
	long objId )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");

	ARM_CorrelManager* correlManager = NULL;
	ARM_Vector* vX;
	ARM_Vector* vY;
	ARM_Matrix*	vZ;

	try
	{
		vX		= CreateARMMatuVectorFromStrVECTOR(X);
		vY		= CreateARMMatuVectorFromStrVECTOR(Y);
		vZ		= CreateARMMatrixFromVECTOR( Z, vX->GetSize(), vY->GetSize() );

		ARM_CorrelManager* prevCorrelManager = NULL;
		if( !GetObjectFromId( &prevCorrelManager, correlManagerId, ARM_CORRELMANAGER ) )
		{
			result.setMsg ("ARM_ERR: correl manager is not of a good type");
			return ARM_KO;
		};

		correlManager = (ARM_CorrelManager*) prevCorrelManager->Clone();

		correlManager->Fill( CCSTringToSTLString( mktTag ), CCSTringToSTLString( intraMktTag ), vX, vY, vZ );
	
		if( !assignObject( correlManager, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}

	catch(Exception& x)
	{
		delete correlManager;

		x.DebugPrint();

		ARM_RESULT();
	}
}



/*! 
 * function to fill a correlation manager
 */
extern long ARMLOCAL_CorrelManager_FillFromMat(
	const CCString& mktTag,
	const CCString& intraMktTag,
	long CorrelMatId,
	long correlManagerId,
	ARM_result& result,
	long objId )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");

	ARM_CorrelManager* correlManager = NULL;

	try
	{
		ARM_CorrelMatrix* correlMat = NULL;
		if( !GetObjectFromId( &correlMat, CorrelMatId, ARM_CORRELMATRIX ) )
		{
			result.setMsg ("ARM_ERR: correl matrix is not of a good type");
			return ARM_KO;
		};

		ARM_CorrelManager* prevCorrelManager = NULL;
		if( !GetObjectFromId( &prevCorrelManager, correlManagerId, ARM_CORRELMANAGER ) )
		{
			result.setMsg ("ARM_ERR: correl manager is not of a good type");
			return ARM_KO;
		};

		correlManager = (ARM_CorrelManager*) prevCorrelManager->Clone();

		correlManager->Fill( CCSTringToSTLString( mktTag ), CCSTringToSTLString( intraMktTag ), correlMat );
	
		if( !assignObject( correlManager, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}

	catch(Exception& x)
	{
		delete correlManager;

		x.DebugPrint();

		ARM_RESULT();
	}
}



/*!
 * Function to compute a correlation with the
 * correlation manager
 */
extern long ARMLOCAL_ComputeCorrelFromCorrelManager(
	const CCString& tmpMktTag,
	const CCString& tmpIntraMktTag,
	double x,
	double y,
	long correlManagerId,
	ARM_result& result )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");

	ARM_CorrelManager* correlManager = NULL;

	try
	{
		ARM_CorrelManager* correlManager = NULL;
		if( !GetObjectFromId( &correlManager, correlManagerId, ARM_CORRELMANAGER ) )
		{
			result.setMsg ("ARM_ERR: correl manager is not of a good type");
			return ARM_KO;
		};

		string mktTag = CCSTringToSTLString( tmpMktTag );
		string intraMktTag = CCSTringToSTLString( tmpIntraMktTag );
		double dResult = correlManager->ComputeCorrelData( mktTag, intraMktTag, x, y);
		result.setDouble(dResult);
		return ARM_OK;
	}


	catch(Exception& x)
	{
		delete correlManager;

		x.DebugPrint();

		ARM_RESULT();
	}
}








/*!
 * Function to compute a correlation curve from 
 * zero coupon and year to year volatility curves
 */

extern long ARMLOCAL_Vol_to_Cor(
		const VECTOR<double>& ZCVolV, 
		const VECTOR<double>& YtYVolV, 
		const VECTOR<double>& MaturityV,
		VECTOR<double>& CorV,
		ARM_result& result )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");
	std::vector<double>& pZCVol		= NULL;	
	std::vector<double>& pYtYVol		= NULL;	
	std::vector<double>& pMaturity	= NULL;	


	try
	{
		pZCVol		= CreateARMGPVectorFromVECTOR( ZCVolV );
		pYtYVol		= CreateARMGPVectorFromVECTOR( YtYVolV );
		pMaturity	= CreateARMGPVectorFromVECTOR( MaturityV );
		std::vector<double> Cor;

		ARM_MarketDataValidator mktValidator;
		mktValidator.Vol_to_Cor( *pZCVol, *pYtYVol, *pMaturity, Cor );

		CorV.resize( Cor.size() );
		for( int i=0; i<Cor.size(); ++i )
		{
			CorV[i] = Cor[i];
		}

		/// to avoid memory leak, delete the ARM_Vector*
		delete pZCVol;
		delete pYtYVol;
		delete pMaturity;

		return ARM_OK;
	}

	catch(Exception& x)
	{
		delete pZCVol;
		delete pYtYVol;
		delete pMaturity;

		x.DebugPrint();

		ARM_RESULT();
	}
}


/*!
 * Function to compute a zero coupon volatility curve from 
 * a correlation and a year to year volatility curves
 */

extern long ARMLOCAL_YtYCor_to_ZC(
		const VECTOR<double>& YtYVolV,
		const VECTOR<double>& CorV,
		const VECTOR<double>& MaturityV,
		VECTOR<double>& ZCVolV,
		ARM_result& result )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");
	std::vector<double>& pCor		= NULL;	
	std::vector<double>& pYtYVol		= NULL;	
	std::vector<double>& pMaturity	= NULL;	




	try
	{
		pCor		= CreateARMGPVectorFromVECTOR( CorV );
		pYtYVol		= CreateARMGPVectorFromVECTOR( YtYVolV );
		pMaturity	= CreateARMGPVectorFromVECTOR( MaturityV );
		std::vector<double> ZCVol;
		

		ARM_MarketDataValidator mktValidator;
		mktValidator.YtYCor_to_ZC( *pYtYVol, *pCor, *pMaturity, ZCVol );

		ZCVolV.resize( ZCVol.size() );
		for( int i=0; i<ZCVol.size(); ++i )
		{
			ZCVolV[i] = ZCVol[i];
		}

		/// to avoid memory leak, delete the std::vector<double>&
		delete pCor;
		delete pYtYVol;
		delete pMaturity;

		return ARM_OK;
	}

	catch(Exception& x)
	{
		delete pCor;
		delete pYtYVol;
		delete pMaturity;

		x.DebugPrint();

		ARM_RESULT();
	}
}


/*!
 * Function to compute a year to year volatility curve from 
 * a correlation and a zero coupon volatility curves
 */

extern long ARMLOCAL_ZCCor_to_YtY(
		const VECTOR<double>& ZCVolV,
		const VECTOR<double>& CorV,
		const VECTOR<double>& MaturityV,
		VECTOR<double>& YtYVolV,
		ARM_result& result )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");
	std::vector<double>& pCor		= NULL;	
	std::vector<double>& pZCVol		= NULL;	
	std::vector<double>& pMaturity	= NULL;	


	try
	{
		pCor		= CreateARMGPVectorFromVECTOR( CorV );
		pZCVol		= CreateARMGPVectorFromVECTOR( ZCVolV );
		pMaturity	= CreateARMGPVectorFromVECTOR( MaturityV );
		std::vector<double> YtYVol;
		
		ARM_MarketDataValidator mktValidator;
		mktValidator.ZCCor_to_YtY( *pZCVol, *pCor, *pMaturity, YtYVol );

		YtYVolV.resize( YtYVol.size() );
		for( int i=0; i<YtYVol.size(); ++i )
		{
			YtYVolV[i] = YtYVol[i];
		}

		/// to avoid memory leak, delete the ARM_Vector*
		delete pCor;
		delete pZCVol;
		delete pMaturity;

		return ARM_OK;
	}

	catch(Exception& x)
	{
		delete pCor;
		delete pZCVol;
		delete pMaturity;

		x.DebugPrint();

		ARM_RESULT();
	}
}


/*!
 * Function to provide the confidence interval for
/* zero coupon & year to year volatilities & correlations
 */



extern long ARMLOCAL_Bounds(
		const VECTOR<double>& ZCVolV,
		const VECTOR<double>& YtYVolV,
		const VECTOR<double>& CorV,
		const VECTOR<double>& MaturityV,
		VECTOR<double>& UBoundV,
		VECTOR<double>& LBoundV,
		const string& choice, 
		string& TBound, 
		ARM_result& result )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");
	std::vector<double>& pCor		= NULL;	
	std::vector<double>& pZCVol		= NULL;	
	std::vector<double>& pMaturity	= NULL;	
	std::vector<double>& pYtYVol		= NULL;	
	


	try
	{
		pCor		= CreateARMGPVectorFromVECTOR( CorV );
		pZCVol		= CreateARMGPVectorFromVECTOR( ZCVolV );
		pYtYVol		= CreateARMGPVectorFromVECTOR( YtYVolV );
		pMaturity	= CreateARMGPVectorFromVECTOR( MaturityV );
		std::vector<double> UBound(1);
		std::vector<double> LBound(1);
		
		
		ARM_MarketDataValidator mktValidator;

		mktValidator.Bounds( pZCVol, pYtYVol, pCor, pMaturity, UBound, LBound, choice, TBound );

		
		UBoundV.resize( UBound.size() );
		int i;
		for( i=0; i<UBound.size(); ++i )
		{
			UBoundV[i] = UBound[i];
		}

		LBoundV.resize( LBound.size() );
		for( i=0; i<LBound.size(); ++i )
		{
			LBoundV[i] = LBound[i];
		}


		ARM_Vector temp(3*UBound.size());
		
		for( i=0; i<LBound.size(); ++i )
		{
			temp[i] = LBound[i];
			temp[i+LBound.size()] = UBound[i];

			if(TBound=="internal") {temp[i+2*LBound.size()]=1.0;} else if(TBound=="external") {temp[i+2*LBound.size()]=0.0;};
		}


		result.setLong(3*UBound.size());

		for ( i=0;i<3*UBound.size();i++)
		{
			result.setArray(temp[i],i);
		}


		/// to avoid memory leak, delete the ARM_Vector*
		delete pCor;
		delete pZCVol;
		delete pYtYVol;
		delete pMaturity;

		return ARM_OK;
	}

	catch(Exception& x)
	{
		delete pCor;
		delete pZCVol;
		delete pYtYVol;
		delete pMaturity;

		x.DebugPrint();

		ARM_RESULT();
	}
}


/*!
 * Function to compute a zero coupon volatility curve from 
 * a correlation and a year to year volatility curves in 
 * the homogeneous case
 */

extern long ARMLOCAL_HmgVol_to_Cor(
		const VECTOR<double>& ZCVolV,
		const VECTOR<double>& YtYVolV,
		const VECTOR<double>& MaturityV,
		VECTOR<double>& CorV,
		const int length,
		ARM_result& result )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");
	std::vector<double>& pZCVol		= NULL;	
	std::vector<double>& pYtYVol		= NULL;	
	std::vector<double>& pMaturity	= NULL;

	try
	{
		pZCVol		= CreateARMGPVectorFromVECTOR( ZCVolV );
		pYtYVol		= CreateARMGPVectorFromVECTOR( YtYVolV );
		pMaturity	= CreateARMGPVectorFromVECTOR( MaturityV );
		std::vector<double> Cor(1);
		ARM_MarketDataValidator mktValidator;
		mktValidator.HmgVol_to_Cor( pZCVol, pYtYVol, pMaturity, length, Cor);

		CorV.resize( Cor.size() );
		for( int i=0; i<Cor.size(); ++i )
		{
			CorV[i] = Cor[i];
		}

		/// to avoid memory leak, delete the ARM_Vector*		
		delete pZCVol;
		delete pYtYVol;
		delete pMaturity;
		return ARM_OK;
	}

	catch(Exception& x)
	{
		delete pZCVol;
		delete pYtYVol;
		delete pMaturity;

		x.DebugPrint();
		ARM_RESULT();
	}
	

}


/*!
 * Function to compute a year to year volatility curve from 
 * a correlation and a zero coupon volatility curves
 * in the homogeneous case
 */

extern long ARMLOCAL_HmgZCCor_to_YtY(
		const VECTOR<double>& ZCVolV,
		const VECTOR<double>& CorV,
		const VECTOR<double>& MaturityV,
		const int length,
		VECTOR<double>& YtYVolV,
		ARM_result& result )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");
	std::vector<double>& pCor		= NULL;	
	std::vector<double>& pZCVol		= NULL;	
	std::vector<double>& pMaturity	= NULL;	


	try
	{
		pCor		= CreateARMGPVectorFromVECTOR( CorV );
		pZCVol		= CreateARMGPVectorFromVECTOR( ZCVolV );
		pMaturity	= CreateARMGPVectorFromVECTOR( MaturityV );
		std::vector<double> YtYVol;
		
		ARM_MarketDataValidator mktValidator;
		mktValidator.HmgZCCor_to_YtY( pZCVol, pCor, pMaturity, length, YtYVol );

		YtYVolV.resize( YtYVol.size() );
		for( int i=0; i<YtYVol.size(); ++i )
		{
			YtYVolV[i] = YtYVol[i];
		}

		/// to avoid memory leak, delete the ARM_Vector*
		delete pCor;
		delete pZCVol;
		delete pMaturity;

		return ARM_OK;
	}

	catch(Exception& x)
	{
		delete pCor;
		delete pZCVol;
		delete pMaturity;

		x.DebugPrint();

		ARM_RESULT();
	}
}

/*!
 * Function to compute a zero coupon volatility curve from 
 * a correlation and a year to year volatility curves
 * in the homogeneous case
 */

extern long ARMLOCAL_HmgYtYCor_to_ZC(
		const VECTOR<double>& YtYVolV, 
		const VECTOR<double>& CorV, 
		const VECTOR<double>& MaturityV,
		const int length,
		VECTOR<double>& ZCVolV,
		ARM_result& result )
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");
	std::vector<double>& pCor		= NULL;	
	std::vector<double>& pYtYVol		= NULL;	
	std::vector<double>& pMaturity	= NULL;	

	try
	{
		pCor		= CreateARMGPVectorFromVECTOR( CorV );
		pYtYVol		= CreateARMGPVectorFromVECTOR( YtYVolV );
		pMaturity	= CreateARMGPVectorFromVECTOR( MaturityV );
		std::vector<double> ZCVol;
		

		ARM_MarketDataValidator mktValidator;
		mktValidator.HmgYtYCor_to_ZC( pYtYVol, pCor, pMaturity, length, ZCVol );

		ZCVolV.resize( ZCVol.size() );
		for( int i=0; i<ZCVol.size(); ++i )
		{
			ZCVolV[i] = ZCVol[i];
		}

		/// to avoid memory leak, delete the ARM_Vector*
		delete pCor;
		delete pYtYVol;
		delete pMaturity;

		return ARM_OK;
	}

	catch(Exception& x)
	{
		delete pCor;
		delete pYtYVol;
		delete pMaturity;

		x.DebugPrint();

		ARM_RESULT();
	}
}



/*!
 * Function to compute the inflation swap's  
 * implied volatility from the year-on-year ones 
 */

extern long ARMLOCAL_VolYoY_to_VolSwp(
		const VECTOR<double>&  DFactor, 
		const VECTOR<double>&  FwdCPI, 
		const VECTOR<double>&  Vol_DF,
		const VECTOR<double>&  Vol_YoY, 
		const VECTOR<double>&  AvgCor, 
		const VECTOR<double>&  Dates, 
		const VECTOR<double>&  Tenors, 
		const double SwpRate, 
		double &Vol_Swp,
		ARM_result& result )
{
	
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");
	std::vector<double>& pDFactor	= NULL;	
	std::vector<double>& pFwdCPI		= NULL;	
	std::vector<double>& pVol_DF		= NULL;	
	std::vector<double>& pVol_YoY	= NULL;	
	std::vector<double>& pAvgCor		= NULL;	
	std::vector<double>& pDates		= NULL;	
	std::vector<double>& pTenors		= NULL;	


	try
	{
		pDFactor	= CreateARMGPVectorFromVECTOR( DFactor );
		pFwdCPI		= CreateARMGPVectorFromVECTOR( FwdCPI );
		pVol_DF		= CreateARMGPVectorFromVECTOR( Vol_DF );
		pVol_YoY	= CreateARMGPVectorFromVECTOR( Vol_YoY );
		pAvgCor		= CreateARMGPVectorFromVECTOR( AvgCor );
		pDates		= CreateARMGPVectorFromVECTOR( Dates );
		pTenors		= CreateARMGPVectorFromVECTOR( Tenors );

		double Vol_SwpI;
		ARM_Convert_Vol  YoY_to_Swp;
		YoY_to_Swp.VolYoY_to_VolSwp( pDFactor, pFwdCPI, pVol_DF, pVol_YoY, pAvgCor, pDates, pTenors, SwpRate, Vol_SwpI);
		Vol_Swp = Vol_SwpI;

		/// to avoid memory leak, delete the ARM_Vector*
		delete pDFactor;
		delete pFwdCPI;
		delete pVol_DF;
		delete pVol_YoY;
		delete pAvgCor;
		delete pDates;
		delete pTenors;
		

		return ARM_OK;
	}

	catch(Exception& x)
	{
		delete pDFactor;
		delete pFwdCPI;
		delete pVol_DF;
		delete pVol_YoY;
		delete pAvgCor;
		delete pDates;
		delete pTenors;

		x.DebugPrint();

		ARM_RESULT();
	}
}






/*!
 * Function to create an inflation curve
 */
extern long ARMLOCAL_SeasonalityManager_Create(
			const VECTOR<CCString>&	monthsList,
			const VECTOR<double>&	seasonSpreadList, 
			const VECTOR<double>&	seasonHorizonList, 
			long SeasonalityCorrectionType,
			ARM_result&				result,
			long					objId
			)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_SeasonalityManager* seasonManager = NULL;
	vector<string> vMonthList;
	std::vector<double>& vSeasonSpreadList  = NULL;
	vector<double> vSeasonHorizonList(0);

	try
	{
		
		int i = 0;
		for( i = 0; i < monthsList.size(); i++)
		{
			vMonthList.push_back(CCSTringToSTLString(monthsList[i]));
		}

		if (seasonHorizonList.size())
		{
			vSeasonHorizonList = seasonHorizonList;
			vSeasonSpreadList = CreateARMGPVectorFromVECTOR(seasonSpreadList);
			seasonManager = new ARM_SeasonalityManager(
				vMonthList, 
				*vSeasonSpreadList,
				vSeasonHorizonList,
				(ARM_SeasonalityManager::CorrectionMode) SeasonalityCorrectionType );
		}
		else
		{
			vSeasonSpreadList = CreateARMGPVectorFromVECTOR(seasonSpreadList);
			seasonManager = new ARM_SeasonalityManager( vMonthList, 
			*vSeasonSpreadList, (ARM_SeasonalityManager::CorrectionMode) SeasonalityCorrectionType );
		}

		delete vSeasonSpreadList;
		
		/// assign object
		if( !assignObject( seasonManager, result, objId ) )
		{
			delete seasonManager;
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		delete vSeasonSpreadList;
		delete seasonManager;

		x.DebugPrint();
		
		ARM_RESULT();
	}
}


ARM_Date searchCPI(const CCString& index, double date, double& CPIVAL)
{
	// on récupère la liste des CPI de reference en fonction de l'index

	int id;
    int size = sizeof(RefIFRF)/sizeof(RefIFRF[0]);

	if (strcmp((const char*)index,"EMU") == 0)
		id = 1;
	else if (strcmp((const char*)index,"EMUT") == 0)
		id = 2;
	else if (strcmp((const char*)index,"IFRF") == 0)
		id = 3;
	else
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Wrong Inflation Index. Check parameters ..." );

		return ARM_KO;
	}

	// On peut faire ce test car la prmière date des tableaux est identique
	if (date < RefIFRF[0].ZCDate)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Invalid date for searchin CPI");
	}

	int i = 1;
	double CPIDate;

	switch (id)
	{
	case 1:
		while ( (i < size) && (date > RefEMU[i].ZCDate))
			i++;

		CPIVAL = RefEMU[i-1].CPIValue;

		CPIDate = RefEMU[i-1].CPIDate;

		break;

	case 2:
		while ( (i < size) && (date > RefEMUT[i].ZCDate))
			i++;

		CPIVAL = RefEMUT[i-1].CPIValue;

		CPIDate = RefEMUT[i-1].CPIDate;

		break;

	case 3:
		while ( (i < size) && (date > RefIFRF[i].ZCDate))
			i++;

		CPIVAL = RefIFRF[i-1].CPIValue;

		CPIDate = RefIFRF[i-1].CPIDate;

		break;
	}

	char sDate[11];

	Local_XLDATE2ARMDATE(CPIDate,sDate);

	ARM_Date tmpDate(sDate);

	return tmpDate;
}


ARM_InfCurv* ARMLOCAL_GetInfZcFromSummitFF_Create(
			const CCString&	index,
			const CCString&	ccy,
			const CCString&	cvname,
			double date,
			ARM_result& result)
{
	ARM_InfCurv* infCurv = NULL;
	FILE *Fp = NULL;
	vector<string> MktTerms;
	vector<double> MktValues;

	try
	{
		char sAsOfDate[11];
		Local_XLDATE2ARMDATE(date,sAsOfDate);
		ARM_Date asOfARMDate( sAsOfDate );

		ARM_Date dateToday;

		if (asOfARMDate > dateToday)
		{
			result.setMsg ("ARM_ERR: Invalid Date");
			return NULL;
		}

		CCString myRepertory;
		CCString myRepertory2;

		CCString FileName, FileName2;

		myRepertory = ((CCString)(armlocal_init->data_folder.c_str()) + YLD_SUMMIT_FILE_LOCATION);
		myRepertory2 = ((CCString)(armlocal_init->data_folder.c_str()) + YLD_SUMMIT_FILE_LOCATION + SUMMIT_FILE_LOCATION_HISTO);

		FileName = myRepertory + (CCString)"YLD_" + ccy + "_" + index + "_" + cvname + ".";
		FileName2 = myRepertory2 + (CCString)"YLD_" + ccy + "_" + index + "_" + cvname + ".";

		char sEch[50];
		char buffer[50];

		if (strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0)
		{
			_ltoa(asOfARMDate.GetYear(),buffer,10);
			FileName += (CCString) (buffer);
			FileName2 += (CCString) (buffer);

			_ltoa(asOfARMDate.GetMonth(),buffer,10);
			if (asOfARMDate.GetMonth() < 10)
			{
				FileName += "0";
				FileName2 += "0";
			}
			FileName += (CCString) buffer;
			FileName2 += (CCString) buffer;

			_ltoa(asOfARMDate.GetDay(),buffer,10);
			if (asOfARMDate.GetDay() < 10)
			{
				FileName += "0";
				FileName2 += "0";
			}
			FileName += (CCString) buffer;
			FileName2 += (CCString) buffer;
		}
		else
		{
			_ltoa(asOfARMDate.GetDay(),buffer,10);
			if (asOfARMDate.GetDay() < 10)
			{
				FileName += "0";
				FileName2 += "0";
			}
			FileName += (CCString) buffer;
			FileName2 += (CCString) buffer;

			_ltoa(asOfARMDate.GetMonth(),buffer,10);
			if (asOfARMDate.GetMonth() < 10)
			{
				FileName += "0";
				FileName2 += "0";
			}
			FileName += (CCString) buffer;
			FileName2 += (CCString) buffer;

			_ltoa(asOfARMDate.GetYear(),buffer,10);
			FileName += (CCString) (buffer + 2);
			FileName2 += (CCString) (buffer + 2);
		}

		FileName += (CCString) ".000";
		FileName2 += (CCString) ".000";

		double val;
		int rc = 0;

		// latest curves
		if ((Fp = fopen(FileName,"r")) == NULL)
		{
			// historical curves
			if ((Fp = fopen(FileName2,"r")) == NULL)
			{
				result.setMsg( CCString("ARM_ERR: Could not open the following files:\n" )
							+ FileName + "\n" + FileName2 + "\nCheck parameters ..." );
				return NULL;
			}
		}

		ARM_Date CPIIndexDate;
		double CPIIndex;
		int j = 0;

		while (rc != EOF)
		{
			rc = fscanf(Fp, "%s",buffer);
			rc = fscanf(Fp, "%s",sEch);

			if(rc != EOF )
			{
				if (strlen(sEch) == 8)
				{
					char sTmpDate[11];
					char sTmpDate1[7];
					strncpy(sTmpDate1,sEch,6);
					sTmpDate1[6] = '\0';
					sprintf(sTmpDate,"%s20%s",sTmpDate1,sEch+6);
					ARM_Date tmpDate(sTmpDate);
					sprintf( buffer, "%f", tmpDate.GetJulian() );
					MktTerms.push_back(CCSTringToSTLString(buffer));

					if (j == 0)
						CPIIndexDate = tmpDate;
				}
				else
				{
					if (strcmp((const char*)sEch,"2D") == 0)
					{
						// cas avant le 10/07/03 : on va chercher le CPI
						// de référence dans les structures en dur
						// on cherche le CPI de reference correspondant
						// a la date de la courbe
						CPIIndexDate = searchCPI(index,date,CPIIndex);
					}
					else
					{
						MktTerms.push_back(CCSTringToSTLString(sEch));
					}
				}
			}
			
			rc = fscanf(Fp, "%lf",&val);
			if(rc != EOF )
			{
				MktValues.push_back(val);
				if (j == 0)
					CPIIndex = val;
			}

			rc = fscanf(Fp, "%s",buffer);

			j++;
		}

		fclose(Fp);


		infCurv = new ARM_InfCurv(asOfARMDate,
								  CCSTringToSTLString(index),
								  CPIIndex,
								  CPIIndexDate,
								  MktTerms,
								  MktValues);
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();

		return NULL;
	}

	return infCurv;
}

/*!
 * Function to create an inflation curve
 * from Summit
 */
long ARMLOCAL_GetInfZcFromSummit_Create(
			const CCString&	index,
			const CCString&	ccy,
			const CCString&	cvname,
			double date,
			long seasonAdjId,
			long seasonAdjModeId,
			ARM_result&				result,
			long					objId)
{
	ARM_InfCurv* zc = NULL;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		char sDate[11];
		Local_XLDATE2ARMDATE(date,sDate);

		if (GetDataRetrieverVersion () >= ETKRETRIEVER)
		{
			CCString xmlResponse;

			ARM_Date myDate(sDate);

			xmlResponse = etoolkit_getXMLMYAndZCFromSummit(index,ccy,cvname,myDate);

			zc = ARMLOCAL_ParseXMLForInfZC(xmlResponse, date, index);

			if (seasonAdjId == K_YES)
			{
				ARM_SeasonalityManager* newSeasonManager = ARMLOCAL_ParseXMLForSeasonMgr(index,ccy,cvname,myDate,seasonAdjModeId);

				zc->SetSeasonalityManager(newSeasonManager);

				if (newSeasonManager)
					delete newSeasonManager;
				newSeasonManager = NULL;
			}
		}
		else
		{
/*			if (seasonAdjId == K_YES)
			{
				result.setMsg ("ARM_ERR: Seasonality Manager not implemented without ETK");
				return ARM_KO;
			}
*/
			ARM_Date asOfDate(sDate);
			zc = ARMLOCAL_GetInfZcFromSummitFF_Create(index,
													  ccy,
													  cvname,
													  date,
													  result);

			if (zc == NULL)
			{
				CCString xmlResponse;

				xmlResponse = etoolkit_getXMLMYAndZCFromSummit(index,ccy,cvname,asOfDate);

				zc = ARMLOCAL_ParseXMLForInfZC(xmlResponse, date, index);
			}

			if (seasonAdjId == K_YES)
			{
				ARM_SeasonalityManager* newSeasonManager = ARMLOCAL_ParseXMLForSeasonMgr(index,ccy,cvname,asOfDate,seasonAdjModeId);

				zc->SetSeasonalityManager(newSeasonManager);

				if (newSeasonManager)
					delete newSeasonManager;
				newSeasonManager = NULL;
			}
		}

		/// assign object
		if( !assignObject( zc, result, objId ) )
		{
			delete zc;
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}




//---------------------------------------------------

long ARMLOCAL_LIVRETACURVE(	double asOfDateDble ,
							long infCurvId,
							long euribCurvId,
							long flagArrondi,
							long infresetManagerId,
							long fixingLivretAId,
							long fixingEuribId,
							long monthForAugustId,
							long monthForFebruaryId,
							ARM_result& result,
							long objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_LivretACurve* myLivretACurve = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	try
	{
		char  charDate[11];
		Local_XLDATE2ARMDATE( asOfDateDble, charDate );
		ARM_Date asOfDate( charDate );

		if ( (monthForAugustId != K_MONTH_DEFAULT) &&
			 ( (monthForAugustId < K_MAY) ||
			   (monthForAugustId > K_JUNE)
			 )
		   )
		{
			result.setMsg ("ARM_ERR: Month for August reset must be between May and June or default");
			return ARM_KO;
		}

		if ( (monthForFebruaryId != K_MONTH_DEFAULT) &&
			 ( (monthForFebruaryId < K_NOVEMBER) ||
			   (monthForFebruaryId > K_DECEMBER)
			 )
		   )
		{
			result.setMsg ("ARM_ERR: Month for February reset must be between November and December or default");
			return ARM_KO;
		}

		ARM_InfCurv* infCurv = NULL;
		if( !GetObjectFromId( &infCurv, infCurvId, ARM_INFCURV) )
		{
			result.setMsg ("ARM_ERR: inflation curve is not of a good type");
			return ARM_KO;
		};

		ARM_ZeroCurve* euribCurv=NULL;
		if( !GetObjectFromId( &euribCurv, euribCurvId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: Euribor curve is not of a good type");
			return ARM_KO;
		}

		
		ARM_ReferenceValue* fixingEurib = NULL;
		ARM_ResetManager* euribresetManager = NULL;
		if(fixingEuribId!=-1)
		{
			if( !GetObjectFromId( &fixingEurib, fixingEuribId, ARM_REFERENCE_VALUE ) )
			{
				if( !GetObjectFromIdwNull( &euribresetManager, fixingEuribId, ARM_RESETMANAGER ) )
				{
					result.setMsg ("ARM_ERR: euribor fixing is not of a good type");
					return ARM_KO;
				}
			}
		}

		ARM_ReferenceValue* fixingLivretA = NULL;
		ARM_ResetManager* livretAresetManager = NULL;
		if(fixingLivretAId!=-1)
		{
			if( !GetObjectFromId( &fixingLivretA, fixingLivretAId, ARM_REFERENCE_VALUE ) )
			{
				if( !GetObjectFromIdwNull( &livretAresetManager, fixingLivretAId, ARM_RESETMANAGER ) )
				{
					result.setMsg ("ARM_ERR: livret A fixing is not of a good type");
					return ARM_KO;
				}
			}
		}
		
		ARM_ResetManager* infresetManager = NULL;
		if(infresetManagerId!=-1)
		{
			if( !GetObjectFromIdwNull( &infresetManager, infresetManagerId, ARM_RESETMANAGER ) )
			{
				result.setMsg ("ARM_ERR: reset manager is not of a good type");
				return ARM_KO;
			};
		}

		if ( (livretAresetManager == NULL) && (euribresetManager == NULL) )
		{
			myLivretACurve = new ARM_LivretACurve(	asOfDate, infCurv, euribCurv, flagArrondi, infresetManager,
													fixingLivretA, fixingEurib, monthForAugustId, monthForFebruaryId);
		}
		else if ( ((livretAresetManager != NULL) && (fixingEurib == NULL))
			|| ((euribresetManager != NULL) && (fixingLivretA == NULL)) )
		{
			myLivretACurve = new ARM_LivretACurve(	asOfDate, infCurv, euribCurv, flagArrondi, infresetManager,
													livretAresetManager, euribresetManager, monthForAugustId, monthForFebruaryId);
		}
		else
		{
			result.setMsg ("ARM_ERR: incompatibility of types between livretA fixing and Euribor Fixing");
			return ARM_KO;
		};

		/// assign object
		if( !assignObject( myLivretACurve, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
		
	catch(Exception& x)
	{
		delete myLivretACurve;
		
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


long ARMLOCAL_LIVRETACURVE_WITHRM(double asOfDateDble,
								  long infCurvId,
								  long euribCurvId,
								  long flagArrondi,
								  long infResetManagerId,
								  long euribResetManagerId,
								  long livretAResetManagerId,
								  long monthForAugustId,
								  long monthForFebruaryId,
								  ARM_result& result,
								  long objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_LivretACurve* myLivretACurve = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	try
	{
		char  charDate[11];
		Local_XLDATE2ARMDATE( asOfDateDble, charDate );
		ARM_Date asOfDate( charDate );

		if ( (monthForAugustId != K_MONTH_DEFAULT) &&
			 ( (monthForAugustId < K_MAY) ||
			   (monthForAugustId > K_JUNE)
			 )
		   )
		{
			result.setMsg ("ARM_ERR: Month for August reset must be between May and June or default");
			return ARM_KO;
		}

		if ( (monthForFebruaryId != K_MONTH_DEFAULT) &&
			 ( (monthForFebruaryId < K_NOVEMBER) ||
			   (monthForFebruaryId > K_DECEMBER)
			 )
		   )
		{
			result.setMsg ("ARM_ERR: Month for February reset must be between November and December or default");
			return ARM_KO;
		}

		ARM_InfCurv* infCurv = NULL;
		if( !GetObjectFromId( &infCurv, infCurvId, ARM_INFCURV) )
		{
			result.setMsg ("ARM_ERR: inflation curve is not of a good type");
			return ARM_KO;
		};

		ARM_ZeroCurve* euribCurv=NULL;
		if( !GetObjectFromId( &euribCurv, euribCurvId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: Euribor curve is not of a good type");
			return ARM_KO;
		}

		
		ARM_ResetManager* euribResetManager = NULL;
		if(euribResetManagerId!=-1)
		{
			if( !GetObjectFromIdwNull( &euribResetManager, euribResetManagerId, ARM_RESETMANAGER ) )
			{
				result.setMsg ("ARM_ERR: Euribor ResetManager is not of a good type");
				return ARM_KO;
			}
		}

		ARM_ResetManager* livretAResetManager = NULL;
		if(livretAResetManagerId!=-1)
		{
			if( !GetObjectFromIdwNull( &livretAResetManager, livretAResetManagerId, ARM_RESETMANAGER ) )
			{
				result.setMsg ("ARM_ERR: livret A ResetManager is not of a good type");
				return ARM_KO;
			}
		}
		
		ARM_ResetManager* infResetManager = NULL;
		if(infResetManagerId!=-1)
		{

			if( !GetObjectFromIdwNull( &infResetManager, infResetManagerId, ARM_RESETMANAGER ) )
			{
				result.setMsg ("ARM_ERR: inflation ResetManager is not of a good type");
				return ARM_KO;
			};
		}

		myLivretACurve = new ARM_LivretACurve(	asOfDate, infCurv, euribCurv, flagArrondi, infResetManager,
												livretAResetManager, euribResetManager, monthForAugustId, monthForFebruaryId);

		/// assign object
		if( !assignObject( myLivretACurve, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
		
	catch(Exception& x)
	{
		delete myLivretACurve;
		
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


long ARMLOCAL_LIVRETACURVEGETRATEDATE(
	long livretACurvId,
	double dateDble ,
	ARM_result& result)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	try
	{
		char  charDate[11];
		Local_XLDATE2ARMDATE( dateDble, charDate );
		ARM_Date dateIn( charDate );

		ARM_LivretACurve* livretACurv = NULL;
		if( !GetObjectFromId( &livretACurv, livretACurvId, ARM_LIVRETACURVE) )
		{
			result.setMsg ("ARM_ERR: Livret A curve is not of a good type");
			return ARM_KO;
		};

		int nb = livretACurv->GetRateDateNb(dateIn);

		double tmp = livretACurv->GetDateTerms()->Elt(nb) + livretACurv->GetAsOfDateJul();
		
		double xlDate = dateDble - (dateIn.GetJulian() - tmp);

		result.setDouble(xlDate);

		return ARM_OK;

	}
		
	catch(Exception& x)
	{	
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


/// From year on year vols to periodic annual swap rate vol
long ARMLOCAL_InfSwoVolCurveFromModel_Create(
	const ARM_Date& asOfDate,
	long InfIRModelId,
	const VECTOR<double>& tenors,
	const VECTOR<double>& expiries,
	long computationMethod,
	ARM_result&	result,
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_VolCurve* volCurve		= NULL;
	ARM_InfBSModel* infIRModel	= NULL;
	std::vector<double>& tenorsVec		= NULL;
	std::vector<double>& expiriesVec		= NULL;

	try
	{
		if( !GetObjectFromId( &infIRModel, InfIRModelId, ARM_INFBSMODEL ) )
		{
			result.setMsg ("ARM_ERR: generic Security is not of a good type");
			return ARM_KO;
		}
		tenorsVec	= CreateARMGPVectorFromVECTOR( tenors );
		expiriesVec = CreateARMGPVectorFromVECTOR( expiries );

		volCurve	= ARM_InfVolComputation_Factory::GenerateSwopVolCurve( 
			asOfDate, infIRModel, tenorsVec, 
			expiriesVec, (ARM_InfVolComputation_Factory::ComputationMethodType) computationMethod);

		/// assign object
		if( !assignObject( volCurve, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete volCurve;
		x.DebugPrint();
		ARM_RESULT();
	}
}
//YK

/// From year on year vols to periodic annual swap rate vol
long ARMLOCAL_InfSwoVolCubeFromModel_Create(
	const ARM_Date& asOfDate,
	long InfIRModelId,
	const VECTOR<double>& tenors,
	const VECTOR<double>& expiries,
	const VECTOR<double>& smiledTenors,
	const VECTOR<double>& strikes,
	long computationMethod,
	ARM_result&	result,
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_VolCube* volCube		= NULL;
	ARM_InfBSModel* infIRModel	= NULL;
	std::vector<double>& tenorsVec		= NULL;
	std::vector<double>& expiriesVec		= NULL;
	std::vector<double>& strikesVec		= NULL;
	std::vector<double>& smiledTenorsVec		= NULL;

	try
	{
		if( !GetObjectFromId( &infIRModel, InfIRModelId, ARM_INFBSMODEL ) )
		{
			result.setMsg ("ARM_ERR: generic Security is not of a good type");
			return ARM_KO;
		}
		tenorsVec	= CreateARMGPVectorFromVECTOR( tenors );
		expiriesVec = CreateARMGPVectorFromVECTOR( expiries );
		smiledTenorsVec = CreateARMGPVectorFromVECTOR( smiledTenors );
		strikesVec = CreateARMGPVectorFromVECTOR( strikes );
		
		volCube	= ARM_InfVolCubeComputation_Factory::GenerateSwopVolCube( 
			asOfDate, infIRModel, tenorsVec, 
			expiriesVec,smiledTenorsVec,strikesVec, (ARM_InfVolComputation_Factory::ComputationMethodType) computationMethod);

		/// assign object
		if( !assignObject( volCube, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete volCube;
		x.DebugPrint();
		ARM_RESULT();
	}
}
//YK

/// From zero coupon vols to periodic OAT swap rate vol
long ARMLOCAL_InfOATSwoVolCurveFromModel_Create(
	const	ARM_Date& asOfDate,
	long	InfIRModelId,
	const	VECTOR<double>& tenors,
	const	VECTOR<double>& expiries,
	double	coupon,
	long	choice,
	long	computationMethod,
	ARM_result&	result,
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_VolCurve* volCurve		= NULL;
	ARM_InfBSModel* infIRModel	= NULL;
	std::vector<double>& tenorsVec	= NULL;
	std::vector<double>& expiriesVec	= NULL;

	try
	{
		if( !GetObjectFromId( &infIRModel, InfIRModelId, ARM_INFBSMODEL ) )
		{
			result.setMsg ("ARM_ERR: generic Security is not of a good type");
			return ARM_KO;
		}
		tenorsVec	= CreateARMGPVectorFromVECTOR( tenors );
		expiriesVec = CreateARMGPVectorFromVECTOR( expiries );
		
		coupon /= CC_NS(ARM_Constants,rateBase);

		volCurve	= ARM_InfVolComputation_Factory::GenerateOATSwopVolCurve(
			asOfDate,
			infIRModel,
			tenorsVec,
			expiriesVec,
			(ARM_InfVolComputation_Factory::ComputationMethodType) computationMethod,
			coupon,
			choice);

		/// assign object
		if( !assignObject( volCurve, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete volCurve;
		x.DebugPrint();
		ARM_RESULT();
	}
}




long ARMLOCAL_infMultiBSMod_Create(	
	const VECTOR<long>& infMultiBSModIdVec,
	ARM_result& result,
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_InfMultiBSModel* myModel = NULL;
	//ARM_InfBSModelPtr modelPtr ;
	try
	{
		vector<ARM_InfBSModelPtr> models;

		for( size_t i=0; i<infMultiBSModIdVec.size(); ++i )
		{
			ARM_InfBSModel* model;
			if( !GetObjectFromId( &model, infMultiBSModIdVec[i], ARM_INFBSMODEL ) )
			{
				result.setMsg ("ARM_ERR: inflation bs model is not of a good type");
				return ARM_KO;
			};
			//modelPtr = ARM_InfBSModelPtr( (ARM_InfBSModel* ) model->Clone() );
			ARM_InfBSModelPtr modelPtr( ( ARM_InfBSModel* ) model->Clone() );
			models.push_back( modelPtr);
		}

		myModel = new ARM_InfMultiBSModel(models );

		/// assign object
		if( !assignObject( myModel, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
		
	catch(Exception& x)
	{
		delete myModel;
		x.DebugPrint();
		ARM_RESULT();
	}
}




/*!
 * Function to set a reset Manager on an inflation curve
 */
extern long ARMLOCAL_InfCurv_SetResetManager(
	const long& infCurvId,
	const long& resetManagerId,
	ARM_result&	result, 
	long		objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_InfCurv* infCurv	= NULL;
	ARM_InfCurv* newCurve	= NULL;

	try
	{
		ARM_ResetManager* resetManager = NULL;
		if( !GetObjectFromId( &resetManager, resetManagerId, ARM_RESETMANAGER ) )
		{
			result.setMsg ("ARM_ERR: reset manager is not of a good type");
			return ARM_KO;
		};

		if( !GetObjectFromId( &infCurv, infCurvId, ARM_INFCURV ) )
		{
			result.setMsg ("ARM_ERR: inf curve is not of a good type");
			return ARM_KO;
		};
		
		newCurve=(ARM_InfCurv*) infCurv->Clone();
		newCurve->SetResetManager(resetManager);
	
		/// assign object
		if( !assignObject( newCurve, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		delete infCurv;
		delete newCurve;
		x.DebugPrint();
		ARM_RESULT();
	}
}



/*!
 * Function to set a reset Manager on an inflation curve
 */
extern long ARMLOCAL_InfCurv_SetSeasonalityManager(
	const long& infCurvId,
	const long& seasonalityManagerId,
	ARM_result&	result, 
	long		objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_InfCurv* infCurv	= NULL;
	ARM_InfCurv* newCurve	= NULL;

	try
	{
		ARM_SeasonalityManager* seasonalityManager = NULL;

		if( !GetObjectFromIdwNull( &seasonalityManager, seasonalityManagerId, ARM_SEASONMANAGER ) )
		{
			result.setMsg ("ARM_ERR: season manager is not of a good type");
			return ARM_KO;
		};

		if( !GetObjectFromId( &infCurv, infCurvId, ARM_INFCURV ) )
		{
			result.setMsg ("ARM_ERR: inf curve is not of a good type");
			return ARM_KO;
		};
		
		newCurve=(ARM_InfCurv*) infCurv->Clone();
		newCurve->SetSeasonalityManager(seasonalityManager);
	
		/// assign object
		if( !assignObject( newCurve, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		delete infCurv;
		delete newCurve;
		x.DebugPrint();
		ARM_RESULT();
	}
}


long ARMLOCAL_GetSeasonMgrFromSummit (const CCString& index,
									  const CCString& ccy,
									  const CCString& cvname,
									  double asOf,
									  long modeId,
									  ARM_result& result,
									  long objId)
{
	long seasonMgrId;

	ARM_SeasonalityManager* newSeasonManager = NULL;
	ARM_SeasonalityManager* prevSeasonManager = NULL;

	char sDate[11];
	Local_XLDATE2ARMDATE(asOf,sDate);

	CCString msg (" ");

	try
	{
//		if (GetETKVersion () >= 1)
//		{
			ARM_Date myDate(sDate);

			newSeasonManager = ARMLOCAL_ParseXMLForSeasonMgr(index,ccy,cvname,myDate,modeId);
/*		}
		else
		{
			result.setMsg ("ARM_ERR: This function is not implemented without ETK");
			return ARM_KO;
		}
*/
		if (newSeasonManager == NULL)
		{
			result.setMsg("Object is Null");
			return ARM_KO;
		}
	}

    catch (Exception& x)
	{
		x.DebugPrint();

		if (newSeasonManager)
			delete newSeasonManager;
		newSeasonManager = NULL;

		ARM_RESULT();
	}

	if(objId == -1)
	{
		CREATE_GLOBAL_OBJECT();

		seasonMgrId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSeasonManager);

		if (seasonMgrId == RET_KO)
		{
			if (newSeasonManager)
				delete newSeasonManager;
			newSeasonManager = NULL;

			result.setMsg ("ARM_ERR: Pb with inserting object");				
			return ARM_KO;
		}

		result.setLong(seasonMgrId);

		return ARM_OK;
	}
	else
	{
		prevSeasonManager = (ARM_SeasonalityManager *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevSeasonManager, ARM_SEASONMANAGER) == 1)
		{
			if (prevSeasonManager)
			{
				delete prevSeasonManager;
				prevSeasonManager = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSeasonManager, objId);

			return ARM_OK;
		}
		else
		{
			result.setMsg ("ARM_ERR: previous object is not of a good type");
			return ARM_KO;
		}
	}
}


long ARMLOCAL_GetResetMgrFromSummit (double C_asof,
									 const CCString& C_index,
									 const CCString& C_source,
									 const CCString& C_ccy,
									 long isInflatIndexId,
									 const CCString& C_term,
									 ARM_result& result,
									 long objId)
{
	
	long resetMgrId;

	ARM_ResetManager* newResetManager = NULL;
	ARM_ResetManager* prevResetManager = NULL;

	char myAsOfDate[11];

	CCString msg (" ");

	try
	{
		Local_XLDATE2ARMDATE(C_asof,myAsOfDate);
		
//		if (GetETKVersion () >= 1)
//		{
			newResetManager = ARMLOCAL_ParseXMLForResetMgr((ARM_Date)myAsOfDate,
														   C_index,
														   C_source,
														   C_ccy,
														   isInflatIndexId,
														   C_term);
/*		}
		else
		{
			result.setMsg ("ARM_ERR: This function is not implemented without ETK");
			return ARM_KO;
		}
*/
		if (newResetManager == NULL)
		{
			result.setMsg("Object is Null");
			
			return ARM_KO;
		}
	}

    catch (Exception& x)
	{
		x.DebugPrint();

		if (newResetManager)
			delete newResetManager;
		newResetManager = NULL;
		
		ARM_RESULT();
	}

	if(objId == -1)
	{
		CREATE_GLOBAL_OBJECT();

		resetMgrId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newResetManager);

		if (resetMgrId == RET_KO)
		{
			if (newResetManager)
				delete newResetManager;
			newResetManager = NULL;

			result.setMsg ("ARM_ERR: Pb with inserting object");				
			
			return ARM_KO;
		}

		result.setLong(resetMgrId);

		
		return ARM_OK;
	}
	else
	{
		prevResetManager = (ARM_ResetManager *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevResetManager, ARM_RESETMANAGER) == 1)
		{
			if (prevResetManager)
			{
				delete prevResetManager;
				prevResetManager = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newResetManager, objId);

			
			return ARM_OK;
		}
		else
		{
			result.setMsg ("ARM_ERR: previous object is not of a good type");
			
			return ARM_KO;
		}
	}
	
}


extern long ARMLOCAL_GP_INFCAPFLOOR (long swapId,
									 long CF,
									 double strike,
									 long strikeId,
									 ARM_result&	result, 
									 long        objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_InfCapFloorRielYield* GP_INFCapFloorObj = NULL;
	
	ARM_ReferenceValue* strikeProfile=NULL;
	ARM_Swap* swapObj = NULL;
    bool isStrike=false;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Caption Calculator" );

		/// Convert curve Id to object if possible
		//coupon
		if(swapId != ARM_NULL_OBJECT)
        {
		    swapObj = dynamic_cast<ARM_Swap*>(LOCAL_PERSISTENT_OBJECTS->GetObject(swapId));

		    if (!swapObj)
		    {
			    result.setMsg ("ARM_ERR: swapObj is not of a good type");
			    return ARM_KO;
		    }
        }

        if(strikeId != ARM_NULL_OBJECT)
        {
		    strikeProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(strikeId));

		    if (!strikeProfile)
		    {
			    result.setMsg ("ARM_ERR: strike is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			strikeProfile = new ARM_ReferenceValue(strike);
            isStrike=true;
        }


	

	    /// Create the Caption calculator
		GP_INFCapFloorObj = new ARM_InfCapFloorRielYield(swapObj,
														CF,
														strikeProfile);

	    /// Free memory
        if(isStrike)
            delete strikeProfile;
		strikeProfile = NULL;
     
       
		/// assign object
		if( !assignObject( GP_INFCapFloorObj, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		/// Free memory
        if(isStrike)
            delete strikeProfile;
        
		delete GP_INFCapFloorObj;
		x.DebugPrint();
		ARM_RESULT();
	}
}



long ARMLOCAL_GP_INFCALLSPREAD (	long			swapId,
									long			CF,
									double			strike,
									long			strikeId,
									ARM_result&		result, 
									long			objId){
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	/// declaration
	ARM_InfCallSpreadYield*		GP_INFDigitalObj	=	NULL;
	ARM_ReferenceValue*			strikeProfile		=	NULL;
	ARM_Swap*					swapObj				=	NULL;
    bool						isStrike			=	false;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Caption Calculator" );

		/// Convert curve Id to object if possible
		//coupon
		if(swapId != ARM_NULL_OBJECT) {
		    
			swapObj = dynamic_cast<ARM_Swap*>(LOCAL_PERSISTENT_OBJECTS->GetObject(swapId));

		    if (!swapObj)  {
			    result.setMsg ("ARM_ERR: swapObj is not of a good type");
			    return ARM_KO;
		    }
        }

        if(strikeId != ARM_NULL_OBJECT)	{

		    strikeProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(strikeId));

		    if (!strikeProfile)   {
			    result.setMsg ("ARM_ERR: strike is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			strikeProfile	= new ARM_ReferenceValue(strike);
            isStrike		= true;
        }

		/// Create the Caption calculator
		GP_INFDigitalObj = new ARM_InfCallSpreadYield( swapObj, CF,	strikeProfile );

	    /// Free memory
        if(isStrike)	delete strikeProfile;
		strikeProfile = NULL;
     
       
		/// assign object
		if( !assignObject( GP_INFDigitalObj, result, objId ) )	return ARM_KO; 
		else	return ARM_OK; 
	}
	
	catch(Exception& x)
	{
		/// Free memory
        if(isStrike)	delete strikeProfile;
        
		delete GP_INFDigitalObj;
		x.DebugPrint();
		ARM_RESULT();
	}
}




long ARMLOCAL_GetReset (long resetMgrId,
						double date,
						ARM_result& result)
{
	double dResult;
	ARM_ResetManager* resetMgr = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char myAsOfDate[11];

	CCString msg (" ");

	try
	{
		Local_XLDATE2ARMDATE(date,myAsOfDate);

		resetMgr = (ARM_ResetManager *) LOCAL_PERSISTENT_OBJECTS->GetObject(resetMgrId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(resetMgr, ARM_RESETMANAGER) == 0)
		{
			result.setMsg ("ARM_ERR: resetMgr is not of a good type");
			return ARM_KO;
		}

		dResult = resetMgr->GetReset(((ARM_Date)myAsOfDate).GetJulian());
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


extern long ARMLOCAL_GP_INFDIGITAL (long payLegId,
									long digitLegId,
									long payOffType,
									double barrier,
									long barrierId,
									long CFId,
									long RecOrPayId,
									ARM_result&	result, 
									long        objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_InfHybridDigital* GP_INFDigitalObj = NULL;
	
	ARM_ReferenceValue* barrierProfile=NULL;
	ARM_SwapLeg* PayLegObj = NULL;
	ARM_SwapLeg* digitLegObj = NULL;
    bool isBarrier=false;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Caption Calculator" );

		/// Convert curve Id to object if possible
		//coupon
		if(payLegId != ARM_NULL_OBJECT)
        {
		    PayLegObj = dynamic_cast<ARM_SwapLeg*>(LOCAL_PERSISTENT_OBJECTS->GetObject(payLegId));

		    if (!PayLegObj)
		    {
			    result.setMsg ("ARM_ERR: PayLegObj is not of a good type");
			    return ARM_KO;
		    }
        }

		if(digitLegId != ARM_NULL_OBJECT)
        {
		    digitLegObj = dynamic_cast<ARM_SwapLeg*>(LOCAL_PERSISTENT_OBJECTS->GetObject(digitLegId));

		    if (!digitLegObj)
		    {
			    result.setMsg ("ARM_ERR: digitLegObj is not of a good type");
			    return ARM_KO;
		    }
        }


        if(barrierId != ARM_NULL_OBJECT)
        {
		    barrierProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(barrierId));

		    if (!barrierProfile)
		    {
			    result.setMsg ("ARM_ERR: barrier is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			barrierProfile = new ARM_ReferenceValue(barrier);
            isBarrier=true;
        }


	

	    /// Create the Caption calculator ARM_InfHybridDigital* GP_INFDigitalObj 
		GP_INFDigitalObj = new ARM_InfHybridDigital(PayLegObj,
													digitLegObj,
													payOffType,													
													barrierProfile,
													CFId,
													RecOrPayId);

	    /// Free memory
        if(isBarrier)
            delete barrierProfile;
		barrierProfile = NULL;
     
       
		/// assign object
		if( !assignObject( GP_INFDigitalObj, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		/// Free memory
        if(isBarrier)
            delete barrierProfile;
        
		delete GP_INFDigitalObj;
		x.DebugPrint();
		ARM_RESULT();
	}
}


extern long ARMLOCAL_HybridInfIrMkt_Create(		const ARM_Date&			asOf, 
												const vector<string>&	keys,
												const vector<long>&		mods,
												ARM_result&				result,	
												long					objId)	{

	if( !GlobalPersistanceOk( result ) )	return ARM_KO;


	CCString msg ("");	

	ARM_MarketData_ManagerRep* infMkt = new ARM_MarketData_ManagerRep(asOf);

	try	{
		if (keys.size() != mods.size() ){
			ARM_THROW( ERR_INVALID_ARGUMENT, " : keys and Mkt Datas should contain the same number of elements " );
		}

		map<string, ARM_Object*> objMap;
		for (int i = 0; i < keys.size(); ++i)
			infMkt->RegisterData(keys[i],LOCAL_PERSISTENT_OBJECTS->GetObject(mods[i]));


		if( !assignObject( infMkt, result, objId ) ){ return ARM_KO; }
		else{	return ARM_OK; }
	}
	
	catch(Exception& x){
		delete infMkt;
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_HybridInfIrPayOff_Create(	const long   & cstCpnCoef,
												const long   & cstOptCoef,
												
												const string & mainCpnName,
												const long	 & mainCpnCoef,
												const string & mainOptName, 
												const long	 & mainOptCoef, 
												
												const string & subCpnName,
												const long	 & subCpnCoef,
												const string & subOptName, 
												const long   & subOptCoef,
												
												const string & supCpnName,
												const long	 & supCpnCoef, 
												const string & supOptName, 
												const long   & supOptCoef,
												ARM_result	 & result,	
												long		   objId)	{

	if( !GlobalPersistanceOk( result ) )	return ARM_KO;


	CCString msg ("");	

	ARM_InfHybridPayOff* hybridInfIrPayOff=NULL;

	try	{
		ARM::ARM_MAP_Curve cpnCurve;
		ARM::ARM_MAP_Curve optCurve;


//==> MainCpn
		string		MainCpnName	= mainCpnName;
		stringToUpper(MainCpnName);	
		ARM_Curve*	MainCpnCoef	= NULL;
		MainCpnCoef = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(mainCpnCoef) );
		
		if ( MainCpnName =="NO" && MainCpnCoef){
			result.setMsg ("ARM_ERR: curve without index name");
			return ARM_KO;
		}

		else if ( MainCpnName !="NO" && !MainCpnCoef){
			result.setMsg ("ARM_ERR: curve without values");
			return ARM_KO;
		}

		else if (MainCpnName !="NO" && MainCpnCoef){
			pair<string,ARM_GP_CurvePtr> p( MainCpnName,ARM_GP_CurvePtr(CreateClone(MainCpnCoef) ) );
			cpnCurve.insert(p);
		}

//==> SubCpn
		string		SubCpnName	= subCpnName;
		stringToUpper(SubCpnName);	
		ARM_Curve*	SubCpnCoef	= NULL;
		SubCpnCoef = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(subCpnCoef) );

		if ( SubCpnName =="NO" && SubCpnCoef){
			result.setMsg ("ARM_ERR: curve without index name");
			return ARM_KO;
		}

		else if ( SubCpnName !="NO" && !SubCpnCoef){
			result.setMsg ("ARM_ERR: curve without values");
			return ARM_KO;
		}

		else if (SubCpnName !="NO" && SubCpnCoef){
			pair<string,ARM_GP_CurvePtr> p( SubCpnName,ARM_GP_CurvePtr(CreateClone(SubCpnCoef) ) );
			cpnCurve.insert(p);
		}


//==> SupCpn
		string		SupCpnName	= supCpnName;
		stringToUpper(SupCpnName);	
		ARM_Curve*	SupCpnCoef	= NULL;
		SupCpnCoef = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(supCpnCoef) );

		if ( SupCpnName =="NO" && SupCpnCoef){
			result.setMsg ("ARM_ERR: curve without index name");
			return ARM_KO;
		}

		else if ( SupCpnName !="NO" && !SupCpnCoef){
			result.setMsg ("ARM_ERR: curve without values");
			return ARM_KO;
		}

		else if (SupCpnName !="NO" && SupCpnCoef){
			pair<string,ARM_GP_CurvePtr> p( SupCpnName,ARM_GP_CurvePtr(CreateClone(SupCpnCoef) ) );
			cpnCurve.insert(p);
		}


//==> CstCpn
		string		CstCpnName	= "NO";	
		ARM_Curve*	CstCpnCoef	= NULL;
		CstCpnCoef = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(cstCpnCoef));
		if (!CstCpnCoef)	{
			result.setMsg ("ARM_ERR: nominal should be a curve");
			return ARM_KO;
		}
		else{
			pair<string,ARM_GP_CurvePtr> p( CstCpnName,ARM_GP_CurvePtr(CreateClone(CstCpnCoef) ) );
			cpnCurve.insert(p);
		}

//==> MainOpt
		string		MainOptName	= mainOptName;
		stringToUpper(MainOptName);	
		ARM_Curve*	MainOptCoef	= NULL;
		MainOptCoef = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(mainOptCoef) );
		
		if ( MainOptName =="NO" && MainOptCoef){
			result.setMsg ("ARM_ERR: curve without index name");
			return ARM_KO;
		}

		else if ( MainOptName !="NO" && !MainOptCoef){
			result.setMsg ("ARM_ERR: curve without values");
			return ARM_KO;
		}

		else if ( MainOptName !="NO" && MainOptCoef){
			pair<string,ARM_GP_CurvePtr> p( MainOptName,ARM_GP_CurvePtr(CreateClone(MainOptCoef) ) );
			optCurve.insert(p);
		}

//==> SubOpt
		string		SubOptName	= subOptName;
		stringToUpper(SubOptName);	
		ARM_Curve*	SubOptCoef	= NULL;
		SubOptCoef = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(subOptCoef) );

		if ( SubOptName =="NO" && SubOptCoef){
			result.setMsg ("ARM_ERR: curve without index name");
			return ARM_KO;
		}

		else if ( SubOptName !="NO" && !SubOptCoef){
			result.setMsg ("ARM_ERR: curve without values");
			return ARM_KO;
		}

		else if ( SubOptName !="NO" && SubOptCoef){
			pair<string,ARM_GP_CurvePtr> p( SubOptName,ARM_GP_CurvePtr(CreateClone(SubOptCoef)) );
			optCurve.insert(p);
		}


//==> SupOpt
		string		SupOptName	= supOptName;
		stringToUpper(SupOptName);	
		ARM_Curve*	SupOptCoef	= NULL;
		SupOptCoef = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(supOptCoef) );

		if ( SupOptName =="NO" && SupOptCoef){
			result.setMsg ("ARM_ERR: curve without index name");
			return ARM_KO;
		}

		else if ( SupOptName !="NO" && !SupOptCoef){
			result.setMsg ("ARM_ERR: curve without values");
			return ARM_KO;
		}

		else if ( SupOptName !="NO" && SupOptCoef){
			pair<string,ARM_GP_CurvePtr> p( SupOptName,ARM_GP_CurvePtr(CreateClone(SupOptCoef) ) );
			optCurve.insert(p);
		}

//==> CstOpt
		string		CstOptName	= "NO";	
		ARM_Curve*	CstOptCoef	= NULL;
		CstOptCoef = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(cstOptCoef));
		if (!CstOptCoef)	{
			result.setMsg ("ARM_ERR: nominal should be a curve");
			return ARM_KO;
		}
		else{
			pair<string,ARM_GP_CurvePtr> p( CstOptName,ARM_GP_CurvePtr(CreateClone(CstOptCoef) ) );
			optCurve.insert(p);
		}


		hybridInfIrPayOff = new ARM_InfHybridPayOff( cpnCurve,  optCurve);
		if ( !assignObject( hybridInfIrPayOff, result, objId ) )		return ARM_KO; 
		else	return ARM_OK; 

	}
	
	catch(Exception& x){
		delete hybridInfIrPayOff;
		x.DebugPrint();
		ARM_RESULT();
	}
}


extern long ARMLOCAL_HybridInfIr_Load(		const long &ins, 
											const long &mkt,
											const long &mod,
											const long &pay,
											ARM_result &result,	
											long	objId)	{

	if( !GlobalPersistanceOk( result ) )	return ARM_KO;

	CCString msg ("");	
	ARM_HybridInfIrLeg* infLeg = NULL;
	ARM_InfPricer* infPricer = NULL;

	try	{

		ARM_Object* tmpIns = LOCAL_PERSISTENT_OBJECTS->GetObject(ins);
		if ( dynamic_cast<ARM_HybridInfIrLeg*> ( tmpIns ) )
			infLeg 	= CreateClone(dynamic_cast<ARM_HybridInfIrLeg*> ( tmpIns ) );
		else{
			ARM_THROW( ERR_INVALID_ARGUMENT, " : the instrument should be a Hybrid Inf Ir Leg" ); }

		ARM_Object* tmpMkt = LOCAL_PERSISTENT_OBJECTS->GetObject(mkt);
		if ( !tmpIns && tmpMkt ){
			ARM_THROW( ERR_INVALID_ARGUMENT, " : The instrument should be built before loading mkt" ); }
		if ( dynamic_cast<ARM_MarketData_ManagerRep*> ( tmpMkt ) )
			infLeg->Init( dynamic_cast<ARM_MarketData_ManagerRep*> ( tmpMkt ) );
		else{
			ARM_THROW( ERR_INVALID_ARGUMENT, " : the market data is not a ARM_MarketData_ManagerRep" ); }
		
		ARM_Object* tmpPay = LOCAL_PERSISTENT_OBJECTS->GetObject(pay);

		ARM_Object* tmpMod = LOCAL_PERSISTENT_OBJECTS->GetObject(mod);

		if ( !tmpIns && tmpPay ){
			ARM_THROW( ERR_INVALID_ARGUMENT, " : The instrument should be built before loading pay" ); }
		if ( dynamic_cast<ARM_InfHybridPayOff*> ( tmpPay ) ){
			infPricer = new ARM::ARM_HybridInfLegPricer(  infLeg, 
															dynamic_cast<ARM_InfHybridPayOff*> ( tmpPay ), 
															dynamic_cast<ARM_InfModel*>  ( tmpMod ) ); 
		}

		if ( infPricer ){
			if( !assignObject( infPricer, result, objId ) ){ return ARM_KO; }
		}
		else
			if( !assignObject( infLeg, result, objId ) ){ return ARM_KO; }
		else{	return ARM_OK; }
	}
	
	catch(Exception& x){
		delete infLeg;
		x.DebugPrint();
		ARM_RESULT();
	}
	return ARM_OK;
}

extern long ARMLOCAL_Inf_GetPrice(
	const long&				pricerId,
	const string&			key,
	ARM_GramFctorArg&		argResult,
	ARM_result&				result )
{
	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	
	try	{

		ARM_Object* obj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId );
		ARM_InfPricer* pricer = dynamic_cast<ARM_InfPricer*>(obj);
		if( !pricer ){
			result.setMsg ("ARM_ERR: the inflation pricer is empty");
			return ARM_KO;
		}

		pricer->Compute();
		
		argResult = pricer->GetFunctor().GetData(key);
	}
	catch(Exception& x){
		x.DebugPrint();
		ARM_RESULT();
	}

	catch (...){
		result.setMsg("ARM_ERR: unrecognized failure in price computation");
        return(ARM_KO);
	}
	return ARM_OK;
}


extern long ARMLOCAL_Inf_GetSchedule(			const long				& legId,
												const string			& key,
												ARM_GramFctorArg		& argResult,
												ARM_result				& result){
	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	
	try	{

		ARM_Object* obj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(legId );
		string tmp = key;
		stringToUpper(tmp);

		ARM_HybridInfIrLeg* leg = dynamic_cast<ARM_HybridInfIrLeg*>(obj);
		if( !leg ){
			result.setMsg ("ARM_ERR: the leg is empty");
			return ARM_KO;
		}

		argResult = leg->GetFunctor().GetData(key);
	}
	catch(Exception& x){
		x.DebugPrint();
		ARM_RESULT();
	}

	catch (...){
		result.setMsg("ARM_ERR: unrecognized failure in price computation");
        return(ARM_KO);
	}
	return ARM_OK;
}


long ARM_HybridInfIrModel_CreateFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	string tmp;

	ARM_GenericParams*	genericParams	= GetGenericParams();
	ARM_InfModel*		model			= NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Hybrid Inf/Ir Model Create" );


		string	modelName		= genericParams->GetParamValue("ModelName").GetString();
		stringToUpper(modelName);

		double	discretisation	=	genericParams->GetParamValue("Discretisation").GetDouble();
		double	domain			=	genericParams->GetParamValue("Domain").GetDouble();
		double	epsilon			=	genericParams->GetParamValue("Epsilon").GetDouble();	
		double	center			=	genericParams->GetParamValue("Center").GetDouble();

		if ( modelName=="BILOG" )
			model = new ARM_InfBiLog( modelName, discretisation, domain, epsilon, center);

		if ( modelName=="NUMBILOG" )
			model = new ARM_InfNumBiLog( modelName, discretisation, domain, epsilon, center);

		if ( modelName=="HK" )
			model = new ARM_InfHK( modelName, discretisation, domain, epsilon, center);


		if ( !assignObject( model, result, objId ) )
			return ARM_KO; 
		else
			return ARM_OK; 
	}

	catch(Exception& x){
		delete	model;
		x.DebugPrint();
		ARM_RESULT();
	}
}


extern long ARMLOCAL_InfSpreadCap_Create(		const string & mainIndex,
												const string & mainType,
												const long	 & mainLeverage,
												const string & subIndex,
												const string & subType,
												const long	 & subLeverage,
												const long	 & strike, 
												const long	 & notional, 
												ARM_result	 & result,	
												long		   objId)	{

	if( !GlobalPersistanceOk( result ) )	return ARM_KO;


	CCString msg ("");	

	ARM_InfHybridPayOff* hybridInfIrPayOff=NULL;

	try	{

		ARM::ARM_MAP_Curve optCurve;
		pair<string,ARM_GP_CurvePtr> p;

		ARM_Curve*	Notional		= CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(notional) ));
		ARM_Curve*	MainLeverage	= CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(mainLeverage) ));
		ARM_Curve*	SubLeverage		= CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(subLeverage) ));
		ARM_Curve*	Strike			= CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(strike) ));

		if ( !Notional ){
			result.setMsg ("ARM_ERR: Notional is not build");
			return ARM_KO;
		}

		*Strike			*= *Notional;

		*MainLeverage	*= *Notional;
		if ( mainType == "INF" ){
			*MainLeverage	*= 100.0;
			*Strike			+= *MainLeverage;
		}

		*SubLeverage	*= *Notional;
		if ( subType == "INF" ){
			*SubLeverage	*= 100.0;
			*Strike			+= *SubLeverage;
		}

		*Strike			*= -1.0;

//==> MainIndex
		string MainIndex = mainIndex;
		stringToUpper(MainIndex);	

		if ( !MainLeverage ){
			result.setMsg ("ARM_ERR: MainLeverage is not build");
			return ARM_KO;
		}
		p = pair<string,ARM_GP_CurvePtr>( MainIndex, ARM_GP_CurvePtr(CreateClone(MainLeverage) ) );
		optCurve.insert(p);

//==> SubIndex
		string SubIndex = subIndex;
		stringToUpper(SubIndex);	

		if ( !SubLeverage ){
			result.setMsg ("ARM_ERR: SubLeverage is not build");
			return ARM_KO;
		}
		p = pair<string,ARM_GP_CurvePtr>( SubIndex,ARM_GP_CurvePtr(CreateClone(SubLeverage) ) );
		optCurve.insert(p);


//==> Strike
		if (!Strike)	{
			result.setMsg ("ARM_ERR: Strike is not build");
			return ARM_KO;
		}
		p = pair<string,ARM_GP_CurvePtr>( "NO",ARM_GP_CurvePtr(CreateClone(Strike) ) );
		optCurve.insert(p);


		hybridInfIrPayOff = new ARM_InfHybridCap( optCurve );
		if ( !assignObject( hybridInfIrPayOff, result, objId ) )		return ARM_KO; 
		else	return ARM_OK; 

	}
	
	catch(Exception& x){
		delete hybridInfIrPayOff;
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_InfSpreadDigital_Create(	const string & mainIndex,
												const string & mainType,
												const long	 & mainLeverage,
												const string & subIndex,
												const string & subType,
												const long	 & subLeverage,
												const long	 & strike, 
												const long	 & notional, 
												ARM_result	 & result,	
												long		   objId)	{

	if( !GlobalPersistanceOk( result ) )	return ARM_KO;


	CCString msg ("");	

	ARM_InfHybridPayOff* hybridInfIrPayOff=NULL;

	try	{
		ARM::ARM_MAP_Curve optCurve;
		pair<string,ARM_GP_CurvePtr> p;

		ARM_Curve*	Notional		= CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(notional) ));
		ARM_Curve*	MainLeverage	= CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(mainLeverage) ));
		ARM_Curve*	SubLeverage		= CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(subLeverage) ));
		ARM_Curve*	Strike			= CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(strike) ));

		if ( !Notional ){
			result.setMsg ("ARM_ERR: Notional is not build");
			return ARM_KO;
		}

		*Strike			*= *Notional;

		*MainLeverage	*= *Notional;
		if ( mainType == "INF" ){
			*MainLeverage	*= 100.0;
			*Strike			+= *MainLeverage;
		}

		*SubLeverage	*= *Notional;
		if ( subType == "INF" ){
			*SubLeverage	*= 100.0;
			*Strike			+= *SubLeverage;
		}

		*Strike			*= -1.0;

//==> MainIndex
		string MainIndex = mainIndex;
		stringToUpper(MainIndex);	

		if ( !MainLeverage ){
			result.setMsg ("ARM_ERR: MainLeverage is not build");
			return ARM_KO;
		}
		p = pair<string,ARM_GP_CurvePtr>( MainIndex, ARM_GP_CurvePtr(CreateClone(MainLeverage) ) );
		optCurve.insert(p);

//==> SubIndex
		string SubIndex = subIndex;
		stringToUpper(SubIndex);	

		if ( !SubLeverage ){
			result.setMsg ("ARM_ERR: SubLeverage is not build");
			return ARM_KO;
		}
		p = pair<string,ARM_GP_CurvePtr>( SubIndex,ARM_GP_CurvePtr(CreateClone(SubLeverage) ) );
		optCurve.insert(p);


//==> Strike
		if (!Strike)	{
			result.setMsg ("ARM_ERR: Strike is not build");
			return ARM_KO;
		}
		p = pair<string,ARM_GP_CurvePtr>( "NO",ARM_GP_CurvePtr(CreateClone(Strike) ) );
		optCurve.insert(p);


		hybridInfIrPayOff = new ARM_InfHybridDigit( optCurve, ARM_GP_CurvePtr(Notional) );

		if ( !assignObject( hybridInfIrPayOff, result, objId ) )		return ARM_KO; 
		else	return ARM_OK; 

	}
	
	catch(Exception& x){
		delete hybridInfIrPayOff;
		x.DebugPrint();
		ARM_RESULT();
	}
}


extern long ARMLOCAL_InfDoubleDigital_Create(	const string & mainIndex,
												const string & mainType,
												const long	 & mainLeverage,
												const long	 & mainStrike,
												const string & subIndex,
												const string & subType,
												const long	 & subLeverage,
												const long	 & substrike, 
												const long	 & notional, 
												const long	 & spread,
												ARM_result	 & result,	
												long		   objId)	{

	if( !GlobalPersistanceOk( result ) )	return ARM_KO;


	CCString msg ("");	

	ARM_InfDoubleDigit* doubleDigit=NULL;

	try	{
		ARM::ARM_MAP_Curve cpnCurve;
		ARM::ARM_MAP_Curve m_optCurve;
		ARM::ARM_MAP_Curve s_optCurve;

		pair<string,ARM_GP_CurvePtr> p;

		ARM_Curve*	Notional		=	CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(notional)		) );
		ARM_Curve*	MainLeverage	=	CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(mainLeverage)	) );
		ARM_Curve*	SubLeverage		=	CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(subLeverage)		) );
		ARM_Curve*	MainStrike		=	CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(mainStrike)		) );
		ARM_Curve*	SubStrike		=	CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(substrike)		) );
		ARM_Curve*	Spread			=	CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(spread)			) );

		ARM_Curve*	OptMainLeverage	=	CreateClone(MainStrike);
		*OptMainLeverage	*= 0.0;
		*OptMainLeverage	+= 1.0;
		ARM_Curve*	OptSubLeverage	=	CreateClone(SubStrike);
		*OptSubLeverage		*= 0.0;
		*OptSubLeverage		+= 1.0;

		if ( !Notional ){
			result.setMsg ("ARM_ERR: Notional is not build");
			return ARM_KO;
		}

		*Spread	*= *Notional;

		if ( mainType == "INF" ){
			*OptMainLeverage	*= 100.0;
			*MainStrike			+= *OptMainLeverage;
			*MainLeverage		*= 100.0;
		}
		*MainLeverage	*= *Notional;

		*SubLeverage	*= *Notional;
		if ( subType == "INF" ){
			*OptSubLeverage	*= 100.0;
			*SubStrike		+= *OptSubLeverage;
			*SubLeverage	*= 100.0;
			*Spread			-= *SubLeverage;
		}

		*MainStrike	*= -1.0;
		*SubStrike	*= -1.0;

//==> MainIndex
		string MainIndex = mainIndex;
		stringToUpper(MainIndex);	

		p = pair<string,ARM_GP_CurvePtr>( MainIndex, ARM_GP_CurvePtr(CreateClone(MainLeverage) ) );
		cpnCurve.insert(p);
		p = pair<string,ARM_GP_CurvePtr>( MainIndex, ARM_GP_CurvePtr(CreateClone(OptMainLeverage) ) );
		m_optCurve.insert(p);
		p = pair<string,ARM_GP_CurvePtr>( "NO", ARM_GP_CurvePtr(CreateClone(MainStrike) ) );
		m_optCurve.insert(p);

//==> SubIndex
		string SubIndex = subIndex;
		stringToUpper(SubIndex);	

		p = pair<string,ARM_GP_CurvePtr>( SubIndex, ARM_GP_CurvePtr(CreateClone(SubLeverage) ) );
		cpnCurve.insert(p);
		p = pair<string,ARM_GP_CurvePtr>( SubIndex, ARM_GP_CurvePtr(CreateClone(OptSubLeverage) ) );
		s_optCurve.insert(p);
		p = pair<string,ARM_GP_CurvePtr>( "NO", ARM_GP_CurvePtr(CreateClone(SubStrike) ) );
		s_optCurve.insert(p);

		p = pair<string,ARM_GP_CurvePtr>( "NO", ARM_GP_CurvePtr(CreateClone(Spread) ) );
		cpnCurve.insert(p);


		doubleDigit = new ARM_InfDoubleDigit( cpnCurve, m_optCurve, s_optCurve );

		if ( !assignObject( doubleDigit, result, objId ) )		return ARM_KO; 
		else	return ARM_OK; 

	}
	
	catch(Exception& x){
		delete doubleDigit;
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_InfCorridor_Create(		const string & mainIndex,
												const string & mainType,
												const long	 & mainLeverage,
												const string & subIndex,
												const string & subType,
												const long	 & subLeverage,
												const long	 & strikeInf, 
												const long	 & strikeSup,
												const long	 & notional, 
												ARM_result	 & result,	
												long		   objId)	{

	if( !GlobalPersistanceOk( result ) )	return ARM_KO;


	CCString msg ("");	

	ARM_InfHybridPayOff* hybridInfIrPayOff=NULL;

	try	{

		ARM::ARM_MAP_Curve cpnCurve;			// borne sup du corridor
		ARM::ARM_MAP_Curve optCurve;			// borne inf du corridor

		pair<string,ARM_GP_CurvePtr> p;

		ARM_Curve*	Notional		= CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(notional) ));
		ARM_Curve*	MainLeverage	= CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(mainLeverage) ));
		ARM_Curve*	SubLeverage		= CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(subLeverage) ));
		ARM_Curve*	StrikeInf		= CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(strikeInf) ));
		ARM_Curve*	StrikeSup		= CreateClone(dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(strikeSup) ));

		if ( !Notional ){
			result.setMsg ("ARM_ERR: Notional is not build");
			return ARM_KO;
		}

		*StrikeInf		*= *Notional;

		*MainLeverage	*= *Notional;
		if ( mainType == "INF" ){
			*MainLeverage	*= 100.0;
			*StrikeInf		+= *MainLeverage;
			*StrikeSup		+= *MainLeverage;
		}

		*SubLeverage	*= *Notional;
		if ( subType == "INF" ){
			*SubLeverage	*= 100.0;
			*StrikeInf		+= *SubLeverage;
			*StrikeSup		+= *SubLeverage;
		}

		*StrikeInf			*= -1.0;
		*StrikeSup			*= -1.0;

//==> MainIndex
		string MainIndex = mainIndex;
		stringToUpper(MainIndex);	

		if ( !MainLeverage ){
			result.setMsg ("ARM_ERR: MainLeverage is not build");
			return ARM_KO;
		}
		p = pair<string,ARM_GP_CurvePtr>( MainIndex, ARM_GP_CurvePtr(CreateClone(MainLeverage) ) );
		optCurve.insert(p);
		cpnCurve.insert(p);

//==> SubIndex
		string SubIndex = subIndex;
		stringToUpper(SubIndex);	

		if ( !SubLeverage ){
			result.setMsg ("ARM_ERR: SubLeverage is not build");
			return ARM_KO;
		}
		p = pair<string,ARM_GP_CurvePtr>( SubIndex,ARM_GP_CurvePtr(CreateClone(SubLeverage) ) );
		optCurve.insert(p);
		cpnCurve.insert(p);

//==> StrikeInf 
		if (!StrikeInf)	{
			result.setMsg ("ARM_ERR: StrikeInf is not build");
			return ARM_KO;
		}
		p = pair<string,ARM_GP_CurvePtr>( "NO",ARM_GP_CurvePtr(CreateClone(StrikeInf) ) );
		cpnCurve.insert(p);

//==> StrikeSup
		if (!StrikeSup)	{
			result.setMsg ("ARM_ERR: StrikeSup is not build");
			return ARM_KO;
		}
		p = pair<string,ARM_GP_CurvePtr>( "NO",ARM_GP_CurvePtr(CreateClone(StrikeSup) ) );
		optCurve.insert(p);

		hybridInfIrPayOff = new ARM_InfCorridor( cpnCurve, optCurve );
		if ( !assignObject( hybridInfIrPayOff, result, objId ) )		return ARM_KO; 
		else	return ARM_OK; 

	}
	
	catch(Exception& x){
		delete hybridInfIrPayOff;
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Inf_GetAdjCorrel(
	const long		& infBsSmiledId,
	ARM_result		& result,	
	long		     objId)
{
	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	
	try	{

		ARM_Object* obj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(infBsSmiledId );
		ARM_InfBSSmiledModel* myModel = dynamic_cast<ARM_InfBSSmiledModel*>(obj);

		if( !myModel ){
			result.setMsg ("ARM_ERR: the inflation pricer is empty");
			return ARM_KO;
		}

		ARM_VolCurve* correl =	myModel->GetAdjConvCorrel();
		
		if ( !assignObject( correl, result, objId ) )		return ARM_KO; 
		else	return ARM_OK; 
	}
	catch(Exception& x){
		x.DebugPrint();
		ARM_RESULT();
	}

	catch (...){
		result.setMsg("ARM_ERR: unrecognized failure in price computation");
        return(ARM_KO);
	}
	return ARM_OK;
}

extern long ARMLOCAL_InfEqHwSV_Laplace(
	const long&				modelId,
	const double&			evalTime,
	const double&			startTime,
	const double&			endTime,
	const double&			xt,
	const double&			vt,
	const double&			k_real,
	const double&			k_imag,
	const bool&				isReal,
	ARM_result&				result ){

	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	
	try	{

		ARM_Object* obj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelId );
		ARM_EQHWSV* mod = dynamic_cast<ARM_EQHWSV*>(obj);
		if( !mod ){
			result.setMsg ("ARM_ERR: model is not of a good type");
			return ARM_KO;
		}

		ARM_PricingStatesPtr states			=	ARM_PricingStatesPtr ( new ARM_PricingStates(1, 2) );
		ARM_GP_Matrix*	modelStates			=	new ARM_GP_Matrix(1,2);
		modelStates->Elt(0,0)=xt;
		modelStates->Elt(0,1)=vt;
		states->SetModelStates(ARM_GP_MatrixPtr(modelStates) );

		std::complex<double> k(k_real,k_imag);

		std::complex<double> res = mod->CptLaplace(	evalTime,
													startTime, 
													endTime,
													states,
													k	);
		
		if ( isReal )
			result.setDouble( std::real(res ) );
		else
			result.setDouble( std::imag(res ) );
	}
	catch(Exception& x){
		x.DebugPrint();
		ARM_RESULT();
	}

	catch (...){
		result.setMsg("ARM_ERR: unrecognized failure in ARMLOCAL_CIRBondPrice");
        return(ARM_KO);
	}
	return ARM_OK;
}

extern long ARMLOCAL_InfEqHwSV_Density(
	const long&				modelId,
	const double&			evalTime,
	const double&			startTime,
	const double&			endTime,
	const double&			xt,
	const double&			vt,
	const double&			x,
	const double&			period,
	const double&			frequency,
	ARM_result&				result ){

	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	
	try	{

		ARM_Object* obj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelId );
		ARM_EQHWSV* mod = dynamic_cast<ARM_EQHWSV*>(obj);
		if( !mod ){
			result.setMsg ("ARM_ERR: model is not of a good type");
			return ARM_KO;
		}

		ARM_PricingStatesPtr states			=	ARM_PricingStatesPtr ( new ARM_PricingStates(1, 2) );
		ARM_GP_Matrix*	modelStates			=	new ARM_GP_Matrix(1,2);
		modelStates->Elt(0,0)=xt;
		modelStates->Elt(0,1)=vt;
		states->SetModelStates(ARM_GP_MatrixPtr(modelStates) );

		
		ARM_VectorPtr res = mod->CptDensity( 	evalTime,
												startTime, 
												endTime,
												x, 
												states,
												period,
												frequency);

		result.setDouble(res->Elt(0) );
	}
	catch(Exception& x){
		x.DebugPrint();
		ARM_RESULT();
	}

	catch (...){
		result.setMsg("ARM_ERR: unrecognized failure in ARMLOCAL_CIRBondPrice");
        return(ARM_KO);
	}
	return ARM_OK;
}