/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_base.cpp,v $
 * Revision 1.1  2004/08/20 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_local_gp_base.cpp,
 *
 *  \brief file for the generic calibration local addins functions
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date August 2004
 */

#include "firstToBeIncluded.h"
#include <GP_Base\gpbase\datestrip.h>
#include <GP_Base\gpbase\gpvector.h>
#include <GP_Base\gpbase\gplinalgconvert.h>
#include <GP_Base\gpbase\warning.h>
#include <GP_Base\gpbase\warningkeeper.h>
#include "ARM_local_wrapper.h"
#include "ARM_local_gp_base.h"
#include "ARM_local_glob.h"
#include <util\fromto.h>		/// to use CreateARMVectorFromVECTOR

//// using the namespace directive to access ARM object!
using ARM::ARM_DateStrip;
using ARM::std::vector<double>;
using ARM::CreateARMVectorFromXLDATEVECTOR;
using ARM::CreateARMGPVectorFromVECTOR;
using ARM::ARM_Warning;
using ARM::ARM_WarningKeeper;


/*!
 * Function to create a date strip
 */
extern long ARMLOCAL_DateStrip_Create(
	double startDate,
	double endDate,
	long resetFreq,				/// reset frequency
	long dayCount,				/// dayCount method used for computing the accrued
	const CCString& resetCalendar,/// calendar used for reset
	/// all this can be defaulted
	long fwdRule,				/// whether fwds are with adjusted dates
	long intRule,				/// whether longerest are K_ADJUSTED
	long stubRule,				/// ability to have K8SHORTSTART etc
	long resetGap,				/// reset gap
	long payFreq,				/// payment frequency
	long payGap,				/// pay gap
	const CCString& payCalendar,/// calendar used for payment
	long resetTiming,			/// whether reset are in arrears or in advance
	long payTiming,				/// whether payment are in arrears or in advance
	long adjFirstdate,			/// adjust the first date to business day
	double refDate,				/// reference date
    long accruedfull,           /// Accrud mode or full mode to generat FwdDates
	long firstDateFwdRule,      /// if adjust 1 date, adjust using this rule
	ARM_result&	result,			/// ARM_Result to return the result
	long objId 					/// objId for persistence reason
)
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;
	
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_DateStrip* dateStrip = NULL;
	
	try
	{
		char myStartDate[20];
		char myEndDate[20];
		char myRefDate[20];

		Local_XLDATE2ARMDATE(startDate, myStartDate);
		Local_XLDATE2ARMDATE(endDate, myEndDate);
		if( refDate == GETDEFAULTVALUE )
			strcpy( myRefDate, GETDEFAULTVALUESTR );
		else
			Local_XLDATE2ARMDATE(refDate, myRefDate);

		int stdSpotDays = GETDEFAULTVALUE;  /// standard spot days default is spotDay of default currency
        int indexTerm   = GETDEFAULTVALUE;

		dateStrip = new ARM_DateStrip( 
			(ARM_Date) myStartDate,
			(ARM_Date) myEndDate,
			resetFreq,
			dayCount,
			resetCalendar,
			fwdRule,
			intRule,
			stubRule,
			resetGap,
			payFreq,
			payGap,
			payCalendar,
			resetTiming,
			payTiming,
			adjFirstdate,
			myRefDate,
			stdSpotDays,
			indexTerm,
            accruedfull,
			firstDateFwdRule);

		/// assign object
		if( !assignObject( dateStrip, result, objId ) )
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
		delete dateStrip;

		x.DebugPrint();
		
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete dateStrip;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


/*!
 * Function to give back data from a datestrip
 */
extern long ARMLOCAL_DateStripGetData(
	long dateStripId,
	long dataType,
	VECTOR<double>& Data,
	ARM_result& result )
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;
	
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		ARM_DateStrip* dateStrip= NULL;
		if( !GetObjectFromId( &dateStrip, dateStripId, ARM_DATESTRIP ) )
		{
			result.setMsg ("ARM_ERR: date strip is not of a good type");
			return ARM_KO;
		};
		
		std::vector<double>& tmpData = dateStrip->GetMemberData( dataType );

		int i;

		/// test whether we have dates in which case, we convert this!
		if(		dataType == K_START_DATES 
			||	dataType == K_END_DATES 
			||  dataType == K_RESET_DATES
			||	dataType == K_PAY_DATES 
			||	dataType == K_FWD_START_DATES
			||	dataType == K_FWD_END_DATES 
			)
		{
			for( i=0; i<tmpData->size(); ++i)
			{
				Data.push_back( JulianToXLDate((*tmpData)[i]) );
			}
		}
		/// otherwise no treatment!
		else
		{
			for( i=0; i<tmpData->size(); ++i)
			{
				Data.push_back( (*tmpData)[i] );
			}
		}

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}






/*!
 * Function to create a date strip
 */
extern long ARMLOCAL_DateStripFromVec_Create(
	const VECTOR<double>&  FlowStartDatesV,		/// Flow start dates and Flow 
												/// end dates are used to compute the
	const VECTOR<double>&  FlowEndDatesV,		/// period of interest 
	const VECTOR<double>&  FwdStartDatesV,		/// fwd start dates
	const VECTOR<double>&  FwdEndDatesV,		/// fwd end dates
	const VECTOR<double>&  ResetDatesV,			/// resetDates
	const VECTOR<double>&  PaymentDatesV,		/// paymentDates
	const VECTOR<double>&  InterestDaysV,		/// numbers of days between 2 periods
	const VECTOR<double>&  InterestTermsV,		/// interest term... conversion of InterestDays 
	ARM_result&	result,			/// ARM_Result to return the result
	long objId 					/// objId for persistence reason
)
{

	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;
	
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_DateStrip* dateStrip = NULL;

	try
	{
		std::vector<double>& FlowStartDates	= CreateARMVectorFromXLDATEVECTOR( FlowStartDatesV );
		std::vector<double>& FlowEndDates		= CreateARMVectorFromXLDATEVECTOR( FlowEndDatesV );
		std::vector<double>& FwdStartDates	= CreateARMVectorFromXLDATEVECTOR( FwdStartDatesV );
		std::vector<double>& FwdEndDates		= CreateARMVectorFromXLDATEVECTOR( FwdEndDatesV );
		std::vector<double>& ResetDates		= CreateARMVectorFromXLDATEVECTOR( ResetDatesV );
		std::vector<double>& PaymentDates		= CreateARMVectorFromXLDATEVECTOR( PaymentDatesV );
		std::vector<double>& InterestDays		= CreateARMGPVectorFromVECTOR( InterestDaysV );
		std::vector<double>& InterestTerms	= CreateARMGPVectorFromVECTOR( InterestTermsV );

		dateStrip = new ARM_DateStrip(
			FlowStartDates,
			FlowEndDates,
			FwdStartDates,
			FwdEndDates,
			ResetDates,
			PaymentDates,
			InterestDays,
			InterestTerms );

		/// because the object clones the vector
		/// we need to delete them!
		delete FlowStartDates;
		delete FlowEndDates;
		delete FwdStartDates;
		delete FwdEndDates;
		delete ResetDates;
		delete PaymentDates;
		delete InterestDays;
		delete InterestTerms;

		/// assign object
		if( !assignObject( dateStrip, result, objId ) )
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
		delete dateStrip;

		x.DebugPrint();
		
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete dateStrip;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}

extern long ARMLOCAL_CreateWarning( const string& message,
	ARM_result&	result,
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_Warning* warning= NULL;
   
	try
	{	
		warning = new ARM_Warning( message );
		
		// assign object
		if( !assignObject( warning, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete warning;
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_SetWarningFlag( 
	const bool val,
	ARM_result&	result )
{
	ARM_Warning::SetPopUp( val );
	result.setString( val? "Popup On" : "Popup Off" );
	return ARM_OK;
}



extern long ARMLOCAL_Warning_Activate( 
	const bool val,
	ARM_result&	result )
{
	ARM_WarningKeeper::KeepWarning = val;
	result.setString( val? "Warning On" : "Warning Off" );
	return ARM_OK;
}

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

