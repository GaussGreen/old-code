/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_base_local.cpp,v $
 * Revision 1.1  2004/02/07 15:08:43  ebenhamou
 * Initial version
 *
 */

#include <ARM\libarm_local\firstToBeIncluded.h>
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include "ARM_xl_gp_base_local.h"
#include <libCCxll\CCxll.h>
#include <ARM\libarm_local\ARM_local_gp_base.h>
#include "ARM_xl_wrapper_local.h"
#include "ARM_xl_trycatch_local.h"
#include "ARM_local_interface.h"
#include <util\fromto.h>
#include "ARM_xl_gp_fctorhelper.h"

#include "util\tech_macro.h"

/*!
 * functor for the DateStrip Generation
 */
class dateStripFunc : public ARMResultLong2LongFunc
{
public:
	/// store the context
	dateStripFunc(
		double startDate,
		double endDate,
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
		int resetTiming,			
		int payTiming,				
		int adjFirstDate,
		double refDate,
        long accruedFull,
		int firstDateFwdRule)
	:
		C_startDate( startDate ),
		C_endDate( endDate ),
		C_resetFreq( resetFreq ),
		C_dayCount( dayCount ),
		C_resetCalendar( resetCalendar ),
		C_fwdRule( fwdRule ),
		C_intRule( intRule ),
		C_stubRule( stubRule ),
		C_resetGap( resetGap ),
		C_payFreq( payFreq ),
		C_payGap( payGap ),
		C_payCalendar( payCalendar ),
		C_resetTiming( resetTiming ),
		C_payTiming( payTiming ),
		C_adjFirstDate( adjFirstDate ),
		C_refDate( refDate ),
        C_accruedFull(accruedFull),
		C_firstDateFwdRule(firstDateFwdRule)
	{};

	/// call the proper function
	long operator()( ARM_result& result, long objId ) 
	{
        char* resetCalendarChar = C_resetCalendar.c_str();
        char* payCalendarChar = C_payCalendar.c_str();
		long success = ARMLOCAL_DateStrip_Create(
			C_startDate,
			C_endDate,
			C_resetFreq,
			C_dayCount,
			resetCalendarChar,
			C_fwdRule,
			C_intRule,
			C_stubRule,
			C_resetGap,
			C_payFreq,
			C_payGap,
			payCalendarChar,
			C_resetTiming,
			C_payTiming,
			C_adjFirstDate,
			C_refDate,
            C_accruedFull,
			C_firstDateFwdRule,
			result,
			objId );
        delete resetCalendarChar;
        delete payCalendarChar;
        return success;
        
	}
private:
	double C_startDate;
	double C_endDate;
	int C_resetFreq;
	int C_dayCount;
	CCString C_resetCalendar;	
	int C_fwdRule;
	int C_intRule;
	int C_stubRule;
	int C_resetGap;				
	int C_payFreq;				
	int C_payGap;					
	CCString C_payCalendar;	
	int C_resetTiming;			
	int C_payTiming;				
	int C_adjFirstDate;
	double C_refDate;
    double C_accruedFull;
	int C_firstDateFwdRule;
};



/*!
 * central function that does the creation of the DateStrip
 * using the LPXLOPER arguments
 *
 */
LPXLOPER Local_DateStripCreateCommon(
	LPXLOPER XL_startDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_resetFreq,
	LPXLOPER XL_dayCount,
	LPXLOPER XL_resetCalendar,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_resetGap,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_resetTiming,
	LPXLOPER XL_payTiming,
	LPXLOPER XL_adjFirstDate,
	LPXLOPER XL_refDate,
    LPXLOPER XL_accruedFull,
	LPXLOPER XL_firstDateFwdRule,
	bool PersistentInXL )
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
			
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();

	double C_startDate;
	double C_endDate;
	long C_resetFreq;
	long C_dayCount;
	CCString C_resetCalendar;
	long C_fwdRule;
	long C_intRule;
	long C_stubRule;
	double C_resetGap;
	long C_payFreq;
	double C_payGap;
	CCString C_payCalendar;
	long C_resetTiming;
	long C_payTiming;
	double C_adjFirstDate;
	double gapDefault = GETDEFAULTVALUE;
	double adjFirstDateDefaultValue = -1.0;
	double C_refDate;
	double refDateDefault = GETDEFAULTVALUE;
    double C_accruedFull;
    double accruedFullDefault = K_ACCRUED;
	CCString C_firstDateFwdRuleStr;
	long C_firstDateFwdRule;

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";
	
	XL_readNumCell(XL_startDate,C_startDate,			" ARM_ERR: start Date: date expected",		C_result );
	XL_readNumCell(XL_endDate,C_endDate,				" ARM_ERR: end Date: date expected",		C_result );
	XL_GETFREQUENCYWD(XL_resetFreq,C_resetFreq, "A",	" ARM_ERR: reset freq: string expected",	C_result );
	XL_GETDAYCOUNTWD( XL_dayCount, C_dayCount,"ACTUAL",	" ARM_ERR: dayCount : string expected",		C_result );
	XL_readStrCellWD(XL_resetCalendar, C_resetCalendar, GETDEFAULTVALUESTR," ARM_ERR: reset Calendar: string expected",C_result);
	XL_GETFWDRULEWD( XL_fwdRule, C_fwdRule, "MF",		" ARM_ERR: fwdRule : string expected",		C_result );
	XL_GETINTRULEWD( XL_intRule, C_intRule, "ADJ",		" ARM_ERR: fwdRule : string expected",		C_result );
	XL_GETCONVRULEWD( XL_stubRule, C_stubRule, "SS",	" ARM_ERR: stub Rule : string expected",	C_result );
	XL_readNumCellWD(XL_resetGap,C_resetGap,gapDefault, " ARM_ERR: reset gap: numeric expected",	C_result);
	XL_GETFREQUENCYWD( XL_payFreq, C_payFreq, "-1",		" ARM_ERR: pay frequency: string expected", C_result );
	XL_readNumCellWD(XL_payGap, C_payGap,gapDefault,	" ARM_ERR: pay gap: numeric expected",		C_result);
	XL_readStrCellWD(XL_payCalendar, C_payCalendar,	GETDEFAULTVALUESTR," ARM_ERR: pay Calendar: string expected",	C_result);
	XL_GETPAYRESETTIMINGWD(XL_resetTiming,C_resetTiming,"ADV","ARM_ERR: reset timing: string expected",	C_result );
	XL_GETPAYRESETTIMINGWD( XL_payTiming, C_payTiming, "ARR","ARM_ERR: pay timing: string expected",		C_result );
	XL_readNumCell(XL_endDate,C_endDate,				" ARM_ERR: end Date: date expected",		C_result );
	XL_readNumCellWD(XL_adjFirstDate,C_adjFirstDate,adjFirstDateDefaultValue," ARM_ERR: adjust First Date: numeric expected",C_result);
	XL_readNumCellWD(XL_refDate, C_refDate,	refDateDefault, " ARM_ERR: refDate: date expected", C_result );
	XL_readNumCellWD(XL_accruedFull, C_accruedFull,	accruedFullDefault, " ARM_ERR: AccruedOrFull: numeric expected", C_result );
	XL_readStrCellWD( XL_firstDateFwdRule, C_firstDateFwdRuleStr, GETDEFAULTVALUESTR, " ARM_ERR: firstDateFwdRule : string expected",	C_result );

	if (C_adjFirstDate == -1)
		C_adjFirstDate = C_intRule;

	if (C_firstDateFwdRuleStr == GETDEFAULTVALUESTR)
		C_firstDateFwdRule = C_fwdRule;
	else
		C_firstDateFwdRule = ARM_ConvFwdRule(C_firstDateFwdRuleStr, C_result);

	/// standard interface to keep track of what has been already defined
	/// in excel
	CCString curClass = LOCAL_DATESTRIP_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	dateStripFunc ourFunc(
			C_startDate,
			C_endDate,
			C_resetFreq,
			C_dayCount,
			C_resetCalendar,
			C_fwdRule,
			C_intRule,
			C_stubRule,
			C_resetGap,
			C_payFreq,
			C_payGap,
			C_payCalendar,
			C_resetTiming,
			C_payTiming,
			C_adjFirstDate,
			C_refDate,
            C_accruedFull,
			C_firstDateFwdRule);

	/// call the general function
	fillXL_Result( curClass, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DateStripCreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/*!
* version that takes into account 
* previous creation of object
*/

__declspec(dllexport) LPXLOPER WINAPI Local_DateStripCreate(
	LPXLOPER XL_startDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_resetFreq,
	LPXLOPER XL_dayCount,
	LPXLOPER XL_resetCalendar,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_sutbuRule,
	LPXLOPER XL_resetGap,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_resetTiming,
	LPXLOPER XL_payTiming,
	LPXLOPER XL_adjFirstDate,
	LPXLOPER XL_refDate,
    LPXLOPER XL_accruedFull,
	LPXLOPER XL_firstDateFwdRule
)
{
	ADD_LOG("Local_DateStripCreate");
	bool PersistentInXL = true;
	return Local_DateStripCreateCommon(
		XL_startDate,
		XL_endDate,	 
		XL_resetFreq,	 
		XL_dayCount,	 
		XL_resetCalendar,	 
		XL_fwdRule,	 
		XL_intRule,	 
		XL_sutbuRule,	 
		XL_resetGap,	 
		XL_payFreq,	 
		XL_payGap,	 
		XL_payCalendar,	 
		XL_resetTiming,	 
		XL_payTiming,	 
		XL_adjFirstDate,
		XL_refDate,
        XL_accruedFull,
		XL_firstDateFwdRule,
		PersistentInXL );
}

/*!
* version for non persistentInXl
*/
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_DateStripCreate(
	LPXLOPER XL_startDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_resetFreq,
	LPXLOPER XL_dayCount,
	LPXLOPER XL_resetCalendar,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_sutbuRule,
	LPXLOPER XL_resetGap,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_resetTiming,
	LPXLOPER XL_payTiming,
	LPXLOPER XL_adjFirstDate,
	LPXLOPER XL_refDate,
    LPXLOPER XL_accruedFull,
	LPXLOPER XL_firstDateFwdRule)
{
	ADD_LOG("Local_PXL_DateStripCreate");
	bool PersistentInXL = false;
	return Local_DateStripCreateCommon(
		XL_startDate,
		XL_endDate,	 
		XL_resetFreq,	 
		XL_dayCount,	 
		XL_resetCalendar,	 
		XL_fwdRule,	 
		XL_intRule,	 
		XL_sutbuRule,	 
		XL_resetGap,	 
		XL_payFreq,	 
		XL_payGap,	 
		XL_payCalendar,	 
		XL_resetTiming,	 
		XL_payTiming,	 
		XL_adjFirstDate,
		XL_refDate,
        XL_accruedFull,
		XL_firstDateFwdRule,
		PersistentInXL);
}




/////////////////////////////////////////
/// Local_DateStripGetData function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_DataFromDateStrip(
    LPXLOPER XL_DateStripId,
    LPXLOPER XL_DataType )
{
	ADD_LOG("Local_DataFromDateStrip");
	/// to remove memory leak put a result holder!
	static XLOPER_Holder XL_resultHolder;
	XLOPER& XL_result = XL_resultHolder.GetResult();

	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	// error
	static int error;
	static char* reason = "";

	/// Get the variables from the XLOper variables
	long C_dateStripId;
	XL_GETOBJID( XL_DateStripId,	C_dateStripId,	" ARM_ERR: date strip: Object expected",	C_result);

	long C_dataType;
	XL_GETMETHOD(XL_DataType,		C_dataType,		" ARM_ERR: data type: string expected",		C_result, ARM_ConvDateStripDataType );

	VECTOR<double> C_DataResult;
	long retCode;
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
		retCode = ARMLOCAL_DateStripGetData( C_dateStripId, C_dataType, C_DataResult, C_result );
	else
		retCode = ARM_KO;
		
	/// feed the LPXLOPER object result 
	if (retCode == ARM_OK)
	{
		/// add these additional lines 
		/// to display blank lines
		const int additionalLinesNb = 150;
		bool fillWithBlank = true;
		XL_writeNumVectorWithOptions( XL_result, C_DataResult, " ARM_ERR: Could not get result data", C_result, additionalLinesNb, fillWithBlank );
	}
	
	/// be ware that ARM_ERR is a macro
	/// hence the bracket are necessary
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DataFromDateStrip" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////
/// Local_DateStripGetData function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_DataSizeFromDateStrip(
    LPXLOPER XL_DateStripId,
    LPXLOPER XL_DataType )
{
	ADD_LOG("Local_DataSizeFromDateStrip");
	/// static for macros
	static XLOPER XL_result;
	ARM_result C_result;
	
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	// error
	static int error;
	static char* reason = "";

	/// Get the variables from the XLOper variables
	/// Get the variables from the XLOper variables
	long C_dateStripId;
	XL_GETOBJID( XL_DateStripId,	C_dateStripId,	" ARM_ERR: date strip: Object expected",	C_result);

	long C_dataType;
	XL_GETMETHOD(XL_DataType,		C_dataType,		" ARM_ERR: data type: string expected",		C_result, ARM_ConvDateStripDataType );

	VECTOR<double> C_DataResult;
	long retCode;
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
		retCode = ARMLOCAL_DateStripGetData( C_dateStripId, C_dataType, C_DataResult, C_result );
	else
		retCode = ARM_KO;
		
	/// feed the LPXLOPER object result 
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_DataResult.size();

	}
	/// be ware that ARM_ERR is a macro
	/// hence the bracket are necessary
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DataSizeFromDateStrip" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



/*!
 * functor for the DateStripVec Generation
 */
class dateStripVecFunc : public ARMResultLong2LongFunc
{
public:
	/// store the context
	dateStripVecFunc(
		const VECTOR<double>&  FlowStartDatesV,		/// Flow start dates 
		const VECTOR<double>&  FlowEndDatesV,		/// Flow end dates 
		const VECTOR<double>&  FwdStartDatesV,		/// Fwd start dates 
		const VECTOR<double>&  FwdEndDatesV,		/// Fwd end dates 
		const VECTOR<double>&  ResetDatesV,		/// resetDates
		const VECTOR<double>&  PaymentDatesV,		/// paymentDates
		const VECTOR<double>&  InterestDaysV,		/// numbers of days between 2 periods
		const VECTOR<double>&  InterestTermsV		/// interest term... conversion of InterestDays 
	)
	:
		C_FlowStartDatesV( FlowStartDatesV ),
		C_FlowEndDatesV( FlowEndDatesV ),
		C_FwdStartDatesV( FwdStartDatesV ),
		C_FwdEndDatesV( FwdEndDatesV ),
		C_ResetDatesV( ResetDatesV ),
		C_PaymentDatesV( PaymentDatesV ),
		C_InterestDaysV( InterestDaysV ),
		C_InterestTermsV( InterestTermsV )
	{};

	/// call the proper function
	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_DateStripFromVec_Create(
			C_FlowStartDatesV,
			C_FlowEndDatesV,
			C_FwdStartDatesV,
			C_FwdEndDatesV,
			C_ResetDatesV,
			C_PaymentDatesV,
			C_InterestDaysV, 
			C_InterestTermsV,
			result,
			objId );		
	}
private:
	VECTOR<double>  C_FlowStartDatesV;
	VECTOR<double>  C_FlowEndDatesV;
	VECTOR<double>  C_FwdStartDatesV;
	VECTOR<double>  C_FwdEndDatesV;
	VECTOR<double>  C_ResetDatesV;
	VECTOR<double>  C_PaymentDatesV;
	VECTOR<double>  C_InterestDaysV;
	VECTOR<double>  C_InterestTermsV;
};



/*!
 * central function that does the creation of the DateStrip
 * using the LPXLOPER arguments
 *
 */
LPXLOPER Local_DateStripVecCreateCommon(
	LPXLOPER XL_FlowStartDatesV,
	LPXLOPER XL_FlowEndDatesV,
	LPXLOPER XL_FwdStartDatesV,
	LPXLOPER XL_FwdEndDatesV,
	LPXLOPER XL_ResetDatesV,
	LPXLOPER XL_PaymentDatesV,
	LPXLOPER XL_InterestDaysV,
	LPXLOPER XL_InterestTermsV,
	bool PersistentInXL )
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";
	

	/// standard interface to keep track of what has been already defined
	/// in excel
	CCString stringId = GetLastCurCellEnvValue ();

	VECTOR<double> C_FlowStartDatesV;
	XL_readNumVector(XL_FlowStartDatesV,		C_FlowStartDatesV,	" ARM_ERR: FlowStartDates values: array of numeric expected",	C_result);
	VECTOR<double> C_FlowEndDatesV;
	XL_readNumVector(XL_FlowEndDatesV,		C_FlowEndDatesV,	" ARM_ERR: FlowEndDates values: array of numeric expected",		C_result);
	VECTOR<double> C_FwdStartDatesV;
	XL_readNumVector(XL_FwdStartDatesV,		C_FwdStartDatesV,	" ARM_ERR: Fwd StartDates values: array of numeric expected",	C_result);
	VECTOR<double> C_FwdEndDatesV;
	XL_readNumVector(XL_FwdEndDatesV,		C_FwdEndDatesV,	" ARM_ERR: FwdEndDates values: array of numeric expected",		C_result);
	VECTOR<double> C_ResetDatesV;
	XL_readNumVector(XL_ResetDatesV,		C_ResetDatesV,		" ARM_ERR: ResetDates values: array of numeric expected",		C_result);
	VECTOR<double> C_PaymentDatesV;
	XL_readNumVector(XL_PaymentDatesV,		C_PaymentDatesV,	" ARM_ERR: PaymentDates values: array of numeric expected",		C_result);
	VECTOR<double> C_InterestDaysV;
	XL_readNumVector(XL_InterestDaysV,		C_InterestDaysV, 	" ARM_ERR: InterestDays values: array of numeric expected",		C_result);
	VECTOR<double> C_InterestTermsV;
	XL_readNumVector(XL_InterestTermsV,		C_InterestTermsV,	" ARM_ERR: InterestTerms values: array of numeric expected",	C_result);

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	dateStripVecFunc ourFunc(
			C_FlowStartDatesV,
			C_FlowEndDatesV,
			C_FwdStartDatesV,
			C_FwdEndDatesV,
			C_ResetDatesV,
			C_PaymentDatesV,
			C_InterestDaysV, 
			C_InterestTermsV );

	/// call the general function
	fillXL_Result( LOCAL_DATESTRIP_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DateStripVecCreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/*!
* version that takes into account 
* previous creation of object
*/

__declspec(dllexport) LPXLOPER WINAPI Local_DateStripVecCreate(
	LPXLOPER XL_FlowStartDatesV,
	LPXLOPER XL_FlowEndDatesV,
	LPXLOPER XL_FwdStartDatesV,
	LPXLOPER XL_FwdEndDatesV,
	LPXLOPER XL_ResetDatesV,
	LPXLOPER XL_PaymentDatesV,
	LPXLOPER XL_InterestDaysV,
	LPXLOPER XL_InterestTermsV
)
{
	ADD_LOG("Local_DateStripVecCreate");
	bool PersistentInXL = true;
	return Local_DateStripVecCreateCommon(
		XL_FlowStartDatesV,
		XL_FlowEndDatesV,
		XL_FwdStartDatesV,
		XL_FwdEndDatesV,
		XL_ResetDatesV,
		XL_PaymentDatesV,
		XL_InterestDaysV,
		XL_InterestTermsV,
		PersistentInXL );
}

/*!
* version for non persistentInXl
*/
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_DateStripVecCreate(
	LPXLOPER XL_FlowStartDatesV,
	LPXLOPER XL_FlowEndDatesV,
	LPXLOPER XL_FwdStartDatesV,
	LPXLOPER XL_FwdEndDatesV,
	LPXLOPER XL_ResetDatesV,
	LPXLOPER XL_PaymentDatesV,
	LPXLOPER XL_InterestDaysV,
	LPXLOPER XL_InterestTermsV 
)
{
	ADD_LOG("Local_PXL_DateStripVecCreate");
 bool PersistentInXL = false;
 return Local_DateStripVecCreateCommon(
		XL_FlowStartDatesV,
		XL_FlowEndDatesV,
		XL_FwdStartDatesV,
		XL_FwdEndDatesV,
		XL_ResetDatesV,
		XL_PaymentDatesV,
		XL_InterestDaysV,
		XL_InterestTermsV,
		PersistentInXL );
}


__declspec(dllexport) LPXLOPER WINAPI Local_Warning_SetPopUpFlag( LPXLOPER XL_Flag )
{
	ADD_LOG("Local_Warning_SetPopUpFlag");
	/// static for macros
	static XLOPER XL_result;
	ARM_result C_result;
	
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
	
		// error
		static int error;
		static char* reason = "";

		double C_FlagDble;
		XL_readNumCell( XL_Flag, C_FlagDble,	" ARM_ERR: Flag: boolean compatible expected",	C_result);
		bool C_Flag= C_FlagDble != 0;
		long retCode = ARMLOCAL_SetWarningFlag( C_Flag, C_result );

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype  = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (C_result.getString() );
			XL_result.xltype |= xlbitDLLFree;
		}
		else { ARM_ERR(); }
		
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Warning_SetPopUpFlag" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_Warning_Activate( LPXLOPER XL_Flag )
{
	ADD_LOG("Local_Warning_Activate");
	/// static for macros
	static XLOPER XL_result;
	ARM_result C_result;
	
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
	
		// error
		static int error;
		static char* reason = "";

		double C_FlagDble;
		XL_readNumCell( XL_Flag, C_FlagDble,	" ARM_ERR: Flag: boolean compatible expected",	C_result);
		bool C_Flag= C_FlagDble != 0;
		long retCode = ARMLOCAL_Warning_Activate( C_Flag, C_result );

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype  = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (C_result.getString() );
			XL_result.xltype |= xlbitDLLFree;
		}
		else { ARM_ERR(); }
		
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Warning_Activate" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
