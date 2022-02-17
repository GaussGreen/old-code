#pragma warning(disable : 4786)

#include <libCCxll\CCxll.h>

#include "ARM_xl_utilities.h"
#include <ARM/libarm/ARM_result.h>
#include <ARM/libarm_local/ARM_local_class.h>

#include <ARM/local_xlarm/ARM_xl_trycatch_local.h>
#include <ARM/local_xlarm/ARM_local_interface.h>
#include <ARM/local_xlarm/ARM_local_interglob.h>
#include <util/fromto.h>
#include <ARM/libarm_local/ARM_local_utilities.h>

#include <GP_Base/gpbase/curve.h>
#include <GP_Base/gpbase/curvetypedef.h>

using namespace ARM; 

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GP_CQSO_Create(
	LPXLOPER XL_startDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_FundingCurrency,
	LPXLOPER XL_UnderlyingCurrency,
	LPXLOPER XL_1stIndex,
	LPXLOPER XL_2ndIndex,
	LPXLOPER XL_ResetDates,
	LPXLOPER XL_Notionals,
	LPXLOPER XL_Margins,
	LPXLOPER XL_Strikes,
	LPXLOPER XL_Leverages1,
	LPXLOPER XL_Leverages2,
	LPXLOPER XL_CpnMin,
	LPXLOPER XL_CpnMax,
	LPXLOPER XL_Fees,
	LPXLOPER XL_ScheduleArguments
	)
{
	bool PersistentInXL = false; 

	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;	

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		CCString C_FundingCurrency;
		CCString C_UnderlyingCurrency;
		CCString C_1stIndex;
		CCString C_2ndIndex;
		
		double C_startDate; 
		double C_endDate; 
		vector<double> C_ResetDates,
			C_Notionals,
			C_Margins,
			C_Strikes,
			C_Leverages1,
			C_Leverages2,
			C_CpnMin,
			C_CpnMax,
			C_Fees;
	
//------------------------------------------------------------------------------
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_FundingCurrency,C_FundingCurrency,
			" ARM_ERR: XL_FundingCurrency string attended",C_result);
		XL_readStrCell(XL_UnderlyingCurrency,C_UnderlyingCurrency,
			" ARM_ERR: XL_FundingCurrency string attended",C_result);
		XL_readStrCell(XL_1stIndex,C_1stIndex,
			" ARM_ERR: XL_FundingCurrency string attended",C_result);
		XL_readStrCell(XL_2ndIndex,C_2ndIndex,
			" ARM_ERR: XL_FundingCurrency string attended",C_result);
		
		XL_readNumCell(XL_startDate,C_startDate,
			" ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate,
			" ARM_ERR: end date: date expected",C_result);
	
		XL_readNumVector(XL_ResetDates,C_ResetDates,
			" ARM_ERR: reset dates : array of numeric expected",C_result);
		XL_readNumVector(XL_Notionals,C_Notionals,
			" ARM_ERR: notionals: array of numeric expected",C_result);
		XL_readNumVector(XL_Margins,C_Margins,
			" ARM_ERR: margins : array of numeric expected",C_result);
		XL_readNumVector(XL_Strikes,C_Strikes,
			" ARM_ERR: strikes : array of numeric expected",C_result);
		XL_readNumVector(XL_Leverages1,C_Leverages1,
			" ARM_ERR: leverages1: array of numeric expected",C_result);
		XL_readNumVector(XL_Leverages2,C_Leverages2,
			" ARM_ERR: leverages2 : array of numeric expected",C_result);
		XL_readNumVector(XL_CpnMin,C_CpnMin,
			" ARM_ERR: CpnMin : array of numeric expected",C_result);
		XL_readNumVector(XL_CpnMax,C_CpnMax,
			" ARM_ERR: CpnMax : array of numeric expected",C_result);
		XL_readNumVector(XL_Fees,C_Fees,
			" ARM_ERR: strikes : array of numeric expected",C_result);
	
		VECTOR<CCString> C_prodParams;
		VECTOR<CCString> C_prodParams_default(7);

		C_prodParams_default[0] = "A";
		C_prodParams_default[1] = "30/360";
		C_prodParams_default[2] = "ADJ";
		C_prodParams_default[3] = "ADV";
		C_prodParams_default[4] = "2";
		C_prodParams_default[5] = "EUR";
		C_prodParams_default[6] = "EUR";
	
		// Schedule
		XL_readStrVectorWD(XL_ScheduleArguments,C_prodParams,C_prodParams_default,
			" ARM_ERR: ScheduleArguments vector : array of string expected",DOUBLE_TYPE,C_result);
		
		int iFrequency = ARM_ConvFrequency (C_prodParams[0],C_result);
		int iDayCounter = ARM_ConvDayCount (C_prodParams[1]);
		int iIsAdjusted = ARM_ConvIntRule(C_prodParams[2]);
		int iResetType = ARM_ConvPayResetRule (C_prodParams[3]);
		int iResetGap = atof(C_prodParams[4]);
		CCString sResetCalendar = C_prodParams[5];
		CCString sPaymtCalendar = C_prodParams[6];

		// Payoff
		size_t i,size= C_ResetDates.size(), tmpSize; 
		
		ARM_GP_T_Vector<double> resetDates(size); 
		for(i=0; i<size; ++i)
			resetDates[i] = C_ResetDates[i]; 
		
		//Conversion date XL -> date ARM !!
		ARM_Date reference(01,01,2000);
		int dateGap = reference.GetJulian()-36526;
		resetDates += dateGap; 

		ARM_GP_T_Vector<double> nominal(size);
		tmpSize = C_Notionals.size();  
		if(tmpSize>size)
			tmpSize = size;
		for (i=0; i<tmpSize; ++i)
			nominal[i] = C_Notionals[i]; 
		for(;i<size; ++i)
			nominal[i] = C_Notionals[tmpSize-1];
					
		ARM_GP_T_Vector<double> strikes(size);
		tmpSize = C_Strikes.size();  
		if(tmpSize>size)
			tmpSize = size;
		for (i=0; i<tmpSize; ++i)
			strikes[i] = C_Strikes[i]; 
		for(;i<size; ++i)
			strikes[i] = C_Strikes[tmpSize-1];

		ARM_GP_T_Vector<double> leveragesShort(size);
		tmpSize = C_Leverages2.size();  
		if(tmpSize>size)
			tmpSize = size;
		for (i=0; i<tmpSize; ++i)
			leveragesShort[i] = C_Leverages2[i]; 
		for(;i<size; ++i)
			leveragesShort[i] = C_Leverages2[tmpSize-1];

		ARM_GP_T_Vector<double> leveragesLong(size);
		tmpSize = C_Leverages1.size();  
		if(tmpSize>size)
			tmpSize = size;
		for (i=0; i<tmpSize; ++i)
			leveragesLong[i] = C_Leverages1[i]; 
		for(;i<size; ++i)
			leveragesLong[i] = C_Leverages1[tmpSize-1];

		ARM_GP_T_Vector<double> cpnMin(size);
		tmpSize = C_CpnMin.size();  
		if(tmpSize>size)
			tmpSize = size;
		for (i=0; i<tmpSize; ++i)
			cpnMin[i] = C_CpnMin[i]; 
		for(;i<size; ++i)
			cpnMin[i] = C_CpnMin[tmpSize-1];

		ARM_GP_T_Vector<double> cpnMax(size);
		tmpSize = C_CpnMax.size();  
		if(tmpSize>size)
			tmpSize = size;
		for (i=0; i<tmpSize; ++i)
			cpnMax[i] = C_CpnMax[i]; 
		for(;i<size; ++i)
			cpnMax[i] = C_CpnMax[tmpSize-1];

		ARM_GP_T_Vector<double> margins(size);
		tmpSize = C_Margins.size();  
		if(tmpSize>size)
			tmpSize = size;
		for (i=0; i<tmpSize; ++i)
			margins[i] = C_Margins[i]; 
		for(;i<size; ++i)
			margins[i] = C_Margins[tmpSize-1];

		ARM_GP_T_Vector<double> fees(size);
		tmpSize = C_Fees.size();  
		if(tmpSize>size)
			tmpSize = size;
		for (i=0; i<tmpSize; ++i)
			fees[i] = C_Fees[i]; 
		for(;i<size; ++i)
			fees[i] = C_Fees[tmpSize-1];

		ARM_Curve curv_nominal(resetDates, nominal, new ARM_StepUpLeftOpenCstExtrapol<double,double>());
		ARM_Curve curv_strikes(resetDates, strikes, new ARM_StepUpLeftOpenCstExtrapol<double,double>());
		ARM_Curve curv_leveragesShort(resetDates, leveragesShort, new ARM_StepUpLeftOpenCstExtrapol<double,double>());
		ARM_Curve curv_leveragesLong(resetDates, leveragesLong, new ARM_StepUpLeftOpenCstExtrapol<double,double>());
		ARM_Curve curv_cpnMin(resetDates, cpnMin, new ARM_StepUpLeftOpenCstExtrapol<double,double>());
		ARM_Curve curv_cpnMax(resetDates, cpnMax, new ARM_StepUpLeftOpenCstExtrapol<double,double>());
		ARM_Curve curv_margins(resetDates, margins, new ARM_StepUpLeftOpenCstExtrapol<double,double>());
		ARM_Curve curv_fees(resetDates, fees, new ARM_StepUpLeftOpenCstExtrapol<double,double>());

		CCString prevClass;
		CCString curClass = LOCAL_GENSEC_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		long objId=-1; 
		long retCode; 

		if(!(!stringId)) {
			prevClass = LocalGetStringObjectClass (stringId);
			objId = LocalGetNumObjectId (stringId);
			if(curClass != prevClass)
				FreeCurCellContent ();
		}
		
		retCode = ARMLOCAL_CQSOCREATE(C_startDate, 
					C_endDate,
					iFrequency,
					iDayCounter,
					iIsAdjusted,
					iResetType,
					iResetGap,
					sResetCalendar,
					sPaymtCalendar,
					C_FundingCurrency,
					C_UnderlyingCurrency,
					C_1stIndex,
					C_2ndIndex,
					curv_nominal, 
					curv_strikes, 
					curv_leveragesShort, 
					curv_leveragesLong, 
					curv_cpnMin, 
					curv_cpnMax, 
					curv_margins,
					curv_fees,
					C_result, 
					objId); 
		
		if(retCode == ARM_OK)
		{
			LocalSetCurCellEnvValue (curClass, objId); 
			stringId = LocalMakeObjectId (objId, curClass);

			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (stringId);
			XL_result.xltype |= xlbitDLLFree;

		}
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GP_CQSO_Create" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
