#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_irindex.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_LIBOR (LPXLOPER XL_liborType,
												   LPXLOPER XL_ccy,
												   LPXLOPER XL_resetFreq,
												   LPXLOPER XL_payFreq,
                                                   LPXLOPER XL_daycount,
												   LPXLOPER XL_intRule)
{
	ADD_LOG("Local_LIBOR ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_liborType;
		long liborTypeId;
		
		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_resetFreq;
		long resetFreqId;

		CCString C_payFreq;
		long payFreqId;

		CCString C_dayCount;
		long dayCount;

		CCString C_intRule;
		long intRuleId;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
		XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
		XL_readStrCellWD(XL_daycount,C_dayCount,"-1"," ARM_ERR: day count: string expected",C_result);
		XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: intRule: string expected",C_result);
			
		if (( C_ccy.GetLen() > 3 )
			&& 
			( !(C_ccy == "DEFAULT"))
		   )
		   ccyIsObject = true;

		if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if (strcmp(C_dayCount,"-1") == 0)
		{
			dayCount = -1;
		}
		else
		{
			dayCount = ARM_ConvDayCount (C_dayCount);
		}

		intRuleId = ARM_ConvIntRule (C_intRule);

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_IRINDEX_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_LIBOR (liborTypeId, ccyIsObject, C_ccy, resetFreqId, payFreqId,dayCount,intRuleId, C_result);

			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();

				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if(curClass == prevClass)
			{
				retCode = ARMLOCAL_LIBOR (liborTypeId, ccyIsObject, C_ccy, resetFreqId, payFreqId,dayCount,intRuleId, C_result, objId);

				if (retCode == ARM_OK)
				{			
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
 			}
			else
			{
				FreeCurCellContent ();
				retCode = ARMLOCAL_LIBOR (liborTypeId, ccyIsObject, C_ccy, resetFreqId, payFreqId,dayCount,intRuleId, C_result);

				if(retCode == ARM_OK)
				{
					objId = C_result.getLong ();
				
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
		}

		if(retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LIBOR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LIBOR (LPXLOPER XL_liborType,
													   LPXLOPER XL_ccy,
													   LPXLOPER XL_resetFreq,
													   LPXLOPER XL_payFreq,
                                                       LPXLOPER XL_daycount,
													   LPXLOPER XL_intRule)
{
	ADD_LOG("Local_PXL_LIBOR ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_liborType;
		long liborTypeId;

		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_resetFreq;
		long resetFreqId;

		CCString C_payFreq;
		long payFreqId;

		CCString C_dayCount;
		long dayCount;

		CCString C_intRule;
		long intRuleId;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
		XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
		XL_readStrCellWD(XL_daycount,C_dayCount,"-1"," ARM_ERR: day count: string expected",C_result);
		XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: intRule: string expected",C_result);
			
		if (( C_ccy.GetLen() > 3 )
			&& 
			( !(C_ccy == "DEFAULT"))
		   )
		   ccyIsObject = true;

		if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if (strcmp(C_dayCount,"-1") == 0)
		{
			dayCount = -1;
		}
		else
		{
			dayCount = ARM_ConvDayCount (C_dayCount);
		}

		intRuleId = ARM_ConvIntRule (C_intRule);

		long retCode;
		long objId;

		CCString curClass = LOCAL_IRINDEX_CLASS;
		CCString stringId;

		retCode = ARMLOCAL_LIBOR (liborTypeId, ccyIsObject, C_ccy, resetFreqId, payFreqId,dayCount, intRuleId, C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
		}

		if(retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_LIBOR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_IRINDEX (LPXLOPER XL_dayCount,
													 LPXLOPER XL_payFreq,
													 LPXLOPER XL_maturity,
													 LPXLOPER XL_compMethod,
													 LPXLOPER XL_fwdRule,
													 LPXLOPER XL_resetTiming,
													 LPXLOPER XL_resetGap,
													 LPXLOPER XL_payTiming,
													 LPXLOPER XL_payGap,
													 LPXLOPER XL_ccy,
													 LPXLOPER XL_indexType,
													 LPXLOPER XL_decompFreq,
													 LPXLOPER XL_intRule,
													 LPXLOPER XL_resetFreq)
{
	ADD_LOG("Local_IRINDEX ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_dayCount;
		long dayCountId;

		CCString C_payFreq;
		long payFreqId;

		CCString C_resetFreq;
		long resetFreqId;

		double C_maturity;
		double C_maturity_default = -1;

		CCString C_compMethod;
		long compMethodId;

		CCString C_fwdRule;
		long fwdRuleId;

		CCString C_resetTiming;
		long resetTimingId;

		double C_resetGap;
		double C_resetGap_default = 10000.;

		CCString C_payTiming;
		long payTimingId;

		double C_payGap;
		double C_payGap_default = 10000.;

		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_indexType;
		long indexTypeId;

		double C_decompFreq;
		double C_decompFreq_default = K_COMP_PROP;
			
		CCString C_intRule;
		long intRuleId;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCellWD(XL_dayCount,C_dayCount,"-1"," ARM_ERR: day count: string expected",C_result);
		XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
		XL_readNumCellWD(XL_maturity,C_maturity,C_maturity_default," ARM_ERR: term: numeric expected",C_result);
		XL_readStrCellWD(XL_compMethod,C_compMethod,"P"," ARM_ERR: computing method: string expected",C_result);
		XL_readStrCellWD(XL_fwdRule,C_fwdRule,"MF"," ARM_ERR: forward rule: string expected",C_result);
		XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing: string expected",C_result);
		XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: reset gap: numeric expected",C_result);
		XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: pay timing: string expected",C_result);
		XL_readNumCellWD(XL_payGap,C_payGap,C_payGap_default," ARM_ERR: pay gap: numeric expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_indexType,C_indexType,"FIXED"," ARM_ERR: index type: string expected",C_result);
		XL_readNumCellWD(XL_decompFreq,C_decompFreq,C_decompFreq_default," ARM_ERR: decomp frequency: string expected",C_result);
		XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: intRule: string expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
		
		if (( C_ccy.GetLen() > 3 )
			&& 
			( !(C_ccy == "DEFAULT"))
		   )
		   ccyIsObject = true;

		indexTypeId = ARM_ConvIrType (C_indexType);

		if (strcmp(C_dayCount,"-1") == 0)
		{
			dayCountId = -1;
		}
		else
		{
			dayCountId = ARM_ConvDayCount (C_dayCount);
		}

		if( strcmp(C_resetFreq , "-1") == 0 )
		{
			C_resetFreq = C_payFreq;
		}

		if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((compMethodId = ARM_ConvCompMeth (C_compMethod, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((fwdRuleId = ARM_ConvFwdRule (C_fwdRule, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

		payTimingId = ARM_ConvPayResetRule (C_payTiming);

		indexTypeId = ARM_ConvIrType (C_indexType);

		intRuleId = ARM_ConvIntRule (C_intRule);

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_IRINDEX_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_IRINDEX (dayCountId, payFreqId,
								   C_maturity, compMethodId,
								   fwdRuleId, resetTimingId,
								   C_resetGap, payTimingId,
								   C_payGap, ccyIsObject, C_ccy,
								   indexTypeId, (long)C_decompFreq,
								   intRuleId,
								   resetFreqId,
								   C_result);

			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();

				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if(curClass == prevClass)
			{
				retCode = ARMLOCAL_IRINDEX (dayCountId, payFreqId,
									   C_maturity, compMethodId,
									   fwdRuleId, resetTimingId,
									   (long)C_resetGap, payTimingId,
									   (long)C_payGap, ccyIsObject, C_ccy,
									   indexTypeId, (long)C_decompFreq,
									   intRuleId,
									   resetFreqId,
									   C_result, objId);

				if(retCode == ARM_OK)
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent ();
				retCode = ARMLOCAL_IRINDEX (dayCountId, payFreqId,
									   C_maturity, compMethodId,
									   fwdRuleId, resetTimingId,
									   (long)C_resetGap, payTimingId,
									   (long)C_payGap, ccyIsObject, C_ccy,
									   indexTypeId, (long)C_decompFreq,
									   intRuleId,
									   resetFreqId,
									   C_result);
 				
				if(retCode == ARM_OK)
				{
					objId = C_result.getLong ();
				
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
		}

		if ( retCode == ARM_OK )
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IRINDEX" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}	


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IRINDEX (LPXLOPER XL_dayCount,
														 LPXLOPER XL_payFreq,
														 LPXLOPER XL_maturity,
														 LPXLOPER XL_compMethod,
														 LPXLOPER XL_fwdRule,
														 LPXLOPER XL_resetTiming,
														 LPXLOPER XL_resetGap,
														 LPXLOPER XL_payTiming,
														 LPXLOPER XL_payGap,
														 LPXLOPER XL_ccy,
														 LPXLOPER XL_indexType,
														 LPXLOPER XL_decompFreq,
														 LPXLOPER XL_intRule,
														 LPXLOPER XL_resetFreq)
{
	ADD_LOG("Local_PXL_IRINDEX ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_dayCount;
		long dayCountId;

   		CCString C_payFreq;
		long payFreqId;

		CCString C_resetFreq;
		long resetFreqId;

		double C_maturity;
		double C_maturity_default = -1;

		CCString C_compMethod;
		long compMethodId;

		CCString C_fwdRule;
		long fwdRuleId;

		CCString C_resetTiming;
		long resetTimingId;

		double C_resetGap;
		double C_resetGap_default = 10000.;

		CCString C_payTiming;
		long payTimingId;

		double C_payGap;
		double C_payGap_default = 10000.;

		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_indexType;
		long indexTypeId;

		double C_decompFreq;
		double C_decompFreq_default = K_COMP_PROP;
			
		CCString C_intRule;
		long intRuleId;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCellWD(XL_dayCount,C_dayCount,"-1"," ARM_ERR: day count: string expected",C_result);
		XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
		XL_readNumCellWD(XL_maturity,C_maturity,C_maturity_default," ARM_ERR: term: numeric expected",C_result);
		XL_readStrCellWD(XL_compMethod,C_compMethod,"P"," ARM_ERR: computing method: string expected",C_result);
		XL_readStrCellWD(XL_fwdRule,C_fwdRule,"MF"," ARM_ERR: forward rule: string expected",C_result);
		XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing: string expected",C_result);
		XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: reset gap: numeric expected",C_result);
		XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: pay timing: string expected",C_result);
		XL_readNumCellWD(XL_payGap,C_payGap,C_payGap_default," ARM_ERR: pay gap: numeric expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_indexType,C_indexType,"FIXED"," ARM_ERR: index type: string expected",C_result);
		XL_readNumCellWD(XL_decompFreq,C_decompFreq,C_decompFreq_default," ARM_ERR: decomp frequency: string expected",C_result);
		XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: intRule: string expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);

		if (( C_ccy.GetLen() > 3 )
			&& 
			( !(C_ccy == "DEFAULT"))
		   )
		   ccyIsObject = true;

		indexTypeId = ARM_ConvIrType (C_indexType);

		if (strcmp(C_dayCount,"-1") == 0)
		{
			dayCountId = -1;
		}
		else
		{
			dayCountId = ARM_ConvDayCount (C_dayCount);
		}
		
		if( strcmp(C_resetFreq , "-1") == 0 )
		{
			C_resetFreq = C_payFreq;
		}
		
		if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((compMethodId = ARM_ConvCompMeth (C_compMethod, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((fwdRuleId = ARM_ConvFwdRule (C_fwdRule, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

		payTimingId = ARM_ConvPayResetRule (C_payTiming);

		indexTypeId = ARM_ConvIrType (C_indexType);

		intRuleId = ARM_ConvIntRule (C_intRule);

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_IRINDEX_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_IRINDEX (dayCountId, payFreqId,
							   C_maturity, compMethodId,
							   fwdRuleId, resetTimingId,
							   C_resetGap, payTimingId,
							   C_payGap, ccyIsObject, C_ccy,
							   indexTypeId, (long)C_decompFreq,
							   intRuleId,
							   resetFreqId,
							   C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
		}

		if ( retCode == ARM_OK )
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_IRINDEX" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_MultiIrindex (LPXLOPER XL_irIndexVec,
														  LPXLOPER XL_weightVec)
{
	ADD_LOG("Local_MultiIrindex ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		VECTOR<CCString> C_irindexVect;
	   	VECTOR<double> C_weightVect;

		// error
		static int error;
		static char* reason = "";
	
		XL_readStrVector(XL_irIndexVec,C_irindexVect," ARM_ERR: irindex vector: array of objets expected",DOUBLE_TYPE, C_result);
		XL_readNumVector(XL_weightVec,C_weightVect," ARM_ERR: weight vector: array of numeric expected", C_result);

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_MULTIINDEX_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if(!stringId)
		{
			retCode =  ARMLOCAL_MultiIrindex(C_irindexVect,
											 C_weightVect, 
											 C_result);
			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();

				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if(curClass == prevClass)
			{
				retCode = ARMLOCAL_MultiIrindex (	C_irindexVect,
													C_weightVect,
													C_result, 
													objId);

				if (retCode == ARM_OK)
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent ();
				retCode = ARMLOCAL_MultiIrindex (C_irindexVect,
												 C_weightVect,
												 C_result, 
												 objId);
			
				if(retCode == ARM_OK)
				{
					objId = C_result.getLong ();
				
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
		}

		if(retCode == ARM_OK)
		{
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MultiIrindex" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MultiIrindex (LPXLOPER XL_irIndexVec,
														  LPXLOPER XL_weightVec)
{
	ADD_LOG("Local_PXL_MultiIrindex ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		VECTOR<CCString> C_irindexVect;
	   	VECTOR<double> C_weightVect;

		// error
		static int error;
		static char* reason = "";
	
		XL_readStrVector(XL_irIndexVec,C_irindexVect," ARM_ERR: irindex vector: array of objets expected",DOUBLE_TYPE, C_result);
		XL_readNumVector(XL_weightVec,C_weightVect," ARM_ERR: weight vector: array of numeric expected", C_result);

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_MULTIINDEX_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		retCode =  ARMLOCAL_MultiIrindex(C_irindexVect,
											 C_weightVect, 
											 C_result);
		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
		}

		if(retCode == ARM_OK)
		{
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_MultiIrindex" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_FixedIndex(LPXLOPER XL_dayCount,
													   LPXLOPER XL_ccy)
{
	ADD_LOG("Local_FixedIndex");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_dayCount;
		long dayCountId;

		CCString C_ccy;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCellWD(XL_dayCount,C_dayCount,"-1"," ARM_ERR: day count: string or numeric expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: ccy: string expected",C_result);
		
		if(C_ccy == "DEFAULT")
		{
			ARM_result ccyres;
			ARMLOCAL_ARM_GetDefaultCurrency (ccyres);
			if(ccyres.getRetCode () != ARM_OK)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			else
			{
				C_ccy = ccyres.getString ();
			}
		}

		dayCountId = ARM_ConvDayCount(C_dayCount);
		
		long retCode;
		long objId;
		CCString prevClass;

		CCString curClass = LOCAL_IRINDEX_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if (!stringId)
		{
			retCode = ARMLOCAL_FixedIndex(dayCountId,C_ccy,C_result);

			if ( retCode == ARM_OK )
			{
				objId = C_result.getLong ();

				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			prevClass = LocalGetStringObjectClass (stringId);

			objId = LocalGetNumObjectId (stringId);

			if ( curClass == prevClass )
			{
				retCode = ARMLOCAL_FixedIndex(dayCountId, C_ccy, C_result, objId);

				if ( retCode == ARM_OK )
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent ();

				retCode = ARMLOCAL_FixedIndex(dayCountId, C_ccy, C_result);

				if ( retCode == ARM_OK )
				{
					objId = C_result.getLong ();

					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
		}

		if ( retCode == ARM_OK )
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FixedIndex" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FixedIndex(LPXLOPER XL_dayCount,
														   LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_FixedIndex");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_dayCount;
		long dayCountId;

		CCString C_ccy;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCellWD(XL_dayCount,C_dayCount,"1"," ARM_ERR: day count: string or numeric expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: ccy: string expected",C_result);
		
		if(C_ccy == "DEFAULT")
		{
			ARM_result ccyres;
			ARMLOCAL_ARM_GetDefaultCurrency (ccyres);
			if(ccyres.getRetCode () != ARM_OK)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			else
			{
				C_ccy = ccyres.getString ();
			}
		}
		
		dayCountId = ARM_ConvDayCount (C_dayCount);

		long retCode;
		long objId;

		CCString curClass = LOCAL_IRINDEX_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_FixedIndex(dayCountId, C_ccy, C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
		}

		if ( retCode == ARM_OK )
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_FixedIndex" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_IRINDEX_MONEY_MARKET (LPXLOPER XL_mmTerm,
																  LPXLOPER XL_ccy)
{
	ADD_LOG("Local_IRINDEX_MONEY_MARKET ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_mmTerm;
		CCString C_ccy;
			
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_mmTerm,C_mmTerm," ARM_ERR: money market term: string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
		
		if(C_ccy == "DEFAULT")
		{
			ARM_result ccyres;
			ARMLOCAL_ARM_GetDefaultCurrency (ccyres);
			if(ccyres.getRetCode () != ARM_OK)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			else
			{
				C_ccy = ccyres.getString ();
			}
		}

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_IRINDEX_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_IRINDEX_MONEY_MARKET (C_mmTerm, C_ccy, C_result);

			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();

				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if(curClass == prevClass)
			{
				retCode = ARMLOCAL_IRINDEX_MONEY_MARKET (C_mmTerm, C_ccy, C_result, objId);

				if(retCode == ARM_OK)
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent ();
				retCode = ARMLOCAL_IRINDEX_MONEY_MARKET (C_mmTerm, C_ccy, C_result);
 				
				if(retCode == ARM_OK)
				{
					objId = C_result.getLong ();

					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
		}

		if(retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IRINDEX_MONEY_MARKET" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IRINDEX_MONEY_MARKET (LPXLOPER XL_mmTerm,
																	  LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_IRINDEX_MONEY_MARKET ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_mmTerm;
		CCString C_ccy;
			
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_mmTerm,C_mmTerm," ARM_ERR: money market term: string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
		
		if(C_ccy == "DEFAULT")
		{
			ARM_result ccyres;
			ARMLOCAL_ARM_GetDefaultCurrency (ccyres);
			if(ccyres.getRetCode () != ARM_OK)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			else
			{
				C_ccy = ccyres.getString ();
			}
		}

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_IRINDEX_CLASS;
		CCString stringId;

		retCode = ARMLOCAL_IRINDEX_MONEY_MARKET (C_mmTerm, C_ccy, C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
		}
		
		if(retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_IRINDEX_MONEY_MARKET" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_CMS(LPXLOPER XL_CMSType,
												LPXLOPER XL_liborType,
												LPXLOPER XL_ccy)
{
	ADD_LOG("Local_CMS");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_CMSType;
		long CMSTypeId;

		CCString C_liborType;
		long liborTypeId;

		CCString C_ccy;
		bool ccyIsObject = false;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_CMSType,C_CMSType," ARM_ERR:  CMS type: string expected",C_result);
		XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		
			
		if (( C_ccy.GetLen() > 3 )
			&& 
			( !(C_ccy == "DEFAULT"))
		   )
		   ccyIsObject = true;

		if((CMSTypeId = ARM_ConvCMIndName (C_CMSType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		

		if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_IRINDEX_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if (!stringId)
		{
			retCode = ARMLOCAL_CMS(CMSTypeId,
								   liborTypeId,
								   ccyIsObject,
								   C_ccy,
								   C_result);

			if ( retCode == ARM_OK )
			{
				objId = C_result.getLong ();

				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			prevClass = LocalGetStringObjectClass (stringId);

			objId = LocalGetNumObjectId (stringId);

			if ( curClass == prevClass )
			{
				retCode = ARMLOCAL_CMS(CMSTypeId,
									   liborTypeId,
									   ccyIsObject,
									   C_ccy,
									   C_result,
									   objId);

			if ( retCode == ARM_OK )
			{
				objId = C_result.getLong ();

				stringId = LocalMakeObjectId (objId, curClass);
			}
			}
			else
			{
				FreeCurCellContent ();

				retCode = ARMLOCAL_CMS(CMSTypeId,
									   liborTypeId,
									   ccyIsObject,
									   C_ccy,
									   C_result);

				if ( retCode == ARM_OK )
				{
					objId = C_result.getLong ();

					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
		}

		if ( retCode == ARM_OK )
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CMS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CMS(LPXLOPER XL_CMSType,
													LPXLOPER XL_liborType,
													LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_CMS");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_CMSType;
		long CMSTypeId;

		CCString C_liborType;
		long liborTypeId;

		CCString C_ccy;
		bool ccyIsObject = false;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_CMSType,C_CMSType," ARM_ERR:  CMS type: string expected",C_result);
		XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: object expected",C_result);
		
		if (( C_ccy.GetLen() > 3 )
			&& 
			( !(C_ccy == "DEFAULT"))
		   )
		   ccyIsObject = true;

		if ((CMSTypeId = ARM_ConvCMIndName (C_CMSType, C_result)) == ARM_DEFAULT_ERR)
		{
		   ARM_ARG_ERR();
		
		   return (LPXLOPER)&XL_result;
		}
		

		if ((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
		   ARM_ARG_ERR();
		
		   return (LPXLOPER)&XL_result;
		}
		
		long retCode;
		long objId;

		
		CCString curClass = LOCAL_IRINDEX_CLASS;
		CCString stringId;

		retCode = ARMLOCAL_CMS(CMSTypeId,
							   liborTypeId,
							   ccyIsObject,
							   C_ccy,
							   C_result);

		if ( retCode == ARM_OK )
		{
		   objId = C_result.getLong ();

		   stringId = LocalMakeObjectId (objId, curClass);
		}

		if ( retCode == ARM_OK )
		{
		   FreeCurCellErr ();
		
		   XL_result.xltype = xltypeStr;
		   XL_result.val.str = XL_StrC2StrPascal (stringId);
		   XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
		   ARM_ERR();
		}

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CMS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_IRINDEX2 (LPXLOPER XL_payFreq,
													 LPXLOPER XL_maturity,
													 LPXLOPER XL_ccy,
													 LPXLOPER XL_indexType,
													 LPXLOPER XL_resetFreq)
{
	ADD_LOG("Local_IRINDEX2 ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
   		CCString C_payFreq;
		long payFreqId;

		CCString C_resetFreq;
		long resetFreqId;

		double C_maturity;
		double C_maturity_default = -1;

		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_indexType;
		long indexTypeId;

		// default values
		long dayCount_default = KACTUAL_360;
    
		long compMethod_default = K_COMP_PROP;

		long fwdRule_default = K_MOD_FOLLOWING;

		long resetTim_default = K_ADVANCE;

		double resetGap_default = 10000.;

		long payTim_default = K_ARREARS;

		double payGap_default = 10000.;

		double decompFreq_default = K_COMP_PROP;
			
		long intRule_default = K_ADJUSTED;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: frequency: string expected",C_result);
		XL_readNumCellWD(XL_maturity,C_maturity,C_maturity_default," ARM_ERR: term: numeric expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_indexType,C_indexType,"FIXED"," ARM_ERR: index type: string expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: frequency: string expected",C_result);

		if (( C_ccy.GetLen() > 3 )
			&& 
			( !(C_ccy == "DEFAULT"))
		   )
		   ccyIsObject = true;

		if( strcmp(C_resetFreq , "-1") == 0 )
		{
			C_resetFreq = C_payFreq;
		}

		if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		indexTypeId = ARM_ConvIrType (C_indexType);

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_IRINDEX_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_IRINDEX (dayCount_default, payFreqId,
										C_maturity, compMethod_default,
										fwdRule_default, resetTim_default,
									   (long) resetGap_default, payTim_default,
									   (long) payGap_default, ccyIsObject, C_ccy,
									   indexTypeId, (long) decompFreq_default,
									   intRule_default,resetFreqId,
									   C_result);

			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();

				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if(curClass == prevClass)
			{
				retCode = ARMLOCAL_IRINDEX (dayCount_default, payFreqId,
											C_maturity, compMethod_default,
											fwdRule_default, resetTim_default,
										   (long) resetGap_default, payTim_default,
										   (long) payGap_default, ccyIsObject, C_ccy,
										   indexTypeId, (long) decompFreq_default,
										   intRule_default,resetFreqId,
										   C_result, objId);

				if(retCode == ARM_OK)
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent ();
				retCode = ARMLOCAL_IRINDEX (dayCount_default, payFreqId,
											C_maturity, compMethod_default,
											fwdRule_default, resetTim_default,
										   (long) resetGap_default, payTim_default,
										   (long) payGap_default, ccyIsObject, C_ccy,
										   indexTypeId, (long) decompFreq_default,
										   intRule_default,resetFreqId,
										   C_result);
 				
				if(retCode == ARM_OK)
				{
					objId = C_result.getLong ();
				
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
		}

		if ( retCode == ARM_OK )
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}
//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IRINDEX2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IRINDEX2 (LPXLOPER XL_payFreq,
												    	 LPXLOPER XL_maturity,
													     LPXLOPER XL_ccy,
													     LPXLOPER XL_indexType,
														 LPXLOPER XL_resetFreq)
{
	ADD_LOG("Local_PXL_IRINDEX2 ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
   		CCString C_payFreq;
		long payFreqId;

		CCString C_resetFreq;
		long resetFreqId;

		double C_maturity;
		double C_maturity_default = -1;

		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_indexType;
		long indexTypeId;

		// default values
		long dayCount_default = KACTUAL_360;
    
		long compMethod_default = K_COMP_PROP;

		long fwdRule_default = K_MOD_FOLLOWING;

		long resetTim_default = K_ADVANCE;

		double resetGap_default = 10000.;

		long payTim_default = K_ARREARS;

		double payGap_default = 10000.;

		double decompFreq_default = K_COMP_PROP;
			
		long intRule_default = K_ADJUSTED;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCellWD(XL_payFreq,C_payFreq, "-1"," ARM_ERR: pay frequency: string expected",C_result);
		XL_readNumCellWD(XL_maturity,C_maturity,C_maturity_default," ARM_ERR: term: numeric expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_indexType,C_indexType,"FIXED"," ARM_ERR: index type: string expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
		
		if (( C_ccy.GetLen() > 3 )
			&& 
			( !(C_ccy == "DEFAULT"))
		   )
		   ccyIsObject = true;

		indexTypeId = ARM_ConvIrType (C_indexType);

		if( strcmp(C_resetFreq , "-1") == 0 )
		{
			C_resetFreq = C_payFreq;
		}

		if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long retCode;
		long objId;

		
		CCString curClass = LOCAL_IRINDEX_CLASS;
		CCString stringId;

		retCode = ARMLOCAL_IRINDEX (dayCount_default, payFreqId,
									C_maturity, compMethod_default,
									fwdRule_default, resetTim_default,
								   (long) resetGap_default, payTim_default,
								   (long) payGap_default, ccyIsObject, C_ccy,
								   indexTypeId, (long) decompFreq_default,
								   intRule_default,resetFreqId,
								   C_result);

		if ( retCode == ARM_OK )
		{
		   objId = C_result.getLong ();

		   stringId = LocalMakeObjectId (objId, curClass);
		}

		if ( retCode == ARM_OK )
		{
		   FreeCurCellErr ();
		
		   XL_result.xltype = xltypeStr;
		   XL_result.val.str = XL_StrC2StrPascal (stringId);
		   XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
		   ARM_ERR();
		}


//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_IRINDEX2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
