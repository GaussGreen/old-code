#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_glob.h>
#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"
#include "ARM_xl_wrapper_local.h"

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_rev.h>
#include <gpbase\env.h>

#include "ARM_xl_trycatch_local.h"

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_REVERSE (LPXLOPER XL_structSwapLeg,
													 LPXLOPER XL_classSwapLeg,
													 LPXLOPER XL_receiveOrPay,
													 LPXLOPER XL_coupon,
													 LPXLOPER XL_exe,
													 LPXLOPER XL_redemp,
													 LPXLOPER XL_classRedemp,
													 LPXLOPER XL_dualDate,
													 LPXLOPER XL_dualStrike,
													 LPXLOPER XL_dualFlag)
{
	ADD_LOG("Local_REVERSE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_structSwapLeg;
	CCString C_classSwapLeg;

	CCString C_receiveOrPay;
	double receiveOrPayId;

	CCString C_coupon;
	CCString C_exe;
	CCString C_redemp;
	CCString C_classRedemp;

    double C_dualDate;
    double C_dualStrike;

    CCString C_dualFlag;
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_structSwapLeg,C_structSwapLeg," ARM_ERR: struct swap leg id: object expected",C_result);
	XL_readStrCell(XL_classSwapLeg,C_classSwapLeg," ARM_ERR: class swap leg id: object expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrCell(XL_coupon,C_coupon," ARM_ERR: coupon id: object expected",C_result);
	XL_readStrCell(XL_exe,C_exe," ARM_ERR: exercise id: object expected",C_result);
	XL_readStrCell(XL_redemp,C_redemp," ARM_ERR: redemption id: object expected",C_result);
	XL_readStrCell(XL_classRedemp,C_classRedemp," ARM_ERR: redemption class id: object expected",C_result);
    XL_readNumCell(XL_dualDate,C_dualDate," ARM_ERR: dualDate: numeric expected",C_result);
    XL_readNumCell(XL_dualStrike,C_dualStrike," ARM_ERR: dualStrike: numeric expected",C_result);
    XL_readStrCell(XL_dualFlag,C_dualFlag," ARM_ERR: dualFlag: string expected",C_result);
   


	if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_REVERSE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_REVERSE (LocalGetNumObjectId (C_structSwapLeg),
							   LocalGetNumObjectId (C_classSwapLeg),
							   receiveOrPayId,
							   LocalGetNumObjectId (C_coupon),
							   LocalGetNumObjectId (C_exe),
							   LocalGetNumObjectId (C_redemp),
							   LocalGetNumObjectId (C_classRedemp),
                               C_dualDate,
                               C_dualStrike,
                               C_dualFlag,
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
			retCode = ARMLOCAL_REVERSE (LocalGetNumObjectId (C_structSwapLeg),
								   LocalGetNumObjectId (C_classSwapLeg),
								   receiveOrPayId,
								   LocalGetNumObjectId (C_coupon),
								   LocalGetNumObjectId (C_exe),
								   LocalGetNumObjectId (C_redemp),
								   LocalGetNumObjectId (C_classRedemp),
                                   C_dualDate,
                                   C_dualStrike,
                                   C_dualFlag,
								   C_result,
								   objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_REVERSE (LocalGetNumObjectId (C_structSwapLeg),
								   LocalGetNumObjectId (C_classSwapLeg),
								   receiveOrPayId,
								   LocalGetNumObjectId (C_coupon),
								   LocalGetNumObjectId (C_exe),
								   LocalGetNumObjectId (C_redemp),
								   LocalGetNumObjectId (C_classRedemp),
                                   C_dualDate,
                                   C_dualStrike,
                                   C_dualFlag,
								   C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_REVERSE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_REVERSE (LPXLOPER XL_structSwapLeg,
														 LPXLOPER XL_classSwapLeg,
														 LPXLOPER XL_receiveOrPay,
														 LPXLOPER XL_coupon,
														 LPXLOPER XL_exe,
														 LPXLOPER XL_redemp,
														 LPXLOPER XL_classRedemp,
														 LPXLOPER XL_dualDate,
														 LPXLOPER XL_dualStrike,
														 LPXLOPER XL_dualFlag)
{
	ADD_LOG("Local_PXL_REVERSE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_structSwapLeg;
	CCString C_classSwapLeg;

	CCString C_receiveOrPay;
	double receiveOrPayId;

	CCString C_coupon;
	CCString C_exe;
	CCString C_redemp;
	CCString C_classRedemp;

    double C_dualDate;
    double C_dualStrike;

    CCString C_dualFlag;
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_structSwapLeg,C_structSwapLeg," ARM_ERR: struct swap leg id: object expected",C_result);
	XL_readStrCell(XL_classSwapLeg,C_classSwapLeg," ARM_ERR: class swap leg id: object expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrCell(XL_coupon,C_coupon," ARM_ERR: coupon id: object expected",C_result);
	XL_readStrCell(XL_exe,C_exe," ARM_ERR: exercise id: object expected",C_result);
	XL_readStrCell(XL_redemp,C_redemp," ARM_ERR: redemption id: object expected",C_result);
	XL_readStrCell(XL_classRedemp,C_classRedemp," ARM_ERR: redemption class id: object expected",C_result);
    XL_readNumCell(XL_dualDate,C_dualDate," ARM_ERR: dualDate: numeric expected",C_result);
    XL_readNumCell(XL_dualStrike,C_dualStrike," ARM_ERR: dualStrike: numeric expected",C_result);
    XL_readStrCell(XL_dualFlag,C_dualFlag," ARM_ERR: dualFlag: string expected",C_result);
   
	if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_REVERSE_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_REVERSE (LocalGetNumObjectId (C_structSwapLeg),
							   LocalGetNumObjectId (C_classSwapLeg),
							   receiveOrPayId,
							   LocalGetNumObjectId (C_coupon),
							   LocalGetNumObjectId (C_exe),
							   LocalGetNumObjectId (C_redemp),
							   LocalGetNumObjectId (C_classRedemp),
                               C_dualDate,
                               C_dualStrike,
                               C_dualFlag,
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

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_REVERSE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_STRUCTREVERSECOUPON(LPXLOPER XL_date,
																LPXLOPER XL_strike,
																LPXLOPER XL_power,
																LPXLOPER XL_callput,
																LPXLOPER XL_xo,
																LPXLOPER XL_floor,
																LPXLOPER XL_cap)
{
	ADD_LOG("Local_STRUCTREVERSECOUPON");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_date;
	VECTOR<double> C_strike;
	VECTOR<double> C_power;
	VECTOR<double> C_callput;
	double C_xo;
	
    VECTOR<double> C_floor;
    VECTOR<double> C_cap;

	// error
	static int error;
	static char* reason = "";

	XL_readNumVector(XL_date,C_date," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_strike,C_strike," ARM_ERR: strike: array of numeric expected",C_result);
	XL_readNumVector(XL_power,C_power," ARM_ERR: power: array of numeric expected",C_result);
	XL_readNumVector(XL_callput,C_callput," ARM_ERR: callput: array of numeric expected",C_result);
	XL_readNumCell(XL_xo,C_xo," ARM_ERR: xo: numeric expected",C_result);
	XL_readNumVector(XL_floor,C_floor," ARM_ERR: floor: array of numeric expected",C_result);
    XL_readNumVector(XL_cap,C_cap," ARM_ERR: cap: array of numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_REVERSECOUPON_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_STRUCTREVERSECOUPON (C_date,
											C_strike,
										    C_power,
										    C_callput,
											C_xo,
                                            C_floor,
                                            C_cap,
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
			retCode = ARMLOCAL_STRUCTREVERSECOUPON (C_date,
											    C_strike,
											    C_power,
											    C_callput,
											    C_xo,
                                                C_floor,
                                                C_cap,
												C_result,
											    objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_STRUCTREVERSECOUPON (C_date,
												C_strike,
											    C_power,
											    C_callput,
												C_xo,
                                                C_floor,
                                                C_cap,
											    C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_STRUCTREVERSECOUPON" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_STRUCTREVERSECOUPON(LPXLOPER XL_date,
																	LPXLOPER XL_strike,
																	LPXLOPER XL_power,
																	LPXLOPER XL_callput,
																	LPXLOPER XL_xo,
																	LPXLOPER XL_floor,
																	LPXLOPER XL_cap)
{
	ADD_LOG("Local_PXL_STRUCTREVERSECOUPON");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_date;
	VECTOR<double> C_strike;
	VECTOR<double> C_power;
	VECTOR<double> C_callput;
	double C_xo;
	
    VECTOR<double> C_floor;
    VECTOR<double> C_cap;

	// error
	static int error;
	static char* reason = "";

	XL_readNumVector(XL_date,C_date," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_strike,C_strike," ARM_ERR: strike: array of numeric expected",C_result);
	XL_readNumVector(XL_power,C_power," ARM_ERR: power: array of numeric expected",C_result);
	XL_readNumVector(XL_callput,C_callput," ARM_ERR: callput: array of numeric expected",C_result);
	XL_readNumCell(XL_xo,C_xo," ARM_ERR: xo: numeric expected",C_result);
	XL_readNumVector(XL_floor,C_floor," ARM_ERR: floor: array of numeric expected",C_result);
    XL_readNumVector(XL_cap,C_cap," ARM_ERR: cap: array of numeric expected",C_result);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_REVERSECOUPON_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_STRUCTREVERSECOUPON (C_date,
										C_strike,
										C_power,
										C_callput,
										C_xo,
                                        C_floor,
                                        C_cap,
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

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_STRUCTREVERSECOUPON" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_REVERSE_CALENDAR (LPXLOPER XL_structSwapLeg,
															  LPXLOPER XL_classSwapLeg,
															  LPXLOPER XL_receiveOrPay,
															  LPXLOPER XL_coupon,
															  LPXLOPER XL_exe,
															  LPXLOPER XL_redemp,
															  LPXLOPER XL_classRedemp,
															  LPXLOPER XL_dualDate,
															  LPXLOPER XL_dualStrike,
															  LPXLOPER XL_dualFlag,
															  LPXLOPER XL_dStartDates,
															  LPXLOPER XL_dEndDates,
															  LPXLOPER XL_dFixingDates,
															  LPXLOPER XL_dPaymentDates,
															  LPXLOPER XL_fStartDates,
															  LPXLOPER XL_fEndDates,
															  LPXLOPER XL_fFixingDates,
															  LPXLOPER XL_fPaymentDates)
{
	ADD_LOG("Local_REVERSE_CALENDAR ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_structSwapLeg;
	CCString C_classSwapLeg;

	CCString C_receiveOrPay;
	double receiveOrPayId;

	CCString C_coupon;
	CCString C_exe;
	CCString C_redemp;
	CCString C_classRedemp;

    double C_dualDate;
    double C_dualStrike;

    CCString C_dualFlag;

	VECTOR<double> C_dStartDates;
	VECTOR<double> C_dEndDates;
	VECTOR<double> C_dFixingDates;
	VECTOR<double> C_dPaymentDates;
	VECTOR<double> C_fStartDates;
	VECTOR<double> C_fEndDates;
	VECTOR<double> C_fFixingDates;
	VECTOR<double> C_fPaymentDates;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_structSwapLeg,C_structSwapLeg," ARM_ERR: struct swap leg id: object expected",C_result);
	XL_readStrCell(XL_classSwapLeg,C_classSwapLeg," ARM_ERR: class swap leg id: object expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrCell(XL_coupon,C_coupon," ARM_ERR: coupon id: object expected",C_result);
	XL_readStrCell(XL_exe,C_exe," ARM_ERR: exercise id: object expected",C_result);
	XL_readStrCell(XL_redemp,C_redemp," ARM_ERR: redemption id: object expected",C_result);
	XL_readStrCell(XL_classRedemp,C_classRedemp," ARM_ERR: redemption class id: object expected",C_result);
    XL_readNumCell(XL_dualDate,C_dualDate," ARM_ERR: dualDate: numeric expected",C_result);
    XL_readNumCell(XL_dualStrike,C_dualStrike," ARM_ERR: dualStrike: numeric expected",C_result);
    XL_readStrCell(XL_dualFlag,C_dualFlag," ARM_ERR: dualFlag: string expected",C_result);
   
	XL_readNumVector(XL_dStartDates,C_dStartDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_dEndDates,C_dEndDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_dFixingDates,C_dFixingDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_dPaymentDates,C_dPaymentDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_fStartDates,C_fStartDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_fEndDates,C_fEndDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_fFixingDates,C_fFixingDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_fPaymentDates,C_fPaymentDates," ARM_ERR: date: array of numeric expected",C_result);


	if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_REVERSE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_REVERSE_CALENDAR (LocalGetNumObjectId (C_structSwapLeg),
							   LocalGetNumObjectId (C_classSwapLeg),
							   receiveOrPayId,
							   LocalGetNumObjectId (C_coupon),
							   LocalGetNumObjectId (C_exe),
							   LocalGetNumObjectId (C_redemp),
							   LocalGetNumObjectId (C_classRedemp),
                               C_dualDate,
                               C_dualStrike,
                               C_dualFlag,
							   C_dStartDates,
							   C_dEndDates,
							   C_dFixingDates,
							   C_dPaymentDates,
							   C_fStartDates,
							   C_fEndDates,
							   C_fFixingDates,
							   C_fPaymentDates,
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
			retCode = ARMLOCAL_REVERSE_CALENDAR (LocalGetNumObjectId (C_structSwapLeg),
								   LocalGetNumObjectId (C_classSwapLeg),
								   receiveOrPayId,
								   LocalGetNumObjectId (C_coupon),
								   LocalGetNumObjectId (C_exe),
								   LocalGetNumObjectId (C_redemp),
								   LocalGetNumObjectId (C_classRedemp),
                                   C_dualDate,
                                   C_dualStrike,
                                   C_dualFlag,
								   C_dStartDates,
								   C_dEndDates,
								   C_dFixingDates,
								   C_dPaymentDates,
								   C_fStartDates,
								   C_fEndDates,
								   C_fFixingDates,
								   C_fPaymentDates,
								   C_result,
								   objId);

			if(retCode == ARM_OK)
			{			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_REVERSE_CALENDAR (LocalGetNumObjectId (C_structSwapLeg),
								   LocalGetNumObjectId (C_classSwapLeg),
								   receiveOrPayId,
								   LocalGetNumObjectId (C_coupon),
								   LocalGetNumObjectId (C_exe),
								   LocalGetNumObjectId (C_redemp),
								   LocalGetNumObjectId (C_classRedemp),
                                   C_dualDate,
                                   C_dualStrike,
                                   C_dualFlag,
								   C_dStartDates,
								   C_dEndDates,
								   C_dFixingDates,
								   C_dPaymentDates,
								   C_fStartDates,
								   C_fEndDates,
								   C_fFixingDates,
								   C_fPaymentDates,
								   C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_REVERSE_CALENDAR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_REVERSE_CALENDAR (LPXLOPER XL_structSwapLeg,
																  LPXLOPER XL_classSwapLeg,
																  LPXLOPER XL_receiveOrPay,
																  LPXLOPER XL_coupon,
																  LPXLOPER XL_exe,
																  LPXLOPER XL_redemp,
																  LPXLOPER XL_classRedemp,
																  LPXLOPER XL_dualDate,
																  LPXLOPER XL_dualStrike,
																  LPXLOPER XL_dualFlag,
																  LPXLOPER XL_dStartDates,
																  LPXLOPER XL_dEndDates,
																  LPXLOPER XL_dFixingDates,
																  LPXLOPER XL_dPaymentDates,
																  LPXLOPER XL_fStartDates,
																  LPXLOPER XL_fEndDates,
																  LPXLOPER XL_fFixingDates,
																  LPXLOPER XL_fPaymentDates)
{
	ADD_LOG("Local_PXL_REVERSE_CALENDAR ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_structSwapLeg;
	CCString C_classSwapLeg;

	CCString C_receiveOrPay;
	double receiveOrPayId;

	CCString C_coupon;
	CCString C_exe;
	CCString C_redemp;
	CCString C_classRedemp;

    double C_dualDate;
    double C_dualStrike;

    CCString C_dualFlag;

	VECTOR<double> C_dStartDates;
	VECTOR<double> C_dEndDates;
	VECTOR<double> C_dFixingDates;
	VECTOR<double> C_dPaymentDates;
	VECTOR<double> C_fStartDates;
	VECTOR<double> C_fEndDates;
	VECTOR<double> C_fFixingDates;
	VECTOR<double> C_fPaymentDates;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_structSwapLeg,C_structSwapLeg," ARM_ERR: struct swap leg id: object expected",C_result);
	XL_readStrCell(XL_classSwapLeg,C_classSwapLeg," ARM_ERR: class swap leg id: object expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrCell(XL_coupon,C_coupon," ARM_ERR: coupon id: object expected",C_result);
	XL_readStrCell(XL_exe,C_exe," ARM_ERR: exercise id: object expected",C_result);
	XL_readStrCell(XL_redemp,C_redemp," ARM_ERR: redemption id: object expected",C_result);
	XL_readStrCell(XL_classRedemp,C_classRedemp," ARM_ERR: redemption class id: object expected",C_result);
    XL_readNumCell(XL_dualDate,C_dualDate," ARM_ERR: dualDate: numeric expected",C_result);
    XL_readNumCell(XL_dualStrike,C_dualStrike," ARM_ERR: dualStrike: numeric expected",C_result);
    XL_readStrCell(XL_dualFlag,C_dualFlag," ARM_ERR: dualFlag: string expected",C_result);
   
	XL_readNumVector(XL_dStartDates,C_dStartDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_dEndDates,C_dEndDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_dFixingDates,C_dFixingDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_dPaymentDates,C_dPaymentDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_fStartDates,C_fStartDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_fEndDates,C_fEndDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_fFixingDates,C_fFixingDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_fPaymentDates,C_fPaymentDates," ARM_ERR: date: array of numeric expected",C_result);


	if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_REVERSE_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_REVERSE_CALENDAR (LocalGetNumObjectId (C_structSwapLeg),
							   LocalGetNumObjectId (C_classSwapLeg),
							   receiveOrPayId,
							   LocalGetNumObjectId (C_coupon),
							   LocalGetNumObjectId (C_exe),
							   LocalGetNumObjectId (C_redemp),
							   LocalGetNumObjectId (C_classRedemp),
                               C_dualDate,
                               C_dualStrike,
                               C_dualFlag,
							   C_dStartDates,
							   C_dEndDates,
							   C_dFixingDates,
							   C_dPaymentDates,
							   C_fStartDates,
							   C_fEndDates,
							   C_fFixingDates,
							   C_fPaymentDates,
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

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_REVERSE_CALENDAR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_REVERSENOTIONAL_CALENDAR (LPXLOPER XL_structSwapLeg,
																	  LPXLOPER XL_classSwapLeg,
																	  LPXLOPER XL_receiveOrPay,
																	  LPXLOPER XL_coupon,
																	  LPXLOPER XL_exe,
																	  LPXLOPER XL_redemp,
																	  LPXLOPER XL_classRedemp,
																	  LPXLOPER XL_notExchFlag,
																	  LPXLOPER XL_dualDate,
																	  LPXLOPER XL_dualStrike,
																	  LPXLOPER XL_dualFlag,
																	  LPXLOPER XL_dStartDates,
																	  LPXLOPER XL_dEndDates,
																	  LPXLOPER XL_dFixingDates,
																	  LPXLOPER XL_dPaymentDates,
																	  LPXLOPER XL_fStartDates,
																	  LPXLOPER XL_fEndDates,
																	  LPXLOPER XL_fFixingDates,
																	  LPXLOPER XL_fPaymentDates)
{
	ADD_LOG("Local_REVERSENOTIONAL_CALENDAR ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_structSwapLeg;
	CCString C_classSwapLeg;

	CCString C_receiveOrPay;
	double receiveOrPayId;

	CCString C_coupon;
	CCString C_exe;
	CCString C_redemp;
	CCString C_classRedemp;

    double C_dualDate;
    double C_dualStrike;
	double C_notExchFlag;

    CCString C_dualFlag;

	VECTOR<double> C_dStartDates;
	VECTOR<double> C_dEndDates;
	VECTOR<double> C_dFixingDates;
	VECTOR<double> C_dPaymentDates;
	VECTOR<double> C_fStartDates;
	VECTOR<double> C_fEndDates;
	VECTOR<double> C_fFixingDates;
	VECTOR<double> C_fPaymentDates;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_structSwapLeg,C_structSwapLeg," ARM_ERR: struct swap leg id: object expected",C_result);
	XL_readStrCell(XL_classSwapLeg,C_classSwapLeg," ARM_ERR: class swap leg id: object expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrCell(XL_coupon,C_coupon," ARM_ERR: coupon id: object expected",C_result);
	XL_readStrCell(XL_exe,C_exe," ARM_ERR: exercise id: object expected",C_result);
	XL_readStrCell(XL_redemp,C_redemp," ARM_ERR: redemption id: object expected",C_result);
	XL_readStrCell(XL_classRedemp,C_classRedemp," ARM_ERR: redemption class id: object expected",C_result);
	XL_readNumCell(XL_notExchFlag,C_notExchFlag," ARM_ERR: notional exch flag: numeric expected",C_result);
    XL_readNumCell(XL_dualDate,C_dualDate," ARM_ERR: dualDate: numeric expected",C_result);
    XL_readNumCell(XL_dualStrike,C_dualStrike," ARM_ERR: dualStrike: numeric expected",C_result);
    XL_readStrCell(XL_dualFlag,C_dualFlag," ARM_ERR: dualFlag: string expected",C_result);
   
	XL_readNumVector(XL_dStartDates,C_dStartDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_dEndDates,C_dEndDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_dFixingDates,C_dFixingDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_dPaymentDates,C_dPaymentDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_fStartDates,C_fStartDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_fEndDates,C_fEndDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_fFixingDates,C_fFixingDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_fPaymentDates,C_fPaymentDates," ARM_ERR: date: array of numeric expected",C_result);


	if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_REVERSE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_REVERSENOTIONAL_CALENDAR (LocalGetNumObjectId (C_structSwapLeg),
													 LocalGetNumObjectId (C_classSwapLeg),
													 receiveOrPayId,
													 LocalGetNumObjectId (C_coupon),
													 LocalGetNumObjectId (C_exe),
													 LocalGetNumObjectId (C_redemp),
													 LocalGetNumObjectId (C_classRedemp),
													 (long)C_notExchFlag,
													 C_dualDate,
													 C_dualStrike,
													 C_dualFlag,
													 C_dStartDates,
													 C_dEndDates,
													 C_dFixingDates,
													 C_dPaymentDates,
													 C_fStartDates,
													 C_fEndDates,
													 C_fFixingDates,
													 C_fPaymentDates,
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
			retCode = ARMLOCAL_REVERSENOTIONAL_CALENDAR (LocalGetNumObjectId (C_structSwapLeg),
														 LocalGetNumObjectId (C_classSwapLeg),
														 receiveOrPayId,
														 LocalGetNumObjectId (C_coupon),
														 LocalGetNumObjectId (C_exe),
														 LocalGetNumObjectId (C_redemp),
														 LocalGetNumObjectId (C_classRedemp),
														 (long)C_notExchFlag,
														 C_dualDate,
														 C_dualStrike,
														 C_dualFlag,
														 C_dStartDates,
														 C_dEndDates,
														 C_dFixingDates,
														 C_dPaymentDates,
														 C_fStartDates,
														 C_fEndDates,
														 C_fFixingDates,
														 C_fPaymentDates,
														 C_result,
														 objId);

			if(retCode == ARM_OK)
			{			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_REVERSENOTIONAL_CALENDAR (LocalGetNumObjectId (C_structSwapLeg),
														 LocalGetNumObjectId (C_classSwapLeg),
														 receiveOrPayId,
														 LocalGetNumObjectId (C_coupon),
														 LocalGetNumObjectId (C_exe),
														 LocalGetNumObjectId (C_redemp),
														 LocalGetNumObjectId (C_classRedemp),
														 (long)C_notExchFlag,
														 C_dualDate,
														 C_dualStrike,
														 C_dualFlag,
														 C_dStartDates,
														 C_dEndDates,
														 C_dFixingDates,
														 C_dPaymentDates,
														 C_fStartDates,
														 C_fEndDates,
														 C_fFixingDates,
														 C_fPaymentDates,
														 C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_REVERSENOTIONAL_CALENDAR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_REVERSENOTIONAL_CALENDAR (LPXLOPER XL_structSwapLeg,
																		  LPXLOPER XL_classSwapLeg,
																		  LPXLOPER XL_receiveOrPay,
																		  LPXLOPER XL_coupon,
																		  LPXLOPER XL_exe,
																		  LPXLOPER XL_redemp,
																		  LPXLOPER XL_classRedemp,
																		  LPXLOPER XL_notExchFlag,
																		  LPXLOPER XL_dualDate,
																		  LPXLOPER XL_dualStrike,
																		  LPXLOPER XL_dualFlag,
																		  LPXLOPER XL_dStartDates,
																		  LPXLOPER XL_dEndDates,
																		  LPXLOPER XL_dFixingDates,
																		  LPXLOPER XL_dPaymentDates,
																		  LPXLOPER XL_fStartDates,
																		  LPXLOPER XL_fEndDates,
																		  LPXLOPER XL_fFixingDates,
																		  LPXLOPER XL_fPaymentDates)
{
	ADD_LOG("Local_PXL_REVERSENOTIONAL_CALENDAR ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_structSwapLeg;
	CCString C_classSwapLeg;

	CCString C_receiveOrPay;
	double receiveOrPayId;

	CCString C_coupon;
	CCString C_exe;
	CCString C_redemp;
	CCString C_classRedemp;

    double C_dualDate;
    double C_dualStrike;
	double C_notExchFlag;

    CCString C_dualFlag;

	VECTOR<double> C_dStartDates;
	VECTOR<double> C_dEndDates;
	VECTOR<double> C_dFixingDates;
	VECTOR<double> C_dPaymentDates;
	VECTOR<double> C_fStartDates;
	VECTOR<double> C_fEndDates;
	VECTOR<double> C_fFixingDates;
	VECTOR<double> C_fPaymentDates;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_structSwapLeg,C_structSwapLeg," ARM_ERR: struct swap leg id: object expected",C_result);
	XL_readStrCell(XL_classSwapLeg,C_classSwapLeg," ARM_ERR: class swap leg id: object expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrCell(XL_coupon,C_coupon," ARM_ERR: coupon id: object expected",C_result);
	XL_readStrCell(XL_exe,C_exe," ARM_ERR: exercise id: object expected",C_result);
	XL_readStrCell(XL_redemp,C_redemp," ARM_ERR: redemption id: object expected",C_result);
	XL_readStrCell(XL_classRedemp,C_classRedemp," ARM_ERR: redemption class id: object expected",C_result);
	XL_readNumCell(XL_notExchFlag,C_notExchFlag," ARM_ERR: notional exch flag: numeric expected",C_result);
    XL_readNumCell(XL_dualDate,C_dualDate," ARM_ERR: dualDate: numeric expected",C_result);
    XL_readNumCell(XL_dualStrike,C_dualStrike," ARM_ERR: dualStrike: numeric expected",C_result);
    XL_readStrCell(XL_dualFlag,C_dualFlag," ARM_ERR: dualFlag: string expected",C_result);
   
	XL_readNumVector(XL_dStartDates,C_dStartDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_dEndDates,C_dEndDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_dFixingDates,C_dFixingDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_dPaymentDates,C_dPaymentDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_fStartDates,C_fStartDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_fEndDates,C_fEndDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_fFixingDates,C_fFixingDates," ARM_ERR: date: array of numeric expected",C_result);
	XL_readNumVector(XL_fPaymentDates,C_fPaymentDates," ARM_ERR: date: array of numeric expected",C_result);


	if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_REVERSE_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_REVERSENOTIONAL_CALENDAR (LocalGetNumObjectId (C_structSwapLeg),
												 LocalGetNumObjectId (C_classSwapLeg),
												 receiveOrPayId,
												 LocalGetNumObjectId (C_coupon),
												 LocalGetNumObjectId (C_exe),
												 LocalGetNumObjectId (C_redemp),
												 LocalGetNumObjectId (C_classRedemp),
												 (long)C_notExchFlag,
												 C_dualDate,
												 C_dualStrike,
												 C_dualFlag,
												 C_dStartDates,
												 C_dEndDates,
												 C_dFixingDates,
												 C_dPaymentDates,
												 C_fStartDates,
												 C_fEndDates,
												 C_fFixingDates,
												 C_fPaymentDates,
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

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_REVERSENOTIONAL_CALENDAR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_POWERREVERSE (LPXLOPER XL_initPeriodLeg,
															  LPXLOPER XL_realFundLeg,
															  LPXLOPER XL_fxUnderLeg,
															  LPXLOPER XL_fxNumLeg,
															  LPXLOPER XL_noticeDates,
															  LPXLOPER XL_cancelDates,
															  LPXLOPER XL_FX0,
															  LPXLOPER XL_fxStep,
															  LPXLOPER XL_capValue,
															  LPXLOPER XL_capStepValue,
															  LPXLOPER XL_floorValue,
															  LPXLOPER XL_floorStepValue,
															  LPXLOPER XL_dualOptionFlag,
															  LPXLOPER XL_dualOptionStrike,
                                                              LPXLOPER XL_redempNoticeDate,
															  LPXLOPER XL_lastLiborFixing,
															  LPXLOPER XL_lastFxSpotFixing)
{
	ADD_LOG("Local_ARM_POWERREVERSE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_initPeriodLeg;
	long initPeriodLegId;
	CCString C_realFundLeg;
	CCString C_fxUnderLeg;
	CCString C_fxNumLeg;

	VECTOR<double> C_dNoticeDates;
	VECTOR<double> C_dCancelDates;
	VECTOR<double> C_dCancelDates_default;

	double   C_FX0_double;
	CCString C_FX0_str;
	long     FX0Type;
    double C_fxStep;

	double C_capValue_double;
	CCString C_capValue_str;
	long capValueType;
	double C_capStepValue;

	double C_floorValue_double;
	CCString C_floorValue_str;
	long floorValueType;
	double C_floorStepValue;

	double C_dualOptionFlag;
	double C_dualOptionStrike;

    double C_redempNoticeDate;
    double C_redempNoticeDateDef = -1.;

    double C_lastLiborFixing;
    double C_lastLiborFixingDef = -1.;

    double C_lastFxSpotFixing;
    double C_lastFxSpotFixingDef = -1.;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCellWD(XL_initPeriodLeg,C_initPeriodLeg,"DEFAULT"," ARM_ERR: init period leg id: object expected",C_result);
	XL_readStrCell(XL_realFundLeg,C_realFundLeg," ARM_ERR: real Fund leg id: object expected",C_result);
	XL_readStrCell(XL_fxUnderLeg,C_fxUnderLeg," ARM_ERR: fx under leg id: object expected",C_result);
	XL_readStrCell(XL_fxNumLeg,C_fxNumLeg," ARM_ERR:  fx Num leg id: object expected",C_result);
    XL_readNumVector(XL_noticeDates,C_dNoticeDates," ARM_ERR: notice Dates: array of num expected",C_result);
    XL_readNumVectorWD(XL_cancelDates,C_dCancelDates,C_dCancelDates_default," ARM_ERR: cancel Dates: array of num expected",C_result);

	XL_readNumCell(XL_dualOptionFlag,C_dualOptionFlag," ARM_ERR: dual option flag: numeric expected",C_result);
    XL_readNumCell(XL_dualOptionStrike,C_dualOptionStrike," ARM_ERR: dual option strike: numeric expected",C_result);
    XL_readNumCellWD(XL_redempNoticeDate,C_redempNoticeDate,C_redempNoticeDateDef," ARM_ERR: notice Date: a date expected",C_result);
    XL_readNumCellWD(XL_lastLiborFixing,C_lastLiborFixing,C_lastLiborFixingDef," ARM_ERR: last Libor Fixing: numeric expected",C_result);
    XL_readNumCellWD(XL_lastFxSpotFixing,C_lastFxSpotFixing,C_lastFxSpotFixingDef," ARM_ERR: last Fx Spot Fixing: numeric expected",C_result);

	XL_readStrOrNumCell(XL_FX0, C_FX0_str, C_FX0_double, FX0Type,
		   " ARM_ERR: FX0: numeric or object ID string expected",C_result);
    
	if ( FX0Type == XL_TYPE_STRING )
	{
		C_FX0_double = (double) LocalGetNumObjectId(C_FX0_str);
		C_fxStep = 0.0;

		FX0Type = 1L;
	}
	else
	{
		XL_readNumCell(XL_fxStep,C_fxStep," ARM_ERR: Fx step: numeric expected",C_result);
		
		FX0Type = 0L;
	}

	XL_readStrOrNumCell(XL_capValue,C_capValue_str,C_capValue_double, capValueType,
			" ARM_ERR: cap Value: numeric or object expected",C_result);

	if ( capValueType == XL_TYPE_STRING )
	{
		C_capValue_double = (double) LocalGetNumObjectId(C_capValue_str);
		C_capStepValue = 0.0;

		capValueType = 1L;
	}
	else
	{
		XL_readNumCell(XL_capStepValue,C_capStepValue," ARM_ERR: Cap Step Value: numeric expected",C_result);
		
		capValueType = 0L;
	}

	XL_readStrOrNumCell(XL_floorValue,C_floorValue_str,C_floorValue_double, floorValueType,
			" ARM_ERR: floor Value: numeric or object expected",C_result);

	if ( floorValueType == XL_TYPE_STRING )
	{
		C_floorValue_double = (double) LocalGetNumObjectId(C_floorValue_str);
		C_floorStepValue = 0.0;

		floorValueType = 1L;
	}
	else
	{
		XL_readNumCell(XL_floorStepValue,C_floorStepValue," ARM_ERR: floor Step Value: numeric expected",C_result);
		
		floorValueType = 0L;
	}

	if(C_initPeriodLeg == "DEFAULT")
	{
		initPeriodLegId = ARM_NULL_OBJECT;
	}
	else
	{
		initPeriodLegId = LocalGetNumObjectId (C_initPeriodLeg);
	}

	if (C_dCancelDates.size() == 0)
		C_dCancelDates = C_dNoticeDates;

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_POWER_REVERSE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_POWERREVERSE (initPeriodLegId,
										 LocalGetNumObjectId (C_realFundLeg),
										 LocalGetNumObjectId (C_fxUnderLeg),
 										 LocalGetNumObjectId (C_fxNumLeg),
										 C_dNoticeDates,
										 C_dCancelDates,
										 FX0Type,
										 C_FX0_double,
										 C_fxStep,
										 capValueType,
										 C_capValue_double,
										 C_capStepValue,
										 floorValueType,
										 C_floorValue_double,
										 C_floorStepValue,
										 (long)C_dualOptionFlag,
										 C_dualOptionStrike,
                                         C_redempNoticeDate,
										 C_lastLiborFixing,
										 C_lastFxSpotFixing,
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
			retCode = ARMLOCAL_POWERREVERSE (initPeriodLegId,
											 LocalGetNumObjectId (C_realFundLeg),
											 LocalGetNumObjectId (C_fxUnderLeg),
 											 LocalGetNumObjectId (C_fxNumLeg),
											 C_dNoticeDates,
											 C_dCancelDates,
											 FX0Type,
											 C_FX0_double,
											 C_fxStep,
											 capValueType,
											 C_capValue_double,
											 C_capStepValue,
											 floorValueType,
											 C_floorValue_double,
											 C_floorStepValue,
											 (long)C_dualOptionFlag,
											 C_dualOptionStrike,
											 C_redempNoticeDate,
											 C_lastLiborFixing,
											 C_lastFxSpotFixing,
											 C_result,
											 objId);

			if(retCode == ARM_OK)
			{			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

			retCode = ARMLOCAL_POWERREVERSE (initPeriodLegId,
											 LocalGetNumObjectId (C_realFundLeg),
											 LocalGetNumObjectId (C_fxUnderLeg),
 											 LocalGetNumObjectId (C_fxNumLeg),
											 C_dNoticeDates,
											 C_dCancelDates,
											 FX0Type,
											 C_FX0_double,
											 C_fxStep,
											 capValueType,
											 C_capValue_double,
											 C_capStepValue,
											 floorValueType,
											 C_floorValue_double,
											 C_floorStepValue,
											 (long)C_dualOptionFlag,
											 C_dualOptionStrike,
                                             C_redempNoticeDate,
											 C_lastLiborFixing,
											 C_lastFxSpotFixing,
											 C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_POWERREVERSE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_POWERREVERSE (LPXLOPER XL_initPeriodLeg,
																  LPXLOPER XL_realFundLeg,
																  LPXLOPER XL_fxUnderLeg,
																  LPXLOPER XL_fxNumLeg,
																  LPXLOPER XL_noticeDates,
																  LPXLOPER XL_cancelDates,
																  LPXLOPER XL_FX0,
																  LPXLOPER XL_fxStep,
																  LPXLOPER XL_capValue,
																  LPXLOPER XL_capStepValue,
																  LPXLOPER XL_floorValue,
																  LPXLOPER XL_floorStepValue,
																  LPXLOPER XL_dualOptionFlag,
																  LPXLOPER XL_dualOptionStrike,
                                                                  LPXLOPER XL_redempNoticeDate,
																  LPXLOPER XL_lastLiborFixing,
																  LPXLOPER XL_lastFxSpotFixing)
{
	ADD_LOG("Local_PXL_ARM_POWERREVERSE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_initPeriodLeg;
	long initPeriodLegId;
	CCString C_realFundLeg;
	CCString C_fxUnderLeg;
	CCString C_fxNumLeg;

	VECTOR<double> C_dNoticeDates;
	VECTOR<double> C_dCancelDates;
	VECTOR<double> C_dCancelDates_default;

	double   C_FX0_double;
	CCString C_FX0_str;
	long     FX0Type;
    double C_fxStep;

	double C_capValue_double;
	CCString C_capValue_str;
	long capValueType;
	double C_capStepValue;

	double C_floorValue_double;
	CCString C_floorValue_str;
	long floorValueType;
	double C_floorStepValue;

	double C_dualOptionFlag;
	double C_dualOptionStrike;

    double C_redempNoticeDate;
    double C_redempNoticeDateDef = -1;

    double C_lastLiborFixing;
    double C_lastLiborFixingDef = -1.;

    double C_lastFxSpotFixing;
    double C_lastFxSpotFixingDef = -1.;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCellWD(XL_initPeriodLeg,C_initPeriodLeg,"DEFAULT"," ARM_ERR: init period leg id: object expected",C_result);
	XL_readStrCell(XL_realFundLeg,C_realFundLeg," ARM_ERR: real Fund leg id: object expected",C_result);
	XL_readStrCell(XL_fxUnderLeg,C_fxUnderLeg," ARM_ERR: fx under leg id: object expected",C_result);
	XL_readStrCell(XL_fxNumLeg,C_fxNumLeg," ARM_ERR:  fx Num leg id: object expected",C_result);
    XL_readNumVector(XL_noticeDates,C_dNoticeDates," ARM_ERR: notice Dates: array of num expected",C_result);
    XL_readNumVectorWD(XL_cancelDates,C_dCancelDates,C_dCancelDates_default," ARM_ERR: cancel Dates: array of num expected",C_result);

	XL_readNumCell(XL_dualOptionFlag,C_dualOptionFlag," ARM_ERR: dual option flag: numeric expected",C_result);
    XL_readNumCell(XL_dualOptionStrike,C_dualOptionStrike," ARM_ERR: dual option strike: numeric expected",C_result);
    XL_readNumCellWD(XL_redempNoticeDate,C_redempNoticeDate,C_redempNoticeDateDef," ARM_ERR: notice Date: a date expected",C_result);
    XL_readNumCellWD(XL_lastLiborFixing,C_lastLiborFixing,C_lastLiborFixingDef," ARM_ERR: last Libor Fixing: numeric expected",C_result);
    XL_readNumCellWD(XL_lastFxSpotFixing,C_lastFxSpotFixing,C_lastFxSpotFixingDef," ARM_ERR: last Fx Spot Fixing: numeric expected",C_result);

	XL_readStrOrNumCell(XL_FX0, C_FX0_str, C_FX0_double, FX0Type,
		   " ARM_ERR: FX0: numeric or object ID string expected",C_result);
    
	if ( FX0Type == XL_TYPE_STRING )
	{
		C_FX0_double = (double) LocalGetNumObjectId(C_FX0_str);
		C_fxStep = 0.0;

		FX0Type = 1L;
	}
	else
	{
		XL_readNumCell(XL_fxStep,C_fxStep," ARM_ERR: Fx step: numeric expected",C_result);
		
		FX0Type = 0L;
	}

	XL_readStrOrNumCell(XL_capValue,C_capValue_str,C_capValue_double, capValueType,
			" ARM_ERR: cap Value: numeric or object expected",C_result);

	if ( capValueType == XL_TYPE_STRING )
	{
		C_capValue_double = (double) LocalGetNumObjectId(C_capValue_str);
		C_capStepValue = 0.0;

		capValueType = 1L;
	}
	else
	{
		XL_readNumCell(XL_capStepValue,C_capStepValue," ARM_ERR: Cap Step Value: numeric expected",C_result);
		
		capValueType = 0L;
	}

	XL_readStrOrNumCell(XL_floorValue,C_floorValue_str,C_floorValue_double, floorValueType,
			" ARM_ERR: floor Value: numeric or object expected",C_result);

	if ( floorValueType == XL_TYPE_STRING )
	{
		C_floorValue_double = (double) LocalGetNumObjectId(C_floorValue_str);
		C_floorStepValue = 0.0;

		floorValueType = 1L;
	}
	else
	{
		XL_readNumCell(XL_floorStepValue,C_floorStepValue," ARM_ERR: floor Step Value: numeric expected",C_result);
		
		floorValueType = 0L;
	}
	
	if(C_initPeriodLeg == "DEFAULT")
	{
		initPeriodLegId = ARM_NULL_OBJECT;
	}
	else
	{
		initPeriodLegId = LocalGetNumObjectId (C_initPeriodLeg);
	}

	if (C_dCancelDates.size() == 0)
		C_dCancelDates = C_dNoticeDates;

    long retCode;
	long objId;
	
	CCString curClass = LOCAL_POWER_REVERSE_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_POWERREVERSE (initPeriodLegId,
									 LocalGetNumObjectId (C_realFundLeg),
									 LocalGetNumObjectId (C_fxUnderLeg),
 									 LocalGetNumObjectId (C_fxNumLeg),
									 C_dNoticeDates,
									 C_dCancelDates,
									 FX0Type,
									 C_FX0_double,
									 C_fxStep,
									 capValueType,
									 C_capValue_double,
									 C_capStepValue,
									 floorValueType,
									 C_floorValue_double,
									 C_floorStepValue,
									 (long)C_dualOptionFlag,
									 C_dualOptionStrike,
                                     C_redempNoticeDate,
									 C_lastLiborFixing,
									 C_lastFxSpotFixing,		 
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

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_POWERREVERSE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////
/// Local_DataFromPowerReverse function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_DataFromPowerReverse(
    LPXLOPER XL_prcsId,
    LPXLOPER XL_DataType )
{
	ADD_LOG("Local_DataFromPowerReverse");
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
	long C_prcsId;
	XL_GETOBJID( XL_prcsId,	C_prcsId,	" ARM_ERR: Power Reverse: Object expected",	C_result);

	long C_dataType;
	XL_GETMETHOD(XL_DataType, C_dataType, " ARM_ERR: data type: string expected", C_result, ARM_ConvDatePowerReverseDataType );

	VECTOR<double> C_DataResult;
	long retCode;
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
		retCode = ARMLOCAL_DatePowerReverseGetData( C_prcsId, C_dataType, C_DataResult, C_result );
	else
		retCode = ARM_KO;
		
	/// feed the LPXLOPER object result 
	if (retCode == ARM_OK)
	{
		/// add these additional lines 
		/// to display blank lines
		const int additionalLinesNb = 100;
		bool fillWithBlank = true;
		FreeCurCellContent ();
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



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GETPRCSDATA (LPXLOPER XL_prcs)
{
	ADD_LOG("Local_ARM_GETPRCSDATA ");
//	ARM_BEGIN();
	
	// return
	static XLOPER XL_result;
	LPXLOPER pxArray;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_prcs;

		VECTOR<CCString> listeString;
		VECTOR<double>* listeDouble;
		int nbColumns;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_prcs,C_prcs," ARM_ERR: prcs: object expected",C_result);

		long retCode;
    
		retCode = ARMLOCAL_GETPRCSDATA (LocalGetNumObjectId(C_prcs),C_result);

		if(retCode == ARM_OK)
		{
			long result;

			listeDouble = new std::vector<double> [12];
			result = LocalExtractInfoFromFilePRCS ((CCString)"123" + (CCString)(CC_NS(ARM,ARM_USERNAME).c_str()), listeString, listeDouble, nbColumns);

			if(result == ARM_OK)
			{
				int nbrows = listeString.size ()-1;
			
				FreeCurCellErr ();
				XL_result.xltype = xltypeMulti;
				XL_result.val.array.columns = nbColumns;
				XL_result.val.array.rows = nbrows; 
				XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbColumns * sizeof (XLOPER));

				for(int i = 0; i < nbrows; i++)
				{
					pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype = xltypeStr;
					pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].val.str = XL_StrC2StrPascal (listeString[i]);
					pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype |= xlbitDLLFree;

					for (int j=1;j<nbColumns;j++)
					{
						pxArray[XL_Coordonnate2Rank (i, j, nbColumns)].xltype = xltypeNum;
						pxArray[XL_Coordonnate2Rank (i, j, nbColumns)].val.num = (listeDouble[j-1])[i];
					}
				}
			}

			delete [] listeDouble;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GETPRCSDATA" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_DELETENEXTCALLFROMPRCS (LPXLOPER XL_prcs,
																		LPXLOPER XL_asOf)
{
	ADD_LOG("Local_ARM_DELETENEXTCALLFROMPRCS ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_prcs;
		double C_asOf;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_prcs,C_prcs," ARM_ERR: prcs id: object expected",C_result);
		XL_readNumCell(XL_asOf,C_asOf," ARM_ERR: asOfDate: date expected",C_result);

		long retCode;
		long objId;
		CCString prevClass;

		CCString curClass = LOCAL_POWER_REVERSE_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if(!stringId)
		{
			retCode = ARMLOCAL_DELETENEXTCALLFROMPRCS (LocalGetNumObjectId (C_prcs),
													   C_asOf,
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
				retCode = ARMLOCAL_DELETENEXTCALLFROMPRCS (LocalGetNumObjectId (C_prcs),
														   C_asOf,
														   C_result,
														   objId);

				if(retCode == ARM_OK)
				{			
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent ();

				retCode = ARMLOCAL_DELETENEXTCALLFROMPRCS (LocalGetNumObjectId (C_prcs),
														   C_asOf,
														   C_result);
			
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_DELETENEXTCALLFROMPRCS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//----------------------------------------------------
// PRCS access function
// if dateType = "C", get vector of Cancellation dates
// if dateType = "N", get vector of Notice dates
//----------------------------------------------------
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GETOPTIONDATES (LPXLOPER XL_prcs,
																LPXLOPER XL_dateType)
{
	ADD_LOG("Local_ARM_GETOPTIONDATES ");
//	ARM_BEGIN();
	
	// return
	static XLOPER XL_result;
	LPXLOPER pxArray;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_prcs;
		CCString C_dateType;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_prcs,C_prcs," ARM_ERR: prcs: object expected",C_result);
		XL_readStrCell(XL_dateType,C_dateType," ARM_ERR: prcs: string expected",C_result);

		long retCode;
    
		retCode = ARMLOCAL_GETOPTIONDATES (LocalGetNumObjectId(C_prcs),C_dateType,C_result);

		if(retCode == ARM_OK)
		{
			int nbColumns = 1;
			int nbRows = (C_result.getStringVector()).size();
			FreeCurCellErr ();
			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = nbColumns;
			XL_result.val.array.rows = nbRows; 
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbRows * nbColumns * sizeof (XLOPER));

			for(int i = 0; i < nbRows; i++)
			{
				pxArray[XL_Coordonnate2Rank (0, i, nbRows)].xltype = xltypeNum;
				pxArray[XL_Coordonnate2Rank (0, i, nbRows)].val.num = Local_ARMDATE2XLDATE ((C_result.getStringVector())[i]);
			}
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GETOPTIONDATES" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GETDUALOPTIONSTRIKE (LPXLOPER XL_prcs)
{
	ADD_LOG("Local_ARM_GETDUALOPTIONSTRIKE ");
//	ARM_BEGIN();
	
	// return
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_prcs;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_prcs,C_prcs," ARM_ERR: prcs: object expected",C_result);

		long retCode;
    
		retCode = ARMLOCAL_GETDUALOPTIONSTRIKE (LocalGetNumObjectId(C_prcs),C_result);

		if(retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GETDUALOPTIONSTRIKE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}