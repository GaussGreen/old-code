#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>

#include "ARM_xl_bond_local.h"

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_bond.h>

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_BOND (LPXLOPER XL_issueDate,
												  LPXLOPER XL_maturityDate,
												  LPXLOPER XL_firstCouponDate,
												  LPXLOPER XL_couponRate,
												  LPXLOPER XL_redemptionPrice,
												  LPXLOPER XL_periodicity,
												  LPXLOPER XL_dayCount,
												  LPXLOPER XL_settleGap,
												  LPXLOPER XL_couponDateFlag,
												  LPXLOPER XL_ccyId)
{
	ADD_LOG("Local_BOND ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_issueDate;
	double C_maturityDate;
	double C_firstCouponDate;
	double C_couponRate;
	double C_redemptionPrice;
	double C_periodicity;

	CCString C_dayCount_str;
	double C_dayCount_double;
	long dayCount_type;

	double C_settleGap;
	double C_settleGap_default = 3;

	double C_couponDateFlag;
	double C_couponDateFlag_default = 3;

	CCString C_ccyId;
	long l_ccyId; 

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_issueDate,C_issueDate," ARM_ERR: issue date: date expected",C_result);
	XL_readNumCell(XL_maturityDate,C_maturityDate," ARM_ERR: maturity date: date expected",C_result);
	XL_readNumCell(XL_firstCouponDate,C_firstCouponDate," ARM_ERR: first coupon date: numeric expected",C_result);
	XL_readNumCell(XL_couponRate,C_couponRate," ARM_ERR: coupon rate: numeric expected",C_result);
	XL_readNumCell(XL_redemptionPrice,C_redemptionPrice," ARM_ERR: redemption price: numeric expected",C_result);
	XL_readNumCell(XL_periodicity,C_periodicity," ARM_ERR: periodicity: numeric expected",C_result);
	XL_readStrOrNumCell(XL_dayCount,C_dayCount_str,C_dayCount_double,dayCount_type," ARM_ERR: day count: string or numeric expected",C_result);
	XL_readNumCellWD(XL_settleGap,C_settleGap,C_settleGap_default," ARM_ERR: settle gap: numeric expected",C_result);
	XL_readNumCellWD(XL_couponDateFlag,C_couponDateFlag,C_couponDateFlag_default," ARM_ERR: coupon date flag: numeric expected",C_result);
	XL_readStrCellWD(XL_ccyId, C_ccyId, "DEFAULT"," ARM_ERR: ccyId: string expected",C_result);

	l_ccyId = LocalGetNumObjectId (C_ccyId);

	if(dayCount_type == XL_TYPE_STRING)
	{
		C_dayCount_double = (double)ARM_ConvDayCount (C_dayCount_str);
	}
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_BOND_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_bond (C_issueDate, C_maturityDate, 
			                C_firstCouponDate, C_couponRate,
							C_redemptionPrice, (long)C_periodicity,
							(long)C_dayCount_double, (long)C_settleGap,
							(long)C_couponDateFlag, l_ccyId , C_result);

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
			if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
			{
				retCode = ARMLOCAL_bond (C_issueDate, C_maturityDate, 
									C_firstCouponDate, C_couponRate,
									C_redemptionPrice, (long)C_periodicity,
									(long)C_dayCount_double, (long)C_settleGap,
									(long)C_couponDateFlag, l_ccyId, C_result, objId);

				LocalSetCurCellEnvValue (curClass, objId); 
				
				stringId = LocalMakeObjectId (objId, curClass);
			}
			else
				retCode = ARM_KO;			
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_bond (C_issueDate, C_maturityDate, 
			                C_firstCouponDate, C_couponRate,
							C_redemptionPrice, (long)C_periodicity,
							(long)C_dayCount_double, (long)C_settleGap,
							(long)C_couponDateFlag, l_ccyId ,C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BOND" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_RISKYBOND (LPXLOPER XL_issueDate,
													   LPXLOPER XL_maturityDate,
												       LPXLOPER XL_firstCouponDate,
												       LPXLOPER XL_couponRate,
												       LPXLOPER XL_redemptionPrice,
												       LPXLOPER XL_periodicity,
												       LPXLOPER XL_dayCount,
												       LPXLOPER XL_settleGap,
												       LPXLOPER XL_couponDateFlag,
												       LPXLOPER XL_ccyId,
													   LPXLOPER XL_sRepo,
												       LPXLOPER XL_ssl,
												       LPXLOPER XL_recoveryRate
													   )
{
	ADD_LOG("Local_RISKYBOND ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_issueDate;
	double C_maturityDate;
	double C_firstCouponDate;
	double C_couponRate;
	double C_redemptionPrice;
	double C_periodicity;
	double C_sRepo;
	double C_sRepo_default = -1.;
	double C_ssl;
	double C_ssl_default = -1.;
	double C_recoveryRate;
	double C_recoveryRate_default = -1.;

	CCString C_dayCount_str;
	double C_dayCount_double;
	long dayCount_type;

	double C_settleGap;
	double C_settleGap_default = 3;

	double C_couponDateFlag;
	double C_couponDateFlag_default = 3;

	CCString C_ccyId;
	long l_ccyId; 

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_issueDate,C_issueDate," ARM_ERR: issue date: date expected",C_result);
	XL_readNumCell(XL_maturityDate,C_maturityDate," ARM_ERR: maturity date: date expected",C_result);
	XL_readNumCell(XL_firstCouponDate,C_firstCouponDate," ARM_ERR: first coupon date: numeric expected",C_result);
	XL_readNumCell(XL_couponRate,C_couponRate," ARM_ERR: coupon rate: numeric expected",C_result);
	XL_readNumCell(XL_redemptionPrice,C_redemptionPrice," ARM_ERR: redemption price: numeric expected",C_result);
	XL_readNumCell(XL_periodicity,C_periodicity," ARM_ERR: periodicity: numeric expected",C_result);
	XL_readStrOrNumCell(XL_dayCount,C_dayCount_str,C_dayCount_double,dayCount_type," ARM_ERR: day count: string or numeric expected",C_result);
	XL_readNumCellWD(XL_settleGap,C_settleGap,C_settleGap_default," ARM_ERR: settle gap: numeric expected",C_result);
	XL_readNumCellWD(XL_couponDateFlag,C_couponDateFlag,C_couponDateFlag_default," ARM_ERR: coupon date flag: numeric expected",C_result);
	XL_readStrCellWD(XL_ccyId, C_ccyId, "DEFAULT"," ARM_ERR: ccyId: string expected",C_result);
	XL_readNumCellWD(XL_sRepo,C_sRepo,C_sRepo_default," ARM_ERR: spread repo: numeric expected",C_result);
	XL_readNumCellWD(XL_ssl,C_ssl,C_ssl_default," ARM_ERR: ssl: numeric expected",C_result);
	XL_readNumCellWD(XL_recoveryRate,C_recoveryRate,C_recoveryRate_default," ARM_ERR: recovery rate: numeric expected",C_result);



	l_ccyId = LocalGetNumObjectId (C_ccyId);

	if(dayCount_type == XL_TYPE_STRING)
	{
		C_dayCount_double = (double)ARM_ConvDayCount (C_dayCount_str);
	}
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_BOND_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_RiskyBond (C_issueDate, C_maturityDate, 
									C_firstCouponDate, C_couponRate,
									C_redemptionPrice, (long)C_periodicity,
									(long)C_dayCount_double, (long)C_settleGap,
									(long)C_couponDateFlag, l_ccyId , C_sRepo,
									C_ssl, C_recoveryRate, C_result);

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
			if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
			{
				retCode = ARMLOCAL_RiskyBond (C_issueDate, C_maturityDate, 
										C_firstCouponDate, C_couponRate,
										C_redemptionPrice, (long)C_periodicity,
										(long)C_dayCount_double, (long)C_settleGap,
										(long)C_couponDateFlag, l_ccyId, 
										C_sRepo, C_ssl, C_recoveryRate, C_result, objId);

				LocalSetCurCellEnvValue (curClass, objId); 
				
				stringId = LocalMakeObjectId (objId, curClass);
			}
			else
				retCode = ARM_KO;			
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_RiskyBond (C_issueDate, C_maturityDate, 
								C_firstCouponDate, C_couponRate,
								C_redemptionPrice, (long)C_periodicity,
								(long)C_dayCount_double, (long)C_settleGap,
								(long)C_couponDateFlag, l_ccyId , C_sRepo, 
								C_ssl, C_recoveryRate, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_RISKYBOND" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_RISKYBONDWITHCF (LPXLOPER XL_asOfDate,
															 LPXLOPER XL_redemptionPrice,
															 LPXLOPER XL_periodicity,
															 LPXLOPER XL_dayCount,
															 LPXLOPER XL_settleGap,
															 LPXLOPER XL_couponDateFlag,
															 LPXLOPER XL_ccyId,
															 LPXLOPER XL_yearTerms,
															 LPXLOPER XL_cashFlows,
															 LPXLOPER XL_sRepo,
															 LPXLOPER XL_ssl,
															 LPXLOPER XL_recoveryRate
															 )
{
	ADD_LOG("Local_RISKYBONDWITHCF ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_asOfDate;
	double C_redemptionPrice;
	double C_periodicity;
	double C_sRepo;
	double C_sRepo_default = -1.;
	double C_ssl;
	double C_ssl_default = -1.;
	double C_recoveryRate;
	double C_recoveryRate_default = -1.;

	CCString C_dayCount_str;
	double C_dayCount_double;
	long dayCount_type;

	double C_settleGap;
	double C_settleGap_default = 3;

	double C_couponDateFlag;
	double C_couponDateFlag_default = 3;

	CCString C_ccyId;
	long l_ccyId; 

	VECTOR<double> C_yearTerms;
	VECTOR<double> C_cashFlows;


	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_asOfDate,C_asOfDate," ARM_ERR: issue date: date expected",C_result);
	XL_readNumCell(XL_redemptionPrice,C_redemptionPrice," ARM_ERR: redemption price: numeric expected",C_result);
	XL_readNumCell(XL_periodicity,C_periodicity," ARM_ERR: periodicity: numeric expected",C_result);
	XL_readStrOrNumCell(XL_dayCount,C_dayCount_str,C_dayCount_double,dayCount_type," ARM_ERR: day count: string or numeric expected",C_result);
	XL_readNumCellWD(XL_settleGap,C_settleGap,C_settleGap_default," ARM_ERR: settle gap: numeric expected",C_result);
	XL_readNumCellWD(XL_couponDateFlag,C_couponDateFlag,C_couponDateFlag_default," ARM_ERR: coupon date flag: numeric expected",C_result);
	XL_readStrCellWD(XL_ccyId, C_ccyId, "DEFAULT"," ARM_ERR: ccyId: string expected",C_result);
	XL_readNumVector(XL_yearTerms, C_yearTerms, " ARM_ERR: Year Terms : numeric vector expected",C_result);
	XL_readNumVector(XL_cashFlows, C_cashFlows, " ARM_ERR: Cash Flows: numeric vector expected",C_result);
	XL_readNumCellWD(XL_sRepo,C_sRepo,C_sRepo_default," ARM_ERR: spread repo: numeric expected",C_result);
	XL_readNumCellWD(XL_ssl,C_ssl,C_ssl_default," ARM_ERR: ssl: numeric expected",C_result);
	XL_readNumCellWD(XL_recoveryRate,C_recoveryRate,C_recoveryRate_default," ARM_ERR: recovery rate: numeric expected",C_result);



	l_ccyId = LocalGetNumObjectId (C_ccyId);

	if(dayCount_type == XL_TYPE_STRING)
	{
		C_dayCount_double = (double)ARM_ConvDayCount (C_dayCount_str);
	}
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_BOND_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_RiskyBondWithCF (C_asOfDate, C_redemptionPrice, (long)C_periodicity,
											(long)C_dayCount_double, (long)C_settleGap,
											(long)C_couponDateFlag, l_ccyId ,
											C_yearTerms, C_cashFlows, C_sRepo,

											C_ssl, C_recoveryRate, C_result);

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
			if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
			{
				retCode = ARMLOCAL_RiskyBondWithCF ( C_asOfDate, C_redemptionPrice, (long)C_periodicity,
													(long)C_dayCount_double, (long)C_settleGap,
													(long)C_couponDateFlag, l_ccyId, 
													C_yearTerms, C_cashFlows, C_sRepo,
													C_ssl, C_recoveryRate, C_result, objId);

				LocalSetCurCellEnvValue (curClass, objId); 
				
				stringId = LocalMakeObjectId (objId, curClass);
			}
			else
				retCode = ARM_KO;			
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_RiskyBondWithCF (C_asOfDate, C_redemptionPrice, (long)C_periodicity,
												(long)C_dayCount_double, (long)C_settleGap,
												(long)C_couponDateFlag, l_ccyId, 
												C_yearTerms, C_cashFlows, C_sRepo,
												C_ssl, C_recoveryRate, C_result, objId);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_RISKYBOND" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_BONDTEC (LPXLOPER XL_issueDate,
													 LPXLOPER XL_maturityDate,
													 LPXLOPER XL_firstCouponDate,
													 LPXLOPER XL_couponRate,
													 LPXLOPER XL_redemptionPrice,
													 LPXLOPER XL_periodicity,
													 LPXLOPER XL_dayCount,
													 LPXLOPER XL_settleGap,
													 LPXLOPER XL_couponDateFlag,
													 LPXLOPER XL_ccyId,
													 LPXLOPER XL_tec,
													 LPXLOPER XL_pfId,
													 LPXLOPER XL_modId)
{
	ADD_LOG("Local_BONDTEC ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_issueDate;
	double C_maturityDate;
	double C_firstCouponDate;
	double C_couponRate;
	double C_redemptionPrice;
	double C_periodicity;

	CCString C_dayCount_str;
	double C_dayCount_double;
	long dayCount_type;

	double C_settleGap;
	double C_settleGap_default = 3;

	double C_couponDateFlag;
	double C_couponDateFlag_default = 3;

	double C_tec;
	double C_tec_default =0;

	CCString C_ccyId;
	long l_ccyId; 

	CCString C_pfId;
	long l_pfId;

	CCString C_modId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_issueDate,C_issueDate," ARM_ERR: issue date: date expected",C_result);
	XL_readNumCell(XL_maturityDate,C_maturityDate," ARM_ERR: maturity date: date expected",C_result);
	XL_readNumCell(XL_firstCouponDate,C_firstCouponDate," ARM_ERR: first coupon date: numeric expected",C_result);
	XL_readNumCell(XL_couponRate,C_couponRate," ARM_ERR: coupon rate: numeric expected",C_result);
	XL_readNumCell(XL_redemptionPrice,C_redemptionPrice," ARM_ERR: redemption price: numeric expected",C_result);
	XL_readNumCell(XL_periodicity,C_periodicity," ARM_ERR: periodicity: numeric expected",C_result);
	XL_readStrOrNumCell(XL_dayCount,C_dayCount_str,C_dayCount_double,dayCount_type," ARM_ERR: day count: string or numeric expected",C_result);
	XL_readNumCellWD(XL_settleGap,C_settleGap,C_settleGap_default," ARM_ERR: settle gap: numeric expected",C_result);
	XL_readNumCellWD(XL_couponDateFlag,C_couponDateFlag,C_couponDateFlag_default," ARM_ERR: coupon date flag: numeric expected",C_result);
	XL_readStrCellWD(XL_ccyId, C_ccyId, "DEFAULT"," ARM_ERR: ccyId: string expected",C_result);
	XL_readNumCellWD(XL_tec, C_tec,C_tec_default," ARM_ERR: tec : numeric expected 0(default) or 1",C_result);
	XL_readStrCellWD(XL_pfId, C_pfId,"DEFAULT"," ARM_ERR: pfTecId : string expected ",C_result);
	XL_readStrCell(XL_modId,C_modId," ARM_ERR: model id: object expected",C_result);

	l_ccyId = LocalGetNumObjectId (C_ccyId);
	l_pfId = LocalGetNumObjectId (C_pfId);

	if(dayCount_type == XL_TYPE_STRING)
	{
		C_dayCount_double = (double)ARM_ConvDayCount (C_dayCount_str);
	}
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_BOND_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_bondTEC (C_issueDate, C_maturityDate, 
			                C_firstCouponDate, C_couponRate,
							C_redemptionPrice, (long)C_periodicity,
							(long)C_dayCount_double, (long)C_settleGap,
							(long)C_couponDateFlag, l_ccyId , C_tec,l_pfId,LocalGetNumObjectId (C_modId), C_result);

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
			if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
			{
				retCode = ARMLOCAL_bondTEC (C_issueDate, C_maturityDate, 
									C_firstCouponDate, C_couponRate,
									C_redemptionPrice, (long)C_periodicity,
									(long)C_dayCount_double, (long)C_settleGap,
									(long)C_couponDateFlag, l_ccyId, C_tec,l_pfId, LocalGetNumObjectId (C_modId), C_result, objId);

				LocalSetCurCellEnvValue (curClass, objId); 
				
				stringId = LocalMakeObjectId (objId, curClass);
			}
			else
				retCode = ARM_KO;			
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_bondTEC (C_issueDate, C_maturityDate, 
			                C_firstCouponDate, C_couponRate,
							C_redemptionPrice, (long)C_periodicity,
							(long)C_dayCount_double, (long)C_settleGap,
							(long)C_couponDateFlag, l_ccyId , C_tec,l_pfId , LocalGetNumObjectId (C_modId), C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BONDTEC" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BOND (LPXLOPER XL_issueDate,
													  LPXLOPER XL_maturityDate,
													  LPXLOPER XL_firstCouponDate,
													  LPXLOPER XL_couponRate,
													  LPXLOPER XL_redemptionPrice,
													  LPXLOPER XL_periodicity,
													  LPXLOPER XL_dayCount,
													  LPXLOPER XL_settleGap,
													  LPXLOPER XL_couponDateFlag,
													  LPXLOPER XL_ccyId)
{
	ADD_LOG("Local_PXL_BOND ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_issueDate;
	double C_maturityDate;
	double C_firstCouponDate;
	double C_couponRate;
	double C_redemptionPrice;
	double C_periodicity;

	CCString C_dayCount_str;
	double C_dayCount_double;
	long dayCount_type;

	double C_settleGap;
	double C_settleGap_default = 3;

	double C_couponDateFlag;
	double C_couponDateFlag_default = 3;

	CCString C_ccyId;
	long l_ccyId; 

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_issueDate,C_issueDate," ARM_ERR: issue date: date expected",C_result);
	XL_readNumCell(XL_maturityDate,C_maturityDate," ARM_ERR: maturity date: date expected",C_result);
	XL_readNumCell(XL_firstCouponDate,C_firstCouponDate," ARM_ERR: first coupon date: numeric expected",C_result);
	XL_readNumCell(XL_couponRate,C_couponRate," ARM_ERR: coupon rate: numeric expected",C_result);
	XL_readNumCell(XL_redemptionPrice,C_redemptionPrice," ARM_ERR: redemption price: numeric expected",C_result);
	XL_readNumCell(XL_periodicity,C_periodicity," ARM_ERR: periodicity: numeric expected",C_result);
	XL_readStrOrNumCell(XL_dayCount,C_dayCount_str,C_dayCount_double,dayCount_type," ARM_ERR: day count: string or numeric expected",C_result);
	XL_readNumCellWD(XL_settleGap,C_settleGap,C_settleGap_default," ARM_ERR: settle gap: numeric expected",C_result);
	XL_readNumCellWD(XL_couponDateFlag,C_couponDateFlag,C_couponDateFlag_default," ARM_ERR: coupon date flag: numeric expected",C_result);
	XL_readStrCellWD(XL_ccyId, C_ccyId, "DEFAULT"," ARM_ERR: ccyId: string expected",C_result);

	l_ccyId=LocalGetNumObjectId (C_ccyId);

	if(dayCount_type == XL_TYPE_STRING)
	{
		C_dayCount_double = (double)ARM_ConvDayCount (C_dayCount_str);
	}
	
	long retCode;
	long objId;
	
	CCString curClass = LOCAL_BOND_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_bond (C_issueDate, C_maturityDate, 
			            C_firstCouponDate, C_couponRate,
						C_redemptionPrice, (long)C_periodicity,
						(long)C_dayCount_double, (long)C_settleGap,
						(long)C_couponDateFlag, l_ccyId, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_BOND" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BONDTEC (LPXLOPER XL_issueDate,
														 LPXLOPER XL_maturityDate,
														 LPXLOPER XL_firstCouponDate,
														 LPXLOPER XL_couponRate,
														 LPXLOPER XL_redemptionPrice,
														 LPXLOPER XL_periodicity,
														 LPXLOPER XL_dayCount,
														 LPXLOPER XL_settleGap,
														 LPXLOPER XL_couponDateFlag,
														 LPXLOPER XL_ccyId,
														 LPXLOPER XL_tec,
														 LPXLOPER XL_pfId,
														 LPXLOPER XL_modId)
{
	ADD_LOG("Local_PXL_BONDTEC ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_issueDate;
	double C_maturityDate;
	double C_firstCouponDate;
	double C_couponRate;
	double C_redemptionPrice;
	double C_periodicity;

	CCString C_dayCount_str;
	double C_dayCount_double;
	long dayCount_type;

	double C_settleGap;
	double C_settleGap_default = 3;

	double C_couponDateFlag;
	double C_couponDateFlag_default = 3;

	double C_tec;
	double C_tec_default =0;

	CCString C_ccyId;
	long l_ccyId; 

	CCString C_pfId;
	long l_pfId;

	CCString C_modId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_issueDate,C_issueDate," ARM_ERR: issue date: date expected",C_result);
	XL_readNumCell(XL_maturityDate,C_maturityDate," ARM_ERR: maturity date: date expected",C_result);
	XL_readNumCell(XL_firstCouponDate,C_firstCouponDate," ARM_ERR: first coupon date: numeric expected",C_result);
	XL_readNumCell(XL_couponRate,C_couponRate," ARM_ERR: coupon rate: numeric expected",C_result);
	XL_readNumCell(XL_redemptionPrice,C_redemptionPrice," ARM_ERR: redemption price: numeric expected",C_result);
	XL_readNumCell(XL_periodicity,C_periodicity," ARM_ERR: periodicity: numeric expected",C_result);
	XL_readStrOrNumCell(XL_dayCount,C_dayCount_str,C_dayCount_double,dayCount_type," ARM_ERR: day count: string or numeric expected",C_result);
	XL_readNumCellWD(XL_settleGap,C_settleGap,C_settleGap_default," ARM_ERR: settle gap: numeric expected",C_result);
	XL_readNumCellWD(XL_couponDateFlag,C_couponDateFlag,C_couponDateFlag_default," ARM_ERR: coupon date flag: numeric expected",C_result);
	XL_readStrCellWD(XL_ccyId, C_ccyId, "DEFAULT"," ARM_ERR: ccyId: string expected",C_result);
	XL_readNumCellWD(XL_tec, C_tec,C_tec_default," ARM_ERR: tec : numeric expected 0(default) or 1",C_result);
	XL_readStrCellWD(XL_pfId, C_pfId,"DEFAULT"," ARM_ERR: pfTecId : string expected ",C_result);
	XL_readStrCell(XL_modId,C_modId," ARM_ERR: model id: object expected",C_result);

	l_ccyId=LocalGetNumObjectId (C_ccyId);
	l_pfId = LocalGetNumObjectId (C_pfId);

	if(dayCount_type == XL_TYPE_STRING)
	{
		C_dayCount_double = (double)ARM_ConvDayCount (C_dayCount_str);
	}
	
	long retCode;
	long objId;
	
	CCString curClass = LOCAL_BOND_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_bondTEC (C_issueDate, C_maturityDate, 
			            C_firstCouponDate, C_couponRate,
						C_redemptionPrice, (long)C_periodicity,
						(long)C_dayCount_double, (long)C_settleGap,
						(long)C_couponDateFlag, l_ccyId, C_tec, l_pfId, LocalGetNumObjectId (C_modId), C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_BONDTEC" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_YTOPRICE (LPXLOPER XL_bondId,
													  LPXLOPER XL_settlement,
													  LPXLOPER XL_yield)
{
	ADD_LOG("Local_YTOPRICE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_bondId;
	double C_settlement;
	double C_yield;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_bondId,C_bondId," ARM_ERR: bond id: object expected",C_result);
	XL_readNumCell(XL_settlement,C_settlement," ARM_ERR: settlement: date expected",C_result);
	XL_readNumCell(XL_yield,C_yield," ARM_ERR: yield: numeric expected",C_result);
	
	long retCode = ARMLOCAL_YTOPRICE (LocalGetNumObjectId (C_bondId), C_settlement, C_yield, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_YTOPRICE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PTOYIELD (LPXLOPER XL_bondId,
													  LPXLOPER XL_settlement,
													  LPXLOPER XL_price)
{
	ADD_LOG("Local_PTOYIELD ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_bondId;
	double C_settlement;
	double C_price;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_bondId,C_bondId," ARM_ERR: bond id: object expected",C_result);
	XL_readNumCell(XL_settlement,C_settlement," ARM_ERR: settlement: date expected",C_result);
	XL_readNumCell(XL_price,C_price," ARM_ERR: price: numeric expected",C_result);
	
	long retCode = ARMLOCAL_PTOYIELD (LocalGetNumObjectId (C_bondId), C_settlement, C_price, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PTOYIELD" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ZEROCOUPON (LPXLOPER XL_maturityDate,
														LPXLOPER XL_dayCount,
														LPXLOPER XL_issueDate,
														LPXLOPER XL_settleGap)
{
	ADD_LOG("Local_ZEROCOUPON ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_maturityDate;
	
	CCString C_dayCount;
	long dayCountId;
	
	double C_issueDate;
	double C_issueDate_default = 0;
	
	double C_settleGap;
	double C_settleGap_default = 0;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_maturityDate,C_maturityDate," ARM_ERR: maturity date: date expected",C_result);
	XL_readStrCellWD(XL_dayCount,C_dayCount,"1"," ARM_ERR: day count: string or numeric expected",C_result);
	XL_readNumCellWD(XL_issueDate,C_issueDate,C_issueDate_default," ARM_ERR: issue date: date expected",C_result);
	XL_readNumCellWD(XL_settleGap,C_settleGap,C_settleGap_default," ARM_ERR: settle gap: numeric expected",C_result);
	
	dayCountId = ARM_ConvDayCount (C_dayCount);

	if(C_issueDate == C_issueDate_default)
	{
		C_issueDate = C_maturityDate - 20 * 365;
	}

	double C_firstCouponDate = C_issueDate;
	double C_couponRate = 0;
	double C_redemptionPrice = 100;
	double C_periodicity = 1;
	double C_couponDateFlag = 3;
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_BOND_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_bond (C_issueDate, C_maturityDate, 
			                C_firstCouponDate, C_couponRate,
							C_redemptionPrice, (long)C_periodicity,
							dayCountId, (long)C_settleGap,
							(long)C_couponDateFlag,-1, C_result);

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
			retCode = ARMLOCAL_bond (C_issueDate, C_maturityDate, 
			                    C_firstCouponDate, C_couponRate,
							    C_redemptionPrice, (long)C_periodicity,
							    dayCountId, (long)C_settleGap,
							    (long)C_couponDateFlag,-1,C_result, objId);
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_bond (C_issueDate, C_maturityDate, 
			                    C_firstCouponDate, C_couponRate,
							    C_redemptionPrice, (long)C_periodicity,
							    dayCountId, (long)C_settleGap,
							    (long)C_couponDateFlag,-1,C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ZEROCOUPON" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZEROCOUPON (LPXLOPER XL_maturityDate,
															LPXLOPER XL_dayCount,
															LPXLOPER XL_issueDate,
															LPXLOPER XL_settleGap)
{
	ADD_LOG("Local_PXL_ZEROCOUPON ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_maturityDate;
	
	CCString C_dayCount;
	long dayCountId;
	
	double C_issueDate;
	double C_issueDate_default = 0;
	
	double C_settleGap;
	double C_settleGap_default = 0;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_maturityDate,C_maturityDate," ARM_ERR: maturity date: date expected",C_result);
	XL_readStrCellWD(XL_dayCount,C_dayCount,"1"," ARM_ERR: day count: string or numeric expected",C_result);
	XL_readNumCellWD(XL_issueDate,C_issueDate,C_issueDate_default," ARM_ERR: issue date: date expected",C_result);
	XL_readNumCellWD(XL_settleGap,C_settleGap,C_settleGap_default," ARM_ERR: settle gap: numeric expected",C_result);
	
	dayCountId = ARM_ConvDayCount (C_dayCount);

	if(C_issueDate == C_issueDate_default)
	{
		C_issueDate = C_maturityDate - 20 * 365;
	}

	double C_firstCouponDate = C_issueDate;
	double C_couponRate = 0;
	double C_redemptionPrice = 100;
	double C_periodicity = 1;
	double C_couponDateFlag = 3;
	
	long retCode;
	long objId;
	
	CCString curClass = LOCAL_BOND_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_bond (C_issueDate, C_maturityDate, 
			            C_firstCouponDate, C_couponRate,
						C_redemptionPrice, (long)C_periodicity,
						dayCountId, (long)C_settleGap,
						(long)C_couponDateFlag,-1, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ZEROCOUPON" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BDFAPRICE (LPXLOPER XL_bondId,
													   LPXLOPER XL_settlement,
													   LPXLOPER XL_actuPrice,
													   LPXLOPER XL_forwardDate,
													   LPXLOPER XL_repoRate)
{
	ADD_LOG("Local_BDFAPRICE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_bondId;
	double C_settlement;
	double C_actuPrice;
	double C_forwardDate;
	double C_repoRate;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_bondId,C_bondId," ARM_ERR: bond id: object expected",C_result);
	XL_readNumCell(XL_settlement,C_settlement," ARM_ERR: settlement: date expected",C_result);
	XL_readNumCell(XL_actuPrice,C_actuPrice," ARM_ERR: actualisation price: numeric expected",C_result);
	XL_readNumCell(XL_forwardDate,C_forwardDate," ARM_ERR: forward date: date expected",C_result);
	XL_readNumCell(XL_repoRate,C_repoRate," ARM_ERR: repo rate: numeric expected",C_result);
				
	long retCode = ARMLOCAL_BDFAPRICE (LocalGetNumObjectId (C_bondId), C_settlement,
							 	  C_actuPrice, C_forwardDate, C_repoRate, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BDFAPRICE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BDREPORATE (LPXLOPER XL_bondId,
														LPXLOPER XL_settlement,
														LPXLOPER XL_actuPrice,
														LPXLOPER XL_forwardDate,
														LPXLOPER XL_forwardPrice)
{
	ADD_LOG("Local_BDREPORATE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_bondId;
	double C_settlement;
	double C_actuPrice;
	double C_forwardDate;
	double C_forwardPrice;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_bondId,C_bondId," ARM_ERR: bond id: object expected",C_result);
	XL_readNumCell(XL_settlement,C_settlement," ARM_ERR: settlement: date expected",C_result);
	XL_readNumCell(XL_actuPrice,C_actuPrice," ARM_ERR: actualisation price: numeric expected",C_result);
	XL_readNumCell(XL_forwardDate,C_forwardDate," ARM_ERR: forward date: date expected",C_result);
	XL_readNumCell(XL_forwardPrice,C_forwardPrice," ARM_ERR: forward price: numeric expected",C_result);
				
	long retCode = ARMLOCAL_BDREPORATE (LocalGetNumObjectId (C_bondId), C_settlement,
							 	   C_actuPrice, C_forwardDate, C_forwardPrice, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BDREPORATE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_YTODURATION (LPXLOPER XL_bondId,
														 LPXLOPER XL_settlement,
														 LPXLOPER XL_actuRate,
														 LPXLOPER XL_flagCpn)
{
	ADD_LOG("Local_YTODURATION ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_bondId;
	double C_settlement;
	double C_actuRate;
	
	double C_flagCpn;
	double C_flagCpn_default = 1;
		
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_bondId,C_bondId," ARM_ERR: bond id: object expected",C_result);
	XL_readNumCell(XL_settlement,C_settlement," ARM_ERR: settlement: date expected",C_result);
	XL_readNumCell(XL_actuRate,C_actuRate," ARM_ERR: actualisation rate: numeric expected",C_result);
	XL_readNumCellWD(XL_flagCpn,C_flagCpn,C_flagCpn_default," ARM_ERR: flag coupon: numeric expected",C_result);
			
	long retCode = ARMLOCAL_YTODURATION (LocalGetNumObjectId (C_bondId),
										 C_settlement,
										 C_actuRate,
										 C_flagCpn,
										 C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_YTODURATION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_YTOCONVEXITY (LPXLOPER XL_bondId,
														  LPXLOPER XL_settlement,
														  LPXLOPER XL_actuRate)
{
	ADD_LOG("Local_YTOCONVEXITY ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_bondId;
	double C_settlement;
	double C_actuRate;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_bondId,C_bondId," ARM_ERR: bond id: object expected",C_result);
	XL_readNumCell(XL_settlement,C_settlement," ARM_ERR: settlement: date expected",C_result);
	XL_readNumCell(XL_actuRate,C_actuRate," ARM_ERR: actualisation rate: numeric expected",C_result);
			
	long retCode = ARMLOCAL_YTOCONVEXITY (LocalGetNumObjectId (C_bondId),
															   C_settlement,
															   C_actuRate,
															   C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_YTOCONVEXITY" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
