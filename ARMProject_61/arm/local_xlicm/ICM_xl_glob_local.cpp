
#pragma warning (disable:4005 4786 )

#include "ICMKernel\glob\icm_mktdatamng.h"
#include <libCCxll\CCxll.h>
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm\ARM_result.h>

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libicm_local\ICM_local_glob.h>
#include <ARM\libicm_local\ICM_local_pwccurve.h>
#include <ARM\libicm_local\ICM_local_pricer.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"
#include "XL_local_xlarm_common.h"

#include "ICMKernel\util\icm_qmatrix.h"


#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>

#include "ExcelTools.h"

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Version()
{
	static XLOPER XL_result;
	try 
	{
		ExcelTools::convert(ICMLOCAL_Version(),&XL_result); 
		return &XL_result ;
	}
	catch(...) 
	{
		ExcelCaller::get().setError("Unknown error"); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return &XL_result; 
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DisplaySV (LPXLOPER XL_pricerId,
														   LPXLOPER XL_typeValues)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	LPXLOPER pxArray;
	try {
	ARM_result C_result;

	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pricerId;
	
	CCString C_typeValues;
	long typeValuesId;

	CCString C_RecOrPay = "R";
	long recId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: pricer id: object expected",C_result);
	XL_readStrCell(XL_typeValues,C_typeValues," ARM_ERR: type Values: string expected",C_result);
	//XL_readStrCellWD(XL_RecOrPay,C_RecOrPay,"R"," ARM_ERR: Rec or Pay: string expected",C_result);

	if((typeValuesId = ARM_ConvTypeValues (C_typeValues, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((recId = ARM_ConvRecOrPay (C_RecOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ARM_Credit_DisplayScheduleValues(LocalGetNumObjectId(C_pricerId),
													  typeValuesId,
													  recId,
													  C_result);

	if(retCode == ARM_OK)
	{
		VECTOR<double> dVal;
		long vecSize;
		retCode = ExtractVectorDoubleFromFile("123",dVal,vecSize);

		if ( (retCode == ARM_OK) && (vecSize != 1) )
		{
			int nbrows = vecSize;
			int nbcolumns = 1;
		
			FreeCurCellErr ();
			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = nbcolumns;
			XL_result.val.array.rows = nbrows; 
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

			for(int i = 0; i < nbrows; i++)
			{
				pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeNum;
				pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.num = dVal[i]; 
			}
		}
		else if ( (retCode == ARM_OK) && (vecSize == 1) )
		{
			int nbrows = 2;
			int nbcolumns = 1;

			FreeCurCellErr ();
			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = nbcolumns;
			XL_result.val.array.rows = nbrows; 
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

			pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].val.num = dVal[0];

			pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].xltype = xltypeErr;
			pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].val.err = xlerrNA;
		}
		else
		{
			ARM_ERR();
		}
	}
	else
	{
		ARM_ERR();
	}

//	ARM_END();
	
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DisplaySD (LPXLOPER XL_instId,
															LPXLOPER XL_typeDates)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	try {

	LPXLOPER pxArray;
	ARM_result C_result;

	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_instId;
	
	CCString C_typeDates;
	long typeDatesId;

	CCString C_RecOrPay="R";
	long recId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_instId,C_instId," ARM_ERR: instrument id: object expected",C_result);
	XL_readStrCell(XL_typeDates,C_typeDates," ARM_ERR: type Values: string expected",C_result);

	if((typeDatesId = ARM_ConvTypeDates (C_typeDates, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((recId = ARM_ConvRecOrPay (C_RecOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ARM_Credit_DisplayScheduleDates(LocalGetNumObjectId(C_instId),
													 typeDatesId,
													 recId,
													 C_result);

	if(retCode == ARM_OK)
	{
		VECTOR<CCString> dDate;
		long vecSize;
		retCode = ExtractVectorDateFromFile("123",dDate,vecSize);

		if ( (retCode == ARM_OK) && (vecSize != 1) )
		{
			int nbrows = vecSize;
			int nbcolumns = 1;

			FreeCurCellErr ();
			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = nbcolumns;
			XL_result.val.array.rows = nbrows; 
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

			for(int i = 0; i < nbrows; i++)
			{
				pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeNum;
				pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.num = Local_ARMDATE2XLDATE(dDate[i]);
			}
		}
		else if ( (retCode == ARM_OK) && (vecSize == 1) )
		{
			int nbrows = 2;
			int nbcolumns = 1;

			FreeCurCellErr ();
			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = nbcolumns;
			XL_result.val.array.rows = nbrows; 
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

			pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].val.num = Local_ARMDATE2XLDATE(dDate[0]);

			pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].xltype = xltypeErr;
			pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].val.err = xlerrNA;
		}
		else
		{
			ARM_ERR();
		}
	}
	else
	{
		ARM_ERR();
	}

//	ARM_END();
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Spread (LPXLOPER XL_pricerId, 
														 LPXLOPER XL_spread)
{

	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable

	double Mtm = 0.0;
	double Mtm_default = 0.0;

	// C variable
	CCString C_pricerId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: pricer id: object expected",C_result);
	XL_readNumCellWD(XL_spread, Mtm, Mtm_default,  " ARM_ERR: MtM: numeric expected",C_result);

	long retCode = ICMLOCAL_Spread (LocalGetNumObjectId (C_pricerId), 
									Mtm, 
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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Spread" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_NPV (LPXLOPER XL_pricerId, 
														 LPXLOPER XL_Option)
{


	static XLOPER XL_result;
	ARM_result C_result;


	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// error
		static int error;
		static char* reason = "";

		std::string pricerId ;
		ExcelTools::convert(XL_pricerId,pricerId); 
 		qCMPMETH measure ;
		ExcelTools::econvert(XL_Option,"NPV",measure); 
		
		double res = ICMLOCAL_NPV (LocalPersistent::get().getObjectId(pricerId), measure); 

 
		ExcelTools::convert(res,&XL_result); 
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_NPV" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_PriceVector (LPXLOPER XL_pricerId, 
														 LPXLOPER XL_Option)
{


	static XLOPER XL_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();


		std::string pricerId ; ExcelTools::convert(XL_pricerId,pricerId); 
		std::string measure; ExcelTools::convert(XL_Option,measure); 
		
		ARM_Vector output; 
		ICMLOCAL_PriceVector(LocalPersistent::get().getObjectId(pricerId), measure,output); 
 
		ExcelTools::convert(output,&XL_result); 
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Price (LPXLOPER XL_pricerId,
														LPXLOPER XL_AsOfDate)
{

	static XLOPER XL_result;
	ARM_result C_result;


	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pricerId;

	double AsOfDate = 0.;
	double AsOfDate_default = -1.;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: security id: object expected",C_result);
	XL_readNumCellWD(XL_AsOfDate,AsOfDate,AsOfDate_default," ARM_ERR: AsOfDate: date expected",C_result);
	
	long retCode = ICMLOCAL_Price (LocalGetNumObjectId (C_pricerId),
								   AsOfDate,
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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Price" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/**
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_PriceVector (LPXLOPER XL_pricerId,
														LPXLOPER XL_AsOfDate)
{

	static XLOPER XL_result;
	ARM_result C_result;


	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pricerId;

	double AsOfDate = 0.;
	double AsOfDate_default = -1.;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: security id: object expected",C_result);
	XL_readNumCellWD(XL_AsOfDate,AsOfDate,AsOfDate_default," ARM_ERR: AsOfDate: date expected",C_result);
	
	ARM_Vector output; 
	ICMLOCAL_PriceVector (LocalGetNumObjectId (C_pricerId),
								   AsOfDate,
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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Price" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
**/
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GetLabel (LPXLOPER XL_CurveId )
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_CurveId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_CurveId,C_CurveId," ARM_ERR: Pay Yield Curve id: object expected",C_result);
	
	long retCode = ICMLOCAL_GetLabel (LocalGetNumObjectId (C_CurveId), C_result);

	if(retCode == ARM_OK)
	{
        FreeCurCellErr ();

        XL_result.xltype = xltypeStr;
        XL_result.val.str = XL_StrC2StrPascal (C_result.getMsg());
        XL_result.xltype |= xlbitDLLFree;
	}
	else
	{
		ARM_ERR();
	}

//	ARM_END();
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_SetLabel (LPXLOPER XL_CurveId,
														   LPXLOPER XL_label)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_CurveId;
	CCString C_label;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_CurveId,C_CurveId," ARM_ERR: Pay Yield Curve id: object expected",C_result);
	XL_readStrCell(XL_label,C_label," ARM_ERR: label : string expected",C_result);
	
	long retCode = ICMLOCAL_SetLabel (LocalGetNumObjectId (C_CurveId), C_label,C_result);

	if(retCode == ARM_OK)
	{
        FreeCurCellErr ();

        XL_result.xltype = xltypeStr;
        XL_result.val.str = XL_StrC2StrPascal (C_result.getMsg());
        XL_result.xltype |= xlbitDLLFree;
	}
	else
	{
		ARM_ERR();
	}

//	ARM_END();
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DTR (LPXLOPER XL_pricerId,
													  LPXLOPER XL_NbDefaults,
													  LPXLOPER XL_SL,
													  LPXLOPER XL_Labels,
													  LPXLOPER XL_RecoveryRates)
{
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	CCString C_pricerId;
    CCString C_SL;
	CCString C_SL_default = "S";
	VECTOR<CCString> C_Labels ; 
	VECTOR<double> C_RecoveryRates;
	double C_NbDefaults = 0 ;

	int l_NbDefaults = 0;
	

    
	long retCode = 0;
	

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: Pricer object expected",C_result);
	XL_readNumCell(XL_NbDefaults,C_NbDefaults," ARM_ERR: Nbpaths numeric expected",C_result);
	XL_readStrCellWD(XL_SL,C_SL,C_SL_default," ARM_ERR: Plot: string ex:S, L expected",C_result);
	XL_readStrVector(XL_Labels,C_Labels," ARM_ERR: Labels: array of strings expected",DOUBLE_TYPE,C_result);
	XL_readNumVector(XL_RecoveryRates,C_RecoveryRates," ARM_ERR: Recovery Rates array: numerics expected",C_result);

	l_NbDefaults = (int) C_NbDefaults;

	if ( C_RecoveryRates.size() != C_Labels.size() || C_RecoveryRates.size() != C_NbDefaults || C_NbDefaults != C_Labels.size() )
	{
		C_result.setMsg("nb of Labels incompatible with nb of Recovery Rates or with Nb of defaults");
		ARM_ERR();
	}
	

	retCode = ICMLOCAL_DTR(LocalGetNumObjectId (C_pricerId),
						  l_NbDefaults,
						  C_SL,
						  C_Labels,
						  C_RecoveryRates,
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

	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Sensitivity (LPXLOPER XL_pricerId,
															  LPXLOPER XL_crvtype,
															  LPXLOPER XL_plot,
															  LPXLOPER XL_label,
															  LPXLOPER XL_Epsilon,
															  LPXLOPER XL_EpsilonGamma)
{
	VECTOR<CCString> VectorOut;
	int NbCol=0;
	VECTOR<CCString> VectorName;

	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pricerId;
    CCString C_plot;
	CCString C_plot_default = "NONE";
    CCString C_label;
	CCString C_label_default = "NONE";

	long retCode = 0;
	// long crvtype = 0;

	double C_Epsilon = 0.;
	double C_Epsilon_default = -999.;
	double C_EpsilonGamma = 0.;
	double C_EpsilonGamma_default = 0.;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: Pricer object expected",C_result);
	XL_readStrCellWD(XL_plot,C_plot,C_plot_default," ARM_ERR: Plot: string ex:2D,1M.. expected",C_result);
	XL_readStrCellWD(XL_label,C_label,C_label_default," ARM_ERR: label: curve label expected",C_result);
	XL_readNumCellWD(XL_Epsilon,C_Epsilon,C_Epsilon_default ," ARM_ERR: numeric expected",C_result);
	XL_readNumCellWD(XL_EpsilonGamma, C_EpsilonGamma, C_EpsilonGamma_default , "ARM_ERR: numeric expected",C_result);

	qSENSITIVITY_TYPE sensiType; 	
	ExcelTools::econvert(XL_crvtype,"NONE",sensiType) ;
	
	long CurveId = LocalGetNumObjectId (C_label);
	ARM_Object* curve = NULL;
	curve = (ARM_Object* )LOCAL_PERSISTENT_OBJECTS->GetObject(CurveId);
	std::string s_label("");
	if (curve)
	{	
		deducelabelforobject(*curve, s_label);
		
	} else {
		s_label = CCSTringToSTLString(C_label);
	}


	retCode = ICMLOCAL_Sensitivity (LocalGetNumObjectId (C_pricerId), 
									sensiType,
									CCSTringToSTLString(C_plot),
									s_label,
									C_Epsilon,
									C_EpsilonGamma,
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

// 
// 	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Sensitivity" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_SurvivalProba (LPXLOPER XL_defcurveId, 
														 LPXLOPER XL_maturity)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable

	double C_maturity = 0.0;

	CCString C_defcurveId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_defcurveId,C_defcurveId," ARM_ERR: curve id: object expected",C_result);
	XL_readNumCell(XL_maturity, C_maturity, " ARM_ERR: maturity : numeric expected",C_result);
	
	long retCode = ICMLOCAL_SurvivalProba (LocalGetNumObjectId(C_defcurveId),C_maturity,C_result);

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
		catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DefaultIntensity (LPXLOPER XL_idCurve,
															LPXLOPER XL_maturity)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable

	double C_maturity = 0.0;
	double C_meth = 0.0;

	CCString C_idCurve;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_idCurve, C_idCurve," ARM_ERR: curve id: object expected",C_result);
	XL_readNumCell(XL_maturity, C_maturity, " ARM_ERR: maturity : numeric expected",C_result);
	//XL_readNumCell(XL_meth, C_meth, " ARM_ERR: maturity : numeric expected",C_result);
	
	long retCode = ICMLOCAL_DefaultIntensity (LocalGetNumObjectId (C_idCurve), 
											  C_maturity, 
											  (long) C_meth,
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
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DefaultProba (LPXLOPER XL_defcurveId, 
														 LPXLOPER XL_maturity)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;


	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable

	double C_maturity = 0.0;

	CCString C_defcurveId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_defcurveId,C_defcurveId," ARM_ERR: curve id: object expected",C_result);
	XL_readNumCell(XL_maturity, C_maturity, " ARM_ERR: maturity : numeric expected",C_result);
	
	long retCode = ICMLOCAL_DefaultProba (LocalGetNumObjectId (C_defcurveId), 
																C_maturity, 
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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_DefaultProba" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DisplayMatrix (LPXLOPER XL_pricerId)
{
	//	ARM_BEGIN();
	VECTOR<CCString> VectorOut;
	int NbCol=0;
	VECTOR<CCString> VectorName;
	

	// return
	static XLOPER XL_result;
	LPXLOPER pxArray;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pricerId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: instrument id: object expected",C_result);

	long retCode = ICMLOCAL_DisplayMatrix(LocalGetNumObjectId(C_pricerId),
										  VectorOut,
										  NbCol,
										  VectorName,
										  C_result);
	if(retCode == ARM_OK)
	{
		int nbrows = VectorOut.size()/NbCol;
		int nbcolumns = NbCol;
		
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows+1; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, (nbrows+2) * nbcolumns * sizeof (XLOPER));

		for(int k = 0; k < nbcolumns; k++)
		{
			pxArray[XL_Coordonnate2Rank (0, k, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (0, k, nbcolumns)].val.str = XL_StrC2StrPascal (VectorName[k]);
			pxArray[XL_Coordonnate2Rank (0, k, nbcolumns)].xltype |= xlbitDLLFree;
		}

		for(int j = 0; j < nbcolumns; j++)
		{
			for(int i = 0; i < nbrows; i++)
			{
			pxArray[XL_Coordonnate2Rank (i+1, j, nbcolumns)].xltype = xltypeNum;

			if (j==0)
			{
				pxArray[XL_Coordonnate2Rank (i+1, j, nbcolumns)].val.num = Local_ARMDATE2XLDATE(VectorOut[i+nbrows*j]); 
			}	
			else
				pxArray[XL_Coordonnate2Rank (i+1, j, nbcolumns)].val.num = atof((const char*)(VectorOut[i+nbrows*j]));
			}
		}
	}
	else
	{
		ARM_ERR();
	}

	VectorName.clear();
	VectorOut.clear();

//	ARM_END();
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_RiskyPV01(LPXLOPER XL_DefCurveId,
														   LPXLOPER XL_Date1,
														   LPXLOPER XL_Date2)
{
	

	// return
	static XLOPER XL_result;
	

	/// to make the infrastructure very robust we put a try catch!
	try 
	{
	ARM_NOCALCIFWIZ();

	std::string defCurve ; ExcelTools::convert(XL_DefCurveId,defCurve); 
	double xlDate ; ExcelTools::convert(XL_Date1,-1.,xlDate); 
	ARM_Date date1 ; if (xlDate!=-1) ExcelTools::convert(XL_Date1,date1); 
	long defCurveId = LocalPersistent::get().getObjectId(defCurve); 
	double res ;  
	if (ExcelTools::isXLDate(XL_Date2)) 
	{
		ARM_Date date2 ; ExcelTools::convert(XL_Date2,date2); 
		if (xlDate==-1) res = ICMLOCAL_RiskyPV01(defCurveId,NULL,date2); 
		else res = ICMLOCAL_RiskyPV01(defCurveId,&date1,date2); 
	}
	else
	{
		std::string tenor; ExcelTools::convert(XL_Date2,tenor); 
		if (xlDate==-1) res = ICMLOCAL_RiskyPV01(defCurveId,NULL,tenor); 
		else res = ICMLOCAL_RiskyPV01(defCurveId,&date1,tenor); 

	}
	ExcelTools::convert(res,&XL_result); 
	return &XL_result; 
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;

}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_RiskyPV01AsSensitivity(LPXLOPER XL_DefCurveId,
														   LPXLOPER XL_Tenor)
{
	

	// return
	static XLOPER XL_result;
	

	/// to make the infrastructure very robust we put a try catch!
	try 
	{
		ARM_NOCALCIFWIZ();
		std::string defCurve ; ExcelTools::convert(XL_DefCurveId,defCurve); 
		std::string tenor ;  ExcelTools::convert(XL_Tenor,tenor); 
		long defCurveId = LocalPersistent::get().getObjectId(defCurve); 
		double res = ICMLOCAL_RiskyPV01AsSensitivity(defCurveId,tenor); 
		ExcelTools::convert(res,&XL_result); 
		return &XL_result; 
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_RiskyDuration(LPXLOPER XL_DefCurveId,
														   LPXLOPER XL_Date)
{
	// return
	static XLOPER XL_result;
	
	/// to make the infrastructure very robust we put a try catch!
	try 
	{
	ARM_NOCALCIFWIZ();

	std::string defCurve ; ExcelTools::convert(XL_DefCurveId,defCurve); 
	long defCurveId = LocalPersistent::get().getObjectId(defCurve); 
	double res ;  
	if (ExcelTools::isXLDate(XL_Date)) 
	{
		ARM_Date date  ; ExcelTools::convert(XL_Date,date); 
		res = ICMLOCAL_RiskyDuration(defCurveId,date); 
	}
	else
	{
		std::string tenor; ExcelTools::convert(XL_Date,tenor); 
		res = ICMLOCAL_RiskyDuration(defCurveId,tenor); 
	}
	ExcelTools::convert(res,&XL_result); 
	return &XL_result; 
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
	
}


/** 
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_RiskyDuration(LPXLOPER XL_DefCurveId,
															   LPXLOPER XL_Date,
															   LPXLOPER XL_Tenor)
{
	int NbCol=0;

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_DefCurveId;
    double C_Date = 0.;
	long retCode = 0;
	CCString C_Tenor;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_DefCurveId,C_DefCurveId," ARM_ERR: DefCurve id: object expected",C_result);
	XL_readNumCell(XL_Date,C_Date," ARM_ERR: numeric expected",C_result);
	XL_readStrCellWD(XL_Tenor,C_Tenor,"NONE","ARM_ERR: Tenor : string expected",C_result);

	retCode = ICMLOCAL_RiskyDuration (LocalGetNumObjectId(C_DefCurveId),C_Date,C_Tenor,C_result);

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

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_RiskyDuration" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

**/ 
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Delivery(LPXLOPER XL_AsOfDate,
														  LPXLOPER XL_Maturity)
{
	VECTOR<double> VectorOut;
	int NbCol=0;

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable
	double C_AsOfDate= 0.;
	CCString C_Maturity;
   
	long retCode = 0;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_AsOfDate,C_AsOfDate," ARM_ERR: Excel date expected",C_result);
	XL_readStrCell(XL_Maturity,C_Maturity," ARM_ERR: String expected",C_result);
	
	retCode = ICMLOCAL_Delivery(C_AsOfDate,C_Maturity,C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		ARM_ERR();
	}
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GetBeta (LPXLOPER XL_pricerId,
														LPXLOPER XL_issuer)
{
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pricerId;
	CCString C_issuer;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_issuer,C_issuer," ARM_ERR: Issuer: string expected",C_result);
	
	long retCode = ICMLOCAL_GetBeta (LocalGetNumObjectId (C_pricerId),
									 C_issuer,	
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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_GetBeta" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CorrMatrix(LPXLOPER XL_Labels,
														   LPXLOPER XL_Coefs,
														   LPXLOPER XL_AsOf,
														   LPXLOPER XL_NAME)
{
	// return
	static XLOPER XL_result ;
	ARM_result C_result;
		
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<std::string> C_Labels; 
	VECTOR<double> C_Coefs; 

	std::string C_Name; 
	
	// error
	static int error;
	static char* reason = "";


	// XL_readStrVector(XL_Labels,C_Labels," ARM_ERR: Labels: array of strings expected",DOUBLE_TYPE,C_result);
	ExcelTools::convert(XL_Labels,C_Labels) ;
	XL_readNumVector(XL_Coefs,C_Coefs," ARM_ERR: Matrix correlation array: numerics expected",C_result);
	ARM_Date AsOf; ExcelTools::convert(XL_AsOf,AsOf) ;

	// XL_readStrCellWD(XL_NAME,C_Name,"USERDEF"," ARM_ERR: pricer id: object expected",C_result);
	ExcelTools::convert(XL_NAME,"USERDEF",C_Name); 
	
	long retCode;
	long objId;
	CCString prevClass;

	double ValueTest = C_Coefs.size()/C_Labels.size();

	if (ValueTest != C_Labels.size())
	{
		C_result.setMsg("Matrix is not symetric or nb of Labels incompatible with Matrix");
		ARM_ERR();
	}

	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;
	CCString stringId = GetLastCurCellEnvValue();
	
	if (!stringId)
	{
		retCode = ICMLOCAL_CORRMATRIX(AsOf,
									  C_Name,	
									  C_Labels,
									  C_Coefs,				  
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
			
		if (curClass == prevClass)
		{
		retCode = ICMLOCAL_CORRMATRIX(AsOf,
									  C_Name,	
									  C_Labels,
									  C_Coefs,				  
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

		retCode = ICMLOCAL_CORRMATRIX(AsOf,
									  C_Name,	
									  C_Labels,
									  C_Coefs,				  
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
		C_Labels.clear(); 
		C_Labels.empty();
		C_Coefs.clear(); 
		C_Coefs.empty();


		ARM_ERR();
	}


	C_Labels.clear(); 
	C_Labels.empty();
	C_Coefs.clear(); 
	C_Coefs.empty();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ..." )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_ExtractCorrelations(LPXLOPER XL_CorrMatrixId,
																	 LPXLOPER XL_labels)
{
	//	ARM_BEGIN();
	int NbCol=0;
	VECTOR<double> VectorOut;
	VECTOR<std::string> C_Labels;

	// return
	static XLOPER XL_result;
	LPXLOPER pxArray;
	try {
	ARM_result C_result;

	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_CorrMatrixId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_CorrMatrixId,C_CorrMatrixId," ARM_ERR: Correlation Matrix id: object expected",C_result);
	// XL_readStrVector(XL_labels, C_Labels ," ARM_ERR: Issuers Labels : string vector expected",DOUBLE_TYPE,C_result);
	ExcelTools::convert(XL_labels, C_Labels ); 

	long retCode = ICMLOCAL_EXTRACTCORRMATRIX(LocalGetNumObjectId(C_CorrMatrixId),
										  C_Labels,
										  VectorOut,
										  C_result);
	if(retCode == ARM_OK)
	{
		NbCol = C_Labels.size();
		int nbrows = VectorOut.size()/NbCol;
		int nbcolumns = NbCol;
		
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows* nbcolumns * sizeof (XLOPER));

		for(int j = 0; j < nbcolumns; j++)
		{
			for(int i = 0; i < nbrows; i++)
			{
			pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.num = VectorOut[i+nbrows*j];
			}
		}
	}
	else
	{
		ARM_ERR();
	}

	C_Labels.clear();
	C_Labels.empty();
	VectorOut.clear();
	VectorOut.empty();

//	ARM_END();
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GetDefProbTranche (LPXLOPER XL_pricerId, 
																	LPXLOPER XL_YF)
{

	
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pricerId;
	double C_YF;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: pricer id: object expected",C_result);
	XL_readNumCell(XL_YF, C_YF, " ARM_ERR: Maturity expected",C_result);

	long retCode = ICMLOCAL_GetDefProbTranche (LocalGetNumObjectId (C_pricerId), 
											C_YF, 
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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_GetDefProbTranche" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GetDuration (LPXLOPER XL_pricerId)
{

	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	// C variable
	CCString C_pricerId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: pricer id: object expected",C_result);

	long retCode = ICMLOCAL_GetDuration(LocalGetNumObjectId (C_pricerId),
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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_GetDuration" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GenSchedule ( LPXLOPER XL_AccStartDate,
															   LPXLOPER XL_AccEndDate,
															   LPXLOPER XL_FixingFreq,
															   LPXLOPER XL_DayCountFrq,
															   LPXLOPER XL_refdate, 	
															   LPXLOPER XL_Currency,
															   LPXLOPER XL_typeDates,
															   LPXLOPER XL_ModFoll,
															   LPXLOPER XL_GapCredit)

{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	LPXLOPER pxArray;

	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	

	double EffectiveDate;
	double EndDate;

	CCString C_FixingFreq;
	int FixingFreq;

	int DayCountFrq;

	double refdate;
	double refdate_default = -1.;

	CCString Currency;
	CCString C_Currency;

	CCString C_DayCountFrq;

	CCString C_typeDates;
	CCString C_typeDates_default = "STARTDATE";
	long typeDatesId;

	CCString C_ModFoll;
	CCString C_ModFoll_default = "MF";
	long ModFollId;

	double CreditGap =0.;
	double CreditGap_default =-3.;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_AccStartDate,EffectiveDate," ARM_ERR: Effective Date date expected",C_result);
	XL_readNumCell(XL_AccEndDate,EndDate," ARM_ERR: EndDate date: date expected",C_result);
	XL_readStrCellWD(XL_FixingFreq,C_FixingFreq,"Q"," ARM_ERR: FixingFreq: string expected",C_result);
	XL_readStrCellWD(XL_DayCountFrq,C_DayCountFrq,"A360"," ARM_ERR: DayCountFrq: string expected",C_result);
	XL_readNumCellWD(XL_refdate,refdate,refdate_default," ARM_ERR: First_period_refdate Rate: date expected",C_result);
	XL_readStrCellWD(XL_Currency,C_Currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_typeDates,C_typeDates,C_typeDates_default," ARM_ERR: type Values: string expected",C_result);
	XL_readStrCellWD(XL_ModFoll,C_ModFoll,C_ModFoll_default," ARM_ERR: type Values: string expected",C_result);
	XL_readNumCellWD(XL_GapCredit,CreditGap,CreditGap_default," ARM_ERR: Effective Date date expected",C_result);

	ARM_result currencyres;
	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);

	if(C_Currency == "DEFAULT")
	{
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_Currency = currencyres.getString ();
		}
	}


	long retCode;

	if((FixingFreq = ARM_ConvFrequency (C_FixingFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	DayCountFrq = ARM_ConvDayCount (C_DayCountFrq);

	if((typeDatesId = ARM_ConvTypeGenDates (C_typeDates, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((ModFollId = ARM_ConvFwdRule (C_ModFoll, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	retCode = ICMLOCAL_GenSchedule(EffectiveDate,
								   EndDate,
								   refdate,
								   FixingFreq,
							       DayCountFrq,
							       C_Currency,
								   typeDatesId,
								   ModFollId,
								   CreditGap,	
							       C_result);
		
	if(retCode == ARM_OK)
	{
		VECTOR<CCString> dDate;
		long vecSize;
		retCode = ExtractVectorDateFromFile("123",dDate,vecSize);

		if ( (retCode == ARM_OK) && (vecSize != 1) )
		{
			int nbrows = vecSize;
			int nbcolumns = 1;

			FreeCurCellErr ();
			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = nbcolumns;
			XL_result.val.array.rows = nbrows; 
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

			for(int i = 0; i < nbrows; i++)
			{
				pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeNum;
				pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.num = Local_ARMDATE2XLDATE(dDate[i]);
			}
		}
		else if ( (retCode == ARM_OK) && (vecSize == 1) )
		{
			int nbrows = 2;
			int nbcolumns = 1;

			FreeCurCellErr ();
			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = nbcolumns;
			XL_result.val.array.rows = nbrows; 
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

			pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].val.num = Local_ARMDATE2XLDATE(dDate[0]);

			pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].xltype = xltypeErr;
			pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].val.err = xlerrNA;
		}
		else
		{
			ARM_ERR();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_GenSchedule" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/**
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GetProbCondNTD (LPXLOPER XL_pricerId, 
																 LPXLOPER XL_Label)
{

	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pricerId;
	CCString C_Label;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: pricer id: object expected",C_result);
	XL_readStrCell(XL_Label, C_Label, " ARM_ERR: : Tenor expected",C_result);

	// 
	// long retCode = ICMLOCAL_GetProbCondNTD (LocalGetNumObjectId (C_pricerId), 
	// 										C_Label, 
	// 										C_result);
	// 

	long retCode = ARM_KO;

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

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_GetProbCondNTD" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
**/ 
/** 
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CorrelationSmile (LPXLOPER XL_pricerId,
																	LPXLOPER XL_smiletype,
																	LPXLOPER XL_mktprice,
																	LPXLOPER XL_Seed,
																	LPXLOPER XL_UpfrontPayment,
																	LPXLOPER XL_datatype)
{
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();


	CCString C_pricerId;
	double C_smiletype;
	double C_smiletype_default = 0. ;
	double C_datatype;
	double C_datatype_default = 0. ;
    double C_mktprice;
	double C_mktprice_default = 0.;
    double C_Seed;
	double C_Seed_default = 0.;

	long retCode = 0;
	
	double C_Upfront;
	double C_Upfront_default = 0.;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: Pricer object expected",C_result);
	XL_readNumCellWD(XL_smiletype,C_smiletype,C_smiletype_default," ARM_ERR: 0 or 1 expected",C_result);
	XL_readNumCellWD(XL_datatype,C_datatype,C_datatype_default," ARM_ERR: 0 or 1 expected",C_result);
	XL_readNumCellWD(XL_mktprice,C_mktprice,C_mktprice_default," ARM_ERR: numeric expected",C_result);
	XL_readNumCellWD(XL_Seed,C_Seed,C_Seed_default," ARM_ERR: numeric expected",C_result);
	XL_readNumCellWD(XL_UpfrontPayment,C_Upfront,C_Upfront_default ," ARM_ERR: numeric expected",C_result);

	retCode = ICMLOCAL_CorrelationSmile (LocalGetNumObjectId (C_pricerId), 
										 C_smiletype,
										 C_mktprice,
										 C_Seed,
										 C_Upfront,
										 C_datatype,
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
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;

}
**/ 

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CloneCorrMatrixBary(LPXLOPER XL_CorrMatrixId,
																	 LPXLOPER XL_Beta,
																	 LPXLOPER XL_UpOrDown)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{

	ARM_NOCALCIFWIZ();

	// C variable
	
	CCString C_CorrMatrixId;
	double C_Beta = 0.;
	double C_UpOrDown = 0.;

	long retCode;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_CorrMatrixId,C_CorrMatrixId," ARM_ERR: FeeLeg: string expected",C_result);
	XL_readNumCell(XL_Beta,C_Beta," ARM_ERR: Beta: numeric expected",C_result);
	XL_readNumCell(XL_UpOrDown,C_UpOrDown," ARM_ERR: UpOrDown: numeric expected",C_result);
	

	int l_UpOrDown = (int) C_UpOrDown;

	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{	
		retCode = ICMLOCAL_CloneCorrMatrixBary(LocalGetNumObjectId (C_CorrMatrixId),
											   C_Beta,
											   C_UpOrDown,
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
		retCode = ICMLOCAL_CloneCorrMatrixBary(LocalGetNumObjectId (C_CorrMatrixId),
											   C_Beta,
											   C_UpOrDown,
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

		retCode = ICMLOCAL_CloneCorrMatrixBary(LocalGetNumObjectId (C_CorrMatrixId),
											   C_Beta,
											   C_UpOrDown,
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

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_CloneCorrMatrixBary" )


	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_BSGreeks (LPXLOPER XL_pricerId,
														   LPXLOPER XL_greektype
														   )
{

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{

	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pricerId;
	CCString C_greektype;
	CCString C_greektype_default = "DELTA";
    
	long retCode = 0;
	long greektype = 0;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: Pricer object expected",C_result);
	XL_readStrCellWD(XL_greektype,C_greektype,C_greektype_default," ARM_ERR: Greek type : string expected",C_result);
	
	if((greektype = ICM_ConvGreekType (C_greektype, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	retCode = ICMLOCAL_BSGreeks (LocalGetNumObjectId (C_pricerId), 
								 greektype,
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
	
	} // try 
	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_CptLeverageLevels" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Greeks (LPXLOPER XL_pricerId,
														   LPXLOPER XL_greektype
														   )
{

	// return
	static XLOPER XL_result;
	ARM_result C_result;
	ARM_XL_TRY_BLOCK_BEGIN
	{

	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pricerId;
	CCString C_greektype;
	CCString C_greektype_default = "DELTA";
    
	long retCode = 0;
	long greektype = 0;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: Pricer object expected",C_result);
	XL_readStrCellWD(XL_greektype,C_greektype,C_greektype_default," ARM_ERR: Greek type : string expected",C_result);
	
	if((greektype = ICM_ConvGreekType (C_greektype, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	retCode = ICMLOCAL_Greeks (LocalGetNumObjectId (C_pricerId), 
								 greektype,
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

	} // try 
	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_CptLeverageLevels" )
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_ImpliedVol(LPXLOPER XL_pricerId,
															LPXLOPER XL_Price
															)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pricerId;
	double C_Price;
	    
	long retCode = 0;
		
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: Pricer object expected",C_result);
	XL_readNumCell(XL_Price, C_Price,  " ARM_ERR: Price: numeric expected",C_result);
	
	retCode = ICMLOCAL_ImpliedVol (LocalGetNumObjectId (C_pricerId), 
								   C_Price,
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

	/// return the result as an LPXLOPER
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_FwdSpreadPricer(LPXLOPER XL_pricerId,
																 LPXLOPER XL_Maturity1,
																 LPXLOPER XL_Maturity2)
{
	

	// return
	static XLOPER XL_result;
	ARM_result C_result;
	try {
	ARM_NOCALCIFWIZ();

	// C variable

	double C_maturity1 = 0.0;
	double C_maturity2 = 0.0;

	CCString C_pricerId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: pricer id: object expected",C_result);
	// XL_readNumCell(XL_Maturity1, C_maturity1, " ARM_ERR: maturity1 : numeric expected",C_result);
	// XL_readNumCell(XL_Maturity2, C_maturity2, " ARM_ERR: maturity2 : numeric expected",C_result);
	ARM_Date matu1,matu2; 
	ExcelTools::convert(XL_Maturity1, matu1) ; 
	ExcelTools::convert(XL_Maturity2, matu2) ; 
	long retCode = ICMLOCAL_FwdSpreadPricer(LocalGetNumObjectId(C_pricerId),
											matu1,
											matu2,
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
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_VirtualCdsSpread(LPXLOPER XL_pricerId, 
																  LPXLOPER XL_Maturity)
{
	
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pricerId;
	
	double C_maturity = 0.0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: pricer id: object expected",C_result);
	XL_readNumCell(XL_Maturity, C_maturity, " ARM_ERR: maturity : numeric expected",C_result);

	long retCode = ICMLOCAL_VirtualCdsSpread (LocalGetNumObjectId (C_pricerId), 
											  C_maturity, 
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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_GetDefProbTranche" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CorrelationStrike(LPXLOPER XL_Labels,
																  LPXLOPER XL_volcurves,
																  LPXLOPER XL_proportions,
																  LPXLOPER XL_smilestrikelow,
																  LPXLOPER XL_smilestrikehight,
																  LPXLOPER XL_Vindex,
																  LPXLOPER XL_AsOf,
																  LPXLOPER XL_NAME,
																  LPXLOPER XL_FullStrikeLow,
																  LPXLOPER XL_FullStrikeUp
																  )
{

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN 
	{
	ARM_NOCALCIFWIZ();
	ICMLOG(" ARM_Credit_CorrelationStrike is Used !! "); 
	// C variable
	VECTOR<CCString> C_volcurves; 
	VECTOR<CCString> C_index; 
	VECTOR<long> I_volcurves; 
	VECTOR<long> I_index; 

	VECTOR<double> C_proportions; 
	VECTOR<double> C_smilestrikelow; 	
	VECTOR<double> C_smilestrikehight; 
	ICM_QMatrix<double> fullStrikeLow,fullStrikeUp,defValue; 

	vector<CCString> DefaultValue;
	DefaultValue.clear();DefaultValue.empty();

	vector<double> DefaultValueDbl;
	DefaultValueDbl.clear();DefaultValueDbl.empty();

	std::string C_Name; 
	ARM_Date AsOf; 
	std::vector<std::string> C_Labels; 

	// error
	static int error;
	static char* reason = "";

	ExcelTools::convert(XL_Labels,C_Labels); 

	XL_readStrVector(XL_volcurves,C_volcurves," ARM_ERR: Labels: object volcurve expected",DOUBLE_TYPE,C_result);
	XL_readStrVectorWD(XL_Vindex,C_index,DefaultValue," ARM_ERR: Labels: object Credit Index expected",DOUBLE_TYPE,C_result);

	XL_readNumVector(XL_proportions,C_proportions," ARM_ERR: Proportions vector: numerics expected",C_result);
	XL_readNumVectorWD(XL_smilestrikelow,C_smilestrikelow,DefaultValueDbl," ARM_ERR: smile strike low: numerics expected",C_result);
	XL_readNumVectorWD(XL_smilestrikehight,C_smilestrikehight,DefaultValueDbl," ARM_ERR: smile strike hight : numerics expected",C_result);

	ExcelTools::convert(XL_AsOf,AsOf); 
	ExcelTools::convert(XL_NAME,"USERDER",C_Name); 

	ExcelTools::convert(XL_FullStrikeLow,defValue,fullStrikeLow); 
	ExcelTools::convert(XL_FullStrikeUp,defValue,fullStrikeUp); 
	
	long retCode = 0;
	long objId = 0;
	CCString prevClass;

	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;
	CCString stringId = GetLastCurCellEnvValue();
	
	if (C_volcurves.size() > 0)	I_volcurves.resize(C_Labels.size());
	if (C_index.size()>0) I_index.resize(C_Labels.size());

	for (int i=0; i<C_Labels.size(); i++)
	{
		I_volcurves[i] = LocalGetNumObjectId(C_volcurves[i]);
		if (C_index.size()>0) I_index[i] = LocalGetNumObjectId(C_index[i]);
	}


	if (!stringId)
	{
		retCode = ICMLOCAL_CORRELATION_STRIKE(AsOf,
											  C_Name,	
											  C_Labels,
											  I_volcurves, 
											  C_proportions,
											  C_smilestrikelow,
											  C_smilestrikehight,
											  fullStrikeLow,
											  fullStrikeUp,
											  I_index,
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
			
		if (curClass == prevClass)
		{
		retCode = ICMLOCAL_CORRELATION_STRIKE(AsOf,
											  C_Name,	
											  C_Labels,
											  I_volcurves, 
											  C_proportions,
											  C_smilestrikelow,
											  C_smilestrikehight,
											  fullStrikeLow,
											  fullStrikeUp,
											  I_index,
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

		retCode = ICMLOCAL_CORRELATION_STRIKE(AsOf,
											  C_Name,	
											  C_Labels,
											  I_volcurves, 
											  C_proportions,
											  C_smilestrikelow,
											  C_smilestrikehight,
											  fullStrikeLow,
											  fullStrikeUp,
											  I_index,
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
		C_Labels.clear(); 
		C_volcurves.clear(); 
		I_volcurves.clear(); 
		C_proportions.clear(); 
		C_smilestrikelow.clear(); 
		C_smilestrikehight.clear(); 

		ARM_ERR();
	}


		C_Labels.clear(); 
		C_volcurves.clear(); 
		I_volcurves.clear(); 
		C_proportions.clear(); 
		C_smilestrikelow.clear(); 
		C_smilestrikehight.clear(); 

	} // ARM TRY
	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_CorrelationStrike" )

	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CorrelationSmileStrike(LPXLOPER XL_Labels,
																  LPXLOPER XL_volcurves,
																  LPXLOPER XL_proportions,
																  LPXLOPER XL_AsOf,
																  LPXLOPER XL_smilestrikelow,
																  LPXLOPER XL_smilestrikehight,
																  LPXLOPER XL_Vindex,
																  LPXLOPER XL_NAME,
																  LPXLOPER XL_FullStrikeLow,
																  LPXLOPER XL_FullStrikeUp
																  )
{

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN 
	{
	ARM_NOCALCIFWIZ();

	// C variable
	vector<string> C_Labels; 
	vector<string> C_volcurves; 
	vector<string> C_index; 
	vector<string> C_index_Empty(0); 
	vector<long> I_volcurves; 
	vector<long> I_index; 

	vector<double> C_proportions; 
	vector<double> C_smilestrikelow; 	
	vector<double> C_smilestrikehight; 
	ICM_QMatrix<double> fullStrikeLow,fullStrikeUp,defMatrix; 

	vector<string> DefaultValue;
	DefaultValue.clear();DefaultValue.empty();

	vector<double> DefaultValueDbl;
	DefaultValueDbl.clear();DefaultValueDbl.empty();

	double C_AsOf = 0.;
	// double C_AsOf_default = -1.;
	string C_Name;

	// error
	static int error;
	static char* reason = "";
	
	ExcelTools::convert(XL_Labels, C_Labels);
	ExcelTools::convert(XL_volcurves,C_volcurves); //," ARM_ERR: Labels: object volcurve expected",DOUBLE_TYPE,C_result);
	ExcelTools::convert(XL_Vindex, C_index_Empty,C_index); //," ARM_ERR: Labels: object Credit Index expected",DOUBLE_TYPE,C_result);

	ExcelTools::convert(XL_proportions,C_proportions); //" ARM_ERR: Proportions vector: numerics expected",C_result);
	ExcelTools::convert(XL_smilestrikelow,DefaultValueDbl,C_smilestrikelow);//,DefaultValueDbl," ARM_ERR: smile strike low: numerics expected",C_result);
	ExcelTools::convert(XL_smilestrikehight,DefaultValueDbl,C_smilestrikehight); //,DefaultValueDbl," ARM_ERR: smile strike hight : numerics expected",C_result);

	ExcelTools::convert(XL_AsOf, -1.,C_AsOf); //C_AsOf_default, " ARM_ERR: AsOf : excel date expected",C_result);
	ExcelTools::convert(XL_NAME,"USERDEF",C_Name) ; //,"USERDEF"," ARM_ERR: pricer id: object expected",C_result);
	ExcelTools::convert(XL_FullStrikeLow,defMatrix,fullStrikeLow); 
	ExcelTools::convert(XL_FullStrikeUp,defMatrix,fullStrikeUp); 
	
	CCString prevClass;

	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;
	CCString stringId = GetLastCurCellEnvValue();
	
	if (C_volcurves.size() > 0)	I_volcurves.resize(C_Labels.size());
	if (C_index.size()>0) I_index.resize(C_Labels.size());

	for (int i=0; i<C_Labels.size(); i++)
	{
		I_volcurves[i] = LocalPersistent::get().getObjectId(C_volcurves[i]);
		if (C_index.size()>0) I_index[i] = LocalPersistent::get().getObjectId(C_index[i]);
	}

	long prevId = ExcelCaller::get().getObjectId();
	long retCode = ICMLOCAL_CORRELATION_SMILE_STRIKE(C_AsOf,
											  C_Name,	
											  C_Labels,
											  I_volcurves, 
											  ARM_Vector(C_proportions),
											  ARM_Vector(C_smilestrikelow),
											  ARM_Vector(C_smilestrikehight),
											  fullStrikeLow,
											  fullStrikeUp,
											  I_index,
											  prevId);
	
		string objectLabel = ExcelCaller::get().setObject(retCode, LOCAL_PRICER_CORRMATRIX_CLASS);
		//ExcelCaller::get().setObject(objectLabel);
		ExcelTools::convert(objectLabel, &XL_result);
	} 
	catch (Exception& e)
	{
		ExcelCaller::get().setError(e.GetErrorString());
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (std::exception& e)
	{
		ExcelCaller::get().setError(e.what() );
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (...)
	{
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	return (LPXLOPER)&XL_result;

}

_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_BetaCorrelation(LPXLOPER XL_Labels,
																LPXLOPER XL_betas,
																LPXLOPER XL_AsOf,
																LPXLOPER XL_NAME,
																LPXLOPER XL_index1,
																LPXLOPER XL_index2)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
		 
	ARM_NOCALCIFWIZ();

	// C variable
	std::vector<std::string> C_Labels; 
	VECTOR<double> C_betas; 
	// double C_AsOf = 0.;
	// double C_AsOf_default = -1.;

	// CCString C_Name;
	std::string C_Name ;

	// error
	static int error;
	static char* reason = "";

	ExcelTools::convert(XL_Labels,C_Labels); 
	ExcelTools::convert(XL_NAME,"",C_Name); // Name optional 
	ARM_Date AsOf ; 
	ExcelTools::convert(XL_AsOf,AsOf); 
	ExcelTools::convert(XL_betas,C_betas); 
	
	std::string index1 ; ExcelTools::convert(XL_index1,"",index1); 
	std::string index2 ; ExcelTools::convert(XL_index2,"",index2); 

	long idIndex1 = LocalPersistent::get().getObjectId(index1); 
	long idIndex2 = LocalPersistent::get().getObjectId(index2); 

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;
	CCString stringId = GetLastCurCellEnvValue();
	
	if (!stringId)
	{
		retCode = ICMLOCAL_BETA_CORRELATION(AsOf,
											C_Name,
											C_Labels,
											C_betas,
											idIndex1,
											idIndex2,
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
			
		if (curClass == prevClass)
		{
		retCode = ICMLOCAL_BETA_CORRELATION(AsOf,
											C_Name,
											C_Labels,
											C_betas,
											idIndex1,
											idIndex2,
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

		retCode = ICMLOCAL_BETA_CORRELATION(AsOf,
											C_Name,
											C_Labels,
											C_betas,
											idIndex1,
											idIndex2,
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
		C_Labels.clear(); 
		C_betas.clear(); 

		ARM_ERR();
	}


		C_Labels.clear(); 
		C_betas.clear(); 
	} 
	catch (Exception& e)
	{
		ExcelCaller::get().setError(e.GetErrorString());
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (std::exception& e)
	{
		ExcelCaller::get().setError(e.what() );
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (...)
	{
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	return (LPXLOPER)&XL_result;
	// return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_FlatCorrel(LPXLOPER XL_Asof,
														   LPXLOPER XL_structName,
														   LPXLOPER XL_Value,
														   LPXLOPER XL_index1,
														   LPXLOPER XL_index2)
{
	// return
	static XLOPER XL_result;
	

	try 
	{
		ARM_NOCALCIFWIZ(); 
		std::string structName; ExcelTools::convert(XL_structName,"",structName); 
		ARM_Date AsOf ; ExcelTools::convert(XL_Asof,AsOf); 
		double correlValue ; ExcelTools::convert(XL_Value,correlValue); 
		std::string index1 ; ExcelTools::convert(XL_index1,"",index1); 
		std::string index2 ; ExcelTools::convert(XL_index2,"",index2); 
		long idIndex1 = LocalPersistent::get().getObjectId(index1); 
		long idIndex2 = LocalPersistent::get().getObjectId(index2); 		
		long prevId= ExcelCaller::get().getObjectId() ; 
		long newId = ICMLOCAL_FLAT_CORRELATION(AsOf,structName,correlValue,idIndex1,idIndex2,prevId); 
		std::string newName = ExcelCaller::get().setObject(newId,LOCAL_PRICER_CORRMATRIX_CLASS); 
		ExcelTools::convert(newName,&XL_result); 
	} 
	catch (Exception& e)
	{
		ExcelCaller::get().setError(e.GetErrorString());
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (std::exception& e)
	{
		ExcelCaller::get().setError(e.what() );
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (...)
	{
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	return &XL_result;
	// return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_SetCorrelation (LPXLOPER XL_mmcId,
																 LPXLOPER XL_correlId)
{

	static XLOPER XL_result;
	ARM_result C_result;


	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_mmcId;
	CCString C_correlId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_mmcId,C_mmcId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_correlId,C_correlId," ARM_ERR: Correl Id: correlation object expected",C_result);
	
	long retCode = ICMLOCAL_SetCorrelation (LocalGetNumObjectId (C_mmcId),
											LocalGetNumObjectId (C_correlId),
											C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = 1.;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Price" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GetExpectedLoss (LPXLOPER XL_pricerId,
																  LPXLOPER XL_YearTerm)
{

	static XLOPER XL_result;
	ARM_result C_result;


	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pricerId;

	double YearTerm = 0.;
		// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: security id: object expected",C_result);
	XL_readNumCell(XL_YearTerm,YearTerm," ARM_ERR: YearTerm: year fraction expected",C_result);
	
	long retCode = ICMLOCAL_GetExpectedLoss (LocalGetNumObjectId (C_pricerId),
											 YearTerm,
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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_GetExpectedLoss" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_AddPeriod (LPXLOPER XL_AsOfDate,
															LPXLOPER XL_Maturity,
															LPXLOPER XL_Currency,
															LPXLOPER XL_AdjRule,
															LPXLOPER XL_AdjCDS)
{

	
	static XLOPER XL_result;

	try {
	double C_AsOfDate =0.;
	CCString C_Maturity;

	CCString C_ccy;

	CCString C_AdjRule;
	bool IsAdj =false;

	CCString C_AdjCDS;
	qCDS_ADJ l_AdjCDS = qCredit_Adjust20;
	
	ARM_result C_result;
	ARM_NOCALCIFWIZ();

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_AsOfDate,C_AsOfDate," ARM_ERR: as of date: date expected",C_result);
	XL_readStrCell(XL_Maturity,C_Maturity," ARM_ERR: Maturity expected ",C_result);
	XL_readStrCellWD(XL_Currency,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_AdjRule,C_AdjRule,"N"," ARM_ERR: Adj Rule: string expected",C_result);
	XL_readStrCellWD(XL_AdjCDS,C_AdjCDS,"STDCDS"," ARM_ERR: Adj CDS: string expected",C_result);

	if(C_ccy == "DEFAULT")
	{
		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_ccy = currencyres.getString ();
		}
	}

	if ((C_AdjRule == "Y") || (C_AdjRule == "Adj"))  IsAdj = true; 

	
	if((l_AdjCDS = ARM_ConvAdjCalCDS (C_AdjCDS, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;


	retCode = ICMLOCAL_Credit_AddPeriod (C_AsOfDate,
											 C_Maturity,
											 C_ccy,
											 IsAdj,
											 l_AdjCDS,
											 C_result) ;

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
	
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;

}



_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GetBaseCorrelFromSummit(LPXLOPER XL_AsOf,
																		LPXLOPER XL_proportions,
																		LPXLOPER XL_CurveId,
																		LPXLOPER XL_IndexName,
																		LPXLOPER XL_CCys,
																		LPXLOPER XL_CvIssuerNames,
																		LPXLOPER XL_smilestrikelow,
																		LPXLOPER XL_smilestrikehight,
																		LPXLOPER XL_NAME)
{
	// return
	static XLOPER XL_result;
	try {
	ARM_result C_result;

	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<CCString> C_CCys; 
	VECTOR<CCString> C_CvIssuerNames; 

	VECTOR<double> C_proportions; 
	VECTOR<double> C_smilestrikelow; 	
	VECTOR<double> C_smilestrikehight; 

	vector<CCString> DefaultValueCCys;
	DefaultValueCCys.push_back("EUR");
	DefaultValueCCys.push_back("USD");

	vector<CCString> DefaultValueCvIssuerNames;
	DefaultValueCvIssuerNames.push_back("ITRAXEU3M10");
	DefaultValueCvIssuerNames.push_back("CDXNA4M10");

	vector<CCString> DefaultValue;
	DefaultValue.clear();

	vector<double> DefaultValueDbl;
	DefaultValueDbl.clear();

	// double C_AsOf = 0.;
	// double C_AsOf_default = -1.;
	CCString C_Name;
	CCString C_IndexName;
	CCString C_CurveId;
	// error
	static int error;
	static char* reason = "";

	// XL_readNumCellWD(XL_AsOf, C_AsOf,C_AsOf_default, " ARM_ERR: AsOf : excel date expected",C_result);
	ARM_Date AsOf ; ExcelTools::convert(XL_AsOf,AsOf); 
	XL_readStrCellWD(XL_CurveId,C_CurveId,"MO"," ARM_ERR: Curve: object expected",C_result);
	XL_readStrCellWD(XL_IndexName,C_IndexName,"TRAX"," ARM_ERR: Index Name: object expected",C_result);
	XL_readStrVectorWD(XL_CCys,C_CCys,DefaultValueCCys," ARM_ERR: CCys: array of strings expected",DOUBLE_TYPE,C_result);
	XL_readStrVectorWD(XL_CvIssuerNames,C_CvIssuerNames,DefaultValueCvIssuerNames," ARM_ERR: CrvsNames: array of strings expected",DOUBLE_TYPE,C_result);
	XL_readNumVector(XL_proportions,C_proportions," ARM_ERR: Proportions vector: numerics expected",C_result);
	XL_readNumVectorWD(XL_smilestrikelow,C_smilestrikelow,DefaultValueDbl," ARM_ERR: smile strike low: numerics expected",C_result);
	XL_readNumVectorWD(XL_smilestrikehight,C_smilestrikehight,DefaultValueDbl," ARM_ERR: smile strike hight : numerics expected",C_result);
	XL_readStrCellWD(XL_NAME,C_Name,"USERDEF"," ARM_ERR: pricer id: object expected",C_result);

	
	long retCode = 0;
	long objId = 0;
	CCString prevClass;

	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;
	CCString stringId = GetLastCurCellEnvValue();
	
	if (!stringId)
	{
		retCode = ICMLOCAL_GetBaseCorrelFromSummit(AsOf,
											  C_IndexName,	
											  C_CurveId,
											  C_CCys,
											  C_CvIssuerNames,
											  C_proportions,
											  C_smilestrikelow,
											  C_smilestrikehight,
											  C_Name,
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
			
		if (curClass == prevClass)
		{
		retCode = ICMLOCAL_GetBaseCorrelFromSummit(AsOf,
											  C_IndexName,	
											  C_CurveId,
											  C_CCys,
											  C_CvIssuerNames,
											  C_proportions,
											  C_smilestrikelow,
											  C_smilestrikehight,
											  C_Name,
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

		retCode = ICMLOCAL_GetBaseCorrelFromSummit(AsOf,
											  C_IndexName,	
											  C_CurveId,
											  C_CCys,
											  C_CvIssuerNames,
											  C_proportions,
											  C_smilestrikelow,
											  C_smilestrikehight,
											  C_Name,
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
}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_IndexCorrelation(LPXLOPER XL_AsOf,
																 LPXLOPER XL_NAME,
																 LPXLOPER XL_CalMethod,	
																 LPXLOPER XL_IndexId,	
																 LPXLOPER XL_VStrikeLow,		
																 LPXLOPER XL_VStrikeHigh,		
																 LPXLOPER XL_VMktBid,
																 LPXLOPER XL_VMktAsk,
																 LPXLOPER XL_VUpfBid,
																 LPXLOPER XL_VUpfAsk,
																 LPXLOPER XL_VInitialCorrel,
																 LPXLOPER XL_RFLBeta0,
																 LPXLOPER XL_CreditParameters,
																 LPXLOPER XL_MMCId)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();
 
	// C variable
	VECTOR<double> C_VStrikeLow; 
	VECTOR<double> C_VStrikeHigh; 	
	VECTOR<double> C_VMktBid; 
	VECTOR<double> C_VMktAsk; 
	VECTOR<double> C_VInitialCorrel; 
	VECTOR<double> C_VUpfBid; 
	VECTOR<double> C_VUpfAsk; 

	vector<double> DefaultValueDbl;
	DefaultValueDbl.clear();DefaultValueDbl.empty();

	double C_AsOf = 0.;
	CCString C_Name;
	CCString C_CalMethod;
	CCString C_IndexId;
	int l_CalMethod;
	double C_RFLBeta0 = 0.;
	double C_RFLBeta0_default = -1.;
	CCString C_Parameters;
	CCString C_MMCId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_AsOf, C_AsOf, " ARM_ERR: AsOf : excel date expected",C_result);
	XL_readStrCell(XL_NAME,C_Name," ARM_ERR: name: string expected",C_result);
	XL_readStrCell(XL_CalMethod,C_CalMethod," ARM_ERR: calmethod: string expected",C_result);
	XL_readStrCell(XL_IndexId,C_IndexId," ARM_ERR: IndexId: string expected",C_result);
	XL_readNumVector(XL_VStrikeLow,C_VStrikeLow," ARM_ERR: StrikeLow: array of numerics expected",C_result);
	XL_readNumVector(XL_VStrikeHigh,C_VStrikeHigh," ARM_ERR: StrikeHigh: array of numerics expected",C_result);
	XL_readNumVector(XL_VMktBid,C_VMktBid," ARM_ERR: MktBid: array of numerics expected",C_result);
	XL_readNumVector(XL_VMktAsk,C_VMktAsk," ARM_ERR: MktAsk: array of numerics expected",C_result);
	XL_readNumVector(XL_VUpfBid,C_VUpfBid," ARM_ERR: UpfBid: array of numerics expected",C_result);
	XL_readNumVector(XL_VUpfAsk,C_VUpfAsk," ARM_ERR: UpfAsk: array of numerics expected",C_result);
	XL_readNumVectorWD(XL_VInitialCorrel,C_VInitialCorrel,DefaultValueDbl," ARM_ERR: InitialCorrel: array of numerics expected",C_result);
	XL_readNumCellWD(XL_RFLBeta0, C_RFLBeta0,C_RFLBeta0_default, " ARM_ERR: AsOf : numeric expected",C_result);
	XL_readStrCellWD(XL_CreditParameters,C_Parameters,"NONE"," ARM_ERR: Parameters: string expected",C_result);
	XL_readStrCellWD(XL_MMCId,C_MMCId,"NONE"," ARM_ERR: Model: string expected",C_result);

	long retCode = 0;
	long objId = 0;
	CCString prevClass;

	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;
	CCString stringId = GetLastCurCellEnvValue();
	
	if((l_CalMethod = ARM_ConvCalibType (C_CalMethod, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if (!stringId)
	{
		retCode = ICMLOCAL_INDEX_CORRELATION(C_AsOf,
								C_Name,
								l_CalMethod,
								LocalGetNumObjectId(C_IndexId),
								C_VStrikeLow,
								C_VStrikeHigh,
								C_VMktBid,
								C_VMktAsk,
								C_VUpfBid,
								C_VUpfAsk,
								C_VInitialCorrel,
								C_RFLBeta0,
								LocalGetNumObjectId(C_Parameters),
								LocalGetNumObjectId(C_MMCId),
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
			
		if (curClass == prevClass)
		{
		retCode = ICMLOCAL_INDEX_CORRELATION(C_AsOf,
								C_Name,
								l_CalMethod,
								LocalGetNumObjectId(C_IndexId),
								C_VStrikeLow,
								C_VStrikeHigh,
								C_VMktBid,
								C_VMktAsk,
								C_VUpfBid,
								C_VUpfAsk,
								C_VInitialCorrel,
								C_RFLBeta0,
								LocalGetNumObjectId(C_Parameters),
								LocalGetNumObjectId(C_MMCId),
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

		retCode = ICMLOCAL_INDEX_CORRELATION(C_AsOf,
								C_Name,
								l_CalMethod,
								LocalGetNumObjectId(C_IndexId),
								C_VStrikeLow,
								C_VStrikeHigh,
								C_VMktBid,
								C_VMktAsk,
								C_VUpfBid,
								C_VUpfAsk,
								C_VInitialCorrel,
								C_RFLBeta0,
								LocalGetNumObjectId(C_Parameters),
								LocalGetNumObjectId(C_MMCId),
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
}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CptCorrel(LPXLOPER XL_correlId,
														   LPXLOPER XL_maturity,
														   LPXLOPER XL_Type)
{
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_correlId;

	// error
	static int error;
	static char* reason = "";

	double c_maturity;
	CCString c_type;

	XL_readStrCell(XL_correlId,C_correlId," ARM_ERR: correl id: object expected",C_result);
	XL_readNumCell(XL_maturity,c_maturity," ARM_ERR: maturity : date expected",C_result);
	XL_readStrCell(XL_Type,c_type," ARM_ERR: type : string expected",C_result);
	
	long retCode = ARM_OK; 
	
	if (c_type == "UP")
	{
	retCode = ICMLOCAL_GetCorrelStrikeUp (LocalGetNumObjectId (C_correlId),
			  				  c_maturity,
							  C_result);
	}
	else
	{
	retCode = ICMLOCAL_GetCorrelStrikeDown (LocalGetNumObjectId (C_correlId),
			  				  c_maturity,
							  C_result);
	}

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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_CptCorrel" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CloneCorrFromModel(LPXLOPER XL_Model)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_Model; 

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_Model,C_Model," ARM_ERR: model id: object expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;
	CCString stringId = GetLastCurCellEnvValue();
	
	if (!stringId)
	{
		retCode = ICMLOCAL_GetCorrelation (LocalGetNumObjectId(C_Model),
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
			
		if (curClass == prevClass)
		{
		retCode = ICMLOCAL_GetCorrelation (LocalGetNumObjectId(C_Model),
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

		retCode = ICMLOCAL_GetCorrelation (LocalGetNumObjectId(C_Model),
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
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CptBaseCorrelation(LPXLOPER XL_AsOf,
																 //LPXLOPER XL_NAME,
																 LPXLOPER XL_CreditParameters,
																 LPXLOPER XL_CalMethod,	
																 LPXLOPER XL_IndexId,	
																 LPXLOPER XL_VStrikeLow,		
																 LPXLOPER XL_VStrikeHigh,		
																 LPXLOPER XL_VMktBid,
																 LPXLOPER XL_VMktAsk,
																 LPXLOPER XL_VUpfBid,
																 LPXLOPER XL_VUpfAsk,
																 LPXLOPER XL_VInitialCorrel,
																 LPXLOPER XL_VDeltaLevrage,
																 LPXLOPER XL_ModelId,
																 LPXLOPER XL_Integrationstep,
																 LPXLOPER XL_LagStartDate,
																 LPXLOPER XL_CreditLag,
																 LPXLOPER XL_VectorPrevIndexId,
																 LPXLOPER XL_MatrixPrevBC,
																 LPXLOPER XL_Step,
																 LPXLOPER XL_CalMeth)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;
	LPXLOPER pxArray;

	try {
	ARM_NOCALCIFWIZ();
 
	// C variable
	VECTOR<double> C_VStrikeLow; C_VStrikeLow.clear();C_VStrikeLow.empty();
	VECTOR<double> C_VStrikeHigh; C_VStrikeHigh.clear();C_VStrikeHigh.empty(); 	
	VECTOR<double> C_VMktBid; C_VMktBid.clear();C_VMktBid.empty();
	VECTOR<double> C_VMktAsk; C_VMktAsk.clear();C_VMktAsk.empty();
	VECTOR<double> C_VInitialCorrel; C_VInitialCorrel.clear();C_VInitialCorrel.empty();
	VECTOR<double> C_VUpfBid; C_VUpfBid.clear();C_VUpfBid.empty();
	VECTOR<double> C_VUpfAsk; C_VUpfAsk.clear();C_VUpfAsk.empty();
	VECTOR<double> C_VDeltaLevrage; C_VDeltaLevrage.clear();C_VDeltaLevrage.empty();

	vector<double> DefaultValueDbl;
	DefaultValueDbl.clear();DefaultValueDbl.empty();

	double C_AsOf = 0.;
//	CCString C_Name;
	CCString C_CalMethod;
	CCString C_IndexId;

	CCString C_ModelId;
	int l_CalMethod;

	double C_Integrationstep = 0.;
	double C_Integrationstep_default = 60.;

	double C_LagStartDate = 0.;
	double C_LagStartDate_default = 1.;

	double C_CreditLag = 0.;
	double C_CreditLag_default = 30.;

	VECTOR<CCString> C_PrevIndexId;
	VECTOR<CCString> C_PrevIndexId_default;

	VECTOR<double> C_PrevMatrixBC; C_PrevMatrixBC.clear();C_PrevMatrixBC.empty();
	VECTOR<double> C_PrevMatrixBC_default; C_PrevMatrixBC_default.clear();C_PrevMatrixBC_default.empty();

	CCString C_CalMeth;
	int l_CalMeth;

	CCString C_Parameters;

	// error
	static int error;
	static char* reason = "";

	double C_Step = 0.;
	double C_Step_default = 0.01;

	XL_readNumCell(XL_AsOf, C_AsOf, " ARM_ERR: AsOf : excel date expected",C_result);
//	XL_readStrCell(XL_NAME,C_Name," ARM_ERR: name: string expected",C_result);
	XL_readStrCell(XL_CalMethod,C_CalMethod," ARM_ERR: calmethod: string expected",C_result);
	XL_readStrCell(XL_IndexId,C_IndexId," ARM_ERR: IndexId: string expected",C_result);
	XL_readNumVector(XL_VStrikeLow,C_VStrikeLow," ARM_ERR: StrikeLow: array of numerics expected",C_result);
	XL_readNumVector(XL_VStrikeHigh,C_VStrikeHigh," ARM_ERR: StrikeHigh: array of numerics expected",C_result);
	XL_readNumVector(XL_VMktBid,C_VMktBid," ARM_ERR: MktBid: array of numerics expected",C_result);
	XL_readNumVector(XL_VMktAsk,C_VMktAsk," ARM_ERR: MktAsk: array of numerics expected",C_result);
	XL_readNumVector(XL_VUpfBid,C_VUpfBid," ARM_ERR: UpfBid: array of numerics expected",C_result);
	XL_readNumVector(XL_VUpfAsk,C_VUpfAsk," ARM_ERR: UpfAsk: array of numerics expected",C_result);
	XL_readNumVectorWD(XL_VInitialCorrel,C_VInitialCorrel,DefaultValueDbl," ARM_ERR: InitialCorrel: array of numerics expected",C_result);
	XL_readNumVectorWD(XL_VDeltaLevrage,C_VDeltaLevrage,DefaultValueDbl," ARM_ERR: DeltaLeverage: numeric array awaited",C_result);
	XL_readStrCellWD(XL_ModelId,C_ModelId,"NONE"," ARM_ERR: IndexId: string expected",C_result);
	XL_readNumCellWD(XL_Integrationstep, C_Integrationstep,C_Integrationstep_default, " ARM_ERR: Intstep : numeric expected",C_result);
	XL_readNumCellWD(XL_LagStartDate, C_LagStartDate,C_LagStartDate_default, " ARM_ERR: LagStart : numeric expected",C_result);
	XL_readNumCellWD(XL_CreditLag, C_CreditLag,C_CreditLag_default, " ARM_ERR: creditlag : numeric expected",C_result);
	XL_readStrVectorWD(XL_VectorPrevIndexId,C_PrevIndexId,C_PrevIndexId_default," ARM_ERR: IndexId: string expected",DOUBLE_TYPE,C_result);
	XL_readNumVectorWD(XL_MatrixPrevBC,C_PrevMatrixBC,C_PrevMatrixBC_default," ARM_ERR: Matrix: array of numerics expected",C_result);
	XL_readNumCellWD(XL_Step, C_Step,C_Step_default, " ARM_ERR: step : numeric expected",C_result);
	XL_readStrCellWD(XL_CalMeth,C_CalMeth,"NEWTON"," ARM_ERR: Newton or Dicho: string expected",C_result);
	XL_readStrCellWD(XL_CreditParameters,C_Parameters,"NONE"," ARM_ERR: Parameters: string expected",C_result);

	long retCode = 0;
	long objId = 0;
	CCString prevClass;

	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;
	CCString stringId = GetLastCurCellEnvValue();
	
	if((l_CalMethod = ARM_ConvCalibType (C_CalMethod, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((l_CalMeth = ARM_ConvCalibMeth (C_CalMeth, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	VECTOR<double> basesCorrelation;
	VECTOR<long> V_PrevIndexId;V_PrevIndexId.resize(C_PrevIndexId.size());

	int kk=0;

	for (kk=0; kk<C_PrevIndexId.size();kk++)
	{V_PrevIndexId[kk]=LocalGetNumObjectId(C_PrevIndexId[kk]);}

	retCode = ICMLOCAL_CPT_BASE_CORRELATION(C_AsOf,
								//C_Name,
								l_CalMethod,
								LocalGetNumObjectId(C_IndexId),
								C_VStrikeLow,
								C_VStrikeHigh,
								C_VMktBid,
								C_VMktAsk,
								C_VUpfBid,
								C_VUpfAsk,
								C_VInitialCorrel,
								C_VDeltaLevrage,
								basesCorrelation,
								LocalGetNumObjectId(C_ModelId),
								(int)C_Integrationstep,
								(int)C_LagStartDate,
								(int)C_CreditLag,
								V_PrevIndexId,
								C_PrevMatrixBC,
								C_Step,
								(qOPTIMIZE_TYPE)l_CalMeth,
								LocalGetNumObjectId(C_Parameters),
								C_result);

	if ( retCode == ARM_OK )
	{			
		int nbrows = basesCorrelation.size();
		int nbcolumns = 1;
		
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		for(int i = 0; i < nbrows; i++)
		{
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.num = basesCorrelation[i]; 
		}

	}
	else
	{
		ARM_ERR();
	}
}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}

																 
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CptLeverageLevels(LPXLOPER XL_pricerId,
																	LPXLOPER XL_Trigger_CorrelationId,
																	LPXLOPER XL_Matrix_Multiples,
																	LPXLOPER XL_Matrix_Flags,
																	LPXLOPER XL_CreditParameters
																	)

{
	static XLOPER XL_result;
	ARM_result C_result;
	LPXLOPER pxArray;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString	C_pricerId;
		CCString	C_CreditParameters;
		CCString	C_Trigger_CorrelationId;

		VECTOR<double> C_Matrix_Multiples;
		VECTOR<double> C_Matrix_Flags;

		ICM_QMatrix<double>* MatrixMultiples;
		vector<double>		VectOfMaturitiesInYF;
		vector<double>		VectOfLosses;

		double YearTerm = 0.;
			// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: security id: object expected",C_result);

		// Credit Manager
		XL_readStrCell(XL_Trigger_CorrelationId, C_Trigger_CorrelationId," ARM_ERR: Trigger Correlation Id: object expected",C_result);

		// Cash Flow Id
		XL_readStrCellWD(XL_CreditParameters, C_CreditParameters,"NONE"," ARM_ERR: Credit Parameters: Object expected",C_result);
		
		// Spreads	
		XL_readNumVector(XL_Matrix_Multiples, C_Matrix_Multiples," ARM_ERR: Matrix_Multiples: matrix of numerics expected",C_result);

		// Spreads	
		XL_readNumVector(XL_Matrix_Flags, C_Matrix_Flags," ARM_ERR: Matrix_Flags: matrix of numerics expected",C_result);

		long retCode = ICMLOCAL_CptLeverageLevels(LocalGetNumObjectId (C_pricerId),
												 LocalGetNumObjectId (C_CreditParameters),
												 LocalGetNumObjectId (C_Trigger_CorrelationId),
												 C_Matrix_Multiples,
												C_Matrix_Flags,
												 VectOfLosses,
												 VectOfMaturitiesInYF,
												 MatrixMultiples,
												 C_result);

		if (retCode == ARM_OK)
		{	
			int	i,j;
			int nbrows		=	(*MatrixMultiples).Getnbrows() + 1;		// Number of Rows + Label
			int nbcolumns	=	(*MatrixMultiples).Getnbcols() + 1;		// Number of Cols + Label
			
			FreeCurCellErr();

			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = nbcolumns;
			XL_result.val.array.rows	= nbrows; 
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

			// Blank, high left corner
			pxArray[XL_Coordonnate2Rank (0,0, nbcolumns)].xltype	=	xltypeStr;
			pxArray[XL_Coordonnate2Rank (0,0, nbcolumns)].val.str	=	XL_StrC2StrPascal ("");
					
			// Multiples Matrice values
			for (i=1; i<nbrows; i++)
			{
				// label for rows
				pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype	=	xltypeNum;
				pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.num	=	VectOfLosses[i-1];

				for (j=1; j<nbcolumns; j++)
				{
					pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype	=	xltypeNum;
					pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.num	=	(*MatrixMultiples)(i-1,j-1); 
				}
			}
		
			for (j=1; j<nbcolumns; j++)
			{
				// En tete du vecteur proba de defaut
				pxArray[XL_Coordonnate2Rank (0, j, nbcolumns)].xltype	=	xltypeNum;
				pxArray[XL_Coordonnate2Rank (0, j, nbcolumns)].val.num	=	VectOfMaturitiesInYF[j-1];
			}
			
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_CptLeverageLevels" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_SetRecovCoef(LPXLOPER XL_SecId,
															LPXLOPER XL_RecovCoef
															)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;


	try {
		ARM_NOCALCIFWIZ();

	// C variable
	CCString C_SecId;
	double C_RecovCoef;
	    
	long retCode = 0;
		
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_SecId,C_SecId," ARM_ERR: Mezzanine object expected",C_result);
	XL_readNumCell(XL_RecovCoef, C_RecovCoef,  " ARM_ERR: Recovery coefficient: numeric expected",C_result);
	
	retCode = ICMLOCAL_SetRecovCoef (LocalGetNumObjectId (C_SecId), 
								   C_RecovCoef,
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

	/// return the result as an LPXLOPER
	}
		catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}



extern void DeduceRefDateAndStub(ARM_Date StartDate,
						  ARM_Date EndDate,
						  ARM_Date StartStubDate,
						  ARM_Date EndStubDate,
						  const int& payFreq,
						  const int & RefDay,
						  const std::string& ccy,
						  char* ReferenceDate,
						  int& stub) ;

__declspec(dllexport) LPXLOPER WINAPI 
ARM_Credit_DeduceRefDateAndStub(LPXLOPER XL_StartDate,
								LPXLOPER XL_EndDate,
								LPXLOPER XL_StartStubDate,
								LPXLOPER XL_EndStubDate,
								LPXLOPER XL_PayFreq,
								LPXLOPER XL_RefDay,
								LPXLOPER XL_Ccy,
								LPXLOPER XL_RefDate,
								LPXLOPER XL_Stub)

{

	static XLOPER XL_result;
	static int error; 
	static char* reason="" ;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		double 	C_StartDate;
		double 	C_EndDate;
		double 	C_StartStubDate;
		double 	C_EndStubDate;
		CCString	C_PayFreq;
		CCString	C_Ccy;
		CCString 	C_RefDate; 
		CCString	C_Stub ; 
		double 	C_RefDay; 

		double  defaultDay=-1; 
		XL_readNumCell(XL_StartDate,C_StartDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_EndDate,C_EndDate," ARM_ERR: end date: date expected",C_result);
		XL_readNumCellWD(XL_StartStubDate,C_StartStubDate,defaultDay," ARM_ERR: start stub date: date expected",C_result);
		XL_readNumCellWD(XL_EndStubDate,C_EndStubDate,defaultDay," ARM_ERR: end stub date: date expected",C_result);
		XL_readStrCellWD(XL_RefDate,C_RefDate,""," ARM_ERR: ref date: date expected",C_result);
		XL_readStrCellWD(XL_PayFreq,C_PayFreq,"Q"," ARM_ERR: PayFrequency: string expected",C_result);
		XL_readStrCellWD(XL_Ccy,C_Ccy,"EUR"," ARM_ERR: currency: string expected",C_result);
		XL_readStrCellWD(XL_Stub,C_Stub,"SS"," ARM_ERR: stub rule: string expected",C_result);
		XL_readNumCell(XL_RefDay,C_RefDay, " ARM_ERR: stub rule: string expected",C_result);
	

		ARM_result C_result ;
		long payFreq ;
		int  stubRuleId;

		if((payFreq = ARM_ConvFrequency (C_PayFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		if((stubRuleId = ARM_ConvStubRule (C_Stub)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		// ARM_Currency ccy((char*)C_Ccy); 

		ARM_Date StartDate(C_StartDate); 
		ARM_Date EndDate(C_EndDate); 
		ARM_Date StartStubDate ;
		ARM_Date EndStubDate ;
		char ReferenceDate[11];

		if (C_StartStubDate!=-1)  Local_XLDATE2ARMDATE(C_StartStubDate,StartStubDate);
		if (C_EndStubDate!=-1) Local_XLDATE2ARMDATE(C_EndStubDate,EndStubDate);
		// if (C_RefDate!=-1) 
		// {
		// 	ARM_Date toto ;
		// 	Local_XLDATE2ARMDATE(C_RefDate,toto);
		// 	toto.JulianToStrDate(ReferenceDate)  ;
		// }
		else strcpy(ReferenceDate,"NULL"); 
		 DeduceRefDateAndStub( StartDate,
						   EndDate,
						   StartStubDate,
						   EndStubDate,
							payFreq,
						  C_RefDay,
						  CCSTringToSTLString(C_Ccy),
						  ReferenceDate,
						  stubRuleId) ;
		
		LPXLOPER pxArray;
		unsigned int nbcolumns = 2 ; 
		unsigned int nbrows = 1 ; 
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		// 
		pxArray[XL_Coordonnate2Rank( 0, 0, nbcolumns)].xltype = xltypeStr | xlbitDLLFree;;
		switch (stubRuleId) 
		{
		case 1 : pxArray[XL_Coordonnate2Rank( 0, 0, nbcolumns)].val.str = XL_StrC2StrPascal("SS"); break ;
		case 2 : pxArray[XL_Coordonnate2Rank( 0, 0, nbcolumns)].val.str = XL_StrC2StrPascal("LS"); break ;
		case 3 : pxArray[XL_Coordonnate2Rank( 0, 0, nbcolumns)].val.str = XL_StrC2StrPascal("SE"); break ;
		case 4 : pxArray[XL_Coordonnate2Rank( 0, 0, nbcolumns)].val.str = XL_StrC2StrPascal("LE"); break ;
		}
		pxArray[XL_Coordonnate2Rank( 0, 1, nbcolumns)].xltype = xltypeStr | xlbitDLLFree;
		pxArray[XL_Coordonnate2Rank( 0, 1, nbcolumns)].val.str = XL_StrC2StrPascal (ReferenceDate);

/*		for(int i = 0; i < nbrows; i++)
		{
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.str = XL_StrC2StrPascal (Matu[i]);
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].val.num = Spread[i]; 
			pxArray[XL_Coordonnate2Rank (i, 2, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (i, 2, nbcolumns)].val.num = Recovery[i]; 
			pxArray[XL_Coordonnate2Rank (i, 3, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (i, 3, nbcolumns)].val.str = XL_StrC2StrPascal (C_Ccy);
			pxArray[XL_Coordonnate2Rank (i, 4, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (i, 4, nbcolumns)].val.str = XL_StrC2StrPascal (C_IndexName);
		}
*/ 		
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_DeduceRefDateAndStub" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_SetInterpolationType(LPXLOPER XL_VolCurveId,
																	  LPXLOPER XL_InterpolType
																	  )

{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;
	try {
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_VolCurveId;

	CCString C_InterpolType;
	CCString C_InterpolType_Default = "LINEAR";

	long InterpolType = 0;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_VolCurveId,C_VolCurveId," ARM_ERR: VolCurve id: object expected",C_result);
	XL_readStrCellWD(XL_InterpolType,C_InterpolType,C_InterpolType_Default," ARM_ERR: Interpolation Type : string expected",C_result);

	if((InterpolType = ARM_ConvInterpMethod (C_InterpolType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	long retCode = ICMLOCAL_SetInterpolationType (LocalGetNumObjectId (C_VolCurveId), 
												  InterpolType,
												  C_result);

	C_result.setDouble(1.);

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
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Sectorial_Correlation(LPXLOPER XL_AsOf,
																	  LPXLOPER XL_StructName,
																		LPXLOPER XL_Sectorial_Correlation_Type,
																		LPXLOPER XL_Labels,
																		LPXLOPER XL_Sector_Membership,
																		LPXLOPER XL_Intra_Sector_Correlation,
																		LPXLOPER XL_Inter_Sector_Correlation,
																		LPXLOPER XL_Sector_Betas,
																		LPXLOPER XL_Sector_Lambdas)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable
	// CCString	C_Sectorial_Correlation_Type;
	// CCString	C_Sectorial_Correlation_Type_Default = "LINEAR";

	// VECTOR<CCString> C_Labels;
	std::vector<std::string> C_Labels; 

	VECTOR<double>	C_Sector_Betas;
	VECTOR<double>	C_Sector_Lambdas;

	VECTOR<double>	C_Sector_Betas_Down;
	VECTOR<double>	C_Sector_Lambdas_Down;

	VECTOR<double>	C_Sector_Membership_db;
	VECTOR<int>		C_Sector_Membership;

    double	C_Intra_Sector_Correlation =0.;
	double	C_Intra_Sector_Correlation_Default=0.;
    double	C_Inter_Sector_Correlation=0.;
    double	C_Inter_Sector_Correlation_Default=0.;

	// long	Sectorial_Correlation_Type = 0;

	vector<double> DefaultValueDbl;
	DefaultValueDbl.clear();DefaultValueDbl.empty();

	CCString C_Name;

	// error
	static int error =0;
	static char* reason = "";

	// XL_readStrCellWD(XL_Sectorial_Correlation_Type, C_Sectorial_Correlation_Type, C_Sectorial_Correlation_Type_Default," ARM_ERR: Sectorial Correlation Type : string expected",C_result);

	// XL_readStrVector(XL_Labels, C_Labels, " ARM_ERR: Labels: array of strings expected",DOUBLE_TYPE, C_result);
	ExcelTools::convert(XL_Labels,C_Labels) ;
	ARM_Date AsOf; ExcelTools::convert(XL_AsOf,AsOf) ;
	std::string StructName; ExcelTools::convert(XL_StructName,StructName) ;
	XL_readNumVector(XL_Sector_Membership, C_Sector_Membership_db," ARM_ERR: Sector Membership: double vector expected", C_result);

	XL_readNumCellWD(XL_Intra_Sector_Correlation, C_Intra_Sector_Correlation, C_Intra_Sector_Correlation_Default, " ARM_ERR: Intra_Sector_Correlation_Default : numeric expected",C_result);
	XL_readNumCellWD(XL_Inter_Sector_Correlation, C_Inter_Sector_Correlation, C_Inter_Sector_Correlation_Default, " ARM_ERR: Intra_Sector_Correlation_Default : numeric expected",C_result);

	XL_readNumVectorWD(XL_Sector_Betas, C_Sector_Betas, DefaultValueDbl, " ARM_ERR: Sector Betas: double vector expected", C_result);
	XL_readNumVectorWD(XL_Sector_Lambdas, C_Sector_Lambdas, DefaultValueDbl, " ARM_ERR: Sector Lambdas: double vector expected", C_result);
	

	qTWO_FACTORS_CORRELATION_TYPE Sectorial_Correlation_Type ; 
	ExcelTools::econvert(XL_Sectorial_Correlation_Type,"LINEAR",Sectorial_Correlation_Type); 

	// CHECK Sectorial Correlation Type
	// if((Sectorial_Correlation_Type = ARM_Conv_Sectorial_Correlation(C_Sectorial_Correlation_Type, C_result)) == ARM_DEFAULT_ERR)
	// {
	// 	ARM_ARG_ERR();
	// 	return (LPXLOPER)&XL_result;
	// }
	
	int	size =0;
	size	=	C_Sector_Membership_db.size();

	C_Sector_Membership.resize(size);

	for (int i=0; i<size; i++)
		C_Sector_Membership[i]	=	(int)	C_Sector_Membership_db[i];

	long retCode =0;
	long objId =0 ;
	CCString prevClass;

	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;
	CCString stringId = GetLastCurCellEnvValue();
	
	if (!stringId)
	{
		retCode = ICMLOCAL_SECTORIAL_CORRELATION(AsOf,
			StructName,
				Sectorial_Correlation_Type,
											C_Labels,
											C_Sector_Membership,
											C_Intra_Sector_Correlation,
											C_Inter_Sector_Correlation,
											C_Sector_Betas,
											C_Sector_Lambdas,
											C_Sector_Betas_Down,
											C_Sector_Lambdas_Down,
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
			
		if (curClass == prevClass)
		{
		retCode = ICMLOCAL_SECTORIAL_CORRELATION(AsOf,
			StructName,
											Sectorial_Correlation_Type,
											C_Labels,
											C_Sector_Membership,
											C_Intra_Sector_Correlation,
											C_Inter_Sector_Correlation,
											C_Sector_Betas,
											C_Sector_Lambdas,
											C_Sector_Betas_Down,
											C_Sector_Lambdas_Down,
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

		retCode = ICMLOCAL_SECTORIAL_CORRELATION(AsOf,
			StructName,
											Sectorial_Correlation_Type,
											C_Labels,
											C_Sector_Membership,
											C_Intra_Sector_Correlation,
											C_Inter_Sector_Correlation,
											C_Sector_Betas,
											C_Sector_Lambdas,
											C_Sector_Betas_Down,
											C_Sector_Lambdas_Down,
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
}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GenTsCorrfromBaseCorr(LPXLOPER XL_irCurveId,
																	LPXLOPER XL_defCurveId,
																	LPXLOPER XL_volcurveId,
																	LPXLOPER XL_creditlag)

{
	static XLOPER XL_result;
	ARM_result C_result;
	LPXLOPER pxArray;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString	C_irCurveId;
		CCString	C_defCurveId;
		CCString	C_volcurveId;

		ICM_QMatrix<double>* C_Matrix;

		double C_creditlag = 0.;

			// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_irCurveId,C_irCurveId," ARM_ERR: ircurve id: object expected",C_result);
		XL_readStrCell(XL_defCurveId, C_defCurveId," ARM_ERR: defcurve Id: object expected",C_result);
		XL_readStrCell(XL_volcurveId, C_volcurveId," ARM_ERR: volcurve: Object expected",C_result);
		XL_readNumCell(XL_creditlag, C_creditlag," ARM_ERR: creditlag: numeric expected",C_result);


		long retCode = ICMLOCAL_GenTsCorrfromBaseCorr (LocalGetNumObjectId(C_irCurveId), 
											LocalGetNumObjectId(C_defCurveId), 
											LocalGetNumObjectId(C_volcurveId), 
											(int)C_creditlag,
											C_Matrix,
											C_result);

		if (retCode == ARM_OK)
		{	
			int	i,j;
			int nbrows		=	(*C_Matrix).Getnbrows();		// Number of Rows + Label
			int nbcolumns	=	(*C_Matrix).Getnbcols();		// Number of Cols + Label
			
			FreeCurCellErr();

			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = nbcolumns;
			XL_result.val.array.rows	= nbrows; 
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

			// Multiples Matrice values
			for (i=0; i<nbrows; i++)
			{
				for (j=0; j<nbcolumns; j++)
				{
					pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype	=	xltypeNum;
					pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.num	=	(*C_Matrix)(i,j); 
				}
			}
			
			if (C_Matrix) delete C_Matrix; C_Matrix=NULL;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_CptLeverageLevels" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DefProbInverse (LPXLOPER XL_defcurveId, 
														 LPXLOPER XL_defproba)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;


	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable

	double C_defprob = 0.0;
	CCString C_defcurveId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_defcurveId,C_defcurveId," ARM_ERR: curve id: object expected",C_result);
	XL_readNumCell(XL_defproba, C_defprob, " ARM_ERR: def proba : numeric expected",C_result);
	
	long retCode = ICMLOCAL_DefProbInverse (LocalGetNumObjectId (C_defcurveId), 
																C_defprob, 
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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_DefProbInverse" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Get_IR_Curve_Moved_In_Time (LPXLOPER XL_DiscountCurveId,
																			LPXLOPER XL_Date
														   )
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;


	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_DiscountCurveId =0.;
    double C_Date = 0.;
	double C_Date_default = 0.;
	
	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_DiscountCurveId, C_DiscountCurveId ," ARM_ERR: Discount Curve : object expected",C_result);
	XL_readNumCellWD(XL_Date, C_Date, C_Date_default," ARM_ERR: numeric expected",C_result);
	
	long retCode;
	long objId;
	CCString prevClass;
	
    CCString curClass = LOCAL_YIELD_CURVE_BASIC_CLASS;
    CCString stringId = GetLastCurCellEnvValue();
/*
    if (!stringId)
    {
        retCode = ICMLOCAL_Get_IR_Curve_Moved_In_Time (LocalGetNumObjectId (C_DiscountCurveId), 
																C_Date, 
																C_result);

        if (retCode == ARM_OK)
        {
            objId = C_result.getLong ();

            LocalSetCurCellEnvValue (curClass, objId); 

            stringId = LocalMakeObjectId (objId, curClass);
        }
    }
    else
    {
        prevClass = LocalGetStringObjectClass(stringId);
        
        objId = LocalGetNumObjectId(stringId);
            
        if ( curClass == prevClass )
        {
			retCode = ICMLOCAL_Get_IR_Curve_Moved_In_Time (LocalGetNumObjectId (C_DiscountCurveId), 
																C_Date, 
																C_result);

			if ( retCode == ARM_OK )
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
        else
        {

			FreeCurCellContent();

			retCode = ICMLOCAL_Get_IR_Curve_Moved_In_Time (LocalGetNumObjectId (C_DiscountCurveId), 
																C_Date, 
																C_result,
																objId);
        
			if ( retCode == ARM_OK )
			{
				objId = C_result.getLong ();

				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
        }
    }
*/
    if (!stringId)
    {
		;
	}
	else
	{
		FreeCurCellContent();
	}

    retCode = ICMLOCAL_Get_IR_Curve_Moved_In_Time (LocalGetNumObjectId (C_DiscountCurveId), 
															C_Date, 
															C_result);

    if (retCode == ARM_OK)
    {
        objId = C_result.getLong ();

        LocalSetCurCellEnvValue (curClass, objId); 

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

//    ARM_END();
    
    return (LPXLOPER)&XL_result;
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Get_IR_Curve_Moved_In_Time" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Math_BivNormale (LPXLOPER XL_x, 
																LPXLOPER XL_y,
																LPXLOPER XL_rho)
{

	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_x;
	double C_y;
	double C_rho;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_x, C_x,  "ARM_ERR: X: numeric expected",C_result);
	XL_readNumCell(XL_y, C_y,  "ARM_ERR: Y: numeric expected",C_result);
	XL_readNumCell(XL_rho, C_rho,  "ARM_ERR: rho: numeric expected",C_result);

	long retCode = ICMLOCAL_Math_Bivariate_normale (C_x,C_y,C_rho,C_result);

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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Spread" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Math_RandUniform (LPXLOPER XL_seed)
{

	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_seed;
	double C_seed_default=-1.;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCellWD(XL_seed, C_seed,C_seed_default,  "ARM_ERR: X: numeric expected",C_result);

	long retCode = ICMLOCAL_Math_random_uniform (C_seed,C_result);

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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Spread" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Math_RandNormal (LPXLOPER XL_a, 
																	LPXLOPER XL_b,
																   LPXLOPER XL_seed)
{

	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_seed;
	double C_seed_default=-1.;

	double C_a;
	double C_b;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_a, C_a, "ARM_ERR: a: numeric expected",C_result);
	XL_readNumCell(XL_b, C_b, "ARM_ERR: b: numeric expected",C_result);
	XL_readNumCellWD(XL_seed, C_seed,C_seed_default,  "ARM_ERR: X: numeric expected",C_result);

	long retCode = ICMLOCAL_Math_random_normal (C_a,C_b,C_seed,C_result);

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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Spread" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_QMatrix(LPXLOPER XL_QMatrix)
{
	// return
	static XLOPER XL_result;

	/// to make the infrastructure very robust we put a try catch!
	try
	{
		long prevId = ExcelCaller::get().getObjectId();
		ICM_QMatrix<double>  lQMatrix;
		ExcelTools::convert(XL_QMatrix,lQMatrix) ; // will throw if empty

		long newId = ICMLOCAL_QMatrix(lQMatrix,prevId);
		
		string objectLabel = ExcelCaller::get().setObject(newId, LOCAL_CREDIT_UTIL);
		//ExcelCaller::get().setObject(objectLabel);
		ExcelTools::convert(objectLabel, &XL_result);
	} 
	catch (Exception& e)
	{
		ExcelCaller::get().setError(e.GetErrorString());
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (std::exception& e)
	{
		ExcelCaller::get().setError(e.what() );
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (...)
	{
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_QuickELoss(LPXLOPER pdefault,
															LPXLOPER recovery,
															LPXLOPER nbnames,
															LPXLOPER strikedw,
															LPXLOPER strikeup,
															LPXLOPER correldw,
															LPXLOPER correlup,
															LPXLOPER lossesno,
															LPXLOPER intstep,
															LPXLOPER lhp)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;


	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double c_pdefault;
	double c_recovery;
	double c_nbnames;
	double c_strikedw;
	double c_strikeup;
	double c_correldw;
	double c_correlup;
	double c_lossesno;
	double c_lossesno_default=-1;
	double c_intstep;
	double c_intstep_default=40;
	CCString C_Lhp;
	CCString C_Lhp_default = "N";

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(pdefault, c_pdefault, " ARM_ERR: pdef : numeric expected",C_result);
	XL_readNumCell(recovery, c_recovery, " ARM_ERR: recovery : numeric expected",C_result);
	XL_readNumCell(nbnames, c_nbnames, " ARM_ERR: nbnames : numeric expected",C_result);
	XL_readNumCell(strikedw, c_strikedw, " ARM_ERR: strikedw : numeric expected",C_result);
	XL_readNumCell(strikeup, c_strikeup, " ARM_ERR: strikeup : numeric expected",C_result);
	XL_readNumCell(correldw, c_correldw, " ARM_ERR: correldw : numeric expected",C_result);
	XL_readNumCell(correlup, c_correlup, " ARM_ERR: correlup : numeric expected",C_result);
	XL_readNumCellWD(lossesno, c_lossesno,c_lossesno_default, " ARM_ERR: lossesno : numeric expected",C_result);
	XL_readNumCellWD(intstep, c_intstep,c_intstep_default, " ARM_ERR: intstep : numeric expected",C_result);
	XL_readStrCellWD(lhp, C_Lhp,C_Lhp_default, " ARM_ERR: lhp : string expected",C_result);
	

	bool b_lhp=false;
	if (C_Lhp=="Y") {b_lhp = true;}

	long retCode = ICMLOCAL_QuickELoss (c_pdefault,c_recovery,
							   (int) c_nbnames,c_strikedw,c_strikeup,c_correldw,
							    c_correlup,(int) c_lossesno,(int) c_intstep, b_lhp, C_result);

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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_DefaultProba" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_QuickCDOPV(LPXLOPER i_discid,
															LPXLOPER i_pdefid,
															LPXLOPER i_startdate,
															LPXLOPER i_enddate,
															LPXLOPER i_frequency,
															LPXLOPER i_rate,
															LPXLOPER i_strikedw,
															LPXLOPER i_strikeup,
															LPXLOPER i_correldw,
															LPXLOPER i_correlup,
															LPXLOPER i_notfeeleg,
															LPXLOPER i_notdefleg,
															LPXLOPER i_recovery,
															LPXLOPER i_nbnames,
															LPXLOPER i_intstep,
															LPXLOPER i_LHP,
															LPXLOPER i_pvtype,
															LPXLOPER i_volcurveTS)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;


	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_discid;
	CCString C_volcurveTS;
	CCString C_pdefid;
	double startdate;
	double enddate;
	double frequency;
	double rate;
	double strikedw;
	double strikeup;
	double correldw;
	double correlup;
	double notfeeleg;
	double notdefleg;
	double recovery;
	double nbnames;
	double intstep;
	double intstep_default=40;
	CCString LHP;
	CCString LHP_default = "N";
	CCString pvtype;
	CCString pvtype_default = "NPV";
	long l_pvtype;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(i_discid, C_discid, " ARM_ERR: ircurve : string expected",C_result);
	XL_readStrCell(i_pdefid, C_pdefid, " ARM_ERR: defcure : string expected",C_result);
	XL_readNumCell(i_startdate, startdate, " ARM_ERR: startdate : numeric expected",C_result);
	XL_readNumCell(i_enddate, enddate, " ARM_ERR: enddate : numeric expected",C_result);
	XL_readNumCell(i_frequency, frequency, " ARM_ERR: frequency : numeric expected",C_result);
	XL_readNumCell(i_rate, rate, " ARM_ERR: rate : numeric expected",C_result);
	XL_readNumCell(i_strikedw, strikedw, " ARM_ERR: strikedw : numeric expected",C_result);
	XL_readNumCell(i_strikeup, strikeup, " ARM_ERR: strikeup : numeric expected",C_result);
	XL_readNumCell(i_correldw, correldw, " ARM_ERR: correldw : numeric expected",C_result);
	XL_readNumCell(i_correlup, correlup, " ARM_ERR: correlup : numeric expected",C_result);
	XL_readNumCell(i_notfeeleg, notfeeleg, " ARM_ERR: notfeeleg : numeric expected",C_result);
	XL_readNumCell(i_notdefleg, notdefleg, " ARM_ERR: notdefleg : numeric expected",C_result);
	XL_readNumCell(i_recovery, recovery, " ARM_ERR: recovery : numeric expected",C_result);
	XL_readNumCell(i_nbnames, nbnames, " ARM_ERR: nbnames : numeric expected",C_result);
	XL_readNumCellWD(i_intstep, intstep,intstep_default, " ARM_ERR: intstep : numeric expected",C_result);
	XL_readStrCellWD(i_LHP, LHP,LHP_default, " ARM_ERR: LHP mode : string expected",C_result);
	XL_readStrCellWD(i_pvtype, pvtype,pvtype_default, " ARM_ERR: pricing type : string expected",C_result);
	XL_readStrCellWD(i_volcurveTS, C_volcurveTS,"NONE", " ARM_ERR: correl : string expected",C_result);

	bool lhp = false;
	if (LHP=="Y") {lhp = true;}

	if((l_pvtype = ICM_ConvCptType (pvtype, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ICMLOCAL_QuickCDOPV( LocalGetNumObjectId(C_discid), 
									    LocalGetNumObjectId(C_pdefid), 
										startdate,
										enddate,
										frequency,
										rate,
										strikedw,
										strikeup,
										correldw,
										correlup,
										notfeeleg,
										notdefleg,
										recovery,
										nbnames,
										intstep,
										lhp,
										l_pvtype,
										LocalGetNumObjectId(C_volcurveTS), 
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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_DefaultProba" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Math_Interpol (LPXLOPER XL_X,
																LPXLOPER XL_Y,
																LPXLOPER XL_VALUE,
																LPXLOPER XL_type,
																LPXLOPER XL_smooth,
																LPXLOPER XL_weights,
																LPXLOPER XL_modespline,
																LPXLOPER XL_withC1condition,
																LPXLOPER XL_leftSlope,
																LPXLOPER XL_rightSlope)
{

	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_x;
	VECTOR<double> C_y;
	double value=0.;
	double type=0.;
	double type_default=0.;
	double smooth=0.;
	double smooth_default=-1.;
	VECTOR<double> C_weights;
	VECTOR<double> C_weights_default;

	double C_modespline=0.;
	double C_withC1condition=0.;
	double C_leftSlope=0.;
	double C_rightSlope=0.;

	double C_modespline_default=1.;
	double C_withC1condition_default=0.;
	double C_leftSlope_default=0.;
	double C_rightSlope_default=0.;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumVector(XL_X, C_x," ARM_ERR: X: double vector expected", C_result);
	XL_readNumVector(XL_Y, C_y," ARM_ERR: Y: double vector expected", C_result);
	XL_readNumCell(XL_VALUE, value, " ARM_ERR: value : numeric expected",C_result);
	XL_readNumCellWD(XL_type, type,type_default, " ARM_ERR: type : numeric expected",C_result);
	XL_readNumCellWD(XL_smooth, smooth,smooth_default, " ARM_ERR: type : numeric expected",C_result);
	XL_readNumVectorWD(XL_weights, C_weights, C_weights_default, " ARM_ERR: weights: double vector expected", C_result);
	XL_readNumCellWD(XL_modespline, C_modespline,C_modespline_default, " ARM_ERR: modespline : numeric expected",C_result);
	XL_readNumCellWD(XL_withC1condition, C_withC1condition,C_withC1condition_default, " ARM_ERR: c1 condition : numeric expected",C_result);
	XL_readNumCellWD(XL_leftSlope, C_leftSlope,C_leftSlope_default, " ARM_ERR: leftSlope : numeric expected",C_result);
	XL_readNumCellWD(XL_rightSlope, C_rightSlope,C_rightSlope_default, " ARM_ERR: rightSlope : numeric expected",C_result);

	long retCode = ICMLOCAL_Math_Interpol(C_x,C_y,value,(int) type,smooth,C_weights,(int)C_modespline,(int)C_withC1condition,
		C_leftSlope,C_rightSlope,C_result);

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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Spread" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Schedule_Info (LPXLOPER 	XL_EffectiveDate,
																LPXLOPER 	XL_MaturityDate,			
																LPXLOPER 	XL_payFrequency,
																LPXLOPER 	XL_ResetFreq ,
																LPXLOPER 	XL_DayCount,
																LPXLOPER 	XL_Stubrule,
																LPXLOPER 	XL_intRule,
																LPXLOPER 	XL_payCalName,
																LPXLOPER 	XL_PayTiming,
																LPXLOPER 	XL_ResetTiming,
																LPXLOPER 	XL_fwdRule,
																LPXLOPER 	XL_IncludeMaturity,
																LPXLOPER 	XL_adj,
																LPXLOPER 	XL_intStartAdj,
																LPXLOPER 	XL_AccDayCount,
																LPXLOPER 	XL_ReferenceDate,
																LPXLOPER 	XL_FirstCpnEffDate,
																LPXLOPER 	XL_AdjCal)
{
	// return
    static XLOPER XL_result;
	ARM_result C_result;
	try
	{
		 long prevId = ExcelCaller::get().getObjectId();
		 
		 double	EffectiveDate = -1;
		 double MaturityDate = -1;
		 string strPayFreq = "";
		 int	payFrequency = 0;
		 string strResetFreq = "";
		 int ResetFreq =0;
		 string Basis = "";
		 int	DayCount =0;
		 std::string strStubrule = "";
		 int	Stubrule =0;
		 string strIntRule = "";
		 int	intrule =0;
		 std::string payCalName = "";

		 string strPayTiming ="";
		 int PayTiming =0;
		 string strResetTiming ="";
		 int ResetTiming =0 ;
		 string strfwdRule = "";
		 int fwdRule=0;
		 string strIncludeMaturity = "N";
		 bool	IncludeMaturity = false;
		 string strAdj = "";
		 int Adj =0;
		 string strStartAdj ="";
		 int	intStartAdj =0;
		 string strAccDayCount = "";
		 int AccDayCount =0;
		 string strReferenceDate="";
		 double ReferenceDate = -1;
		 string strFirstCpnEffDate = "";
		 double FirstCpnEffDate = -1;
		 string strAdjCal="";
		 int AdjCal = 0;
		 
		 ExcelTools::convert( XL_EffectiveDate,EffectiveDate);
		 // if(EffectiveDate == -1) { ARM_ARG_ERR();	return (LPXLOPER)&XL_result;}
		 ExcelTools::convert( XL_MaturityDate,MaturityDate);
		 // if(MaturityDate == -1) { ARM_ARG_ERR();	return (LPXLOPER)&XL_result;}
		 ExcelTools::convert( XL_payFrequency,strPayFreq);	
			if (strPayFreq.empty() ) strPayFreq = "Q";
			if (strPayFreq== "-1" )
			   payFrequency = K_DEF_FREQ;
			else
				if ((payFrequency = ARM_ConvFrequency(strPayFreq.c_str(), C_result)) == ARM_DEFAULT_ERR)	{	ARM_ARG_ERR();	return (LPXLOPER)&XL_result;	}
			
		ExcelTools::convert( XL_ResetFreq , "M",strResetFreq);
		//	if (strResetFreq.empty() ) strResetFreq = "M";
			if((ResetFreq = ARM_ConvFrequency (strResetFreq.c_str(), C_result)) == ARM_DEFAULT_ERR)	{	ARM_ARG_ERR();	return (LPXLOPER)&XL_result;	}

		ExcelTools::convert( XL_DayCount,"A360",Basis);
		//	if (Basis.empty() ) Basis = "A360";
			if((DayCount = ARM_ConvDayCount (Basis.c_str())) == ARM_DEFAULT_ERR)	{	ARM_ARG_ERR();	return (LPXLOPER)&XL_result;	}

		ExcelTools::convert( XL_Stubrule,"LS",strStubrule);
		//	if (strStubrule.empty()) strStubrule = "LS";
			if((Stubrule = ARM_ConvStubRule (strStubrule.c_str())) == ARM_DEFAULT_ERR)	{	ARM_ARG_ERR();	return (LPXLOPER)&XL_result;	}
	
		ExcelTools::convert( XL_intRule,"ADJ", strIntRule);
		//	if (strIntRule.empty()) strIntRule = "ADJ";
			intrule = ARM_ConvIntRule (strIntRule.c_str());

		ExcelTools::convert( XL_payCalName, "EUR", payCalName);
		//	if (payCalName.empty() ) payCalName = "EUR";
		
		ExcelTools::convert( XL_PayTiming,"ARR",strPayTiming);
		//	if (strPayTiming.empty()) strPayTiming = "ARR";
			PayTiming = ARM_ConvPayResetRule (strPayTiming.c_str());
			
		ExcelTools::convert( XL_ResetTiming,"ADV",strResetTiming);
		//	if (strResetTiming.empty()) strResetTiming = "ADV";
			ResetTiming = ARM_ConvPayResetRule (strResetTiming.c_str());

		ExcelTools::convert( XL_fwdRule,"MF", strfwdRule);
		//	if(strfwdRule.empty()) strfwdRule = "MF";
			if((fwdRule =  ARM_ConvFwdRule (strfwdRule.c_str(), C_result)) == ARM_DEFAULT_ERR)	{	ARM_ARG_ERR();	return (LPXLOPER)&XL_result; }

		ExcelTools::convert( XL_IncludeMaturity,"N",strIncludeMaturity);
		if (strIncludeMaturity == "Y" ) IncludeMaturity = true;

		ExcelTools::convert( XL_adj,"ADJ", strAdj);
		//	if (strAdj.empty()) strAdj = "ADJ";
			if((Adj = ARM_ConvIntRule (strAdj.c_str())) == ARM_DEFAULT_ERR){	ARM_ARG_ERR();	return (LPXLOPER)&XL_result;	}
		
		ExcelTools::convert( XL_intStartAdj, "", strStartAdj);
		if (strStartAdj.empty()) intStartAdj = K_ADJUSTED;
		else if( intStartAdj = ARM_ConvStartAdjRule (strStartAdj.c_str()) == ARM_DEFAULT_ERR)	{	ARM_ARG_ERR();	return (LPXLOPER)&XL_result;	}

		ExcelTools::convert( XL_AccDayCount, "A360", strAccDayCount);
		//	if (strAccDayCount.empty()) strAccDayCount = "A360";
			if((AccDayCount = ARM_ConvDayCount (strAccDayCount.c_str())) == ARM_DEFAULT_ERR){	ARM_ARG_ERR();	return (LPXLOPER)&XL_result;	}
		
		ExcelTools::convert( XL_ReferenceDate,double(-1), ReferenceDate);
		//if (strReferenceDate.empty() ) ReferenceDate = -1;

		ExcelTools::convert( XL_FirstCpnEffDate,double(-1),FirstCpnEffDate);
		//	if (strFirstCpnEffDate.empty() ) FirstCpnEffDate = -1;

		ExcelTools::convert( XL_AdjCal,"STDCDS", strAdjCal);
		//	if (strAdjCal.empty() ) strAdjCal = "STDCDS";
			if((AdjCal = ARM_ConvAdjCalCDS (CCString(strAdjCal.c_str()), C_result)) == ARM_DEFAULT_ERR){	ARM_ARG_ERR();	return (LPXLOPER)&XL_result;}


		long newId = ICMLOCAL_SCHEDULE_INFO(EffectiveDate,
											MaturityDate,			
											payFrequency,
											ResetFreq ,
											DayCount,
											Stubrule,
											intrule,
											payCalName,
											PayTiming,
											ResetTiming,
											fwdRule,
											IncludeMaturity,
											Adj,
											intStartAdj,
											AccDayCount,
											ReferenceDate,
											FirstCpnEffDate,
											AdjCal,
											prevId);
	
		
		string objectLabel = ExcelCaller::get().setObject(newId, LOCAL_SCHEDULE_INFO_CLASS);
		//ExcelCaller::get().setObject(objectLabel);
		ExcelTools::convert(objectLabel, &XL_result);	
	} 
	catch (Exception& e)
	{
		ExcelCaller::get().setError(e.GetErrorString());
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (std::exception& e)
	{
		ExcelCaller::get().setError(e.what() );
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (...)
	{
		ExcelCaller::get().setError("unknown error");
		ExcelTools::convert("ARM_ERR", &XL_result);
	}

 

    return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Math_CF_SpreadOption(LPXLOPER XL_yearterm,
																LPXLOPER XL_strike,
																LPXLOPER XL_correlation,
																LPXLOPER XL_coefs,
																LPXLOPER XL_spots,
																LPXLOPER XL_vol,
																LPXLOPER XL_intstep)
{

	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double yearterm=0.;
	double strike=0.;
	double correlation=0.;
	VECTOR<double> C_coefs;
	VECTOR<double> C_spots;
	VECTOR<double> C_vol;

	double intstep=0.;
	double intstep_default=40.;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_yearterm, yearterm, " ARM_ERR: yearterm : numeric expected",C_result);
	XL_readNumCell(XL_strike, strike, " ARM_ERR: strike : numeric expected",C_result);
	XL_readNumCell(XL_correlation, correlation, " ARM_ERR: correlation : numeric expected",C_result);
	XL_readNumVector(XL_coefs, C_coefs," ARM_ERR: coefs: double vector expected", C_result);
	XL_readNumVector(XL_spots, C_spots," ARM_ERR: spots: double vector expected", C_result);
	XL_readNumVector(XL_vol, C_vol," ARM_ERR: vol: double vector expected", C_result);
	XL_readNumCellWD(XL_intstep, intstep,intstep_default, " ARM_ERR: intstep : numeric expected",C_result);

	long retCode = ICMLOCAL_Math_CF_SpreadOption(yearterm,strike,correlation,C_coefs,C_spots,C_vol,(int)intstep,C_result);

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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Spread" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Calibrator (LPXLOPER XL_security_vector,
															 LPXLOPER XL_price_vectorBid,
															 LPXLOPER XL_price_vectorAsk,
															 LPXLOPER XL_parameters_vector,
															 LPXLOPER XL_tsparams_vector,
															 LPXLOPER XL_model,
															 LPXLOPER XL_pricertype,
															 LPXLOPER XL_parameters,
															 LPXLOPER XL_parameters_inf,
															 LPXLOPER XL_parameters_sup,
															 LPXLOPER XL_pricingtype,
															 LPXLOPER XL_Optimparameters)
{

	static XLOPER XL_result;
	ARM_result C_result;
	int i=0;


	/// to make the infrastructure very robust we put a try catch!
	try
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<CCString> C_security_vector; 
	VECTOR<double> C_price_vectorBid; 
	VECTOR<double> C_price_vectorAsk; 
	VECTOR<CCString> C_parameters_vector;
	VECTOR<double> C_tsparams_vector; 
	CCString C_model;
	CCString C_parameters;
	CCString C_parameters_inf;
	CCString C_parameters_sup;
	CCString C_paramforts;
	CCString C_pricingtype;
	CCString C_Optimparameters;

	// error
	static int error;
	static char* reason = "";

	long newId = 0;

	XL_readStrVector(XL_security_vector,C_security_vector," ARM_ERR: vector expected",DOUBLE_TYPE,C_result);	
	XL_readNumVector(XL_price_vectorBid,C_price_vectorBid," ARM_ERR: vector expected",C_result);	
	XL_readNumVector(XL_price_vectorAsk,C_price_vectorAsk," ARM_ERR: vector expected",C_result);	
	XL_readStrVector(XL_parameters_vector,C_parameters_vector," ARM_ERR: vector expected",DOUBLE_TYPE,C_result);	
	XL_readNumVector(XL_tsparams_vector,C_tsparams_vector," ARM_ERR: vector expected",C_result);	
	XL_readStrCell(XL_model,C_model," ARM_ERR: object expected",C_result);
	XL_readStrCell(XL_parameters,C_parameters," ARM_ERR: object expected",C_result);
	XL_readStrCell(XL_parameters_inf,C_parameters_inf," ARM_ERR: object expected",C_result);
	XL_readStrCell(XL_parameters_sup,C_parameters_sup," ARM_ERR: object expected",C_result);
	XL_readStrCell(XL_pricingtype,C_pricingtype," ARM_ERR: object expected",C_result);
	XL_readStrCellWD(XL_Optimparameters,C_Optimparameters,""," ARM_ERR: object expected",C_result);

	VECTOR<long> l_security_vector;
	
	for (i=0;i<C_security_vector.size();i++)
	{l_security_vector.push_back(LocalGetNumObjectId (C_security_vector[i]));}

	std::string  c_PricerType;ExcelTools::convert(XL_pricertype,"UNKNOWN",c_PricerType); 
	long pricerType = ARM_ConvCreditPricerType(c_PricerType);

	int pricingtype = ICM_ConvCptType(C_pricingtype,C_result);
	
	long prevId= ExcelCaller::get().getObjectId() ; 

	newId = ICMLOCAL_Calibrator (l_security_vector,
										C_price_vectorBid,
										C_price_vectorAsk,
										C_parameters_vector,
										C_tsparams_vector,
										LocalGetNumObjectId(C_model),
										pricerType,
										LocalGetNumObjectId(C_parameters),
										LocalGetNumObjectId(C_parameters_inf),
										LocalGetNumObjectId(C_parameters_sup),
										pricingtype,
										LocalGetNumObjectId(C_Optimparameters),
										prevId);

	std::string newName = ExcelCaller::get().setObject(newId,LOCAL_CREDIT_CASHFLOWS_CLASS); 
	ExcelTools::convert(newName,&XL_result); 

	}
	catch (Exception& e)
	{
		ExcelCaller::get().setError(e.GetErrorString());
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (std::exception& e)
	{
		ExcelCaller::get().setError(e.what() );
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (...)
	{
		ExcelTools::convert("ARM_ERR", &XL_result);
	}

	return &XL_result;
}
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Random_Generator (LPXLOPER XL_RandomType,
															 LPXLOPER XL_paramId)
{

	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	try
	{
	ARM_NOCALCIFWIZ();

	// C variable

	std::string RandomType; 
	ExcelTools::convert(XL_RandomType,"",RandomType); 

	std::string ParameterId;
	std::string ParameterIdDefault="";
	long ParamLong = -1;
	long newId = 0;
	
	// error
	static int error;
	static char* reason = "";
	ExcelTools::convert(XL_paramId, "",ParameterId);
	if(ParameterIdDefault != ParameterId)
		ParamLong = LocalGetNumObjectId (CCString(ParameterId.c_str()));
	
	long prevId= ExcelCaller::get().getObjectId() ; 
	qRAN_GEN rType;
	bool res=false;
	string list;
	ICM_EnumsCnv::cnv(RandomType,rType,res, list);
	if(!res){
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_Credit_Random_Generator: type of random not defined list is" << list) ;
	}
	
	newId = ICMLOCAL_RandomGenerator (rType,
										ParamLong, 
										prevId);

	std::string newName = ExcelCaller::get().setObject(newId,LOCAL_CREDIT_RANDOM_GENERATOR);
	ExcelTools::convert(newName,&XL_result); 

	}
	catch (Exception& e)
	{
		ExcelCaller::get().setError(e.GetErrorString());
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (std::exception& e)
	{
		ExcelCaller::get().setError(e.what() );
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (...)
	{
		ExcelTools::convert("ARM_ERR", &XL_result);
	}

	return &XL_result;
}
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GenerateOneRandom (LPXLOPER XL_RandomId)
{

	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	try
	{
	ARM_NOCALCIFWIZ();

	// C variable

	CCString C_RandomId;
	double RandomNb;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_RandomId,C_RandomId," ARM_ERR: object expected ",C_result);

	ICMLOCAL_GenerateOneRandom (LocalGetNumObjectId (C_RandomId), 
								RandomNb);

	ExcelTools::convert(RandomNb,&XL_result); 

	}
	catch (Exception& e)
	{
		ExcelCaller::get().setError(e.GetErrorString());
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (std::exception& e)
	{
		ExcelCaller::get().setError(e.what() );
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (...)
	{
		ExcelTools::convert("ARM_ERR", &XL_result);
	}

	return &XL_result;
}
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GenerateRandoms (LPXLOPER XL_RandomId,
																  LPXLOPER XL_DimVector)
{

	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	try
	{
	ARM_NOCALCIFWIZ();

	// C variable

	CCString C_RandomId;
	int DimVector;

	// error
	static int error;
	static char* reason = "";
	
	ExcelTools::convert(XL_DimVector,DimVector); 

	XL_readStrCell(XL_RandomId,C_RandomId," ARM_ERR: object expected ",C_result);
	//XL_readNumCell(XL_DimVector,DimVector," ARM_ERR: Dim Vector integer expected",C_result);

	ARM_Vector RandomVector(DimVector,0.);

	ICMLOCAL_GenerateRandoms (LocalGetNumObjectId (C_RandomId), 
								RandomVector);

	ExcelTools::convert(RandomVector,&XL_result); 

	}
	catch (Exception& e)
	{
		ExcelCaller::get().setError(e.GetErrorString());
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (std::exception& e)
	{
		ExcelCaller::get().setError(e.what() );
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (...)
	{
		ExcelTools::convert("ARM_ERR", &XL_result);
	}

	return &XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_ResetRandom (LPXLOPER XL_RandomId)
{

	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	try
	{
	ARM_NOCALCIFWIZ();

	// C variable

	CCString C_RandomId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_RandomId,C_RandomId," ARM_ERR: object expected ",C_result);

	ICMLOCAL_ResetRandom (LocalGetNumObjectId (C_RandomId));

	ExcelTools::convert("Reset done",&XL_result); 

	}
	catch (Exception& e)
	{
		ExcelCaller::get().setError(e.GetErrorString());
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (std::exception& e)
	{
		ExcelCaller::get().setError(e.what() );
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (...)
	{
		ExcelTools::convert("ARM_ERR", &XL_result);
	}

	return &XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_SetMatuLabel (LPXLOPER XL_CurveId,
														   LPXLOPER XL_labels)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_CurveId;
	VECTOR<CCString> C_labels;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_CurveId,C_CurveId," ARM_ERR: Pay Yield Curve id: object expected",C_result);
	XL_readStrVector(XL_labels,C_labels," ARM_ERR: vector expected",DOUBLE_TYPE,C_result);	
	
	long retCode = ICMLOCAL_SetMatuLabel (LocalGetNumObjectId (C_CurveId), C_labels,C_result);

	if(retCode == ARM_OK)
	{
        FreeCurCellErr ();

        XL_result.xltype = xltypeNum;
        XL_result.val.num =  C_result.getDouble ();
        XL_result.xltype |= xlbitDLLFree;
	}
	else
	{
		ARM_ERR();
	}

//	ARM_END();
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	
	return (LPXLOPER)&XL_result;
}