
#pragma warning(disable :4005 4786 )

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include <ARM\libicm_local\ICM_local_Cair.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>
#include "ExcelTools.h"


_declspec(dllexport) LPXLOPER WINAPI ARM_CreditManager_MarketData(
														LPXLOPER XL_CreditManagerId,
														LPXLOPER XL_ValDate, 
														LPXLOPER XL_Currency,
														LPXLOPER XL_SummitIRCurve,
														LPXLOPER XL_MarketDataParameters,
														LPXLOPER XL_ImposedIRValues
														)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{

		ARM_NOCALCIFWIZ();

	// C variable
	double C_ValDate =0.;
	CCString C_ccy;

	CCString C_ircurve;
	long ircurveId = 0;

	CCString	C_MarketDataParameters;
	CCString	C_ImposedIRValues;
	CCString	C_CreditManagerId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_CreditManagerId,C_CreditManagerId," ARM_ERR: Credit Manager Id: object expected",C_result);
	XL_readNumCell(XL_ValDate,C_ValDate," ARM_ERR: as of date: date expected",C_result);
	XL_readStrCellWD(XL_Currency,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCell(XL_SummitIRCurve,C_ircurve," ARM_ERR: Curve Id expected ",C_result);
	XL_readStrCellWD(XL_MarketDataParameters, C_MarketDataParameters,"NONE"," ARM_ERR: Market Data Parameters Id: Object expected",C_result);
	XL_readStrCellWD(XL_ImposedIRValues, C_ImposedIRValues,"NONE"," ARM_ERR: Imposed IR Values Id: Object expected",C_result);

	long retCode;

	retCode = ICMLOCAL_SetCreditManager_MarketData(
									LocalGetNumObjectId (C_CreditManagerId),
									C_ValDate,
									C_ccy,
									LocalGetNumObjectId (C_ircurve),
									LocalGetNumObjectId (C_MarketDataParameters),
									LocalGetNumObjectId (C_ImposedIRValues),
									C_result);

	if (retCode == ARM_OK)
	{			
		FreeCurCellErr();

		XL_result.xltype = xltypeStr;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_CreditManager_MarketData")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI ARM_CreditManager_CreditData(
														LPXLOPER XL_CreditManagerId,
														LPXLOPER XL_Labels, 
														LPXLOPER XL_Description,
														LPXLOPER XL_CreditSpreads,
														LPXLOPER XL_Maturities,
														LPXLOPER XL_CreditDataParameters,
														LPXLOPER XL_HedgesCDSMaturity,
														LPXLOPER XL_CreditDataCDOSquareParameters,
														LPXLOPER XL_CreditDataCDOSquareData
														)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<CCString> C_Labels; 
	VECTOR<CCString> C_Maturities; 

	VECTOR<double> C_Spreads;

	CCString	C_CreditManagerId;
	CCString	C_Description;
	CCString	C_CreditDataParameters;
	CCString	C_HedgesCDSMaturity;
	CCString	C_CreditDataCDOSquareParameters;
	CCString	C_CreditDataCDOSquareData;

	vector<double> MatIssuersSpread;

	// error
	static int error;
	static char* reason = "";

	// Credit Manager
	XL_readStrCell(XL_CreditManagerId,C_CreditManagerId," ARM_ERR: Credit Manager Id: object expected",C_result);

	// Labels
	XL_readStrVector(XL_Labels,C_Labels," ARM_ERR: Labels: array of strings expected",DOUBLE_TYPE,C_result);

	// Cash Flow Id
	XL_readStrCellWD(XL_Description, C_Description,"NONE"," ARM_ERR: Description Id: Object expected",C_result);

	// Maturities
	XL_readStrVector(XL_Maturities,C_Maturities," ARM_ERR: Labels: array of strings expected",DOUBLE_TYPE,C_result);

	// Cash Flow Id
	XL_readStrCellWD(XL_CreditDataParameters, C_CreditDataParameters,"NONE"," ARM_ERR: Credit Data Parameters: Object expected",C_result);

	// Spreads	
	XL_readNumVector(XL_CreditSpreads, C_Spreads," ARM_ERR: Spreads : matrix of numerics expected",C_result);

	// Hedges CDS Maturity
	XL_readStrCellWD(XL_HedgesCDSMaturity, C_HedgesCDSMaturity,"TO_MATURITY"," ARM_ERR: Hedges CDS Maturity: string expected",C_result);

	// CDO^2 Parameters Cash Flow Id
	XL_readStrCellWD(XL_CreditDataCDOSquareParameters, C_CreditDataCDOSquareParameters,"NONE"," ARM_ERR: Credit Data CDO Square Parameters: Object expected",C_result);

	// CDO^2 Data Cash Flow Id
	XL_readStrCellWD(XL_CreditDataCDOSquareData, C_CreditDataCDOSquareData,"NONE"," ARM_ERR: Credit Data CDO Square Data: Object expected",C_result);

	long retCode;

	retCode = ICMLOCAL_SetCreditManager_CreditData(
									LocalGetNumObjectId (C_CreditManagerId),
									C_Labels,
									LocalGetNumObjectId (C_Description),
									C_Spreads,
									C_Maturities,
									LocalGetNumObjectId (C_CreditDataParameters),
									C_HedgesCDSMaturity,
									LocalGetNumObjectId (C_CreditDataCDOSquareParameters),
									LocalGetNumObjectId (C_CreditDataCDOSquareData),
									C_result);

	if (retCode == ARM_OK)
	{			
		FreeCurCellErr();

		XL_result.xltype = xltypeStr;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_CreditManager_CreditData")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;

}

_declspec(dllexport) LPXLOPER WINAPI ARM_CreditManager_CreditModel(
														LPXLOPER XL_CreditManagerId,
														LPXLOPER XL_CreditModelParameters, 
														LPXLOPER XL_CorrelationValue,
														LPXLOPER XL_BetaVector,
														LPXLOPER XL_BaseCorrelationStrikesArray,
														LPXLOPER XL_BaseCorrelationValuesArray,
														LPXLOPER XL_CorrelationMatrixId,
														LPXLOPER XL_CorrelationId,
														LPXLOPER XL_FactorLoadingParameters
														)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<CCString> C_Labels; 
	VECTOR<CCString> C_Maturities; 

	VECTOR<double> C_BetaVector;

	double	CorrelationValue;

//	double	DefaultCorrelationValue = 0.5;

	CCString	C_CreditManagerId;
	CCString	C_CreditModelParameters;
	CCString	C_CorrelationId;
	CCString	C_CorrelationMatrixId;
	CCString	C_CreditFactorLoadingParameters;

	VECTOR<double> C_BaseCorrelationStrikesArray;
	VECTOR<double> C_BaseCorrelationValuesArray;

	double	C_Correlation = 0;
	double	C_Correlation_default = 0.5;
	double	size = 0.;

	// error
	static int error;
	static char* reason = "";

	// Credit Manager
	XL_readStrCell(XL_CreditManagerId,C_CreditManagerId," ARM_ERR: Credit Manager Id: object expected",C_result);

	// Cash Flow Id
	XL_readStrCellWD(XL_CreditModelParameters, C_CreditModelParameters,"NONE"," ARM_ERR: Parameters Id: Object expected",C_result);

	// Correlation Single Value
	XL_readNumCell(XL_CorrelationValue, CorrelationValue," ARM_ERR: Correlation : numeric expected",C_result);

	// Beta Vector	
	XL_readNumVector(XL_BetaVector, C_BetaVector," ARM_ERR: Beta Vector: vector of numerics expected",C_result);

	// Base Correlation Strikes Array	
	XL_readNumVector(XL_BaseCorrelationStrikesArray, C_BaseCorrelationStrikesArray," ARM_ERR: Base Correlation Strikes Array: vector of numerics expected",C_result);

	// Base Correlation Values Array	
	XL_readNumVector(XL_BaseCorrelationValuesArray, C_BaseCorrelationValuesArray," ARM_ERR: Base Correlation Values Array: vector of numerics expected",C_result);

	// Correlation Matrix Object
    XL_readStrCell(XL_CorrelationMatrixId, C_CorrelationMatrixId ," ARM_ERR: Correlation Matrix: object expected",C_result);

	// Correlation Object
    XL_readStrCell(XL_CorrelationId, C_CorrelationId ," ARM_ERR: Correlation: object expected",C_result);

	// Factor Loading Cash Flow Id
	XL_readStrCellWD(XL_FactorLoadingParameters, C_CreditFactorLoadingParameters,"NONE"," ARM_ERR: Factor Loading - Parameters Id: Object expected",C_result);

	long retCode;

	retCode = ICMLOCAL_SetCreditManager_CreditModel(
									LocalGetNumObjectId (C_CreditManagerId),
									LocalGetNumObjectId (C_CreditModelParameters),
									CorrelationValue,
									C_BetaVector,
									C_BaseCorrelationStrikesArray,
									C_BaseCorrelationValuesArray,
									LocalGetNumObjectId(C_CorrelationMatrixId),
									LocalGetNumObjectId(C_CorrelationId),
									LocalGetNumObjectId (C_CreditFactorLoadingParameters),
									C_result);

	if (retCode == ARM_OK)
	{			
		FreeCurCellErr();

		XL_result.xltype = xltypeStr;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_CreditManager_CreditModel")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI ARM_CreditManager_CreditProduct(
														LPXLOPER XL_CreditManagerId,
														LPXLOPER XL_CreditProductDefaultLeg, 
														LPXLOPER XL_CreditProductPremiumLeg,
														LPXLOPER XL_CreditProductPricingParameters
														)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

	// C variable
	CCString	C_CreditManagerId;
	CCString	C_CreditProductDefaultLegId;
	CCString	C_CreditProductPremiumLegId;
	CCString	C_CreditProductPricingParametersId;

	// error
	static int error;
	static char* reason = "";

	// Credit Manager
	XL_readStrCell(XL_CreditManagerId,C_CreditManagerId," ARM_ERR: Credit Manager Id: object expected",C_result);

	// Credit Product Default Leg
    XL_readStrCell(XL_CreditProductDefaultLeg, C_CreditProductDefaultLegId ," ARM_ERR: Credit Product Default Leg Id: object expected",C_result);

	// Credit Product Premium Leg
    XL_readStrCell(XL_CreditProductPremiumLeg, C_CreditProductPremiumLegId ," ARM_ERR: Credit Product Premium Leg Id: object expected",C_result);

	// Credit Product Pricing Parameters
    XL_readStrCell(XL_CreditProductPricingParameters, C_CreditProductPricingParametersId ," ARM_ERR: Credit Product Pricing Parameters Id: object expected",C_result);

	long retCode;

	retCode = ICMLOCAL_SetCreditManager_CreditProduct(
									LocalGetNumObjectId (C_CreditManagerId),
									LocalGetNumObjectId (C_CreditProductDefaultLegId),
									LocalGetNumObjectId (C_CreditProductPremiumLegId),
									LocalGetNumObjectId (C_CreditProductPricingParametersId),
									C_result);

	if (retCode == ARM_OK)
	{			
		FreeCurCellErr();

		XL_result.xltype = xltypeStr;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_CreditManager_CreditProduct")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI ARM_CreditManager_GetData(
														LPXLOPER XL_CreditManagerId,
														LPXLOPER XL_DataLabel
														)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

	// C variable
	CCString	C_CreditManagerId;
	CCString	C_DataLabel;

	// error
	static int error;
	static char* reason = "";

	// Credit Manager
	XL_readStrCell(XL_CreditManagerId,C_CreditManagerId," ARM_ERR: Credit Manager Id: object expected",C_result);

	// Data Label
	XL_readStrCellWD(XL_DataLabel, C_DataLabel, "NPV"," ARM_ERR: data label: string expected", C_result);

	long retCode;

	retCode = ICMLOCAL_GetCreditManager_DataFromLabel(
									LocalGetNumObjectId (C_CreditManagerId),
									C_DataLabel,
									C_result);

	if (retCode == ARM_OK)
	{
		// just a double, waiting for some arrays
		FreeCurCellErr();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_CreditManager_GetData")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI ARM_CreditManager_CreditCalibrator(
														LPXLOPER XL_CreditManagerId,
														LPXLOPER XL_Description,
														LPXLOPER XL_Maturities
														)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<CCString> C_Maturities; 

	CCString	C_CreditManagerId;
	CCString	C_Description;

	vector<double> Lags;
	vector<double> P_Default;

	// error
	static int error;
	static char* reason = "";

	// Credit Manager
	XL_readStrCell(XL_CreditManagerId,C_CreditManagerId," ARM_ERR: Credit Manager Id: object expected",C_result);

	// Cash Flow Id
	XL_readStrCellWD(XL_Description, C_Description,"NONE"," ARM_ERR: Description Id: Object expected",C_result);

	// Maturities
	XL_readStrVector(XL_Maturities,C_Maturities," ARM_ERR: Labels: array of strings expected",DOUBLE_TYPE,C_result);


	long retCode;

	retCode = ICMLOCAL_SetCreditManager_CreditCalibrator(
									LocalGetNumObjectId (C_CreditManagerId),
									LocalGetNumObjectId (C_Description),
									C_Maturities,
									C_result);

	if (retCode == ARM_OK)
	{			
		FreeCurCellErr();

		XL_result.xltype = xltypeStr;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_CreditManager_CreditData")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;

}


_declspec(dllexport) LPXLOPER WINAPI ARM_CreditManager_CorrelationCalibrator(
														LPXLOPER XL_CreditManagerId,
														LPXLOPER XL_DataDescriptionMatrix,
														LPXLOPER XL_ParametersVector,
														LPXLOPER XL_TranchesId
														)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

	// C variable
	CCString	C_CreditManagerId;
	CCString	C_DescriptionId;
	CCString	C_ParametersId;

	VECTOR<CCString>	V_TranchesId;
	VECTOR<long>		V_L_TranchesId;

	// error
	static int error;
	static char* reason = "";

	// Credit Manager
	XL_readStrCell(XL_CreditManagerId,C_CreditManagerId," ARM_ERR: Credit Manager Id: object expected",C_result);

	// Data Matrix Id
	XL_readStrCellWD(XL_DataDescriptionMatrix, C_DescriptionId,"NONE"," ARM_ERR: Description Id: Object expected",C_result);

	// Data Parameters Id
	XL_readStrCellWD(XL_ParametersVector, C_ParametersId,"NONE"," ARM_ERR: Parameters Id: Object expected",C_result);

	// Tranches Ids
    XL_readStrVector(XL_TranchesId, V_TranchesId ," ARM_ERR: Tranches Id: object vector expected",DOUBLE_TYPE,C_result);

	long retCode;

	int	i, size;
	size	= V_TranchesId.size();

    for (i=0; i<size; i++)
        V_L_TranchesId.push_back (LocalGetNumObjectId (V_TranchesId[i]));

	retCode = ICMLOCAL_SetCreditManager_CorrelationCalibrator(
									LocalGetNumObjectId (C_CreditManagerId),
									LocalGetNumObjectId (C_DescriptionId),
									LocalGetNumObjectId (C_ParametersId),
									V_L_TranchesId,
									C_result);

	if (retCode == ARM_OK)
	{			
		FreeCurCellErr();

		XL_result.xltype = xltypeStr;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_CreditManager_CreditData")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;

}


_declspec(dllexport) LPXLOPER WINAPI ARM_Get_CreditManager()
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {

	ARM_NOCALCIFWIZ();

	// error
	static int error;
	static char* reason = "";

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_PRICER_CLASS;
	CCString stringId = GetLastCurCellEnvValue();

	if (!stringId)
	{
		retCode = ICMLOCAL_CreditManager(C_result);

		if (retCode == ARM_OK)
		{
			objId = C_result.getLong();

			LocalSetCurCellEnvValue(curClass, objId); 

			stringId = LocalMakeObjectId(objId, curClass);
		}
	}
	else
	{
		prevClass = LocalGetStringObjectClass (stringId);
		
		objId = LocalGetNumObjectId(stringId);
			
		if (curClass == prevClass)
		{
			retCode = ICMLOCAL_CreditManager(C_result,
									 objId);

			if (retCode == ARM_OK)
			{			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent();

			retCode = ICMLOCAL_CreditManager(C_result);

			if (retCode == ARM_OK)
			{
				objId = C_result.getLong();
			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId(objId, curClass);
			}
		}
	}

	if (retCode == ARM_OK)
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


_declspec(dllexport) LPXLOPER WINAPI ARM_CreditManager_GetDataMatrix(
														LPXLOPER XL_CreditManagerId,
														LPXLOPER XL_DataLabel
														)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

	// C variable
	CCString	C_CreditManagerId;
	CCString	C_DataLabel;

	// error
	static int error;
	static char* reason = "";

	// Credit Manager
	XL_readStrCell(XL_CreditManagerId,C_CreditManagerId," ARM_ERR: Credit Manager Id: object expected",C_result);

	// Data Label
	XL_readStrCellWD(XL_DataLabel, C_DataLabel, "NPV"," ARM_ERR: data label: string expected", C_result);

	long retCode;

	int	i,j;
	int	OutputNbRows;
	int	OutputNbCols;
	vector<double*>		OutputMatrix;
	VECTOR<CCString>	OutputLabels;

	LPXLOPER pxArray;

	retCode = ICMLOCAL_GetCreditManager_DataMatrixFromLabel(
									LocalGetNumObjectId (C_CreditManagerId),
									C_DataLabel,
									OutputMatrix,
									OutputLabels,
									OutputNbRows,
									OutputNbCols,
									C_result);


	if (retCode == ARM_OK)
	{
		int nbrows		=	OutputNbRows+1;		// add a string... label
		int nbcolumns	=	OutputNbCols;		// Outputs
		
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		// First of all, leave it blank
		for (i=0;i<nbrows;i++)
		{
			for (j=0;j<nbcolumns;j++)
			{
				pxArray[XL_Coordonnate2Rank (i,j, nbcolumns)].xltype = xltypeStr;
				pxArray[XL_Coordonnate2Rank (i,j, nbcolumns)].val.str = XL_StrC2StrPascal ("");
			}
		}
				
		// Label
		for (j=0;j<nbcolumns;j++)
		{			
			pxArray[XL_Coordonnate2Rank (0, j, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (0, j, nbcolumns)].val.str = XL_StrC2StrPascal(OutputLabels[j]); 
		}
	
		// All outputs
		for (i=1;i<nbrows;i++)
			for (j=0;j<nbcolumns;j++)
			{
				pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype = xltypeNum;
				pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.num = OutputMatrix[j][i-1];
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_CreditManager_GetData")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}