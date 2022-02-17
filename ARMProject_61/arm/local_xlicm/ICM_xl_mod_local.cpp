
#pragma warning(disable :4005 4786 )

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include <ARM\libicm_local\ICM_local_mod.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"
#include "ExcelTools.h"
#include "ARM_xl_trycatch_local.h"




__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_ModelMultiCurves (LPXLOPER XL_DefCurvesID,
																   LPXLOPER XL_DiscountCurveID,
																   LPXLOPER XL_CorrelationId,
																   LPXLOPER XL_RecoveryRates,
																   LPXLOPER XL_VolCurve,
																   LPXLOPER XL_CloneOrNot,
																   LPXLOPER XL_InfCurv,
																   LPXLOPER XL_CpnIrCurv)
{
//    ARM_BEGIN();
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

    // C variable
	int i = 0;
	double C_NbOfCurve =0.;

	VECTOR<CCString> V_DefCurvesID;
	VECTOR<long>	V_LDefCurvesID;
	VECTOR<double> V_RecoveryRates;

	double C_Correlation = 0;
	double C_Correlation_default = 0.5;
	double size = 0.;

	CCString C_DiscountCurveID =0.;
	CCString C_CorrID =0.;
	CCString C_CorrelationId;
	CCString C_VolCurve;
	CCString C_InfCurv;
	CCString C_CpnIrCurv;

    // error
    static int error;
    static char* reason = "";

    XL_readStrVector(XL_DefCurvesID, V_DefCurvesID ," ARM_ERR: Default Curves : object vector expected",DOUBLE_TYPE,C_result);
    XL_readStrCell(XL_DiscountCurveID, C_DiscountCurveID ," ARM_ERR: Discount Curve : object expected",C_result);

	C_NbOfCurve = V_DefCurvesID.size();

    for (i = 0; i < C_NbOfCurve; i++)
    {
        V_LDefCurvesID.push_back (LocalGetNumObjectId (V_DefCurvesID[i]));
    }

	if (XL_RecoveryRates->xltype != 128) 
	{
	XL_readNumVector(XL_RecoveryRates, V_RecoveryRates," ARM_ERR: Recovery Rates : vector of numerics expected",C_result);
	}
	else
	    for (i = 0; i < C_NbOfCurve; i++)
			V_RecoveryRates.push_back(-999.);

    XL_readStrCellWD(XL_CorrelationId, C_CorrelationId ,""," ARM_ERR: Correlation : object expected",C_result);

	long retCode;
    long objId;
    CCString prevClass;

	CCString C_CloneOrNot;
	bool b_CloneOrNot = true; 
    
	if (V_DefCurvesID.size() != V_RecoveryRates.size())
	{
		C_result.setMsg("Incompatibilité de taille entre les vecteurs & matrices en input du modele");
		ARM_ERR();
	}
	XL_readStrCellWD(XL_VolCurve, C_VolCurve,"" ," ARM_ERR: Volatility Curve : object expected",C_result);
	XL_readStrCellWD(XL_CloneOrNot,C_CloneOrNot,"Y"," ARM_ERR: Y or N expected",C_result);
	XL_readStrCellWD(XL_InfCurv,C_InfCurv,"NONE"," ARM_ERR: Cpn Inflation curve : object expected",C_result);
	XL_readStrCellWD(XL_CpnIrCurv,C_CpnIrCurv,"NONE"," ARM_ERR: Cpn Zero curve : object expected",C_result);

	if (!(C_CloneOrNot=="Y")) b_CloneOrNot = false;

    CCString curClass = LOCAL_MULTICURVESMODEL_CLASS;
    CCString stringId = GetLastCurCellEnvValue();

    if (!stringId)
    {
        retCode = ICMLOCAL_ModelMultiCurves( C_NbOfCurve, 
									   V_LDefCurvesID,
									   LocalGetNumObjectId(C_DiscountCurveID),
									   V_RecoveryRates ,
									   LocalGetNumObjectId(C_CorrelationId),
									   LocalGetNumObjectId(C_VolCurve),
									   b_CloneOrNot,
									   LocalGetNumObjectId(C_InfCurv),
									   LocalGetNumObjectId(C_CpnIrCurv),
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
        prevClass = LocalGetStringObjectClass(stringId);
        
        objId = LocalGetNumObjectId(stringId);
            
        if ( curClass == prevClass )
        {
        retCode = ICMLOCAL_ModelMultiCurves( C_NbOfCurve, 
									   V_LDefCurvesID,
									   LocalGetNumObjectId(C_DiscountCurveID),
									   V_RecoveryRates ,
									   LocalGetNumObjectId(C_CorrelationId),
									   LocalGetNumObjectId(C_VolCurve),
									   b_CloneOrNot,
									   LocalGetNumObjectId(C_InfCurv),
									   LocalGetNumObjectId(C_CpnIrCurv),
									   C_result,
									   objId);


			if ( retCode == ARM_OK )
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
        else
        {
            FreeCurCellContent();

			retCode = ICMLOCAL_ModelMultiCurves( C_NbOfCurve, 
									   V_LDefCurvesID,
									   LocalGetNumObjectId(C_DiscountCurveID),
									   V_RecoveryRates ,
									   LocalGetNumObjectId(C_CorrelationId),
									   LocalGetNumObjectId(C_VolCurve),
									   b_CloneOrNot,
									   LocalGetNumObjectId(C_InfCurv),
									   LocalGetNumObjectId(C_CpnIrCurv),
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

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_FRN" )


    
    return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_ModelMultiCvMktDataMng (LPXLOPER XL_DefCurvesID,
																		LPXLOPER XL_DiscountCurveID,
																		LPXLOPER XL_RecoveryRates,
																		LPXLOPER XL_Correlation,
																		LPXLOPER XL_MktDataMng,
																		LPXLOPER XL_VolCurve,
																		LPXLOPER XL_CloneOrNot)
{
//    ARM_BEGIN();

    // return
    static XLOPER XL_result;
	if (ExcelCaller::get().isCalledByWizard()) 
		{
			ExcelTools::convert("",&XL_result);
			return &XL_result; 
		}
	try {
		// C variable
		int i = 0;
		double C_NbOfCurve =0.;

		VECTOR<string> V_DefCurvesID;
		
		VECTOR<double> V_RecoveryRates,vDefault;

		double C_Correlation = 0;
		double C_Correlation_default = 0.5;
		double size = 0.;

		string C_DiscountCurveID("");
		string C_MktDatMng("");
		string C_VolCurve("");
		string C_CloneOrNot("");
		string C_CorrelId("");

		ExcelTools::convert(XL_DefCurvesID, V_DefCurvesID);
		ExcelTools::convert(XL_RecoveryRates, vDefault,V_RecoveryRates);
		ExcelTools::convert(XL_Correlation, "",C_CorrelId);
		ExcelTools::convert(XL_DiscountCurveID, "",C_DiscountCurveID);
		ExcelTools::convert(XL_MktDataMng, "",C_MktDatMng);
		ExcelTools::convert(XL_VolCurve, "",C_VolCurve);
		ExcelTools::convert(XL_CloneOrNot,"Y",C_CloneOrNot);

		C_NbOfCurve = V_DefCurvesID.size();
		VECTOR<long>	V_DefCurvesLong;

		for (i = 0; i < C_NbOfCurve; i++)
		{
			V_DefCurvesLong.push_back (LocalGetNumObjectId (CCString(V_DefCurvesID[i].c_str())));
		}

		if (V_DefCurvesID.size() != V_RecoveryRates.size())
		{
			ICMTHROW(ERR_INVALID_ARGUMENT, "ARM_Credit_ModelMultiCvMktDataMng : Size of DefcurveV et recoRateV != ");
		}
		
		bool b_CloneOrNot = true; 
		if ((C_CloneOrNot=="N")) b_CloneOrNot = false;

		long prevId = ExcelCaller::get().getObjectId();

		long retCode = ICMLOCAL_ModelMultiCurves( C_NbOfCurve, 
										   V_DefCurvesLong,
										   LocalGetNumObjectId(CCString(C_DiscountCurveID.c_str())),
										   V_RecoveryRates ,
										   LocalGetNumObjectId(CCString(C_CorrelId.c_str())),
										   LocalGetNumObjectId(CCString(C_MktDatMng.c_str())),
										   LocalGetNumObjectId(CCString(C_VolCurve.c_str())),
										   b_CloneOrNot, prevId);
		string objLabel = ExcelCaller::get().setObject(retCode, LOCAL_MULTICURVESMODEL_CLASS);
		//ExcelCaller::get().setObject(objLabel);
		ExcelTools::convert(objLabel, &XL_result);
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

/*
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_SetPropMixCopule (LPXLOPER XL_ModelId,
																   LPXLOPER XL_PropIndep,
																   LPXLOPER XL_PropFullCorrel)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_ModelId;

	double PropIndep = 0.;
	double PropFullCorrel = 0.;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_ModelId,C_ModelId," ARM_ERR: security id: object expected",C_result);
	XL_readNumCell(XL_PropIndep,PropIndep," ARM_ERR: PropIndep: numeric expected",C_result);
	XL_readNumCell(XL_PropFullCorrel,PropFullCorrel," ARM_ERR: PropFullCorrel: date expected",C_result);
	
	long retCode = ICMLOCAL_SetPropMixCopule (LocalGetNumObjectId (C_ModelId),
											  PropIndep,
											  PropFullCorrel,
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

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_FRN" )


	
	return (LPXLOPER)&XL_result;
}

*/
_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DefProbModelNew(LPXLOPER XL_DefProb ,
															 LPXLOPER XL_Disc,
															 LPXLOPER XL_VolCurve)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_DefProb; 
	CCString C_Disc;
	CCString C_VolCurve;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_DefProb,C_DefProb," ARM_ERR: curve id: DefProbCurve object expected",C_result);
	XL_readStrCell(XL_Disc,C_Disc," ARM_ERR: curve id: DiscCurve object expected",C_result);
	XL_readStrCellWD(XL_VolCurve,C_VolCurve,""," ARM_ERR: curve id: Volatility object expected",C_result);
	
	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_DEFPROBMODEL_CLASS;
	CCString stringId = GetLastCurCellEnvValue();
	
	if (!stringId)
	{
		retCode = ICMLOCAL_DefProbModel (LocalGetNumObjectId (C_DefProb),
										 LocalGetNumObjectId (C_Disc),
										 LocalGetNumObjectId (C_VolCurve),
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
			retCode = ICMLOCAL_DefProbModel (LocalGetNumObjectId (C_DefProb),
											 LocalGetNumObjectId (C_Disc),
											 LocalGetNumObjectId (C_VolCurve),
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

			retCode = ICMLOCAL_DefProbModel (LocalGetNumObjectId (C_DefProb),
											 LocalGetNumObjectId (C_Disc),
											 LocalGetNumObjectId (C_VolCurve),
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
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ..." )
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_MetaModel (LPXLOPER XL_ModelID,
															LPXLOPER XL_CashFlowID)
{
//    ARM_BEGIN();

    // return
    static XLOPER XL_result;
    ARM_result C_result;
	ARM_XL_TRY_BLOCK_BEGIN
	{

	ARM_NOCALCIFWIZ();

    // C variable
	int i = 0;
	double C_NbOfCurve =0.;

	VECTOR<CCString> V_DefCurvesID;
	VECTOR<long>	V_LDefCurvesID;
	double size = 0.;

	VECTOR<CCString> V_CashFlowId;
	VECTOR<int> V_LCashFlowId;

    // error
    static int error;
    static char* reason = "";

    XL_readStrVector(XL_ModelID, V_DefCurvesID ," ARM_ERR: Models : object vector expected",DOUBLE_TYPE,C_result);
	XL_readStrVector(XL_CashFlowID,V_CashFlowId," ARM_ERR: PricerType: string expected",DOUBLE_TYPE,C_result);

	C_NbOfCurve = V_DefCurvesID.size();

    for (i = 0; i < C_NbOfCurve; i++)
    {
        V_LDefCurvesID.push_back (LocalGetNumObjectId (V_DefCurvesID[i]));
		V_LCashFlowId.push_back (ARM_ConvCreditPricerType (V_CashFlowId[i], C_result));
    }

	long retCode;
    long objId;
    CCString prevClass;
    
    CCString curClass = LOCAL_MULTICURVESMODEL_CLASS;
    CCString stringId = GetLastCurCellEnvValue();




    if (!stringId)
    {
        retCode = ICMLOCAL_MetaModel(C_NbOfCurve, 
									 V_LDefCurvesID,
									 V_LCashFlowId,
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
        prevClass = LocalGetStringObjectClass(stringId);
        
        objId = LocalGetNumObjectId(stringId);
            
        if ( curClass == prevClass )
        {
        retCode = ICMLOCAL_MetaModel(C_NbOfCurve, 
									 V_LDefCurvesID,
									 V_LCashFlowId,
									 C_result,
									 objId);


			if ( retCode == ARM_OK )
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
        else
        {
            FreeCurCellContent();

			retCode = ICMLOCAL_MetaModel(C_NbOfCurve, 
										V_LDefCurvesID,
										V_LCashFlowId,
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

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ..." )


    
    return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_MarketDataMng (LPXLOPER XL_MktData)
{
//    ARM_BEGIN();

    // return
    static XLOPER XL_result;
	try
	{
		long prevId = ExcelCaller::get().getObjectId();

		vector<string> lLabelsStr;
		ExcelTools::convert(XL_MktData,lLabelsStr) ;
		VECTOR<long> lLabelsId;
		lLabelsId.resize(lLabelsStr.size());
		for (int i = 0; i < lLabelsStr.size(); i++)  
			lLabelsId[i] = LocalGetNumObjectId ((CCString)lLabelsStr[i].c_str());

		long newId = ICMLOCAL_MarketDataMng(lLabelsId,prevId);
		
		string objectLabel = ExcelCaller::get().setObject(newId, LOCAL_MARKETDATAMANAGER_CLASS);
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

/*
// MarketDataMng from Laurent Jacquel. Used in Options pricing
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_MarketDataMngMMC(LPXLOPER XL_MktData)
{
//    ARM_BEGIN();

    // return
    static XLOPER XL_result;
	try
	{
		long prevId = ExcelCaller::get().getObjectId();

		vector<string> lLabelsStr;
		ExcelTools::convert(XL_MktData,lLabelsStr) ;
		VECTOR<long> lLabelsId;
		lLabelsId.resize(lLabelsStr.size());
		for (int i = 0; i < lLabelsStr.size(); i++)  
			lLabelsId[i] = LocalGetNumObjectId ((CCString)lLabelsStr[i].c_str());

		long newId = ICMLOCAL_MarkerDataMng(lLabelsId,prevId);
		
		string objectLabel = ExcelCaller::get().makeObjectLabel(newId, LOCAL_MARKETDATAMANAGER_CLASS);
		ExcelCaller::get().setObject(objectLabel);
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
*/
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Customized_MultiCurves(
														LPXLOPER XL_ModelMultiCurvesId,
														LPXLOPER XL_DiscountCurveId,
														LPXLOPER XL_Labels,
														LPXLOPER XL_Spreads,
														LPXLOPER XL_Spreads_Maturities,
														LPXLOPER XL_Data_Description, 
														LPXLOPER XL_Market_Parameters,
														LPXLOPER XL_CDO_Square_Data,
														LPXLOPER XL_CDO_Square_Parameters,
														LPXLOPER XL_Correlation
														)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();
	// C variable
	VECTOR<CCString> C_Labels; 
	VECTOR<CCString> C_Maturities; 

	VECTOR<double> C_Spreads;

	CCString	C_ModelMultiCurvesId;
	CCString	C_DiscountCurveId	=	0.;
	CCString	C_CorrelationId;
	
	CCString	C_Data_Description;
	CCString	C_Market_Parameters;
	CCString	C_CDO_Square_Parameters;
	CCString	C_CDO_Square_Data;

	CCString	C_ModelMultiCurves_Default = "";

	vector<double> MatIssuersSpread;

	// error
	static int error;
	static char* reason = "";

	// Model Multi Curves Id - Optional
	XL_readStrCellWD(XL_ModelMultiCurvesId, C_ModelMultiCurvesId, "" ," ARM_ERR: Credit Manager Id : object expected",C_result);
	
	// Discount Curve
    XL_readStrCell(XL_DiscountCurveId, C_DiscountCurveId," ARM_ERR: Discount Curve : object expected",C_result);

	// Labels
	XL_readStrVector(XL_Labels, C_Labels," ARM_ERR: Labels: array of strings expected",DOUBLE_TYPE,C_result);

	// Spreads	
	XL_readNumVector(XL_Spreads, C_Spreads," ARM_ERR: Spreads : matrix of numerics expected",C_result);

	// Spreads Maturities
	XL_readStrVector(XL_Spreads_Maturities, C_Maturities," ARM_ERR: Maturities: array of strings expected",DOUBLE_TYPE,C_result);

	// Cash Flow Id
	XL_readStrCellWD(XL_Data_Description, C_Data_Description,"NONE"," ARM_ERR: DData escription Id: Object expected",C_result);

	// Cash Flow Id
	XL_readStrCellWD(XL_Market_Parameters, C_Market_Parameters,"NONE"," ARM_ERR: Market Parameters: Object expected",C_result);

	// CDO^2 Parameters Cash Flow Id
	XL_readStrCellWD(XL_CDO_Square_Parameters, C_CDO_Square_Parameters,"NONE"," ARM_ERR: Credit Data CDO Square Parameters: Object expected",C_result);

	// CDO^2 Data Cash Flow Id
	XL_readStrCellWD(XL_CDO_Square_Data, C_CDO_Square_Data,"NONE"," ARM_ERR: Credit Data CDO Square Data: Object expected",C_result);

	// Correlation
	XL_readStrCell(XL_Correlation, C_CorrelationId," ARM_ERR: Correlation Id: object expected",C_result);


	// ----------------------------------------------------------------------

	CCString curClass = LOCAL_CUSTOMIZED_CREDIT_MULTI_CURVES_CLASS;
	CCString stringId = GetLastCurCellEnvValue();

	long	retCode;
    long	objId;
    CCString	prevClass;

	if (!stringId)
	{
		retCode = ICMLOCAL_Customized_Credit_MultiCurves(
										LocalGetNumObjectId(C_DiscountCurveId),
										C_Labels,
										C_Spreads,
										C_Maturities,
										LocalGetNumObjectId(C_Data_Description),
										LocalGetNumObjectId(C_Market_Parameters),
										LocalGetNumObjectId(C_CDO_Square_Parameters),
										LocalGetNumObjectId(C_CDO_Square_Data),
										LocalGetNumObjectId(C_CorrelationId),
										LocalGetNumObjectId(C_ModelMultiCurvesId),
										C_result);

		if (retCode == ARM_OK)
		{
			objId = C_result.getLong();

			LocalSetCurCellEnvValue(curClass, objId); 

			stringId = LocalMakeObjectId(objId, curClass);
		}
	}
	else
	{
		prevClass = LocalGetStringObjectClass(stringId);
    
		objId = LocalGetNumObjectId(stringId);
        
		if (curClass == prevClass)
		{
			retCode = ICMLOCAL_Customized_Credit_MultiCurves( 
										LocalGetNumObjectId(C_DiscountCurveId),
										C_Labels,
										C_Spreads,
										C_Maturities,
										LocalGetNumObjectId(C_Data_Description),
										LocalGetNumObjectId(C_Market_Parameters),
										LocalGetNumObjectId(C_CDO_Square_Parameters),
										LocalGetNumObjectId(C_CDO_Square_Data),
										LocalGetNumObjectId(C_CorrelationId),
										LocalGetNumObjectId(C_ModelMultiCurvesId),
								   C_result,
								   objId);


			if (retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue(curClass, objId); 

				stringId = LocalMakeObjectId(objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent();

			retCode = ICMLOCAL_Customized_Credit_MultiCurves(
										LocalGetNumObjectId(C_DiscountCurveId),
										C_Labels,
										C_Spreads,
										C_Maturities,
										LocalGetNumObjectId(C_Data_Description),
										LocalGetNumObjectId(C_Market_Parameters),
										LocalGetNumObjectId(C_CDO_Square_Parameters),
										LocalGetNumObjectId(C_CDO_Square_Data),
										LocalGetNumObjectId(C_CorrelationId),
										LocalGetNumObjectId(C_ModelMultiCurvesId),
									   C_result);
    
				if (retCode == ARM_OK)
				{
					objId = C_result.getLong();

					LocalSetCurCellEnvValue(curClass, objId); 

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
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ... " )

    
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_SetVolCurve(LPXLOPER XL_ModelId,
															  LPXLOPER XL_VolCurve)
{
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable
	string C_ModelId;
	string C_VolCurve;
	
	// error
	static int error;
	static char* reason = "";

	ExcelTools::convert(XL_ModelId, C_ModelId);
	ExcelTools::convert(XL_VolCurve, C_VolCurve);
	long retCode = ICMLOCAL_SetVolCurve(LocalGetNumObjectId (CCString(C_ModelId.c_str())),
										  LocalGetNumObjectId (CCString(C_VolCurve.c_str())),
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