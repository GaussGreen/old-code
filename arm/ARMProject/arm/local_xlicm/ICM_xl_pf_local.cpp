
#pragma warning(disable :4005 4786 )

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include <ARM\libicm_local\ICM_local_pf.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include "ExcelTools.h"




__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Portfolio (LPXLOPER XL_SecuritiesID,
															LPXLOPER XL_CashFlowID)
{
//    ARM_BEGIN();

    // return
    static XLOPER XL_result;
    ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

    // C variable
	int i = 0;
	double C_NbOfCurve =0.;

	VECTOR<CCString> V_SecuritiesID;
	VECTOR<long>	V_LSecuritiesID;

	long CashFlowId = 0;
	CCString C_CashFlowId;

	double size = 0.;

    // error
    static int error;
    static char* reason = "";

    XL_readStrVector(XL_SecuritiesID, V_SecuritiesID ," ARM_ERR: Securities : object vector expected",DOUBLE_TYPE,C_result);
	XL_readStrCellWD(XL_CashFlowID,C_CashFlowId,"NONE"," ARM_ERR: frequency: string expected",C_result);

	C_NbOfCurve = V_SecuritiesID.size();

    for (i = 0; i < C_NbOfCurve; i++)
    {
        V_LSecuritiesID.push_back (LocalGetNumObjectId (V_SecuritiesID[i]));
    }

	CashFlowId = LocalGetNumObjectId (C_CashFlowId);

	long retCode;
    long objId;
    CCString prevClass;
    
    CCString curClass = LOCAL_CREDIT_PORTFOLIO_CLASS;
    CCString stringId = GetLastCurCellEnvValue();

    if (!stringId)
    {
        retCode = ICMLOCAL_Portfolio(V_LSecuritiesID,
									 CashFlowId,
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
        retCode = ICMLOCAL_Portfolio(V_LSecuritiesID,
									CashFlowId,
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

			retCode = ICMLOCAL_Portfolio(V_LSecuritiesID,
										CashFlowId,
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

//    ARM_END();
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


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_LegBasket (LPXLOPER XL_SecuritiesID,
															LPXLOPER XL_CashFlowID)
{
//    ARM_BEGIN();

    // return
    static XLOPER XL_result;
    ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

    // C variable
	int i = 0;
	double C_NbOfCurve =0.;

	VECTOR<CCString> V_SecuritiesID;
	VECTOR<long>	V_LSecuritiesID;

	long CashFlowId = 0;
	CCString C_CashFlowId;

	double size = 0.;

    // error
    static int error;
    static char* reason = "";

    XL_readStrVector(XL_SecuritiesID, V_SecuritiesID ," ARM_ERR: Securities : object vector expected",DOUBLE_TYPE,C_result);
	XL_readStrCellWD(XL_CashFlowID,C_CashFlowId,"NONE"," ARM_ERR: frequency: string expected",C_result);

	C_NbOfCurve = V_SecuritiesID.size();

    for (i = 0; i < C_NbOfCurve; i++)
    {
        V_LSecuritiesID.push_back (LocalGetNumObjectId (V_SecuritiesID[i]));
    }

	CashFlowId = LocalGetNumObjectId (C_CashFlowId);

	long retCode;
    long objId;
    CCString prevClass;
    
    CCString curClass = LOCAL_CREDIT_PORTFOLIO_CLASS;
    CCString stringId = GetLastCurCellEnvValue();

    if (!stringId)
    {
        retCode = ICMLOCAL_Leg_Basket(V_LSecuritiesID,
									 CashFlowId,
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
        retCode = ICMLOCAL_Leg_Basket(V_LSecuritiesID,
									CashFlowId,
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

			retCode = ICMLOCAL_Leg_Basket(V_LSecuritiesID,
										CashFlowId,
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

//    ARM_END();
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

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Collateral (LPXLOPER XL_Labels,
															 LPXLOPER XL_Notionals
															 // LPXLOPER XL_IsInDefault
															 )
{

    // return
    static XLOPER XL_result;
    

	try {
	ARM_NOCALCIFWIZ();

    // C variable

    // error

	ARM_Vector V_Notionals; ExcelTools::convert(XL_Notionals,V_Notionals); 
	std::vector<string> Labels;   ExcelTools::convert(XL_Labels,Labels); 
	// std::vector<string> def_isindefault(Labels.size(),"N"); 
	// std::vector<string> IsIndefault_ ; ExcelTools::convert(XL_IsInDefault,def_isindefault,IsIndefault_); 
	
	// for (i = 0; i < IsIndefault_.size(); i++)
    // {
    //    if (IsIndefault_[i]=="Y") Vb_IsInDefault[i]=true;
	//    else Vb_IsInDefault[i]=false;
    // }


	long prevId = ExcelCaller::get().getObjectId() ;
	long newId = ICMLOCAL_Collateral(Labels,
									  V_Notionals,
									  // Vb_IsInDefault,
									  prevId);

		std::string label = ExcelCaller::get().setObject(newId,LOCAL_CREDIT_PORTFOLIO_CLASS); 
		ExcelTools::convert(label,&XL_result); 
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

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_VariableCollateral (LPXLOPER XL_Labels,
															 LPXLOPER XL_NotionalsId)
{

    // return
    static XLOPER XL_result;
    
	try {
		ARM_NOCALCIFWIZ();

		std::vector<std::string> notionalsStringIds; ExcelTools::convert(XL_NotionalsId,notionalsStringIds); 
		std::vector<string> Labels;   ExcelTools::convert(XL_Labels,Labels); 

		ICM_Vector<long> notionalIds(notionalsStringIds.size()) ; 
		for(unsigned int i=0;i<notionalsStringIds.size();i++) 
			notionalIds[i]=LocalPersistent::get().getObjectId(notionalsStringIds[i]); 

		long prevId = ExcelCaller::get().getObjectId() ;
		long newId = ICMLOCAL_VariableCollateral(Labels,
										  notionalIds,
										  prevId);

		std::string label = ExcelCaller::get().setObject(newId,LOCAL_CREDIT_PORTFOLIO_CLASS); 
		ExcelTools::convert(label,&XL_result); 
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

