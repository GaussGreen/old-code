
#pragma warning(disable :4005 4786 )

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_class.h>

#include <ARM\libicm_local\ICM_local_cf.h>

#include "ARM_local_interface.h"
#include "ARM_local_interglob.h"
#include "ICMKernel/util/icm_matrix.h"

#include "ExcelTools.h"



__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CashFlows (LPXLOPER XL_CashFlows)
{
//    ARM_BEGIN();

    // return
    static XLOPER XL_result;
    ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

    // C variable
	int i = 0,j =0;
	double nbcf =0.;
	long nbrows=0,nbcolumns=0;

	VECTOR<CCString> C_matrix;
	C_matrix.clear();

    // error
    static int error;
    static char* reason = "";

    XL_readStrVectorAndSize(XL_CashFlows,nbrows,nbcolumns,C_matrix," ARM_ERR: Cash Flows : object matrix expected",DOUBLE_TYPE,C_result);

	double size = C_matrix.size();

	if(size == 0)
	{
		C_result.setMsg ("ARM_ERR: check your Cash Flows Matrix");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
    long objId;
    CCString prevClass;
    
    CCString curClass = LOCAL_CREDIT_CASHFLOWS_CLASS;
    CCString stringId = GetLastCurCellEnvValue();

    if (!stringId)
    {
        retCode = ICMLOCAL_CashFlows(C_matrix,
									 (int)nbrows,
									 (int)nbcolumns,	
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
        retCode = ICMLOCAL_CashFlows(C_matrix,
									 (int)nbrows,
									 (int)nbcolumns,	
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

			retCode = ICMLOCAL_CashFlows(C_matrix,
										(int)nbrows,
										(int)nbcolumns,	
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
C_matrix.clear();
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



__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Parameters (LPXLOPER XL_CashFlows)
{
//    ARM_BEGIN();

    // return
    static XLOPER XL_result;
    ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

    // C variable
	int i = 0,j =0;
	double nbcf =0.;
	long nbrows=0,nbcolumns=0;

	VECTOR<CCString> C_matrix;
	C_matrix.clear();

    // error
    static int error;
    static char* reason = "";

    XL_readStrVectorAndSize(XL_CashFlows,nbrows,nbcolumns,C_matrix," ARM_ERR: Cash Flows : object matrix expected",DOUBLE_TYPE,C_result);

	double size = C_matrix.size();

	if(size == 0)
	{
		C_result.setMsg ("ARM_ERR: check your Cash Flows Matrix");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
    long objId;
    CCString prevClass;
    
    CCString curClass = LOCAL_CREDIT_CASHFLOWS_CLASS;
    CCString stringId = GetLastCurCellEnvValue();

    if (!stringId)
    {
        retCode = ICMLOCAL_Parameters(C_matrix,
									 (int)nbrows,
									 (int)nbcolumns,	
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
        retCode = ICMLOCAL_Parameters(C_matrix,
									 (int)nbrows,
									 (int)nbcolumns,	
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

			retCode = ICMLOCAL_Parameters(C_matrix,
										(int)nbrows,
										(int)nbcolumns,	
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

	C_matrix.clear();
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


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_PropertyList(LPXLOPER xlAttrNames,
															  LPXLOPER xlAttrValues,
															  LPXLOPER xlAttrTypes)
{

	static XLOPER XL_result;
	try {
		ARM_NOCALCIFWIZ();
		ICM_PropertyList pl ;
		ExcelTools::convert(xlAttrNames,xlAttrValues,xlAttrTypes,pl); 
		long prevId = ExcelCaller::get().getObjectId(); 
		long id = LocalPersistent::get().adopt(pl.Clone(),prevId); 
		std::string objName = ExcelCaller::get().setObject(id,"LPLST"); 
		ExcelTools::convert(objName,&XL_result); 
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





