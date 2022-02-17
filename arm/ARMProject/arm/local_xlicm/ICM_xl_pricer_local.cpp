
#pragma warning(disable :4005 4786 )

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include <ARM\libicm_local\ICM_local_pricer.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"
#include "ExcelTools.h"


_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Pricer(LPXLOPER XL_Security, 
													   LPXLOPER XL_Model,
													   LPXLOPER XL_PricerType,
													   LPXLOPER XL_NbPaths,
													   LPXLOPER XL_Params,
													   LPXLOPER XL_AsOf)
{
	// return
	static XLOPER XL_result;

	try {
		ARM_NOCALCIFWIZ();

		std::string  C_Security; ExcelTools::convert(XL_Security,C_Security); 
		std::string  C_Model;ExcelTools::convert(XL_Model,C_Model); 
		std::string  C_PricerType;ExcelTools::convert(XL_PricerType,"UNKNOWN",C_PricerType); 
		int nbPaths ; ExcelTools::convert(XL_NbPaths,-999,nbPaths); 
		std::string  C_Params;ExcelTools::convert(XL_Params,"",C_Params); 
		
		double dAsOf ; ExcelTools::convert(XL_AsOf,-1.,dAsOf); 
		long pricerType = ARM_ConvCreditPricerType(C_PricerType); 

		long prevId = ExcelCaller::get().getObjectId(); 
		long newId = ICMLOCAL_Pricer(
			dAsOf==-1?(ARM_Date*)0:&ARM_Date(XLDateToJulian(dAsOf)),
			LocalPersistent::get().getObjectId(C_Security),
			LocalPersistent::get().getObjectId(C_Model),
			pricerType,
			nbPaths,
			LocalPersistent::get().getObjectId(C_Params),
			prevId); 

		FreeCurCellErr() ;
		std::string objName = ExcelCaller::get().setObject(newId,LOCAL_PRICER_CLASS); 
		//ExcelCaller::get().setObject(objName); 
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

/** 

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Debug_Function (LPXLOPER XL_pricerId, 
																 LPXLOPER XL_Double,
																 LPXLOPER XL_Data)
{
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();


	CCString C_pricerId;
	VECTOR<double> V_Data;
	double C_Double;

	long retCode = 0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pricerId,C_pricerId," ARM_ERR: Pricer object expected",C_result);
	XL_readNumCell(XL_Double, C_Double," ARM_ERR: Double: Numeric expected",C_result);


	if (XL_Data->xltype != 128) 
	{
		XL_readNumVector(XL_Data, V_Data," ARM_ERR: Data : vector of numerics expected",C_result);
	}
	
	retCode = ICMLOCAL_Debug_Function (LocalGetNumObjectId (C_pricerId), 
										C_Double,
										V_Data,
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
		V_Data.clear(); 
		V_Data.empty();

	}

	V_Data.clear(); 
	V_Data.empty();
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

_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_PricerDefaultCdsNew(LPXLOPER XL_Security ,
																LPXLOPER XL_Model,
																LPXLOPER XL_AsOf)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_Security; 
	CCString C_Model;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_Security,C_Security," ARM_ERR: Security Id: Object expected",C_result);
	XL_readStrCell(XL_Model,C_Model," ARM_ERR: model id: Object expected",C_result);
	double dAsOf; ExcelTools::convert(XL_AsOf,-1.,dAsOf); 
	
	long prevId = ExcelCaller::get().getObjectId(); 
	long newId = ICMLOCAL_PricerDefaultCdsNew(
		dAsOf==-1?(ARM_Date*)0:&ARM_Date(XLDateToJulian(dAsOf)),
		LocalGetNumObjectId (C_Security),
		LocalGetNumObjectId (C_Model),
		prevId); 
	std::string newName = ExcelCaller::get().setObject(newId,LOCAL_PRICER_DEF_CLASS); 
	ExcelTools::convert(newName,&XL_result); 

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




_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_SetVolatility(LPXLOPER XL_PricerId,
															  LPXLOPER XL_VolCurve)
{
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pricerId;
	CCString C_VolCurve;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_PricerId,C_pricerId," ARM_ERR: pricer id: object expected",C_result);
	XL_readStrCellWD(XL_VolCurve,C_VolCurve,""," ARM_ERR: curve id: Volatility object expected",C_result);
	

	long retCode = ICMLOCAL_SetVolatility(LocalGetNumObjectId (C_pricerId),
										  LocalGetNumObjectId (C_VolCurve),
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

// Is it only a temporary code?
_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Pricer_GetData(
														LPXLOPER XL_PricerId,
														LPXLOPER XL_DataLabel
														)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	/// to make the infrastructure very robust we put a try catch!
	ARM_NOCALCIFWIZ();

	// C variable
	CCString	C_PricerId;
	CCString	C_DataLabel;

	// error
	static int error;
	static char* reason = "";

	// Credit Manager
	XL_readStrCell(XL_PricerId,C_PricerId," ARM_ERR: Pricer Id: object expected",C_result);

	// Data Label
	XL_readStrCellWD(XL_DataLabel, C_DataLabel, "NPV"," ARM_ERR: data label: string expected", C_result);

	long retCode;

	retCode = ICMLOCAL_GetPricer_DataFromLabel(
									LocalGetNumObjectId (C_PricerId),
									C_DataLabel,
									C_result);

	// try to discriminate
	if (retCode == ARM_OK)
	{
		// just a double, waiting for some arrays
		FreeCurCellErr();

		// discrimate between STRING or NUMERICS DATA from LABEL
		string	MyString(C_DataLabel);

		if (MyString.find("STR_") == string::npos)
		{
			// NUMERICS
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
		}
		else
		{
			// STRING
			XL_result.xltype	=	xltypeStr;
			XL_result.val.str	=	XL_StrC2StrPascal(C_result.getMsg());
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

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
