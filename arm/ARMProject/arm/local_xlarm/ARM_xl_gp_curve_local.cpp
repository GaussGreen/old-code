/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_curvve_local.cpp,v $
 * Revision 1.1  2004/02/07 15:08:43  ebenhamou
 * Initial version
 *
 */

#include <ARM\libarm_local\firstToBeIncluded.h>
#include "ARM_xl_gp_curve_local.h"
#include <ARM\libarm_local\ARM_local_gp_curve.h>
#include <libCCxll\CCxll.h>
#include "ARM_xl_wrapper_local.h"
#include "ARM_xl_trycatch_local.h"
#include "ARM_xl_gp_fctorhelper.h"
#include "ARM_local_interface.h"
#include <GP_Base\gpbase\vectormanip.h>
#include <GP_Base\gpbase\stringmanip.h>

#include "util\tech_macro.h"

////////////////////////////////////////////////
///  THIS INTERFACE IS MORE FACTORISED THAN THE
///  REST OF ARM .... IN ORDER TO RECODE
///  ONLY THE PART TO DEFINE A FUNCTION CALL
///  WE PROVIDE FOR EACH A FUNCTOR
///  WE ALSO USE THE SAME API FOR PERSISTENT 
///  AND NON PERSISTENT ADDIN
///  BY COMMON DECISION.... WE KEEP THIS METHOD
///  TO THIS FILES AND DO NOT SPREAD ELSEWHERE
////////////////////////////////////////////////

///----------------------------------------------
///----------------------------------------------
///             curve construction
/// Inputs :
///     vector of abscisses
///     vector of ordinates
///		nb of rows and cols
///		sort abscisses
///		interpolator name
///----------------------------------------------
///----------------------------------------------


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_GenericCurve_Create_Common(
	LPXLOPER XL_newAbscisses,
	LPXLOPER XL_ordinates,
	LPXLOPER XL_interpolatorType,
	LPXLOPER XL_sortAbscisses,
	LPXLOPER XL_alwaysMulti,
	bool PersistentInXL )
{	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		VECTOR<double> C_newAbscisses;
		XL_readNumVector(XL_newAbscisses,C_newAbscisses," ARM_ERR: abscisses: array of numeric expected",C_result);

		VECTOR<double> C_ordinates;
		long C_rowsNb, C_colsNb;
		XL_readNumVectorAndSize( XL_ordinates, C_rowsNb, C_colsNb, C_ordinates, " ARM_ERR: ordinates vector or matrix expected", C_result);

		CCString C_interpolatorType;
		const char interpolatorTypeDefault[7]("LINEAR");
		XL_readStrCellWD( XL_interpolatorType,	C_interpolatorType,	interpolatorTypeDefault, " ARM_ERR: interpolator type string expected",	C_result);
		string interpolatorName = CCSTringToSTLString( C_interpolatorType );
		ARM::stringToUpper( interpolatorName );

		double C_sortAbscisses;
		double sortAbscissesDefault = 0;
		XL_readNumCellWD( XL_sortAbscisses, C_sortAbscisses, sortAbscissesDefault,	" ARM_ERR: sort abscisses: boolean compatible expected",	C_result);
		bool sortAbscisses = C_sortAbscisses != 0;

		double C_alwaysMulti;
		double alwaysMultiDefault = 0;
		XL_readNumCellWD( XL_alwaysMulti, C_alwaysMulti, alwaysMultiDefault,	" ARM_ERR: always multi: boolean compatible expected",	C_result);
		bool alwaysMulti = C_alwaysMulti != 0;

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc7Args< vector<double>,vector<double>,long,long,bool,string,bool >  ourFunc(
			C_newAbscisses,
			C_ordinates,
			C_rowsNb,
			C_colsNb,
			sortAbscisses,
			interpolatorName,
			alwaysMulti,
			ARMLOCAL_GenericCurve_Create );
		
		/// call the general function
		fillXL_Result( LOCAL_GENERICCURVE_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GenericCurve_Create" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenericCurve_Create(
	LPXLOPER XL_newAbscisses,
	LPXLOPER XL_ordinates,
	LPXLOPER XL_interpolatorType,
	LPXLOPER XL_sortAbscisses,
	LPXLOPER XL_alwaysMulti)
{
	ADD_LOG("Local_GenericCurve_Create");
	bool PersistentInXL = true;
	return Local_GenericCurve_Create_Common(
		XL_newAbscisses,
		XL_ordinates,
		XL_interpolatorType,
		XL_sortAbscisses,
		XL_alwaysMulti,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenericCurve_Create(
	LPXLOPER XL_newAbscisses,
	LPXLOPER XL_ordinates,
	LPXLOPER XL_interpolatorType,
	LPXLOPER XL_sortAbscisses,
	LPXLOPER XL_alwaysMulti )
{
	ADD_LOG("Local_PXL_GenericCurve_Create");
	bool PersistentInXL = false;
	return Local_GenericCurve_Create_Common(
		XL_newAbscisses,
		XL_ordinates,
		XL_interpolatorType,
		XL_sortAbscisses,
		XL_alwaysMulti,
		PersistentInXL );
}



__declspec(dllexport) LPXLOPER WINAPI Local_GenericCurve_Interpolate(
	LPXLOPER XL_CurveId,
	LPXLOPER XL_abscisse )
{
	ADD_LOG("Local_GenericCurve_Interpolate");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;	

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	

		CCString C_genericCurveString;
		XL_readStrCell( XL_CurveId, C_genericCurveString, " ARM_ERR: Curve Id: Object expected",C_result);
		long C_genericCurveId = LocalGetNumObjectId(C_genericCurveString);
		
		double C_Abscisse;
		XL_readNumCell( XL_abscisse, C_Abscisse, " ARM_ERR: abscisse: numeric expected",	C_result);
		
		vector<double> vecResult;
		long retCode = ARMLOCAL_GenericCurve_Interpolate( C_genericCurveId, C_Abscisse, vecResult, C_result );

		if( retCode == ARM_KO )
		{
			ARM_ERR();
		}
		else
		{
			if(vecResult.size() == 1 )
			{
				XL_result.xltype  = xltypeNum;
				XL_result.val.num = vecResult[0];
			}
			else
			{
				/// add these additional lines with blank lines
				const int additionalLinesNb = 100;
				bool fillWithBlank = true;
				XL_writeNumVectorWithOptions( XL_result, vecResult, " ARM_ERR: Could not get argResult array data", C_result, additionalLinesNb, fillWithBlank );
			}
			FreeCurCellErr ();
		}
	}	
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GenericCurve_Interpolate" )

	return (LPXLOPER)&XL_result;
}





///----------------------------------------------
///----------------------------------------------
///             cpt curve 
/// Inputs :
///     id of the previous curve
///		new vector of abscisses
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////

class CptCurveFunc : public ARMResultLong2LongFunc
{
public:
	CptCurveFunc(	long curveId, const vector<double>& newAbscisses )
	:		C_curveId(curveId), C_newAbscisses(newAbscisses)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_GenericCurve_CptCurve(
			C_curveId,
			C_newAbscisses,
            result,
            objId);			
	}
private:
	long C_curveId;
	vector<double> C_newAbscisses;
};



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_GenericCurve_CptCurve_Common(
	LPXLOPER XL_CurveId,
	LPXLOPER XL_newAbscisses,
	bool PersistentInXL )
{	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		CCString C_genericCurveString;
		XL_readStrCell( XL_CurveId, C_genericCurveString, " ARM_ERR: Curve Id: Object expected",C_result);
		long C_genericCurveId = LocalGetNumObjectId(C_genericCurveString);

		VECTOR<double> C_newAbscisses;
		XL_readNumVector(XL_newAbscisses,C_newAbscisses," ARM_ERR: abscisses: array of numeric expected",C_result);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		CptCurveFunc ourFunc(C_genericCurveId,C_newAbscisses);
		
		/// call the general function
		fillXL_Result( LOCAL_GENERICCURVE_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GenericCurve_CptCurve_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenericCurve_CptCurve(
	LPXLOPER XL_CurveId,
	LPXLOPER XL_newAbscisses )
{
	ADD_LOG("Local_GenericCurve_CptCurve");
	bool PersistentInXL = true;
	return Local_GenericCurve_CptCurve_Common(
		XL_CurveId,
		XL_newAbscisses,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenericCurve_CptCurve(
	LPXLOPER XL_CurveId,
	LPXLOPER XL_newAbscisses )
{
	ADD_LOG("Local_PXL_GenericCurve_CptCurve");
	bool PersistentInXL = false;
	return Local_GenericCurve_CptCurve_Common(
		XL_CurveId,
		XL_newAbscisses,
		PersistentInXL );
}



///----------------------------------------------
///----------------------------------------------
///             insert point to a curve
/// Inputs :
///     id of the previous curve
///		new abscisse
///		new ordinate
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////

class InsertCurveFunc : public ARMResultLong2LongFunc
{
public:
	InsertCurveFunc (			long curveId, 
		double					newAbscisse,
		const vector<double>&	newOrdinate )
	:		C_curveId(curveId), 
			C_newAbscisse(newAbscisse),
			C_newOrdinate(newOrdinate)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_GenericCurve_Insert(
			C_curveId,
			C_newAbscisse,
			C_newOrdinate,
            result,
            objId);			
	}
private:
	long C_curveId;
	double C_newAbscisse;
	vector<double> C_newOrdinate;
};

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_GenericCurve_Insert_Common(
	LPXLOPER XL_CurveId,
	LPXLOPER XL_newAbscisse,
	LPXLOPER XL_newOrdinate,
	bool PersistentInXL )
{	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		CCString C_genericCurveString;
		XL_readStrCell( XL_CurveId, C_genericCurveString, " ARM_ERR: Curve Id: Object expected",C_result);
		long C_genericCurveId = LocalGetNumObjectId(C_genericCurveString);

		double C_newAbscisse;
		XL_readNumCell( XL_newAbscisse, C_newAbscisse, " ARM_ERR: new abscisse: numeric expected",	C_result);

		VECTOR<double> C_newOrdinate;
		XL_readNumVectorOrNum(XL_newOrdinate,C_newOrdinate," ARM_ERR: new ordinate: numeric or array of numeric expected",C_result);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		InsertCurveFunc ourFunc(C_genericCurveId,C_newAbscisse,C_newOrdinate);
		
		/// call the general function
		fillXL_Result( LOCAL_GENERICCURVE_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GenericCurve_Insert_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenericCurve_Insert(
	LPXLOPER XL_CurveId,
	LPXLOPER XL_newAbscisse,
	LPXLOPER XL_newOrdinate )
{
	ADD_LOG("Local_GenericCurve_Insert");
	bool PersistentInXL = true;
	return Local_GenericCurve_Insert_Common(
		XL_CurveId,
		XL_newAbscisse,
		XL_newOrdinate,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenericCurve_Insert(
	LPXLOPER XL_CurveId,
	LPXLOPER XL_newAbscisse,
	LPXLOPER XL_newOrdinate )
{
	ADD_LOG("Local_PXL_GenericCurve_Insert");
	bool PersistentInXL = false;
	return Local_GenericCurve_Insert_Common(
		XL_CurveId,
		XL_newAbscisse,
		XL_newOrdinate,
		PersistentInXL );
}


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_FlatSurface_Create_Common(
	LPXLOPER XL_value,
	bool PersistentInXL )
{	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		double C_value;
		XL_readNumCell(	 XL_value,	C_value, " ARM_ERR: value : double expected",				C_result);

		/// a function with a context
		exportFunc1Arg< double  >  ourFunc(C_value, ARMLOCAL_FlatSurface_Create );

		/// call the general function
		fillXL_Result( LOCAL_FLATSURFACE_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CstManager_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_FlatSurface_Create(
	LPXLOPER XL_Value)
{
	ADD_LOG("Local_FlatSurface_Create");
	bool PersistentInXL = true;
	return Local_FlatSurface_Create_Common(
		XL_Value,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FlatSurface_Create(
	LPXLOPER XL_Value)
{
	ADD_LOG("Local_PXL_FlatSurface_Create");
	bool PersistentInXL = false;
	return Local_FlatSurface_Create_Common(
		XL_Value,
		PersistentInXL );
}



/////////////////////////////////////////////////////////////
/// function to create a linear surface
/////////////////////////////////////////////////////////////
LPXLOPER Local_LinSurface_Create_Common(
	LPXLOPER XL_X1,
	LPXLOPER XL_X2,
	LPXLOPER XL_X3,
    LPXLOPER XL_interpolatorType,
	bool PersistentInXL )
{	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		vector<double> C_X1;
		XL_readNumVector(XL_X1,C_X1," ARM_ERR: X1: array of numeric expected",C_result);

		vector<double> C_X2;
		XL_readNumVector(XL_X2,C_X2," ARM_ERR: X2: array of numeric expected",C_result);

		VECTOR<double> C_X3;
		long C_rowsNb, C_colsNb;
		XL_readNumVectorAndSize( XL_X3, C_rowsNb, C_colsNb, C_X3, " ARM_ERR: ordinates vector or matrix expected", C_result);

        CCString C_interpolatorType;

		const char interpolatorTypeDefault[14]("LINEAR_COLUMN");
		XL_readStrCellWD( XL_interpolatorType,	C_interpolatorType,	interpolatorTypeDefault, " ARM_ERR: interpolator type string expected",	C_result);
		string interpolatorName = CCSTringToSTLString( C_interpolatorType );
		ARM::stringToUpper( interpolatorName );

		/// a function with a context
		exportFunc6Args< vector<double >, vector<double>, vector<double>, long, long, string >  ourFunc(C_X1, C_X2, C_X3, C_rowsNb, C_colsNb, interpolatorName,ARMLOCAL_LinSurface_Create );

		/// call the general function
		fillXL_Result( LOCAL_LINSURFACE_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LinSurface_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_LinSurface_Create(
	LPXLOPER XL_X1,
	LPXLOPER XL_X2,    
	LPXLOPER XL_X3,
    LPXLOPER XL_interpolType)
{
	ADD_LOG("Local_LinSurface_Create");
	bool PersistentInXL = true;
	return Local_LinSurface_Create_Common(
		XL_X1,
		XL_X2,
		XL_X3,
        XL_interpolType,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LinSurface_Create(
	LPXLOPER XL_X1,
	LPXLOPER XL_X2,
	LPXLOPER XL_X3,
    LPXLOPER XL_interpolType)
{
	ADD_LOG("Local_PXL_LinSurface_Create");
	bool PersistentInXL = false;
	return Local_LinSurface_Create_Common(
		XL_X1,
		XL_X2,
		XL_X3,
        XL_interpolType,
		PersistentInXL );
}




/////////////////////////////////////////////////////////////
/// central function that interpolates a surface
/////////////////////////////////////////////////////////////
LPXLOPER Local_Surface_Interpolate_Common(
	LPXLOPER XL_surfaceId,
	LPXLOPER XL_X1,
	LPXLOPER XL_X2,
	bool PersistentInXL )
{	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		CCString XL_surfaceIdString;
		XL_readStrCell( XL_surfaceId, XL_surfaceIdString, " ARM_ERR: Curve Id: Object expected",C_result);
		long C_surfaceId = LocalGetNumObjectId(XL_surfaceIdString);

		double C_X1;
		XL_readNumCell(	 XL_X1,	C_X1, " ARM_ERR: X1: double expected",				C_result);

		double C_X2;
		XL_readNumCell(	 XL_X2,	C_X2, " ARM_ERR: X2: double expected",				C_result);

		/// a function with a context
		long retCode = ARMLOCAL_Surface_Interpolate(C_surfaceId, C_X1, C_X2, C_result, -1 );

		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Surface_Interpolate_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Surface_Interpolate(
	LPXLOPER XL_SurfaceId,
	LPXLOPER XL_X2,
	LPXLOPER XL_X3)
{
	ADD_LOG("Local_Surface_Interpolate");
	bool PersistentInXL =true;
	return Local_Surface_Interpolate_Common(
		XL_SurfaceId,
		XL_X2,
		XL_X3,
		PersistentInXL );
}




/////////////////////////////////////////////////////////////
/// function to create a linear surface
/////////////////////////////////////////////////////////////
LPXLOPER Local_Surface_Insert_Common(
	LPXLOPER XL_SurfaceId,
	LPXLOPER XL_X1,
	LPXLOPER XL_X2,
	LPXLOPER XL_X3,
	bool PersistentInXL )
{	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		CCString XL_surfaceIdString;
		XL_readStrCell( XL_SurfaceId, XL_surfaceIdString, " ARM_ERR: Curve Id: Object expected",C_result);
		long C_surfaceId = LocalGetNumObjectId(XL_surfaceIdString);

		double C_X1;
		XL_readNumCell(	 XL_X1,	C_X1, " ARM_ERR: X1: double expected",				C_result);

		double C_X2;
		XL_readNumCell(	 XL_X2,	C_X2, " ARM_ERR: X2: double expected",				C_result);

		double C_X3;
		XL_readNumCell(	 XL_X3,	C_X3, " ARM_ERR: X3: double expected",				C_result);

		/// a function with a context
		exportFunc4Args< long, double, double, double >  ourFunc(C_surfaceId, C_X1, C_X2, C_X3, ARMLOCAL_Surface_Insert );

		/// call the general function
		fillXL_Result( LOCAL_LINSURFACE_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LinSurface_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Surface_Insert(
	LPXLOPER XL_SurfaceId,
	LPXLOPER XL_X1,
	LPXLOPER XL_X2,
	LPXLOPER XL_X3 )
{
	ADD_LOG("Local_Surface_Insert");
	bool PersistentInXL = true;
	return Local_Surface_Insert_Common(
		XL_SurfaceId,
		XL_X1,
		XL_X2,
		XL_X3,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Surface_Insert(
	LPXLOPER XL_SurfaceId,
	LPXLOPER XL_X1,
	LPXLOPER XL_X2,
	LPXLOPER XL_X3 )
{
	ADD_LOG("Local_PXL_Surface_Insert");
	bool PersistentInXL = false;
	return Local_Surface_Insert_Common(
		XL_SurfaceId,
		XL_X1,
		XL_X2,
		XL_X3,
		PersistentInXL );
}



/////////////////////////////////////////////////////////////
/// function to covert vol from summit to linear surface
/////////////////////////////////////////////////////////////
LPXLOPER Local_FromVolSummitToSurfce_Create_Common(
	LPXLOPER XL_VolId,
    LPXLOPER XL_InterpolatorType,
	bool PersistentInXL )
{	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		CCString XL_volIdString;
		XL_readStrCell( XL_VolId, XL_volIdString, " ARM_ERR: Vol Id: Object expected",C_result);
		long C_volId = LocalGetNumObjectId(XL_volIdString);

        CCString C_interpolatorType;
		const char interpolatorTypeDefault[14]("LINEAR_COLUMN");
		XL_readStrCellWD( XL_InterpolatorType,	C_interpolatorType,	interpolatorTypeDefault, " ARM_ERR: interpolator type string expected",	C_result);
		string interpolatorName = CCSTringToSTLString( C_interpolatorType );
		ARM::stringToUpper( interpolatorName );

		/// a function with a context
		exportFunc2Args< long, string >  ourFunc(C_volId, interpolatorName,ARMLOCAL_FromVolSummitToSurfce_Create );

		/// call the general function
		fillXL_Result( LOCAL_LINSURFACE_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LinSurface_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_FromVolSummitToSurfce_Create(
	LPXLOPER XL_VolId,
    LPXLOPER XL_InterpolType)
{
	ADD_LOG("Local_FromVolSummitToSurfce_Create");
	bool PersistentInXL = true;
	return Local_FromVolSummitToSurfce_Create_Common(
		XL_VolId,
        XL_InterpolType,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FromVolSummitToSurfce_Create(
	LPXLOPER XL_VolId,
    LPXLOPER XL_InterpolType)
{
	ADD_LOG("Local_PXL_FromVolSummitToSurfce_Create");
	bool PersistentInXL = false;
	return Local_FromVolSummitToSurfce_Create_Common(
		XL_VolId,
        XL_InterpolType,
		PersistentInXL );
}

/////////////////////////////////////////////////////////////
/// function to covert vol from summit to ARM_Curve
/// Useful for CorrelFXFX curves
/////////////////////////////////////////////////////////////
LPXLOPER Local_FromVolSummitToCurve_Create_Common(
	LPXLOPER XL_VolId,
    LPXLOPER XL_calcMethod,
	bool PersistentInXL )
{	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		CCString XL_volIdString;
		XL_readStrCell( XL_VolId, XL_volIdString, " ARM_ERR: Vol Id: Object expected",C_result);
		long C_volId = LocalGetNumObjectId(XL_volIdString);

		CCString C_calcMethod;
		long calcModId;

		XL_readStrCellWD(XL_calcMethod,C_calcMethod,"LIN"," ARM_ERR: calculation method: string expected",C_result);

		if((calcModId = ARM_ConvCalculationMethod(C_calcMethod, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// a function with a context
		exportFunc2Args< long, long >  ourFunc(C_volId, calcModId, ARMLOCAL_FromVolSummitToCurve_Create);

		/// call the general function
		fillXL_Result( LOCAL_GENERICCURVE_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FromVolSummitToCurve_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_FromVolSummitToCurve_Create(
	LPXLOPER XL_VolId,
    LPXLOPER XL_calcMethod)
{
	ADD_LOG("Local_FromVolSummitToCurve_Create");
	bool PersistentInXL = true;
	return Local_FromVolSummitToCurve_Create_Common(
		XL_VolId,
        XL_calcMethod,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FromVolSummitToCurve_Create(
	LPXLOPER XL_VolId,
    LPXLOPER XL_calcMethod)
{
	ADD_LOG("Local_PXL_FromVolSummitToCurve_Create");
	bool PersistentInXL = false;
	return Local_FromVolSummitToCurve_Create_Common(
		XL_VolId,
        XL_calcMethod,
		PersistentInXL );
}

/////////////////////////////////////////////////////////////
/// function to covert vol from summit to linear surface
/////////////////////////////////////////////////////////////
LPXLOPER Local_CurveMatrix_Create_Common(
	LPXLOPER XL_CorrelsId,
    LPXLOPER XL_NbRows,
	LPXLOPER XL_NbCols,
	bool PersistentInXL )
{	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		VECTOR<CCString> XL_CorrelsIdString;
		XL_readStrVector( XL_CorrelsId, XL_CorrelsIdString, " ARM_ERR: Correls Id: Object expected",DOUBLE_TYPE,C_result);
		VECTOR<long> C_CorrelsId(XL_CorrelsIdString.size());
		
		for (size_t i = 0; i < XL_CorrelsIdString.size(); ++i)
			C_CorrelsId[i] = LocalGetNumObjectId(XL_CorrelsIdString[i]);

        double C_NbRowsDbl;
		XL_readNumCell( XL_NbRows,	C_NbRowsDbl, " ARM_ERR: NbRows numeric expected",	C_result);
		long C_NbRows = C_NbRowsDbl;

		double C_NbColsDbl;
		XL_readNumCell( XL_NbCols,	C_NbColsDbl, " ARM_ERR: NbCols numeric expected",	C_result);
		long C_NbCols = C_NbColsDbl;
		

		/// a function with a context
		exportFunc3Args< VECTOR<long>, long, long >  ourFunc(C_CorrelsId, C_NbRows, C_NbCols, ARMLOCAL_CurveMatrix_Create );

		/// call the general function
		fillXL_Result( LOCAL_LINSURFACE_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LinSurface_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CurveMatrix_Create(
	LPXLOPER XL_CorrelsId,
    LPXLOPER XL_NbRows,
	LPXLOPER XL_NbCols)
{
	ADD_LOG("Local_CurveMatrix_Create");
	bool PersistentInXL = true;
	return Local_CurveMatrix_Create_Common(
		XL_CorrelsId,
        XL_NbRows,
		XL_NbCols,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CurveMatrix_Create(
	LPXLOPER XL_CorrelsId,
    LPXLOPER XL_NbRows,
	LPXLOPER XL_NbCols)
{
	ADD_LOG("Local_PXL_CurveMatrix_Create");
	bool PersistentInXL = false;
	return Local_CurveMatrix_Create_Common(
		XL_CorrelsId,
        XL_NbRows,
		XL_NbCols,
		PersistentInXL );
}


