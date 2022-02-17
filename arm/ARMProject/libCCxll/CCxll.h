/* -------------------------------------------------------------------------
 
   File: %M% 
   Path: %P% 
   Description: Interface de lecture/ecriture des donnees Excel (type LPXLOPER)
   Created: 00/01/03
   Author: Charles-Emmanuel MUSY
   Modified: %E% %U% 
   Last maintained by: Charles-Emmanuel MUSY
   Revision: %I% 
 
   -------------------------------------------------------------------------
 
   Note:
 
   ------------------------------------------------------------------------- */

#ifndef CCXLL_H
#define CCXLL_H

#include <CCcommon.h>
SCCS_ID (XL_ccxll_h_SccsId, "%W%, modified %E%");

#include <windows.h>

#ifdef __cplusplus
extern "C" {
#endif
#include "xlcall.h"
#include <framewrk.h>
#ifdef __cplusplus
}
#endif

#include <CCString.h>

#ifdef STL_WIN32
#include <vector>
#define VECTOR std::vector
#else	// STL_WIN32
#pragma warning(disable:4786)
#include <vector.h>
#define VECTOR  vector
#endif	// STL_WIN32

#define XL_NO_ERROR			-1	
#define XL_ERROR			1

#define LONG_TYPE			0
#define DOUBLE_TYPE			1

#define XL_TYPE_DOUBLE			0
#define XL_TYPE_LONG			1
#define XL_TYPE_STRING			2
#define XL_TYPE_DOUBLE_ARRAY	3
#define XL_TYPE_LONG_ARRAY		4
#define XL_TYPE_STRING_ARRAY	5

#define XL_TYPE_DEFAULT_VALUE   -999.
#define CLASS_NAME_BUFFER		50

CCEXTERN_FUNCTION (BOOL GetHwnd, (HWND*));
extern CCString XL_getCaller ();
extern CCString XL_getActiveCell ();
extern void XL_getActiveCellContent (LPXLOPER xRef);
CCEXTERN_FUNCTION (BOOL XL_IsCalledByFuncWiz, (void));
extern int IsVBA();

CCEXTERN_FUNCTION (int XL_getNumCell, (LPXLOPER, double*, double* defaultValue = NULL));
CCEXTERN_FUNCTION (int XL_getNumVector, (LPXLOPER xlTab, VECTOR<double>&, VECTOR<double>* = NULL));
CCEXTERN_FUNCTION (int XL_getNumVectorWithHole, (LPXLOPER xlTab, VECTOR<double>& vec, VECTOR<double>* vecDefault = NULL));
CCEXTERN_FUNCTION (int XL_getNumVectorAndSize, (LPXLOPER xlTab, long&, long&, VECTOR<double>&));
CCEXTERN_FUNCTION (int XL_getNumVectorVector, (LPXLOPER xlTab, VECTOR< VECTOR<double> >&));
CCEXTERN_FUNCTION (int XL_getNumArray, (LPXLOPER, double**, long*));

extern int XL_getStrCell (LPXLOPER, int*, char**, char*, char*, char*, CCString&, const char* defaultValue = NULL, int long_or_float = LONG_TYPE);
CCEXTERN_FUNCTION (char** XL_getStrArray, (LPXLOPER, long*, int*, char**, char*, char*, char*, int long_or_float = LONG_TYPE));
CCEXTERN_FUNCTION (int XL_getStrVector,	  (LPXLOPER, char**, char*, char*, char*, VECTOR<CCString>&, VECTOR<CCString>* = NULL, int long_or_float = LONG_TYPE));
CCEXTERN_FUNCTION (int XL_getStrVectorWD, (LPXLOPER, char**, char*, char*, char*, VECTOR<CCString>&, VECTOR<CCString>* = NULL, int long_or_float = LONG_TYPE, const char* = ""));
CCEXTERN_FUNCTION (int XL_getStrVectorAndSize,	(LPXLOPER, char**, char*, char*, char*, long&, long&, VECTOR<CCString>&, VECTOR<CCString>* = NULL, int long_or_float = LONG_TYPE));
CCEXTERN_FUNCTION (int XL_getStrVectorAndSizeWD,(LPXLOPER, char**, char*, char*, char*, long&, long&, VECTOR<CCString>&, VECTOR<CCString>* = NULL, int long_or_float = LONG_TYPE, const char* = ""));
CCEXTERN_FUNCTION (int XL_getStrVectorSizeAndType, (LPXLOPER, char**, char*, char*, char*, long&, long&, VECTOR<CCString>&, VECTOR<CCString>*, VECTOR<long>&, VECTOR<long>*, int long_or_float = DOUBLE_TYPE));

extern CCString XL_StrPascal2StrC (LPSTR);
CCEXTERN_FUNCTION (LPSTR XL_StrC2StrPascal, (const CCString&));

CCEXTERN_FUNCTION (int XL_Coordonnate2Rank, (int, int, int));

// Set functions
CCEXTERN_FUNCTION (int XL_setNumCell, ( XLOPER&, double, int additionalLinesNb = 0) );
CCEXTERN_FUNCTION (int XL_setNumVector, ( XLOPER& , const VECTOR<double>&, int additionalLinesNb = 0, bool fillWithBlank = true, bool filterSpecificValue=false, double specificValue =-1 ) );
CCEXTERN_FUNCTION (int XL_setStrVector, ( XLOPER& XL_result, const VECTOR<CCString>& vecStr, int additionalLinesNb = 0,bool fillWithBlank = true ) );
CCEXTERN_FUNCTION (int XL_setStrLineVector, ( XLOPER& XL_result, const VECTOR<CCString>& vecStr, int additionalLinesNb = 0,bool fillWithBlank = true ) );
CCEXTERN_FUNCTION (int XL_setStrAndNumVector, ( XLOPER& XL_result, const VECTOR<CCString>& , const VECTOR<double>&, int additionalLinesNb = 0,bool fillWithBlank = true) );
CCEXTERN_FUNCTION (int XL_setStrMatrixSizeAndType, ( XLOPER&, const VECTOR<CCString>&, const VECTOR<long>&, long , long) );
CCEXTERN_FUNCTION (int XL_setStrMatrixSizeAndTypeWithOptions, ( XLOPER&, const VECTOR<CCString>&, const VECTOR<long>&, long , long ,int additionalLinesNb = 0,bool fillWithBlank = true ) );
CCEXTERN_FUNCTION (int XL_setNumMatrixSize, ( XLOPER&, const VECTOR<double>&, long , long) );
CCEXTERN_FUNCTION (int XL_setNumMatrixSizeWithOptions, ( XLOPER&, const VECTOR<double>&, long , long,int additionalLinesNb = 0,bool fillWithBlank = true) );

#define XL_ERROR_VALUE_INCONSISTENCY	" #ERR 100 valeur inconsistante"
#define XL_ERROR_VALUE_MISSING			" #ERR 103 la valeur n'est pas renseignee"
#define XL_ERROR_VALUE_IS_ERR			" #ERR 104 la valeur passee est un code d'erreur"

#define ARM_ARG_ERR()				XL_result.xltype = xltypeStr;\
									XL_result.val.str = "\007ARM_ERR";\
									SetCurCellErrValue (C_result.getMsg ());

#define ARM_ERR_AND_EXIT( error_msg, result) \
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER)&XL_result;\
	}


#define XL_readNumCell(XL_var,C_var,error_msg,result)	\
	if((error = XL_getNumCell (XL_var, &C_var)) != XL_NO_ERROR)\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER)&XL_result;\
	}

										
#define XL_readNumCellWD(XL_var,C_var,dvalue,error_msg,result)	\
	if((error = XL_getNumCell (XL_var, &C_var, &dvalue)) != XL_NO_ERROR)\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER)&XL_result;\
	}

#define XL_readNumArray(XL_array,C_array,C_size,error_msg,result)	\
	if((error = XL_getNumArray (XL_array, &C_array, &C_size)) != XL_NO_ERROR)\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER)&XL_result;\
	}

#define XL_readNumVector(XL_array,vec,error_msg,result)	\
	if((error = XL_getNumVector (XL_array, vec)) != XL_NO_ERROR)\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER)&XL_result;\
	}

#define XL_readNumVectorWithHole(XL_array,vec,error_msg,result)	\
	if((error = XL_getNumVectorWithHole (XL_array, vec)) != XL_NO_ERROR)\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER)&XL_result;\
	}

#define XL_readNumVectorOrNum(XL_array,vec,error_msg,result)	\
	if((error = XL_getNumVector(XL_array, vec)) != XL_NO_ERROR)\
	{\
		double C_var;\
		XL_readNumCell(XL_array,C_var,error_msg,result); \
		vec.push_back(C_var); \
	}


#define XL_readNumVectorWD(XL_array,vec,vecDefault,error_msg,result)	\
	if((error = XL_getNumVector (XL_array, vec, &vecDefault)) != XL_NO_ERROR)\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER)&XL_result;\
	}

#define XL_readNumVectorWithHoleWD(XL_array,vec,vecDefault,error_msg,result)	\
	if((error = XL_getNumVectorWithHole (XL_array, vec, &vecDefault)) != XL_NO_ERROR)\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER)&XL_result;\
	}

#define XL_readNumVectorAndSize(XL_array,nbrows,nbcolumns,vec,error_msg,result)		\
	if((error = XL_getNumVectorAndSize (XL_array,nbrows,nbcolumns, vec)) != XL_NO_ERROR)\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER)&XL_result;\
	}

#define XL_readNumVectorVector(XL_array,vecvec,error_msg,result)		\
	if((error = XL_getNumVectorVector (XL_array,vecvec)) != XL_NO_ERROR)\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER)&XL_result;\
	}

#define XL_readStrOrNumVector(XL_var,C_Str,vec,type,error_msg,result) \
	if((error = XL_getNumVector(XL_var, vec)) != XL_NO_ERROR) \
	{\
		XL_getStrCell (XL_var, &error, &reason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, C_Str);\
		if ( error != XL_NO_ERROR )\
		{\
	  	    if ( error != xlerrNull )\
			{\
				CCString local_msg (error_msg);\
				result.setMsg(local_msg);\
				ARM_ARG_ERR();\
				return((LPXLOPER) &XL_result);\
		  	}\
		}\
		else\
		{\
			type = XL_TYPE_STRING;\
		}\
	}\
	else\
	{\
		type = XL_TYPE_DOUBLE_ARRAY;\
		C_Str.Set ("");\
	}

#define XL_readStrOrNumVectorWD(XL_var,C_Str,vec,vecDefault,type,error_msg,result) \
	if((error = XL_getNumVector (XL_var, vec, &vecDefault)) != XL_NO_ERROR)\
	{\
		XL_getStrCell (XL_var, &error, &reason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, C_Str);\
		if ( error != XL_NO_ERROR )\
		{\
	  	    if ( error != xlerrNull )\
			{\
				CCString local_msg (error_msg);\
				result.setMsg(local_msg);\
				ARM_ARG_ERR();\
				return((LPXLOPER) &XL_result);\
		  	}\
		}\
		else\
		{\
			type = XL_TYPE_STRING;\
		}\
	}\
	else\
	{\
		type = XL_TYPE_DOUBLE_ARRAY;\
		C_Str.Set ("");\
	}

#define XL_readNumOrStrVector(XL_var,C_Num,vec,type,error_msg,result)		\
	if((error = XL_getNumCell (XL_var, &C_Num)) != XL_NO_ERROR)\
	{\
		error = XL_getStrVector (XL_var, &reason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, vec, NULL, DOUBLE_TYPE); \
		if ( error != XL_NO_ERROR )\
		{\
	  	    if ( error != xlerrNull )\
			{\
				CCString local_msg (error_msg);\
				result.setMsg(local_msg);\
				ARM_ARG_ERR();\
				return((LPXLOPER) &XL_result);\
		  	}\
		}\
		else\
		{\
			type = XL_TYPE_STRING_ARRAY;\
		}\
	}\
	else\
	{\
		type = XL_TYPE_DOUBLE;\
	}


#define XL_readStrCell(XL_Str,C_Str,error_msg,result)		\
	XL_getStrCell (XL_Str, &error, &reason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, C_Str);\
	if(error != XL_NO_ERROR)\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER)&XL_result;\
	}

#define XL_readStrCellWD(XL_Str,C_Str,dvalue,error_msg,result)		\
	XL_getStrCell (XL_Str, &error, &reason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, C_Str, dvalue);\
    if(error != XL_NO_ERROR)\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER)&XL_result;\
	}

#define XL_readStrArray(XL_array,C_array,C_size,error_msg,convert,result)		\
	C_array = XL_getStrArray (XL_array, &C_size, &error, &reason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, convert);\
	if(error != XL_NO_ERROR)\
	{\
		if(error != xlerrNull)\
		{\
			CCString local_msg (error_msg);\
			result.setMsg (local_msg);\
			ARM_ARG_ERR();\
			return (LPXLOPER)&XL_result;\
		}\
	}\

#define XL_readStrVector(XL_array,vec,error_msg,convert,result)		\
	if((error = XL_getStrVector (XL_array, &reason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, vec, NULL, convert)) != XL_NO_ERROR)\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER)&XL_result;\
	}

#define XL_readStrVectorWD(XL_array,vec,vecDefault,error_msg,convert,result)		\
	if((error = XL_getStrVector (XL_array, &reason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, vec, &vecDefault, convert)) != XL_NO_ERROR)\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER)&XL_result;\
	}

#define XL_readStrVectorAndSize(XL_array,nbrows,nbcolumns,vec,error_msg,convert,result)		\
	if((error = XL_getStrVectorAndSize (XL_array, &reason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, nbrows, nbcolumns, vec, NULL, convert )) != XL_NO_ERROR)\
	{\
			CCString local_msg (error_msg);\
			result.setMsg (local_msg);\
			ARM_ARG_ERR();\
			return (LPXLOPER)&XL_result;\
	}

#define XL_readStrVectorAndSizeWD(XL_array,nbrows,nbcolumns,vec,vecDefault,error_msg,convert,result)		\
	if((error = XL_getStrVectorAndSizeWD(XL_array, &reason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, nbrows, nbcolumns, vec, &vecDefault, convert, "" )) != XL_NO_ERROR)\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER)&XL_result;\
	}


#define XL_readStrVectorSizeAndType(XL_array,nbrows,nbcolumns,vec,type,error_msg,convert,result)		\
	if((error = XL_getStrVectorSizeAndType(XL_array,&reason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, nbrows, nbcolumns, vec, NULL, type, NULL, convert )) != XL_NO_ERROR)\
	{\
			CCString local_msg(error_msg);\
			local_msg = local_msg + CCString(reason); \
			result.setMsg (local_msg);\
			ARM_ARG_ERR();\
			return (LPXLOPER)&XL_result;\
	}

#define XL_readStrVectorSizeAndTypeWD(XL_array,nbrows,nbcolumns,vec,vecDefault,type,typeDefault,error_msg,convert,result)		\
	if((error = XL_getStrVectorSizeAndType(XL_array,&reason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, nbrows, nbcolumns, vec, vecDefault, type, typeDefault, convert )) != XL_NO_ERROR)\
	{\
			CCString local_msg (error_msg);\
			result.setMsg (local_msg);\
			ARM_ARG_ERR();\
			return (LPXLOPER)&XL_result;\
	}


#define XL_readStrOrNumCell(XL_var,C_Str,C_var,type,error_msg,result)	\
	if((error = XL_getNumCell (XL_var, &C_var)) != XL_NO_ERROR)\
	{\
		XL_getStrCell (XL_var, &error, &reason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, C_Str);\
		if(error != XL_NO_ERROR)\
		{\
	  		if(error != xlerrNull)\
			{\
				CCString local_msg (error_msg);\
				result.setMsg (local_msg);\
				ARM_ARG_ERR();\
				return (LPXLOPER)&XL_result;\
		  	}\
		}\
		else\
		{\
			type = XL_TYPE_STRING;\
			C_var = 0.0;\
		}\
	}\
	else\
	{\
		type = XL_TYPE_DOUBLE;\
		C_Str.Set ("");\
	}

#define XL_readStrOrNumCellWD(XL_var,C_Str,C_var,dvalue,type,error_msg,result)		\
	if((error = XL_getNumCell (XL_var, &C_var, &dvalue)) != XL_NO_ERROR)\
	{\
		XL_getStrCell (XL_var, &error, &reason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, C_Str);\
		if(error != XL_NO_ERROR)\
		{\
	  		if(error != xlerrNull)\
			{\
				CCString local_msg (error_msg);\
				result.setMsg (local_msg);\
				ARM_ARG_ERR();\
				return (LPXLOPER)&XL_result;\
		  	}\
		}\
		else\
		{\
			type = XL_TYPE_STRING;\
			C_var = 0.0;\
		}\
	}\
	else\
	{\
		type = XL_TYPE_DOUBLE;\
		C_Str.Set ("");\
	}


/// function to specify the size
#define XL_writeNumCell( XL_var, valuedble, error_msg, result )	\
	if(  (error=XL_setNumCell( XL_var, valuedble )) != XL_NO_ERROR )\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER) &XL_var;\
	}

/// Various functions to write vector or matrix
/// options is to add additional lines
/// and to fill these additional lines with blank string
/// therefore 
///		- additionalLinesNb is an integer
///		- fillWithBlank is a boolean (true to get blank string)

/// write a column of numeric
#define XL_writeNumVectorWithOptions( XL_var, vec, error_msg, result, additionalLinesNb, fillWithBlank )	\
	if( (XL_setNumVector( XL_var, vec, additionalLinesNb,fillWithBlank)) != XL_NO_ERROR )\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER) &XL_var;\
	}

/// function to filter specific value on top
#define XL_writeNumVectorWithOptionsAndFilter( XL_var, vec, error_msg, result, additionalLinesNb, fillWithBlank, filterSpecificValue, specificValue )	\
	if( (XL_setNumVector( XL_var, vec, additionalLinesNb,fillWithBlank,filterSpecificValue,specificValue)) != XL_NO_ERROR )\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER) &XL_var;\
	}

/// very simple interface
#define XL_writeNumVector( XL_var, vec, error_msg, result ) \
	XL_writeNumVectorWithOptions( XL_var, vec, error_msg, result, 0, false )


/// write one column of string
#define XL_writeStrVectorWithOptions( XL_var, vecStr, error_msg, result, additionalLinesNb, fillWithBlank )	\
	if( (XL_setStrVector( XL_var, vecStr, additionalLinesNb,fillWithBlank)) != XL_NO_ERROR )\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER) &XL_var;\
	}
#define XL_writeStrVector( XL_var, vecStr, error_msg, result )	\
	XL_writeStrVectorWithOptions( XL_var, vecStr, error_msg, result, 0, false )	

/// write line column of string
#define XL_writeStrLineVectorWithOptions( XL_var, vecStr, error_msg, result, additionalLinesNb, fillWithBlank )	\
	if( (XL_setStrLineVector( XL_var, vecStr, additionalLinesNb,fillWithBlank)) != XL_NO_ERROR )\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER) &XL_var;\
	}


/// write on two columns with first column being string the other one being numeric
#define XL_writeStrAndNumVectorWithOptions( XL_var, vecString, vecDble, error_msg, result, additionalLinesNb, fillWithBlank )	\
	if((error = XL_setStrAndNumVector( XL_var, vecString, vecDble, additionalLinesNb, fillWithBlank) ) != XL_NO_ERROR )\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER) &XL_var;\
	}
#define XL_writeStrAndNumVector( XL_var, vecString, vecDble, error_msg, result ) \
	XL_writeStrAndNumVector( XL_var, vecString, vecDble, error_msg, result, 0, false )	


#define XL_writeStrMatrixSizeAndType( XL_var, values, types, nbRows, nbCols, error_msg, result ) \
	if( (XL_setStrMatrixSizeAndType( XL_var, values, types, nbRows, nbCols)) != XL_NO_ERROR )\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER) &XL_var;\
	}

#define XL_writeNumMatrixSize( XL_var, values, nbRows, nbCols, error_msg, result ) \
	if( (XL_setNumMatrixSize( XL_var, values, nbRows, nbCols)) != XL_NO_ERROR )\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER) &XL_var;\
	}
#define XL_writeNumMatrixSizeWithOptions( XL_var, values, nbRows, nbCols, error_msg, result,additionalLinesNb, fillWithBlank ) \
	if( (XL_setNumMatrixSizeWithOptions( XL_var, values, nbRows, nbCols,additionalLinesNb, fillWithBlank)) != XL_NO_ERROR )\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER) &XL_var;\
	}



#define XL_writeStrMatrixSizeAndTypeWithOptions( XL_var, values, types, nbRows, nbCols, error_msg, result,additionalLinesNb, fillWithBlank ) \
	if( (XL_setStrMatrixSizeAndTypeWithOptions( XL_var, values, types, nbRows, nbCols,additionalLinesNb, fillWithBlank)) != XL_NO_ERROR )\
	{\
		CCString local_msg (error_msg);\
		result.setMsg (local_msg);\
		ARM_ARG_ERR();\
		return (LPXLOPER) &XL_var;\
	}

/// class to hold an xloper and make sure it releases the memory!
class XLOPER_Holder
{
public:
	static void FreeStringVec( XLOPER& XL_result );
	static int FreeXLOPER_Multi( XLOPER& XL_result );
	XLOPER_Holder();
	~XLOPER_Holder();
	XLOPER& GetResult() { return itsResult;}
private:
	XLOPER itsResult;
};


#ifdef CCxll_cc

#include <stdio.h>

#include <winuser.h>

#include <CCmessage.h>

#define XL_CALLER_SIZE				256
#define CLASS_NAME_BUFFER			50

typedef struct _EnumStruct
{
	HWND			hwnd;
  	unsigned short	wLoword;
} EnumStruct;

typedef struct _EnumStructWiz
{
	BOOL bFuncWiz;
	short hwndXLMain;
} EnumStructWiz, FAR* LPEnumStructWiz;

#endif	// CCxll_cc

#endif	// CCXLL_H

// EOF %M%
