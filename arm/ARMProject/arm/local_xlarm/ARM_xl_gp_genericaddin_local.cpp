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
#include "ARM\libarm_local\ARM_local_gp_genericaddin.h"
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"
#include <libCCxll\CCxll.h>
#include "ARM_xl_wrapper_local.h"
#include "ARM_xl_trycatch_local.h"
#include "ARM_xl_gp_fctorhelper.h"
#include "ARM_local_interface.h"
#include <GP_Base\gpbase\vectormanip.h>
#include <GP_Base\gpbase\stringmanip.h>

#include "gpbase/ostringstream.h"

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_GenericAddin_Common(
	LPXLOPER XL_FunctionName,
	LPXLOPER XL_ParamNames,
	LPXLOPER XL_ParamValues,
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

		CCString C_FunctionNameXL;
		XL_readStrCell( XL_FunctionName,	C_FunctionNameXL, " ARM_ERR: Function Name string expected",C_result);
		string C_FunctionName = CCSTringToSTLString( C_FunctionNameXL );

		VECTOR<CCString> C_ParamNames;
		XL_readStrVector(XL_ParamNames,C_ParamNames," ARM_ERR: Param Names: array of string expected", XL_TYPE_STRING,C_result);

		long nbRows, nbCols;
		VECTOR<CCString> ParamValues;
		VECTOR<long> ParamValuesType;
		XL_readStrVectorSizeAndType(XL_ParamValues,nbRows,nbCols,ParamValues,ParamValuesType," ARM_ERR: ParamValues : array of string or double expected",DOUBLE_TYPE,C_result);
		if (nbCols != 1) 
			CC_NS(std,swap)(nbRows,nbCols);

		if (nbCols != 1)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ParamValues is a vector");

		if (C_ParamNames.size() != nbRows)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ParamNames and ParamsValues should have the same size.");
		

		ARM_GenericAddinDesc addinDesc = ARM_GenericAddinsDesc::TheAddinsTable()->GetAddin(C_FunctionName);

		ARM_GenericAddinFunctor* ourFunc = addinDesc.GetFunctor();

		ARM_GenericParams genParams(addinDesc.GetFunctionName());

		for (int i = 0; i < C_ParamNames.size(); ++i)
		{
			string nameSTL = CCSTringToSTLString(C_ParamNames[i]);

			ARM_GenericParamDesc paramDesc = addinDesc.GetParams().GetParam(nameSTL);

			/// input can be :
			/// - a number
			if (((ParamValuesType[i] == XL_TYPE_LONG)||(ParamValuesType[i] == XL_TYPE_DOUBLE)) && ((paramDesc.GetType() == GA_DOUBLE) || (paramDesc.GetType() == GA_DOUBLE_VECTOR)) )
			{
				genParams.SetParamValue(nameSTL,ARM_GenericParamValue(atof(ParamValues[i])));
			}
			/// - a vector of number
			else if ((ParamValuesType[i] == XL_TYPE_STRING) && (paramDesc.GetType() == GA_DOUBLE_VECTOR))
			{
				/// GA_DOUBLE_VECTOR can be a double or a vector of double
				long objectId = LocalGetNumObjectId(ParamValues[i]);
				if (objectId == ARM_KO)
					genParams.SetParamValue(nameSTL,ARM_GenericParamValue(CCSTringToSTLString(ParamValues[i])));
				else
					genParams.SetParamValue(nameSTL,ARM_GenericParamValue(objectId));
			}
			/// - a string
			else if ((ParamValuesType[i] == XL_TYPE_STRING) && (paramDesc.GetType() == GA_STRING))
			{
				genParams.SetParamValue(nameSTL,ARM_GenericParamValue(CCSTringToSTLString(ParamValues[i])));
			}
			/// - a vector of string
			else if ((ParamValuesType[i] == XL_TYPE_STRING) && (paramDesc.GetType() == GA_STRING_VECTOR))
			{
				/// GA_STRING_VECTOR can be a string or a vector of string
				long objectId = LocalGetNumObjectId(ParamValues[i]);
				if (objectId == ARM_KO)
					genParams.SetParamValue(nameSTL,ARM_GenericParamValue(CCSTringToSTLString(ParamValues[i])));
				else
					genParams.SetParamValue(nameSTL,ARM_GenericParamValue(objectId));
			}
			/// - an object
			else if ((ParamValuesType[i] == XL_TYPE_STRING) && (paramDesc.GetType() == GA_OBJECT))
			{
				long objectId = LocalGetNumObjectId(ParamValues[i]);

				genParams.SetParamValue(nameSTL,ARM_GenericParamValue(objectId));
			}
			else if(ParamValuesType[i] == xltypeNil)
			{
			}
			else
			{
				CC_Ostringstream os;

				os << nameSTL << " should be ";

				if (paramDesc.GetType() == GA_DOUBLE)
					os << "a double";
				else if (paramDesc.GetType() == GA_DOUBLE_VECTOR)
					os << "a vector object of double";
				else if (paramDesc.GetType() == GA_STRING)
					os << "a string";
				else if (paramDesc.GetType() == GA_STRING_VECTOR)
					os << "a vector object of string";
				else
					os << "an object";

				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str());;
			}
		}

		addinDesc.fillDefaultParams(&genParams);

		ourFunc->SetGenericParams(&genParams);

		if (addinDesc.IsRetObj())
		{
			if (ourFunc->ClassName())
			{
				/// call the general function
				fillXL_Result( CCString( ourFunc->ClassName()), *ourFunc, C_result, XL_result, PersistentInXL );
			}
			else
			{
				/// call the general function
				fillXL_Result_withName( *ourFunc, C_result, XL_result, PersistentInXL );
			}
		}
		else
		{
			double objId=-1;
			long retCode = (*ourFunc)(C_result,objId);

			if (retCode == ARM_OK)
			{
				int nbRows = ourFunc->GetNbRows();
				int nbCols = ourFunc->GetNbCols();

				const vector <ARM_GenericAddinFunctor::RetStruct>& retValues = ourFunc->GetRetValues();

				VECTOR<CCString> C_OutResult(nbRows*nbCols, "");
				VECTOR<long> types(nbRows*nbCols, xltypeStr);
				
				for (int i = 0; i < nbRows*nbCols; ++i)
				{
					if (retValues[i].type == ARM_GenericAddinFunctor::DOUBLE)
					{
						char outstr[30];
						sprintf(outstr,"%lf",retValues[i].dblVal);
						C_OutResult[i] = outstr;
						types[i] = xltypeNum;
					}
					else
					{
						C_OutResult[i] = retValues[i].strVal.c_str();
					}
				}

				XL_writeStrMatrixSizeAndType(XL_result, C_OutResult, types, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result);
			}
			else
			{
				ARM_ERR();
			}
		}
		
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GenericAddin" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

LPXLOPER Local_GenericAddin(
	LPXLOPER XL_FunctionName,
	LPXLOPER XL_ParamNames,
	LPXLOPER XL_ParamValues)
{
	bool PersistentInXL = true;

	return Local_GenericAddin_Common(
	XL_FunctionName,
	XL_ParamNames,
	XL_ParamValues,
	PersistentInXL);
}

LPXLOPER PXL_Local_GenericAddin(
	LPXLOPER XL_FunctionName,
	LPXLOPER XL_ParamNames,
	LPXLOPER XL_ParamValues)
{
	bool PersistentInXL = false;

	return Local_GenericAddin_Common(
	XL_FunctionName,
	XL_ParamNames,
	XL_ParamValues,
	PersistentInXL);
}

LPXLOPER Local_GenericAddin_HelperCommon(
	LPXLOPER XL_FunctionName,
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

		CCString C_FunctionNameXL;
		CCString DefFunctionName = "";
		XL_readStrCellWD( XL_FunctionName,	C_FunctionNameXL, DefFunctionName, " ARM_ERR: Function Name string expected",C_result);
		string C_FunctionName = CCSTringToSTLString( C_FunctionNameXL );
		
		exportFunc1Arg< string> ourFunc (C_FunctionName, ARMLOCAL_GenericAddin_Helper);
		

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GenericAddin_Helper" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

LPXLOPER Local_GenericAddin_Helper(
	LPXLOPER XL_FunctionName)
{
	bool PersistentInXL = true;

	return Local_GenericAddin_HelperCommon(
	XL_FunctionName,
	PersistentInXL);
}

LPXLOPER PXL_Local_GenericAddin_Helper(
	LPXLOPER XL_FunctionName)
{
	bool PersistentInXL = false;

	return Local_GenericAddin_HelperCommon(
	XL_FunctionName,
	PersistentInXL);
}

LPXLOPER Local_GenericAddin_ParamNames(
	LPXLOPER XL_FunctionName,
	LPXLOPER XL_WithDefault,
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

		CCString C_FunctionNameXL;
		XL_readStrCell( XL_FunctionName,	C_FunctionNameXL, " ARM_ERR: Function Name string expected",C_result);
		string C_FunctionName = CCSTringToSTLString( C_FunctionNameXL );

		CCString C_WithDefault;
		CCString defWithDefault = "N";
		XL_readStrCellWD( XL_WithDefault,	C_WithDefault, defWithDefault,  " ARM_ERR: With Default string expected",C_result);
		int withDefaultInt;
		if((withDefaultInt = ARM_ConvYesOrNo (C_WithDefault, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		bool withDefault = (withDefaultInt?1:0);
		
		vector<string> paramNames;
		long retCode = ARMLOCAL_GenericAddin_ParamNames(C_FunctionName, withDefault, paramNames, C_result);

		VECTOR<CCString> results;

		int i;
		for (i = 0; i < paramNames.size(); ++i)
			results.push_back(paramNames[i].c_str());

		/// feed the LPXLOPER object result
		if (retCode == ARM_OK)
		{
			XL_writeStrVector( XL_result, results, " ARM_ERR: Could not get result data", C_result);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GenericAddin_ParamNames" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}