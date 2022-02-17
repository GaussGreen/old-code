/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_wrapper_local.h,v $
 * Revision 1.1  2004/01/13 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_xl_wrapper_local.h
 *
 *  \brief file various function to factorise the interface
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date September 2003
 */


#ifndef ARM_XL_WRAPPER_LOCAL_H
#define ARM_XL_WRAPPER_LOCAL_H

#include <ARM\libarm\ARM_result.h>
#include "ARM\libarm_local\ARM_local_gp_genericaddin.h"
#include <libCCxll\CCxll.h>

/*
 * In this header
 * various small utility to faciliate 
 * the code of excel interface
 */

/// forward declaration
#include <functional>
using std::binary_function;

/*!
 * General function for all excel interface 
 * intended to avoid the monstruous copy and
 * paste in excel interface
 *
 * the class T has to be a functor that 
 * supports long operator()( ARM_result, long )
 * derived from binary_function to make it easily usable
 * for STL
 *
 */
/*class ARMResultLong2LongFunc : public binary_function< ARM_result, long, long >
{
public:
	/// pure virtual function to
	/// force redefinition
	virtual	long operator()( ARM_result& result, long objId ) = 0;
};*/


/*!
 * General function to fill result
 * handle persistence
 * the functor is previously constructed with a context
 * see for instance ARM_xl_infcurv_local
 */
void fillXL_Result( const CCString& curClass,
				   ARMResultLong2LongFunc& CreateObjFctor, 
				   ARM_result& C_result,
				   XLOPER& XL_result,
				   bool PersistentInXL,
				   bool GetNameFromResult=false);

/*!
 * General macro when there is an error
 */
void fillXL_Result_withName(ARMResultLong2LongFunc& CreateObjFctor, 
							ARM_result& C_result,
							XLOPER& XL_result,
							bool PersistentInXL );
/*!
 * General macro when there is an error
 */
#define ARM_ARG_ERR_AND_EXIT() \
	{\
		ARM_ARG_ERR();\
		return (LPXLOPER)&XL_result;\
	}


/*!
 * various methods for converting string to long type
 * because of the define interface for excel
 * we are forced to use define
 * not very robust code but required in order to use 
 * function XL_readStrCellWD wtih parameters:
 *		-reason, error, XL_Result
 *
 * the type of the arguments are the following
 *		- LPXLOPER XL_Arg,
 *		- ARM_result& C_result,
 *		- const char* typeOfArg
 *		- long (*GetResultFunc) ( const CCString&, ARM_result&),
 *		- const CCString& defaultValue
 *		- returnValue
 */

/*!
 * By convention all these macros are in capital letter
 * they aim at providing the same interface as 
 * the other excel macros
 *
 * we enclose these macros in bracket to avoid splitting the instructions
 */


#define XL_GETMETHODWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result, function )\
	{\
		CCString C_ArgStr_##XL_Arg;\
		XL_readStrCellWD( XL_Arg, C_ArgStr_##XL_Arg, defaultValue, errorMsg, C_result);\
		returnValue = function( C_ArgStr_##XL_Arg, C_result );\
		if( returnValue == ARM_DEFAULT_ERR )\
			ARM_ARG_ERR_AND_EXIT();\
	}
#define XL_GETMETHODWD_NOCRESULT( XL_Arg, returnValue, defaultValue, errorMsg, function )\
	{\
		CCString C_ArgStr_##XL_Arg;\
		XL_readStrCellWD( XL_Arg, C_ArgStr_##XL_Arg, defaultValue, errorMsg, C_result);\
		returnValue = function( C_ArgStr_##XL_Arg);\
		if( returnValue == ARM_DEFAULT_ERR )\
			ARM_ARG_ERR_AND_EXIT();\
	}

/*!
 * XL_Arg is the LPXLOPER argument that contains the input
 * return value is the variable that is used to get the value
 * defaultValue is a proper default in case of missing data
 * errorMsg is a char* of the type " ARM_ERR: blah blah: expected"
 * C_result is an ARM_result variable used to store errors
 */
#define XL_GETFREQUENCYWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result, ARM_ConvFrequency );

#define XL_GETDAYCOUNTWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD_NOCRESULT( XL_Arg, returnValue, defaultValue, errorMsg, ARM_ConvDayCount );

#define XL_GETFWDRULEWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result, ARM_ConvFwdRule );

#define XL_GETINTRULEWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD_NOCRESULT( XL_Arg, returnValue, defaultValue, errorMsg, ARM_ConvIntRule );

#define XL_GETPAYRESETTIMINGWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD_NOCRESULT( XL_Arg, returnValue, defaultValue, errorMsg, ARM_ConvPayResetRule );

#define XL_GETCONVRULEWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD_NOCRESULT( XL_Arg, returnValue, defaultValue, errorMsg, ARM_ConvStubRule );

#define XL_GETNXCHANGEWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result, ARM_NotionalExchange );

#define XL_GETNXTYPEWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result, ARM_NotionalType );

#define XL_GETRCVORPAYWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result, ARM_ConvRecOrPay );

#define XL_GETCPIDLYINTERPWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result, ARM_ConvCPIDailyInterpMethod );

#define XL_GETVOLTYPEWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result, ARM_ConvVolType );

#define XL_GETSTRIKETYPEWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result, ARM_StrikeCode );

#define XL_GETDIMTYPEWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result, ARM_GetDimType );

#define XL_GETCAPORFLOORWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result, ARM_ConvCapOrFloor );

#define XL_GETCONVSHAPETYPEWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result, ARM_ConvShapeType );

#define XL_GETCONVCALLORPUTWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result, ARM_ConvCallOrPut );

#define XL_GETCONVSABRFLAGWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result ) \
	XL_GETMETHODWD( XL_Arg, returnValue, defaultValue, errorMsg, C_result, ARM_ConvGP_CFSABR_ImplicitVol_Formula_Extended_Flag );

 
 
 /*!
 * Version with no default
 */
#define XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, function )\
	{\
		CCString C_ArgStr_##XL_Arg;\
		XL_readStrCell( XL_Arg, C_ArgStr_##XL_Arg, errorMsg, C_result);\
		returnValue = function( C_ArgStr_##XL_Arg, C_result );\
		if( returnValue == ARM_DEFAULT_ERR )\
			ARM_ARG_ERR_AND_EXIT();\
	}

#define XL_GETFREQUENCY( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvFrequency );

#define XL_GETDAYCOUNT( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvDayCount );

#define XL_GETFWDRULE( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvFwdRule );

#define XL_GETINTRULE( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvIntRule );

#define XL_GETPAYRESETTIMING( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvPayResetRule );

#define XL_GETCONVRULE( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvStubRule );

#define XL_GETNXCHANGE( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_NotionalExchange );

#define XL_GETNXTYPE( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_NotionalType );

#define XL_GETRCVORPAY( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvRecOrPay );

#define XL_GETCPIDLYINTERP( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvCPIDailyInterpMethod );

#define XL_GETVOLTYPE( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvVolType );

#define XL_GETSTRIKETYPE( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_StrikeCode );

#define XL_GETDIMTYPE( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_GetDimType );

#define XL_GETCAPORFLOOR( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvCapOrFloor );
 
#define XL_GETINFSWAPTYPE( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvInfSwapType );

#define XL_GETCONVSHAPETYPE( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvShapeType );
 
#define XL_GETCONVCALLORPUT( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvCallOrPut );

#define XL_GETCONVSPREADDERIVEDSTRUCT( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvGP_CFSpreadDigitalOption_Formula_DerivedStruct );

#define XL_GETCONVSABRFLAG( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvGP_CFSABR_ImplicitVol_Formula_Extended_Flag );

#define XL_GETCONVOPTIMIZATIONFLAG( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvGP_CFOptimization_ObjectiveFuntion_Algorithm_Flag );

#define XL_GETCONVHESTONINTERPOMETH( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ArgConvGP_CFHeston_Vector_InterpolationMethod_Flag );

#define XL_GETCONVBARRIEREINOUT( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvGP_CFBS_EuropeanBarriere_Formula_InOut_Flag );

#define XL_GETCONVBARRIEREUPDOWN( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ConvGP_CFBS_EuropeanBarriere_Formula_UpDown_Flag );

#define XL_GETCONVBARRIERETYPE( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ArgConvGP_CFBS_EuropeanBarriere_Formula_OptionType_Flag );

#define XL_GETCONVPARTIALBARRIERETYPE( XL_Arg, returnValue, errorMsg, C_result ) \
	XL_GETMETHOD( XL_Arg, returnValue, errorMsg, C_result, ARM_ArgConvGP_CFBS_PartialTime_Barriere_End_Formula_OptionType_Flag );
 
/*!
 * a little bit specific hence the dedicated function
 * ARM_NULL_OBJECT is defined by #define ARM_NULL_OBJECT -11111
 */
#define XL_GETOBJIDWD(  XL_Arg, returnValue, defaultValue, errorMsg, C_result )\
	{\
		CCString C_ArgStr_##XL_Arg;\
		XL_readStrCellWD(XL_Arg, C_ArgStr_##XL_Arg,	defaultValue, errorMsg, C_result);\
		if (  C_ArgStr_##XL_Arg == defaultValue )\
			returnValue = ARM_NULL_OBJECT;\
		else\
			returnValue = LocalGetNumObjectId( C_ArgStr_##XL_Arg );\
	}

/*
 * Macro to get objId directly with creation of the CCString
 */
#define XL_GETOBJID( XL_Arg, returnValue, errorMsg, C_result ) \
	{\
		CCString ObjId_##Xl_arg;\
		XL_readStrCell( XL_Arg, ObjId_##Xl_arg, errorMsg,C_result);\
		returnValue = LocalGetNumObjectId( ObjId_##Xl_arg );\
	}

/*
 * Macro for general failure message
 */
#define ARM_GENERAL_ERROR_MSSG( mssg ) \
	{ \
		Exception x(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, mssg ); \
		x.DebugPrint(); \
	}	


#endif /* ARM_XL_WRAPPER_LOCAL_H */