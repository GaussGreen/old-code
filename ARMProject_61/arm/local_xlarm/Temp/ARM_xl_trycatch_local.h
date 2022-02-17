/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_trycatch_local.h,v $
 * Revision 1.1  2004/01/13 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_xl_trycatch_local.h
 *
 *  \brief file various macros to handle exception in the interface
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date September 2003
 */


#ifndef ARM_XL_TRYCATCH_LOCAL_H
#define ARM_XL_TRYCATCH_LOCAL_H

#include "expt.h"
#include <CCString.h>

/// macro equivalent to ARM_Result: the difference being that in the local_xlarm
///     the ARM_result variable is called C_resutl and the CCString msg is not defined!

// JLA: don't use limited buffer that will crash. 
//
// #define XL_SET_MSSG_TO_RESULT()			C_result.setRetCode(ARM_KO);    \
// 										char* var = new char[300];      \
// 										x.GetErrorMessage(var);         \
// 										CCString msg( var );            \
// 										C_result.setMsg(msg);           \
//                                      delete var;
// 

#define XL_SET_MSSG_TO_RESULT()			{ C_result.setRetCode(ARM_KO);    \
										C_result.setMsgString(x.GetErrorString()) ; }


#define ARM_ARG_ERR_WITH_MSG( msg )									\
	Exception x(__LINE__, __FILE__,	ERR_INVALID_ARGUMENT, msg );	\
	x.DebugPrint();													\
	XL_SET_MSSG_TO_RESULT()	;										\
	ARM_ARG_ERR();										


/// debug version!
#if defined(_DEBUG)
	#define ARM_XL_TRY_IS_ACTIVE

	#define ARM_XL_TRY_BLOCK_BEGIN	try

	#define ARM_XL_TRY_BLOCK_END         

	#define ARM_XL_CATCH_ARM_EXPT								\
		catch(Exception& x)										\
		{														\
			x.DebugPrint();										\
			C_result.setRetCode(ARM_KO);						\
			C_result.setMsgString(x.GetErrorString()) ;			\
			ARM_ARG_ERR();										\
		}														\
		catch(std::exception& x)								\
		{														\
			C_result.setRetCode(ARM_KO);						\
			C_result.setMsgString(x.what())	;					\
			ARM_ARG_ERR();										\
		}														

	#define ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( msg )				\
		catch( ... )											\
		{														\
			Exception x(__LINE__, __FILE__,	ERR_INVALID_ARGUMENT, msg ); \
			x.DebugPrint();										\
			C_result.setRetCode(ARM_KO);						\
			C_result.setMsgString(x.GetErrorString()) ;			\
			ARM_ARG_ERR();										\
		}

#endif

/// release version
#if !defined(_DEBUG)
	#define ARM_XL_TRY_IS_ACTIVE

	#define ARM_XL_TRY_BLOCK_BEGIN	try

	#define ARM_XL_TRY_BLOCK_END         

	#define ARM_XL_CATCH_ARM_EXPT								\
		catch(Exception& x)										\
		{														\
			x.DebugPrint();										\
			C_result.setRetCode(ARM_KO);						\
			C_result.setMsgString(x.GetErrorString());			\
			ARM_ARG_ERR();										\
		}														\
		catch(std::exception& x)								\
		{														\
			C_result.setRetCode(ARM_KO);						\
			C_result.setMsgString(x.what())	;					\
			ARM_ARG_ERR();										\
		}														

	#define ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( msg )				\
		catch( ... )											\
		{														\
			Exception x(__LINE__, __FILE__,	ERR_INVALID_ARGUMENT, msg ); \
			x.DebugPrint();										\
			C_result.setRetCode(ARM_KO);						\
			C_result.setMsgString(x.GetErrorString());			\
			ARM_ARG_ERR();										\
		}

#endif		





#endif /* ARM_XL_TRYCATCH_LOCAL_H */
