/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: gramargs.h,v $
 * Revision 1.2  2003/13/08 19:46:06  ebenhamou
 * Remove unecessary type
 *
 * Revision 1.1  2003/10/08 16:46:06  ebenhamou
 * Initial revision
 *
 *
 */



/*! \file gramargs.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_GRAMARGS_H
#define _INGPINFRA_GRAMARGS_H

#include "gpbase/port.h"
#include <string>

CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////////////////////////////////////////
/// grammar function argument types (GFAT).  These are used both 
/// in declaring grammar functions and
/// in doing type inference and conversions 
/// during cash flow table parsing and translation.
///////////////////////////////////////////////////////////////////////////

enum GramFuncArgType 
{
	/// null terminator in arg lists.
	GFAT_TERMINATOR = 0,

	/// A marker for tableau entries to say the end!
	GFAT_UNKNOWN_TYPE,

	/// Types returned from parsing
	GFAT_DOUBLE_TYPE,			    /// double!
	GFAT_VECTOR_TYPE,			    /// Vector if MC, collection of states if BI. can be almost everything!
	GFAT_MATRIX_TYPE,			    /// Matrix for intermediate payoffs
    GFAT_CURVE_TYPE ,               /// Curve fore a collection of strikes, notionels,...!
	GFAT_STRING_TYPE,			    /// string
	GFAT_STRINGVEC_TYPE,			/// string vector
	GFAT_STRINGVECTRANS_TYPE,		/// string vector transposed

	GFAT_DATE_TYPE,				    /// date either as double or date expr
	GFAT_MATURITY_TYPE,			    /// maturity		
	GFAT_DATEORMATU_TYPE,		    /// date or maturity
	GFAT_DATESTRIP_TYPE,		    /// dates strip


	GFAT_MULTITOKENSTRING_TYPE,	    /// multi tokens string (stopped either by comma or by eot (end of text!)
	GFAT_FUTUREREF_TYPE,		    /// exclusively for PV!
	GFAT_VECTOR_OR_FUTUREREF_TYPE,  /// PV or vector of states
	GFAT_MODEL_TYPE,			    /// model

	/// combined type for dates
	GFAT_DATE_OR_DOUBLE_TYPE,			/// can be a double or a date!
	GFAT_DATE_OR_VECTOR_TYPE,			/// can be a double or a date!
	GFAT_DATE_VECTOR_OR_DOUBLE_TYPE,	/// can be a date, a vector or a double!
	GFAT_VECTOR_OR_CURVE_TYPE,			/// can be a vector or a curve!
	GFAT_STRING_OR_CURVE_TYPE,			/// can be a string or a curve!
};


/// conversion to string
/// used to print an argument
CC_NS( std, string ) GramFuncArgTypeToString( GramFuncArgType t );

CC_END_NAMESPACE()


#endif
/*-----------------------------------------------------------------------------------------*/
/*---- End Of File ----*/
