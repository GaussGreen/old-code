/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: stringconvert.h,v $
 * Revision 1.1  2003/10/08 16:45:06  ebenhamou
 * Initial revision
 *
 *
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file stringconvert.h
 *
 *  \brief files to convert a string to various type and vice versa
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPBASE_STRINGCONVERT_H
#define _INGPBASE_STRINGCONVERT_H

/// use our macro for namespace
#include "port.h"
#include <string>
CC_USING_NS(std,string)


CC_BEGIN_NAMESPACE( ARM )

double StringMaturityToYearTerm( const string& maturity );
string YearTermToStringMaturity( double yearTerm );
int FromIndexAndTermToIndexType( const string& term, const string& index );
string FromIndexTypeToTermAndType( int indexType, int& );

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
