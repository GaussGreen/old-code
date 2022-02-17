/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: rateconversion.h,v $
 * Revision 1.1  2003/10/08 16:45:06  ebenhamou
 * Initial revision
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file rateconversion.h
 *
 *  \brief files to convert a rate
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPBASE_RATECONVERSION_H
#define _INGPBASE_RATECONVERSION_H

/// use our macro for namespace
#include "port.h"

CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////////////////////////////////////
/// function to convert from a decompounding method to another decompounding method
///////////////////////////////////////////////////////////////////////
double ConvertRateToRate(double rate, double term, int fromCompMeth, int toCompMeth);

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
