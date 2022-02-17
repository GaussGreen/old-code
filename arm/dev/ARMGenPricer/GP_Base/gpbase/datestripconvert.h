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

/*! \file datestripconvert.h
 *
 *  \brief files to convert a swapleg schudele to datestrip
 *
 *	\author  E.Ezzine
 *	\version 1.0
 *	\date June 2006
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPBASE_DATESTRIPCONVERT_H
#define _INGPBASE_DATESTRIPCONVERT_H

/// use our macro for namespace
#include "port.h"
#include "swapleg.h"

CC_BEGIN_NAMESPACE( ARM )
class ARM_DateStrip;

// Function to convert a SwapLeg Schedule into a DateStrip
ARM_DateStrip* SwapLegToDateStrip( const ARM_SwapLeg& swapleg);

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
