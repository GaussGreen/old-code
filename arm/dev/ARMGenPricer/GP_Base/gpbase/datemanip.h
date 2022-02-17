/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file datemanip.h
 *	\brief manipulation on dates
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPBASE_DATEMANIP_H
#define _INGPBASE_DATEMANIP_H

#include "port.h"
#include "env.h"

CC_BEGIN_NAMESPACE( ARM )

const double XLOriginDateInJulianDate = 2415019.0;
extern double ConvertXLDateToJulian( double d );
extern double ConvertJulianToXLDate( double d );

CC_END_NAMESPACE()

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
