/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file datemanip.cpp
 *	\brief manipulation on dates
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date November 2004
 */

#include "gpbase/datemanip.h"

CC_BEGIN_NAMESPACE( ARM )

double ConvertXLDateToJulian( double d )
{
	return d+XLOriginDateInJulianDate;
}

extern double ConvertJulianToXLDate( double d )
{
	return d-XLOriginDateInJulianDate;
}

CC_END_NAMESPACE()

/*----------------------------------------------------------------------*/
/*---- End of file ----*/
