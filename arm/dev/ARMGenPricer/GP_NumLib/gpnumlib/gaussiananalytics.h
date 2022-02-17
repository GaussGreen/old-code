/*!
 *
 * Copyright (c) CDC IXIS CM November 2004 Paris
 *
 *	\file gaussiananalytics.h
 *
 *  \brief General file for the composite generator
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */

#ifndef _INGPNUMLIB_GAUSSIANANALYTICS_H
#define _INGPNUMLIB_GAUSSIANANALYTICS_H

#include "gpbase/port.h"
#include "gpbase/env.h"

CC_BEGIN_NAMESPACE( ARM )


class ARM_GaussianAnalytics
{
public:
	static double  dNormal(double d);
	static double  cdfNormal(double d);
	static double cdfNormal2(double d);
	/// inverse Functions
	static double  cdfNormal_Inv(double d);
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
