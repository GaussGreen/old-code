/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file primenb.h
 *  \brief General file about prime nb
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPNUMLIB_PRIMENB_H
#define _INGPNUMLIB_PRIMENB_H

#include "gpbase/port.h"
#include "gpbase/env.h"

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////
/// \class ARM_PrimeNb
//////////////////////////////
struct ARM_PrimeNb
{
	static int Prime ( int n );
	static int Prime_ge ( int n );
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
