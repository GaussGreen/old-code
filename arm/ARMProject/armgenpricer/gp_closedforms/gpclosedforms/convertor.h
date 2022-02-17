/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions to convert
 *
 *	\file convertor.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date November 2004
 */
 
#ifndef _GP_CF_CONVERTOR_H
#define _GP_CF_CONVERTOR_H

#include "gpbase/port.h"

CC_BEGIN_NAMESPACE(ARM)

double ConvertFXOptionPriceToStrike( double targetPrice, double fwd, double totalVol, int callput, double initValue );

CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
