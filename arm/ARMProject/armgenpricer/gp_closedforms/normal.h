/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file normal.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_NORMAL_H
#define _GP_CF_NORMAL_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"


CC_BEGIN_NAMESPACE(ARM)

double NormalCDF(double x);

double NormalPDF(double x);

double NormalCDFInverse(double x);

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/