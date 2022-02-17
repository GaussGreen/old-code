/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file eigenvalues.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date March 2006
 */
 
#ifndef _GP_CF_EIGENVALUES_H
#define _GP_CF_EIGENVALUES_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"


CC_BEGIN_NAMESPACE(ARM)


double eigenvalues(double rho12,double rho13, double rho14,double rho23,double rho24,double rho34);


CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/