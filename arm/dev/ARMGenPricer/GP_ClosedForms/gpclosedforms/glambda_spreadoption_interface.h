/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file glambda_spreadoption_interface.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date sep 2005
 */
 
#ifndef _GP_CF_GLAMBDA_SPREADOPTION_INTERFACE_H
#define _GP_CF_GLAMBDA_SPREADOPTION_INTERFACE_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

CC_BEGIN_NAMESPACE(ARM)


Export_Student_GLambda_Power_Digital_SpreadOption(double l1a,double l2a,double l3a,double l4a,double l5a,double l6a,
												double l1b,double l2b,double l3b,double l4b,double l5b,double l6b,
												double copula_corr,double copula_degre,
												double a20,double b20,double k20,int n);



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


