/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file asian_lognormal.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_ASIAN_LOGNORMAL_H
#define _GP_CF_ASIAN_LOGNORMAL_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"


CC_BEGIN_NAMESPACE(ARM)



///////////////////////////////////////////////////////////////////////
///  
///			   Pricing Fonctions
///
///////////////////////////////////////////////////////////////////////

double GemanYorAsianVanillaCall(double S,double k,double T,double r,double v,double alpha,int n);


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

