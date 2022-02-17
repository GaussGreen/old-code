/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_lognormal.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_SPREADOPTION_NORMAL_H
#define _GP_CF_SPREADOPTION_NORMAL_H


#include <glob/firsttoinc.h>
#include "gpbase/port.h"


CC_BEGIN_NAMESPACE(ARM)

///////////////////////////////////////////////////////////////////////
///  
///			 Fonctions de base 
///
///////////////////////////////////////////////////////////////////////


/// callput =  K_CALL  for call
/// callput =  K_PUT  for put



double SpreadDigitalOption_N(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int otype);

double Vega1SpreadDigitalOption_N(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int otype);

double Vega2SpreadDigitalOption_N(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int otype);





CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
