/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file shifted_lognormal_interface.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
///////////////////////////////////////////////////////////////////////////////
///
///					Process :
///			dS= (S+beta)dt+ sigma*dW
///			
///
/////////////////////////////////////////////////////////////////////////////:

 
#ifndef _GP_CF_SHIFTED_LOGNORMAL_INTERFACE_H
#define _GP_CF_SHIFTED_LOGNORMAL_INTERFACE_H


#include <glob/firsttoinc.h>
#include "gpbase/port.h"

CC_BEGIN_NAMESPACE(ARM)

/////////////////////////////////////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
/////////////////////////////////////////////////////////////////////////////////////////////////////
double Export_shifted_lognormal_VanillaOption(
									double F,
									double K,
									double V0,
									double t,
									double alpha,
									int callorput
									);

///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_shifted_lognormal_VanillaOption(int i,
									double F,
									double K,
									double V0,
									double t,
									double alpha,
									int callorput
									);

///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_shifted_lognormal_VanillaOption(int i,int j,
									double F,
									double K,
									double V0,
									double t,
									double alpha,
									int callorput
									);





CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
