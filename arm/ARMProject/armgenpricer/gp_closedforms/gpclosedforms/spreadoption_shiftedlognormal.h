/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_shiftedlognormal.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_SPREADOPTION_SHIFTEDLOGNORMAL_H
#define _GP_CF_SPREADOPTION_SHIFTEDLOGNORMAL_H

#include "firsttoinc.h"
#include "gpbase/port.h"


#include "basic_distributions.h" /// for ArgumentList_Checking_Result


CC_BEGIN_NAMESPACE(ARM)


double ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula_value2(const ArgumentList& a);


/// all the functionalities are implemented through templated pricing  and the templated function PowerSpreadOption_pricing



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

