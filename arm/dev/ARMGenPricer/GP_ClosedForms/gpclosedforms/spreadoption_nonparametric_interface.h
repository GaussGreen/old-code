/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 03/15/2007
 *
 *	\file spreadoption_nonparametric_interface.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date March 2007
 */
 
#ifndef _GP_CF_SPREADOPTION_NONPARAMETRIC_INTERFACE_H
#define _GP_CF_SPREADOPTION_NONPARAMETRIC_INTERFACE_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "gpclosedforms/basic_distributions.h" /// for ArgumentList_Checking_Result
#include "gpclosedforms/powerspreadoption.h"		/// for instanciation of the templates
#include "gpbase/gpvector.h"
#include "gpclosedforms/smile_nonparametric.h"


using namespace std;

CC_BEGIN_NAMESPACE(ARM)


double Export_Nonparametric_CompleteSpreadoption(
		ARM_GP_Vector* strike_Vec1,
		ARM_GP_Vector* vol_Vec1,
		ARM_GP_Vector* strike_Vec2,
		ARM_GP_Vector* vol_Vec2,double S1,double S2,
		double index_begin1,double index_end1,int flag_begin1,int flag_end1,
		double index_begin2,double index_end2,int flag_begin2,int flag_end2,
		double correlation,double maturity,double a1,double b1,double k1,double a2,double b2,double k2,
		int nbsteps,int algorithm,int smiletype
	 );


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

