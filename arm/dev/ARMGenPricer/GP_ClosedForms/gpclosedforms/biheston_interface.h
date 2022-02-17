/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file biheston.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date Avril 2007
 */
///////////////////////////////////////////////////////////////////////////////
///
///					Process :
///			dS= S dW1
///			dV=lambda*(Vinf-V)dt +nu*V1^(1/2) dW2 
///			dW1.dW2=rho1*dt
///
///			dS2= S2 dW3
///			dV2=lambda2*(Vinf2-V2)dt +nu2*V2^(1/2) dW4 
///			dW3.dW4=rho2*dt
///
///			dW1.dW3=rhos*dt
///			dW1.dW4=rhoc12*dt
///			dW2.dW3=rhoc21*dt
///			dW2.dW4=rho2v*dt

/////////////////////////////////////////////////////////////////////////////:

 
#ifndef _GP_CF_BIHESTON_INTERFACE_H
#define _GP_CF_BIHESTON_INTERFACE_H


#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/typedef.h"

CC_BEGIN_NAMESPACE(ARM)


double Export_BiShiftedHeston_VanillaOption(
											double		C_F1,
											double		C_V1,
											double		C_Vinfini1,
											double		C_lambda1,
											double		C_nu1,
											double		C_rho1,
											double		C_gamma1,
											double		C_F2,
											double		C_V2,
											double		C_Vinfini2,
											double		C_lambda2,
											double		C_nu2,
											double		C_rho2,
											double		C_gamma2,
											double		rhos,
											double		rhov,
											double		rhoc12,
											double		rhoc21,
											double		C_k,
											double		C_T,
											double		C_CallPut,
											double		C_LambdaB,
											double		C_Flag,
											double		C_nbfirst,
											double		C_nb,
											double		C_NbStage,
											double		C_NbOscill,
											double		C_prec
											);





CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
