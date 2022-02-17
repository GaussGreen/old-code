/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file mepi_interface.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2006
 */
///////////////////////////////////////////////////////////////////////////////
///
///			
///			
///			
///			
///
///			
///
/////////////////////////////////////////////////////////////////////////////:

 
#ifndef _GP_CF_MEPI_INTERFACE_H
#define _GP_CF_MEPI_INTERFACE_H


#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/typedef.h"
#include "gpbase/gplinalgconvert.h"
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE(ARM)

/////////////////////////////////////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

vector<double>* Export_Mepi_EDP_VanillaOption(
							  double P0,
							  double f0,
							  double T,
							  double K,
							  double R,
							  double Emin,
							  double Lmax,
							  double gamma0,
							  double gamma1,  
							  double sig,
							  double lambda,
							  double sigJ,
							  double r,
							  double s,
							  double mu,
							  double fees,
							  double voldrift,
							  double volvol,
							  int CallPut,
							  ARM_GP_Vector* A_params);

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
