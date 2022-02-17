/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file gaussian_integrals.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_CEV_INTERFACE_H
#define _GP_CF_CEV_INTERFACE_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

CC_BEGIN_NAMESPACE(ARM)

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   Forward value  CEV VanillaOption
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


double Export_CEV_VanillaOption(double f, double K, double T,double drift, double sig, double beta,int callput, int nbsteps);

///  
///			1st Derivatives  CEV VanillaOption
///

double Export_CEV_VanillaOption(int i,double f, double K, double T,double drift, double sig, double beta,int callput, int nbsteps);	



///  
///			2nd Derivatives  CEV VanillaOption
///

double Export_CEV_VanillaOption(int i,int j,double f, double K, double T,double drift, double sig, double beta,int callput, int nbsteps);




/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   Forward value  CEV DoubleBarrierOption
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


double Export_CEV_DoubleBarrierOption(double f, double K, double T,double barrdown,double barrup, double drift, double sig, double beta,int callput, int nbsteps);


///  
///			1st Derivatives  CEV DoubleBarrierOption
///

double Export_CEV_DoubleBarrierOption(int i,double f, double K, double T,double barrdown,double barrup,double drift, double sig, double beta,int callput, int nbsteps);


///  
///			2nd Derivatives  CEV DoubleBarrierOption
///

double Export_CEV_DoubleBarrierOption(int i,int j,double f, double K, double T,double barrdown,double barrup,double drift, double sig, double beta,int callput, int nbsteps);


/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   Forward value  CEV BarrierOption
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


double Export_CEV_BarrierOption(double f, double K, double T,double barrier, double drift, double sig, double beta, int optype,int callput, int nbsteps);


///  
///			1st Derivatives  CEV DoubleBarrierOption
///

double Export_CEV_BarrierOption(int i,double f, double K, double T,double barrier, double drift, double sig, double beta, int optype,int callput, int nbsteps);


///  
///			2nd Derivatives  CEV DoubleBarrierOption
///

double Export_CEV_BarrierOption(int i,int j,double f, double K, double T,double barrier, double drift, double sig, double beta, int optype,int callput, int nbsteps);






CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


