/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file asian_lognormal_interface.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
///////////////////////////////////////////////////////////////////////////////
///
///		Asian Option in a Lognormal world, Geman Yor Formula		
///
/////////////////////////////////////////////////////////////////////////////:

 
#ifndef _GP_CF_ASIAN_LOGNORMAL_INTERFACE_H
#define _GP_CF_ASIAN_LOGNORMAL_INTERFACE_H


#include "firsttoinc.h"
#include "gpbase/port.h"



CC_BEGIN_NAMESPACE(ARM)

/////////////////////////////////////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
/////////////////////////////////////////////////////////////////////////////////////////////////////
double Export_Lognormal_Asian_VanillaOption(
											double S,
											double k,
											double T,
											double r,
											double v,
											double alpha,
											int callput,
											int n
											);

///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Lognormal_Asian_VanillaOption(int i,
											double S,
											double k,
											double T,
											double r,
											double v,
											double alpha,
											int callput,
											int n
											);

///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Lognormal_Asian_VanillaOption(int i,int j,
											double S,
											double k,
											double T,
											double r,
											double v,
											double alpha,
											int callput,
											int n
											);





CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
