/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file merton.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_MERTON_INTERFACE_H
#define _GP_CF_MERTON_INTERFACE_H


#include <glob/firsttoinc.h>
#include "gpbase/port.h"



CC_BEGIN_NAMESPACE(ARM)

/////////////////////////////////////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
/////////////////////////////////////////////////////////////////////////////////////////////////////
double Export_Merton_JumpDiffusion(double F,
						double K,
						double t,
						double sigma,
						double lambda,
						double muJ, 
						double sigmaJ,
						int callorput,int nb);

///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Merton_JumpDiffusion(int i,
						double F,
						double K,
						double t,
						double sigma,
						double lambda,
						double muJ, 
						double sigmaJ,
						int callorput,int nb);

///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Merton_JumpDiffusion(int i,int j,
							double F,
							double K,
							double t,
							double sigma,
							double lambda,
							double muJ, 
							double sigmaJ,
							int callorput,int nb);



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
