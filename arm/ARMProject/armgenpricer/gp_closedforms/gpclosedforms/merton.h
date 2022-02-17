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
 
#ifndef _GP_CF_MERTON_H
#define _GP_CF_MERTON_H


#include "firsttoinc.h"
#include "gpbase/port.h"



CC_BEGIN_NAMESPACE(ARM)


/////////////////////////////////////////////////////////////////////////////////////////////////////
///  
///			Auxiliary Functions 
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


double merton(double S,double K,double t,double sigma, double lambda, double muJ, double sigmaJ,int nb);

double MertonOption(double S, double K, int CallPut, double t, double sigma, double lambda1, double U1, double lambda2, double U2, int nb = 20);


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
