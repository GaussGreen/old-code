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
 
#ifndef _GP_CF_TRISPREADOPTION_LOGNORMAL_H
#define _GP_CF_TRISPREADOPTION_LOGNORMAL_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"


CC_BEGIN_NAMESPACE(ARM)




//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
///
///   TriSpreadDigitalCall :  
///	  which pays 1 if a0+a1*S1+a2*S2+a3*S3>0
///
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

double TriSpreadDigitalCall(	   double S1,double S2,double S3,
								   double v1,double v2,double v3,
								   double m1,double m2,double m3,
								   double r12,double r13,double r23,
								   double a0,double a1,double a2,double a3,
								   double T,int n);

// include handling of put !
double TriSpreadDigitalOption(double S1,double S2,double S3,
								   double v1,double v2,double v3,
								   double m1,double m2,double m3,
								   double r12,double r13,double r23,
								   double a0,double a1,double a2,double a3,
								   double T,int callput,int n);

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
///
///   TriSpreadCall 
///	  which pays (a0+a1*S1+a2*S2+a3*S3)^+
///
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

double TriSpreadCall(			   double S1,double S2,double S3,
								   double v1,double v2,double v3,
								   double m1,double m2,double m3,
								   double r12,double r13,double r23,
								   double a0,double a1,double a2,double a3,
								   double T,int n);


///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///
///
///  double digital option where the condition are : 
///			a0 + a1*S1+a2*S2+a3*S3 >0
///		and b0+b1*S1>0
///
///
//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
double TriSpreadDigitalCall2(	   double S1,double S2,double S3,
								   double v1,double v2,double v3,
								   double m1,double m2,double m3,
								   double r12,double r13,double r23,
								   double a0,double a2,double a3,
								   double b0,double b1,
								   double T,int n);

// include handling of put !

double TriSpreadDigitalOption2(	   double S1,double S2,double S3,
								   double v1,double v2,double v3,
								   double m1,double m2,double m3,
								   double r12,double r13,double r23,
								   double a0,double a2,double a3,
								   double b0,double b1,
								   double T,int callput,int n);




CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

