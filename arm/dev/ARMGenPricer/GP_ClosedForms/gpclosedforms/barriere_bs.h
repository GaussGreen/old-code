/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file barriere_bs.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_BARRIERE_BS_H
#define _GP_CF_BARRIERE_BS_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"


CC_BEGIN_NAMESPACE(ARM)

double BS_SingleBarrierOption(double S,double X, double H,double K,double T, double sig,double r,double b,int callput,int optiontype);
double BS_DoubleBarrierCall(double S,double X, double L, double U,double T, double sig,double r,double b,double delta1,double delta2,int nbterms);
double BS_DoubleBarrierPut(double S,double X, double L, double U,double T, double sig,double r,double b,double delta1,double delta2,int nbterms);
double BS_DoubleBarrierOption(double S,double X, double L, double U,double T, double sig,double r,double b,double delta1,double delta2,int callput,int nbterms);
double BS_PartialTime_Start_SingleBarrierCall(double S,double X, double H,double K,double t1,double T, double sig,double r,double b,int optiontype);
double BS_PartialTime_Start_SingleBarrierOption(double S,double X, double H,double K,double t1,double T, double sig,double r,double b,int callput,int optiontype);
double BS_PartialTime_End_SingleBarrierCall(double S,double X, double H,double K,double t1,double T, double sig,double r,double b,int optiontype);
double BS_PartialTime_End_SingleBarrierOption(double S,double X, double H,double K,double t1,double T, double sig,double r,double b,int callput,int optiontype);
double TwoAsset_Single_Barrier_Option(double S1,double S2,double X,double H,double T,double sig1, double sig2, double rho,
									double r,double b1,double b2, int callput, int optiontype);






CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

