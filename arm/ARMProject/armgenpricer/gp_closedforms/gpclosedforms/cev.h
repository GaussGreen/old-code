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
 
#ifndef _GP_CF_CEV_H
#define _GP_CF_CEV_H

#include "firsttoinc.h"
#include "gpbase/port.h"


CC_BEGIN_NAMESPACE(ARM)


double CEV_Call(double f, double K, double T,double drift, double sig, double alpha,int nbsteps);
double CEV_DoubleBarrierCalldeltakm(double k,double m,double a,double b);
double CEV_DoubleBarrierCalldeltakm2(double k,double m,double a,double b);
double CEV_DoubleBarrierCallComputekn(double beta,double mu, double delta, double L, double U, int n);
double CEV_DoubleBarrierCallNormalisation (double K, double beta,double mu, double delta, double L, double U, double kn);
double CEV_DoubleBarrierCallComputecn (double K, double beta,double mu, double delta, double L, double U, double kn, double Nn);
double CEV_DoubleBarrierCallComputeEigenValue (double beta,double mu, double delta, double L, double U, double kn);
double CEV_DoubleBarrierCallComputeEigenFunction (double S, double beta,double mu, double delta, double L, double U, double kn, double Nn);
double CEV_DoubleBarrierCall(double S,double K,double T,double r,double beta,double mu,double delta,double L,double U,int  n);
double CEV_SingleBarrierUpandOutCallComputekn(double beta,double mu,double delta,double U,double i);
double CEV_SingleBarrierUpandOutCallNormalisation(double K,double beta,double mu,double delta,double U,double kn);
double CEV_SingleBarrierUpandOutCallComputecn(double K,double beta,double mu,double delta,double U,double kn,double Nn);
double CEV_SingleBarrierUpandOutCallComputeEigenValue(double beta,double mu,double delta,double U,double kn);
double CEV_SingleBarrierUpandOutCallComputeEigenFunction(double S,double beta,double mu,double delta,double U,double kn,double Nn);
double CEV_SingleBarrierUpandOutCall(double S,double K,double T,double r,double beta,double mu,double delta,double U,int  n);


double CEV_SingleBarrierDownandOutPutComputekn(double beta,double mu,double delta,double L,double n);
double CEV_SingleBarrierDownandOutPutNormalisation(double K,double beta,double mu,double delta,double L,double kn);
double CEV_SingleBarrierDownandOutPutComputecn(double K,double beta,double mu,double delta,double L,double kn,double Nn);
double CEV_SingleBarrierDownandOutPutComputeEigenValue(double beta,double mu,double delta,double L,double kn);
double CEV_SingleBarrierDownandOutPutComputeEigenFunction(double S,double beta,double mu,double delta,double L,double kn,double Nn);

double CEV_SingleBarrierDownandOutPut(double S,double K,double T,double r,double beta,double mu,double delta,double L,int  n);




CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


