/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file extendedsabrformula.h
 *
 *  \brief
 *
 *	\author  O.Croissant & E.Ezzine
 *	\version 1.0
 *	\date April 2004
 */
 
#ifndef _GP_CF_EXTENDEDSABRFORMULA_H
#define _GP_CF_EXTENDEDSABRFORMULA_H

CC_BEGIN_NAMESPACE(ARM)


double SABR_ComputeImpliedVol( double f,
       double K,
       double T,
       double alpha,
       double beta,
       double rho, 
       double nu,
       int SABRflag );

double SABR_ComputeImpliedVolAroundATM( double f,
       double K,
       double T,
       double alpha,
       double beta,
       double rho, 
       double nu,
       int SABRflag );


double SABR_ComputeAlphaFromSigmaATM(double f,
		                             double K,
		                             double tex,
		                             double sigma, 
		                             double beta,
		                             double rho, 
		                             double nu,
		                             int SABRflag);



void SABR_DetermineAllDerivatives( double f,
           double K,
           double T,
           double alpha,
           double beta,
           double rho, 
           double nu,
           int flag,
		   double* implvol,
           double* der_alpha, 
           double* der_beta, 
           double* der_rho, 
           double* der_nu, 
           double* der_f,
		   bool der_alphaflag, 
           bool der_betaflag, 
           bool der_rhoflag, 
           bool der_nuflag, 
           bool der_fflag);

void SABR_DetermineGenDerivatives( double f,
	double K,
	double T, 
	double alpha,
	double beta, 
	double rho, 
	double nu,
	int flag,
	double* implvol,
	double* der_alpha, 
	double* der_beta, 
	double* der_rho, 
	double* der_nu, 
	double* der_f);

void SABR_AroundAtTheMoney_DetermineDerivatives( double f,
          double K,
          double T, 
          double alpha, 
          double beta, 
          double rho, 
          double nu,
          int flag,
          double* implvol,
          double* der_alpha,
          double* der_beta, 
          double* der_rho, 
          double* der_nu, 
          double* der_f,
          bool der_alphaflag,
          bool der_betaflag,
          bool der_rhoflag, 
          bool der_nuflag, 
          bool der_fflag);

double SABR_ComputePartialDerivative( double f,
       double K,
       double T,
       double alpha,
       double beta,
       double rho, 
       double nu,
       const string& ParamStr,
       int SABRflag );

double SABR_ComputePartialDerivativeAroundATM( double f,
       double K,
       double T,
       double alpha,
       double beta,
       double rho, 
       double nu,
       const string& ParamStr,
       int SABRflag );


CC_END_NAMESPACE()
#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

