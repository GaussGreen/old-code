/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file eigenvalues.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date March 2006
 */
 
#ifndef _GP_CF_EIGENVALUES_H
#define _GP_CF_EIGENVALUES_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"


CC_BEGIN_NAMESPACE(ARM)


void eigenvalues4(double rho12,double rho13, double rho14,double rho23,double rho24,double rho34,double* e1,double* e2,double * e3,double* e4);

void eigenvalues3(double rho12,double rho13, double rho23,double* e1,double* e2,double * e3);

void eigenvalues4_AllDerivatives(double rho12,double rho13, double rho14,
								 double rho23,double rho24,double rho34,
								 double* e1,double* e2,double * e3,double* e4,
								 double* e1_der_rho12,double* e2_der_rho12,double * e3_der_rho12,double* e4_der_rho12,
								 double* e1_der_rho13,double* e2_der_rho13,double * e3_der_rho13,double* e4_der_rho13,
								 double* e1_der_rho14,double* e2_der_rho14,double * e3_der_rho14,double* e4_der_rho14,
								 double* e1_der_rho23,double* e2_der_rho23,double * e3_der_rho23,double* e4_der_rho23,
								 double* e1_der_rho24,double* e2_der_rho24,double * e3_der_rho24,double* e4_der_rho24,
								 double* e1_der_rho34,double* e2_der_rho34,double * e3_der_rho34,double* e4_der_rho34,
								 int rho12_flag,int rho13_flag, int rho14_flag,
								 int rho23_flag,int rho24_flag, int rho34_flag
								 );


void SymmetricPolynomOfEigenValues4(double rho12,double rho13, double rho14,
				  double rho23,double rho24,double rho34,
				  double* p1,double* p2,double * p3,double* p4);


void SymmetricPolynomOfEigenValues4_AllDerivatives(double rho12,double rho13, double rho14,
								 double rho23,double rho24,double rho34,
								 double* e1,double* e2,double * e3,double* e4,
								 double* e1_der_rho12,double* e2_der_rho12,double * e3_der_rho12,double* e4_der_rho12,
								 double* e1_der_rho13,double* e2_der_rho13,double * e3_der_rho13,double* e4_der_rho13,
								 double* e1_der_rho14,double* e2_der_rho14,double * e3_der_rho14,double* e4_der_rho14,
								 double* e1_der_rho23,double* e2_der_rho23,double * e3_der_rho23,double* e4_der_rho23,
								 double* e1_der_rho24,double* e2_der_rho24,double * e3_der_rho24,double* e4_der_rho24,
								 double* e1_der_rho34,double* e2_der_rho34,double * e3_der_rho34,double* e4_der_rho34,
								 int rho12_flag,int rho13_flag, int rho14_flag,
								 int rho23_flag,int rho24_flag, int rho34_flag
								 );

void TriSABR_Eigenvalues(
			double rho1,double rho2,double rho3,
			double rhos12, double rhos23, double rhos13,
			double rhov12, double rhov23, double rhov13,
			double rhoc12, double rhoc13,
			double rhoc21, double rhoc23,
			double rhoc31, double rhoc32,
			double* e1,double* e2,double* e3,double* e4,double* e5,double* e6);




CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/