
/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 03/15/2006
 *
 *  basic functions for the closed form framework 
 *
 *	\file eigenvalues.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date March 2006
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>

#include "gpnumlib/gaussiananalytics.h"


/// gpbase

#include "gpbase/gplinalgtypedef.h"

#include "gpbase/stringconvert.h"

#include "gpbase/gpvector.h"

#include "gpbase/gplinalgconvert.h"

#include "gpbase/gpmatrix.h"

#include "gpbase/gpmatrixlinalg.h"

#include "gpbase/numericconstant.h"



#include "expt.h"   // for the exceptions

#include "gpclosedforms/eigenvalues.h"




CC_BEGIN_NAMESPACE(ARM)

void eigenvalues3(double rho12,double rho13,double rho23,double* e1,double* e2,double * e3)
{
	int colNb=3;
    ARM_GP_Matrix* outMatrix = NULL;
	ARM_GP_Vector eigenValues(colNb);
	ARM_GP_Matrix* TestedMatrix = new ARM_GP_Matrix(colNb,colNb) ; 
	
	(*TestedMatrix)(0,0)=1.;
	(*TestedMatrix)(0,1)=rho12;
	(*TestedMatrix)(0,2)=rho13;

	(*TestedMatrix)(1,0)=rho12;
	(*TestedMatrix)(1,1)=1.;
	(*TestedMatrix)(1,2)=rho23;

	(*TestedMatrix)(2,0)=rho13;
	(*TestedMatrix)(2,1)=rho23;
	(*TestedMatrix)(2,2)=1.;

	outMatrix = ACPTransformation( TestedMatrix, eigenValues );

	(*e1)=eigenValues[0];
	(*e2)=eigenValues[1];
	(*e3)=eigenValues[2];
	
	delete TestedMatrix;
	delete outMatrix;

	return ;
}

void eigenvalues4(double rho12,double rho13, double rho14,
				  double rho23,double rho24,double rho34,
				  double* e1,double* e2,double * e3,double* e4)
{
	int colNb=4;
    ARM_GP_Matrix* outMatrix = NULL;
	ARM_GP_Vector eigenValues(colNb);
	ARM_GP_Matrix* TestedMatrix = new ARM_GP_Matrix(colNb,colNb) ; 
	
	(*TestedMatrix)(0,0)=1.;
	(*TestedMatrix)(0,1)=rho12;
	(*TestedMatrix)(0,2)=rho13;
	(*TestedMatrix)(0,3)=rho14;

	(*TestedMatrix)(1,0)=rho12;
	(*TestedMatrix)(1,1)=1.;
	(*TestedMatrix)(1,2)=rho23;
	(*TestedMatrix)(1,3)=rho24;

	(*TestedMatrix)(2,0)=rho13;
	(*TestedMatrix)(2,1)=rho23;
	(*TestedMatrix)(2,2)=1.;
	(*TestedMatrix)(2,3)=rho34;

	(*TestedMatrix)(3,0)=rho14;
	(*TestedMatrix)(3,1)=rho24;
	(*TestedMatrix)(3,2)=rho34;
	(*TestedMatrix)(3,3)=1.;


	outMatrix = ACPTransformation( TestedMatrix, eigenValues );

	(*e1)=eigenValues[0];
	(*e2)=eigenValues[1];
	(*e3)=eigenValues[2];
	(*e4)=eigenValues[3];
	
	delete TestedMatrix;
	delete outMatrix;

	return ;
}

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
								 )
{
double e1_base,e2_base,e3_base,e4_base,e1_shifted,e2_shifted,e3_shifted,e4_shifted;
double shift=0.00000001;
eigenvalues4(rho12,rho13,rho14,rho23,rho24,rho34,&e1_base,&e2_base, &e3_base,&e4_base);
if(rho12_flag>0)
{
	eigenvalues4(rho12+shift,rho13,rho14,rho23,rho24,rho34,&e1_shifted,&e2_shifted, &e3_shifted,&e4_shifted);
	(*e1_der_rho12)=(e1_shifted-e1_base)/shift;
	(*e2_der_rho12)=(e2_shifted-e2_base)/shift;
	(*e3_der_rho12)=(e3_shifted-e3_base)/shift;
	(*e4_der_rho12)=(e4_shifted-e4_base)/shift;
}
if(rho13_flag>0)
{
	eigenvalues4(rho12,rho13+shift,rho14,rho23,rho24,rho34,&e1_shifted,&e2_shifted, &e3_shifted,&e4_shifted);
	(*e1_der_rho13)=(e1_shifted-e1_base)/shift;
	(*e2_der_rho13)=(e2_shifted-e2_base)/shift;
	(*e3_der_rho13)=(e3_shifted-e3_base)/shift;
	(*e4_der_rho13)=(e4_shifted-e4_base)/shift;
}
if(rho14_flag>0)
{
	eigenvalues4(rho12,rho13,rho14+shift,rho23,rho24,rho34,&e1_shifted,&e2_shifted, &e3_shifted,&e4_shifted);
	(*e1_der_rho14)=(e1_shifted-e1_base)/shift;
	(*e2_der_rho14)=(e2_shifted-e2_base)/shift;
	(*e3_der_rho14)=(e3_shifted-e3_base)/shift;
	(*e4_der_rho14)=(e4_shifted-e4_base)/shift;
}
if(rho23_flag>0)
{
	eigenvalues4(rho12,rho13,rho14,rho23+shift,rho24,rho34,&e1_shifted,&e2_shifted, &e3_shifted,&e4_shifted);
	(*e1_der_rho23)=(e1_shifted-e1_base)/shift;
	(*e2_der_rho23)=(e2_shifted-e2_base)/shift;
	(*e3_der_rho23)=(e3_shifted-e3_base)/shift;
	(*e4_der_rho23)=(e4_shifted-e4_base)/shift;
}
if(rho24_flag>0)
{
	eigenvalues4(rho12,rho13,rho14,rho23,rho24+shift,rho34,&e1_shifted,&e2_shifted, &e3_shifted,&e4_shifted);
	(*e1_der_rho24)=(e1_shifted-e1_base)/shift;
	(*e2_der_rho24)=(e2_shifted-e2_base)/shift;
	(*e3_der_rho24)=(e3_shifted-e3_base)/shift;
	(*e4_der_rho24)=(e4_shifted-e4_base)/shift;
}
if(rho34_flag>0)
{
	eigenvalues4(rho12,rho13,rho14,rho23,rho24,rho34+shift,&e1_shifted,&e2_shifted, &e3_shifted,&e4_shifted);
	(*e1_der_rho34)=(e1_shifted-e1_base)/shift;
	(*e2_der_rho34)=(e2_shifted-e2_base)/shift;
	(*e3_der_rho34)=(e3_shifted-e3_base)/shift;
	(*e4_der_rho34)=(e4_shifted-e4_base)/shift;
}
	(*e1)=e1_base;
	(*e2)=e2_base;
	(*e3)=e3_base;
	(*e4)=e4_base;
	
}







//// Symmetric functions of eigen values



void SymmetricPolynomOfEigenValues4(double rho12,double rho13, double rho14,
				  double rho23,double rho24,double rho34,
				  double* p1,double* p2,double * p3,double* p4)
{
	
	double e1,e2,e3,e4;
	eigenvalues4(rho12,rho13,rho14,rho23,rho24,rho34,&e1,&e2,&e3,&e4);

	(*p1)=e1+e2+e3+e4;
	(*p2)=e1*e2+e1*e3+e1*e4+e2*e3+e2*e4+e3*e4;
	(*p3)=e1*e1*e3+e1*e2*e4+e1*e3*e4+e2*e3*e4;
	(*p4)=e1*e2*e3*e4;

	return ;
}

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
								 )
{
	double e1_base,e2_base,e3_base,e4_base,e1_shifted,e2_shifted,e3_shifted,e4_shifted;
	double shift=0.00000001;
	SymmetricPolynomOfEigenValues4(rho12,rho13,rho14,rho23,rho24,rho34,&e1_base,&e2_base, &e3_base,&e4_base);
	if(rho12_flag>0)
	{
		SymmetricPolynomOfEigenValues4(rho12+shift,rho13,rho14,rho23,rho24,rho34,&e1_shifted,&e2_shifted, &e3_shifted,&e4_shifted);
		(*e1_der_rho12)=(e1_shifted-e1_base)/shift;
		(*e2_der_rho12)=(e2_shifted-e2_base)/shift;
		(*e3_der_rho12)=(e3_shifted-e3_base)/shift;
		(*e4_der_rho12)=(e4_shifted-e4_base)/shift;
	}
	if(rho13_flag>0)
	{
		SymmetricPolynomOfEigenValues4(rho12,rho13+shift,rho14,rho23,rho24,rho34,&e1_shifted,&e2_shifted, &e3_shifted,&e4_shifted);
		(*e1_der_rho13)=(e1_shifted-e1_base)/shift;
		(*e2_der_rho13)=(e2_shifted-e2_base)/shift;
		(*e3_der_rho13)=(e3_shifted-e3_base)/shift;
		(*e4_der_rho13)=(e4_shifted-e4_base)/shift;
	}
	if(rho14_flag>0)
	{
		SymmetricPolynomOfEigenValues4(rho12,rho13,rho14+shift,rho23,rho24,rho34,&e1_shifted,&e2_shifted, &e3_shifted,&e4_shifted);
		(*e1_der_rho14)=(e1_shifted-e1_base)/shift;
		(*e2_der_rho14)=(e2_shifted-e2_base)/shift;
		(*e3_der_rho14)=(e3_shifted-e3_base)/shift;
		(*e4_der_rho14)=(e4_shifted-e4_base)/shift;
	}
	if(rho23_flag>0)
	{
		SymmetricPolynomOfEigenValues4(rho12,rho13,rho14,rho23+shift,rho24,rho34,&e1_shifted,&e2_shifted, &e3_shifted,&e4_shifted);
		(*e1_der_rho23)=(e1_shifted-e1_base)/shift;
		(*e2_der_rho23)=(e2_shifted-e2_base)/shift;
		(*e3_der_rho23)=(e3_shifted-e3_base)/shift;
		(*e4_der_rho23)=(e4_shifted-e4_base)/shift;
	}
	if(rho24_flag>0)
	{
		SymmetricPolynomOfEigenValues4(rho12,rho13,rho14,rho23,rho24+shift,rho34,&e1_shifted,&e2_shifted, &e3_shifted,&e4_shifted);
		(*e1_der_rho24)=(e1_shifted-e1_base)/shift;
		(*e2_der_rho24)=(e2_shifted-e2_base)/shift;
		(*e3_der_rho24)=(e3_shifted-e3_base)/shift;
		(*e4_der_rho24)=(e4_shifted-e4_base)/shift;
	}
	if(rho34_flag>0)
	{
		SymmetricPolynomOfEigenValues4(rho12,rho13,rho14,rho23,rho24,rho34+shift,&e1_shifted,&e2_shifted, &e3_shifted,&e4_shifted);
		(*e1_der_rho34)=(e1_shifted-e1_base)/shift;
		(*e2_der_rho34)=(e2_shifted-e2_base)/shift;
		(*e3_der_rho34)=(e3_shifted-e3_base)/shift;
		(*e4_der_rho34)=(e4_shifted-e4_base)/shift;
	}
	
	(*e1)=e1_base;
	(*e2)=e2_base;
	(*e3)=e3_base;
	(*e4)=e4_base;
}

///////////////////////////////////////////////////////////////////////////////////////
///
///					Computation of the eigenvalues of  triSABR correlation matrix
///
///////////////////////////////////////////////////////////////////////////////////////

/// 
void TriSABR_Eigenvalues(
			double rho1,double rho2,double rho3,
			double rhos12, double rhos23, double rhos13,
			double rhov12, double rhov23, double rhov13,
			double rhoc12, double rhoc13,
			double rhoc21, double rhoc23,
			double rhoc31, double rhoc32,
			double* e1,double* e2,double* e3,double* e4,double* e5,double* e6)


{
	int colNb=6;
    ARM_GP_Matrix* outMatrix = NULL;
	ARM_GP_Vector eigenValues(colNb);
	ARM_GP_Matrix* TestedMatrix = new ARM_GP_Matrix(colNb,colNb) ; 
	
	(*TestedMatrix)(0,0)=1.;
	(*TestedMatrix)(0,1)=rhos12;
	(*TestedMatrix)(0,2)=rhos13;
	(*TestedMatrix)(0,3)=rho1;
	(*TestedMatrix)(0,4)=rhoc12;
	(*TestedMatrix)(0,5)=rhoc13;

	(*TestedMatrix)(1,0)=rhos12;
	(*TestedMatrix)(1,1)=1.;
	(*TestedMatrix)(1,2)=rhos23;
	(*TestedMatrix)(1,3)=rhoc21;
	(*TestedMatrix)(1,4)=rho2;
	(*TestedMatrix)(1,5)=rhoc23;

	(*TestedMatrix)(2,0)=rhos13;
	(*TestedMatrix)(2,1)=rhos23;
	(*TestedMatrix)(2,2)=1.;
	(*TestedMatrix)(2,3)=rhoc31;
	(*TestedMatrix)(2,4)=rhoc32;
	(*TestedMatrix)(2,5)=rho3;

	(*TestedMatrix)(3,0)=rho1;
	(*TestedMatrix)(3,1)=rhoc21;
	(*TestedMatrix)(3,2)=rhoc31;
	(*TestedMatrix)(3,3)=1.;
	(*TestedMatrix)(3,4)=rhov12;
	(*TestedMatrix)(3,5)=rhov13;

	(*TestedMatrix)(4,0)=rhoc12;
	(*TestedMatrix)(4,1)=rho2;
	(*TestedMatrix)(4,2)=rhoc32;
	(*TestedMatrix)(4,3)=rhov12;
	(*TestedMatrix)(4,4)=1.;
	(*TestedMatrix)(4,5)=rhov23;

	(*TestedMatrix)(5,0)=rhoc13;
	(*TestedMatrix)(5,1)=rhoc23;
	(*TestedMatrix)(5,2)=rho3;
	(*TestedMatrix)(5,3)=rhov13;
	(*TestedMatrix)(5,4)=rhov23;
	(*TestedMatrix)(5,5)=1.;


	outMatrix = ACPTransformation( TestedMatrix, eigenValues );

	(*e1)=eigenValues[0];
	(*e2)=eigenValues[1];
	(*e3)=eigenValues[2];
	(*e4)=eigenValues[3];
	(*e5)=eigenValues[4];
	(*e6)=eigenValues[5];
	
	delete TestedMatrix;
	delete outMatrix;

	return ;
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/