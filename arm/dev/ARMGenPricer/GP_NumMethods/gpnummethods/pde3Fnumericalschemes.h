/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pde1Fnumericalschemes.h
 *	\author  K Belkheir
 *	\version 1.0
 *	\date July 2006
 */

#ifndef _INGPINFRA_PDE3FNUMSCHEME_H
#define _INGPINFRA_PDE3FNUMSCHEME_H

#include "pdenumericalschemes.h"
#include "gpbase\assignop.h"

CC_BEGIN_NAMESPACE( ARM )


///////////////////////////////////////////////////////////////////////////////
/// This class defines an algorithm to solve a 
///	3 factors PDE using the Craig and Sneyd 
/// algorithm.
/// Ref: An Alternating-Direction Implicit Scheme For
/// Parabolic Equations with Mixed Derivatives
/// I.J.D Craig & A.D. Sneyd 1988.
/// The general pattern of the PDE is the following 
/// one;
/// du/dt =
///		q_xx*du/dx^2 + q_yy*du/dy^2 + q_zz*du/dz^2 +
///		q_xy*du/dxdy + q_yz*du/dydz + q_zx*du/dzdx +
///		p_x*du/dx+ p_y*du/dy + p_z*du/dz +
///		o*u
///	with a general terminal condition u(T,x,y,z) = f(x,y,z)
/// The 10 parameters q_xx, q_yy, q_zz, q_xy, q_yz, q_zx, p_x, p_y, p_z, o
/// can be time and state dependant.
///
/// The parameter theta1, theta2 and theta3 define the implicitness of the 
/// scheme for respectivly:
/// _ theta1 -> the pure diffusion term (q_xx, q_yy, q_zz)
///	_ theta2 -> the convection/drift term (p_x, p_y, p_z)
/// _ theta3 -> the actualisation term (o)
///
/// The weightind parameter lamba is used in the predictor/corrector step.
///////////////////////////////////////////////////////////////////////////////

class ARM_PDE3FCraigSneydNumericalScheme : public ARM_PDE3FNumericalScheme
{
 public:
	 enum BoundConditionType
	{
		None = 0,
		Dirichlet,
		VonNeumann,
		Belkheir
	};
		
	ARM_PDE3FCraigSneydNumericalScheme(
		size_t NX = 100, 
		size_t NY = 2, 
		size_t NZ = 2,
		double theta1 = 0.5, 
		double theta2 = 0.5, 
		double theta3 = 0.5, 
		BoundConditionType CL = VonNeumann, 
		double lambda = 0.0,
		GridType gridType = StdDev,
		const ARM_GP_Matrix&  gridData=ARM_GP_Matrix(0),
		const ARM_GP_Matrix&  schedulerData=ARM_GP_Matrix(0));
    
	ARM_PDE3FCraigSneydNumericalScheme(const ARM_PDE3FCraigSneydNumericalScheme & rhs);
	ASSIGN_OPERATOR(ARM_PDE3FCraigSneydNumericalScheme)
	
	//destructeur
    ~ARM_PDE3FCraigSneydNumericalScheme();

	virtual void Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, int toTimeIdx);
	void Induct1( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, int toTimeIdx);
	void Induct2( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, int toTimeIdx);
	void Induct3( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, int toTimeIdx);

	virtual void BuildTransitionMatrixes( const ARM_PricingModel& model, const ARM_GP_VectorPtr& timeSteps, const ARM_PricingStatesPtr& states, double* timeStepStart = NULL, double* timeStepEnd = NULL) {};
	
	/// Standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new ARM_PDE3FCraigSneydNumericalScheme(*this); }

	virtual NumericalSchemeTypes getNumericalSchemeType() const { return ARM_PDENumericalScheme::CS3F; }

	//Accessors
	const double GetLambda() const { return itsLambda; }

	
  private:
	//Common methods
	virtual ARM_PricingStatesPtr Init( const ARM_PricingModel& model );
	virtual ARM_PricingStatesPtr Init1( const ARM_PricingModel& model );
	ARM_PricingStatesPtr Init2( const ARM_PricingModel& model );
	ARM_PricingStatesPtr Init3( const ARM_PricingModel& model );
	void InitReconstructAndGeometry();
	void InitABC(int toTimeIdx);
	void InitABC2(int toTimeIdx);
	void UpdatePayoffs(ARM_GP_Vector& vec, int toTimeIdx);
	void UpdatePayoffs2(ARM_GP_Vector& vec, int toTimeIdx);
	void UpdatePayoffs3(ARM_GP_Vector& vec, int toTimeIdx);
	void UpdateIntermediatePayoffs(ARM_GP_Matrix& vec, int toTimeIdx);
	void UpdateIntermediatePayoffs2(ARM_PayoffStates& vec, int toTimeIdx);
	void UpdateIntermediatePayoffs3(ARM_GP_Matrix& vec, int toTimeIdx);
	//Specific method;
	void Reconstruct(ARM_GP_Vector &U_tilde, ARM_GP_Vector &oldU, ARM_GP_Vector &newU );
	void Reconstruct2(ARM_GP_Vector &U);

	double itsTheta1;
	double itsTheta2;
	double itsTheta3;
	double itsLambda;
	BoundConditionType itsBC;

	// Geometry coefficient

	ARM_GP_VectorPtr itsRxx;
	ARM_GP_VectorPtr itsRyy;
	ARM_GP_VectorPtr itsRzz;
	ARM_GP_VectorPtr itsRxy;
	ARM_GP_VectorPtr itsRyz;
	ARM_GP_VectorPtr itsRzx;
	ARM_GP_VectorPtr itsSx;
	ARM_GP_VectorPtr itsSy;
	ARM_GP_VectorPtr itsSz;

	ARM_GP_VectorPtr itsPrevO;
	ARM_GP_VectorPtr itsPrevPx;
	ARM_GP_VectorPtr itsPrevPy;
	ARM_GP_VectorPtr itsPrevPz;
	ARM_GP_VectorPtr itsPrevQxx;
	ARM_GP_VectorPtr itsPrevQyy;
	ARM_GP_VectorPtr itsPrevQzz;
	ARM_GP_VectorPtr itsPrevQxy;
	ARM_GP_VectorPtr itsPrevQyz;
	ARM_GP_VectorPtr itsPrevQzx;

	ARM_GP_VectorPtr itsO;
	ARM_GP_VectorPtr itsPx;
	ARM_GP_VectorPtr itsPy;
	ARM_GP_VectorPtr itsPz;
	ARM_GP_VectorPtr itsQxx;
	ARM_GP_VectorPtr itsQyy;
	ARM_GP_VectorPtr itsQzz;
	ARM_GP_VectorPtr itsQxy;
	ARM_GP_VectorPtr itsQyz;
	ARM_GP_VectorPtr itsQzx;

	ARM_GP_VectorPtr itsA1;
	ARM_GP_VectorPtr itsB1;
	ARM_GP_VectorPtr itsC1;

	ARM_GP_VectorPtr itsA2;
	ARM_GP_VectorPtr itsB2;
	ARM_GP_VectorPtr itsC2;

	ARM_GP_VectorPtr itsA3;
	ARM_GP_VectorPtr itsB3;
	ARM_GP_VectorPtr itsC3;

	ARM_GP_T_Vector<bool> itsIsBound;
	ARM_GP_T_Vector<size_t> itsReconstructFrom;
	ARM_GP_T_Vector<size_t> itsReconstructTo;

	ARM_GP_VectorPtr itsD;
	ARM_GP_VectorPtr itsGam;
	ARM_GP_VectorPtr itsUtilde;
	ARM_GP_VectorPtr itsUtilde_temp;
	ARM_GP_VectorPtr itsU_temp;
	ARM_GP_VectorPtr itsU_new;
};

CC_END_NAMESPACE()

#endif