/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pde1Fnumericalschemes.h
 *	\author  A Schauly
 *	\version 1.0
 *	\date July 2005
 */

#ifndef _INGPINFRA_PDE2FNUMSCHEME_H
#define _INGPINFRA_PDE2FNUMSCHEME_H

#include "pdenumericalschemes.h"

CC_BEGIN_NAMESPACE( ARM )


class ARM_PDE2FExplicitNumericalScheme : public ARM_PDE2FNumericalScheme
{
private:
	/// Matrixes are PentaDiagonal
	/// Predictor matrixes
	ARM_VectorPtrVector itsXPDiagElts;
	ARM_VectorPtrVector itsYPDiagElts;
	ARM_VectorPtrVector itsXPOperatorUElts;
	ARM_VectorPtrVector itsXPOperatorLElts;
	ARM_VectorPtrVector itsYPOperatorUElts;
	ARM_VectorPtrVector itsYPOperatorDElts;

	/// Corrector matrixes
	ARM_VectorPtrVector itsXCDiagElts;
	ARM_VectorPtrVector itsYCDiagElts;
	ARM_VectorPtrVector itsXCOperatorUElts;
	ARM_VectorPtrVector itsXCOperatorLElts;
	ARM_VectorPtrVector itsYCOperatorUElts;
	ARM_VectorPtrVector itsYCOperatorDElts;

	/// CrossCoeffs used by both
	ARM_VectorPtr itsCrossCoeffs;

	ARM_GP_VectorPtr itsPayoffs;
	ARM_GP_VectorPtr itsPayoffs2;

	void UpdatePayoffs(ARM_GP_Vector& vec, int toTimeIdx);
	void UpdatePayoffs(ARM_GP_Matrix& vec, int toTimeIdx);

	void UpdateVectorWithVectors( ARM_GP_Vector& prevVec, ARM_GP_Vector& nextVec, ARM_GP_VectorPtr& XDiag, ARM_GP_VectorPtr& YDiag, ARM_GP_VectorPtr& UU, ARM_GP_VectorPtr& LL, ARM_GP_VectorPtr& U, ARM_GP_VectorPtr& L, double CorrelTerm );

	/// New Boundaries when matrix size changes
	void TruncateMatrixesAndStates( ARM_PricingStatesPtr& states, int toTimeIdx );
	/// Vector used by truncation.
	/// ARM_VectorPtr NullVector;


public:
	/// Reimplemented Especially for Explicit Numerical Scheme
	virtual void Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, int toTimeIdx);

	/// With truncation
	/// virtual void BuildTransitionMatrixes( const ARM_GP_VectorPtr& timeSteps,const ARM_GP_MatrixPtr& relativeDrifts,const ARM_GP_MatrixPtr& absoluteDrifts,const ARM_GP_MatrixPtr& vols,const ARM_GP_MatrixPtr& d1Vols,const ARM_GP_MatrixPtr& correls, const ARM_GP_T_Vector<size_t>& TruncationSizeVector, const ARM_PricingStatesPtr& states);
	/// Without truncation
	virtual void BuildTransitionMatrixes( const ARM_PricingModel& model, const ARM_GP_VectorPtr& timeSteps, const ARM_PricingStatesPtr& states, double* timeStepStart = NULL, double* timeStepEnd = NULL );

	//// Constructors
	ARM_PDE2FExplicitNumericalScheme( size_t DiscretizationPointsNb, double StdDevNb, PDE_Truncator* Truncator ) : ARM_PDE2FNumericalScheme( DiscretizationPointsNb, StdDevNb, Truncator ) {}
	ARM_PDE2FExplicitNumericalScheme() : ARM_PDE2FNumericalScheme() {}
	ARM_PDE2FExplicitNumericalScheme( const ARM_PDE2FExplicitNumericalScheme& rhs ) : ARM_PDE2FNumericalScheme( rhs ) {}

	/// Standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new ARM_PDE2FExplicitNumericalScheme(*this); }

	virtual NumericalSchemeTypes getNumericalSchemeType() const { return ARM_PDENumericalScheme::Explicit2F; }
};


class ARM_PDE2FADINumericalScheme : public ARM_PDE2FNumericalScheme
{
private: 
	/// The three tridiagonal matrixes used
	ARM_VectorPtrVector itsRBYUOperators;
	ARM_VectorPtrVector itsLAXUOperators;
	ARM_VectorPtrVector itsLBYUOperators;

	ARM_VectorPtrVector itsRBYDOperators;
	ARM_VectorPtrVector itsLAXDOperators;
	ARM_VectorPtrVector itsLBYDOperators;

	ARM_VectorPtrVector itsRBYDiagOperators;
	ARM_VectorPtrVector itsLAXDiagOperators;
	ARM_VectorPtrVector itsLBYDiagOperators;

	ARM_GP_VectorPtr itsCorrelFactors;

	/// Storing matrixes inverses
	ARM_GP_VectorPtr itsALNormalizationTerms;
	ARM_GP_VectorPtr itsALOtherCoeffs;
	ARM_GP_VectorPtr itsBLNormalizationTerms;
	ARM_GP_VectorPtr itsBLOtherCoeffs;

	ARM_GP_VectorPtr itsTimeSteps;

	/// For internal computation
	double** vector1;
	double** vector2;

	/// Code factorization for payoffs upadte
	void UpdatePayoffs(std::vector<double>& vec,const std::vector<double>& PredictorCorrector, int toTimeIdx);
	void UpdatePayoffs(ARM_GP_Matrix& vec, int payoffStateIdx, int toTimeIdx);


	/// Same as usual
	void PrecomputeMatrixesForInversion( const ARM_GP_VectorPtr& LowerTerm, const ARM_GP_VectorPtr& DiagTerm, const ARM_GP_VectorPtr& UpperTerm, 
		ARM_GP_VectorPtr& NormalizationTerms, ARM_GP_VectorPtr& OtherCoeff);

	/// Storing n+1 states for Predictor/Corrector stuff
	ARM_MatrixPtrVector nextPayoffStates;
	ARM_VectorPtrVector nextPayoffStatesPayoffs;
	ARM_VectorPtrVector nextPayoffs;

	/// Storing 3*V^n - V^n+1
	ARM_MatrixPtrVector diffPayoffStates;
	ARM_VectorPtrVector diffPayoffStatesPayoffs;
	ARM_VectorPtrVector diffPayoffs;

	/// Filling nextPayoffStuff and diffPayoffStruff
	void PrecomputeCrossedTerms( const ARM_PricingStatesPtr& states );

public:
	//// Constructors
	ARM_PDE2FADINumericalScheme() : ARM_PDE2FNumericalScheme(), vector1(NULL), vector2(NULL) { }
	ARM_PDE2FADINumericalScheme( size_t DiscretizationPointsNb, double StdDevNb, PDE_Truncator* Truncator ) : ARM_PDE2FNumericalScheme( DiscretizationPointsNb, StdDevNb, Truncator ), vector1(NULL), vector2(NULL) {}
	ARM_PDE2FADINumericalScheme( const ARM_PDE2FADINumericalScheme & rhs ) : ARM_PDE2FNumericalScheme( rhs ), vector1(NULL), vector2(NULL) { }
	~ARM_PDE2FADINumericalScheme();
	virtual void Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, int toTimeIdx);

	/// With truncation
	/// virtual void BuildTransitionMatrixes( const ARM_GP_VectorPtr& timeSteps,const ARM_GP_MatrixPtr& relativeDrifts,const ARM_GP_MatrixPtr& absoluteDrifts,const ARM_GP_MatrixPtr& vols,const ARM_GP_MatrixPtr& d1Vols,const ARM_GP_MatrixPtr& correls, const ARM_GP_T_Vector<size_t>& TruncationSizeVector, const ARM_PricingStatesPtr& states);
	/// Without truncation
	virtual void BuildTransitionMatrixes( const ARM_PricingModel& model, const ARM_GP_VectorPtr& timeSteps, const ARM_PricingStatesPtr& states, double* timeStepStart = NULL, double* timeStepEnd = NULL );

	/// Standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new ARM_PDE2FADINumericalScheme(*this); }

	virtual NumericalSchemeTypes getNumericalSchemeType() const { return ARM_PDENumericalScheme::ADI2F; }
};

CC_END_NAMESPACE()

#endif