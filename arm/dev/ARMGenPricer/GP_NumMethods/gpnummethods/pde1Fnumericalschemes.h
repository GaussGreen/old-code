/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pde1Fnumericalschemes.h
 *	\author  A Schauly
 *	\version 1.0
 *	\date July 2005
 */

#ifndef _INGPINFRA_PDE1FNUMSCHEME_H
#define _INGPINFRA_PDE1FNUMSCHEME_H

#include "pdenumericalschemes.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_PDE1FExplicitNumericalScheme : public ARM_PDE1FNumericalScheme
{
private:
	/// Matrixes are Tridiagonal
	ARM_VectorPtrVector itsDiagonalElts;
	ARM_VectorPtrVector itsUpperDiagonalElts;
	ARM_VectorPtrVector itsLowerDiagonalElts;

	/// New limit conditions (for truncation)
	/// ARM_VectorPtrVector itsULLimitConditions;
	/// ARM_VectorPtrVector itsDRLimitConditions;

protected:

public:
	/// Reimplemented Especially for Explicit Numerical Scheme
	virtual void Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, int toTimeIdx);

	/// With truncation
	/// virtual void BuildTransitionMatrixes( const ARM_GP_VectorPtr& timeSteps,const ARM_GP_MatrixPtr& relativeDrifts,const ARM_GP_MatrixPtr& absoluteDrifts,const ARM_GP_MatrixPtr& vols,const ARM_GP_MatrixPtr& d1Vols,const ARM_GP_MatrixPtr& correls, const ARM_GP_T_Vector<size_t>& TruncationSizeVector, const ARM_PricingStatesPtr& states);
	/// Without truncation
	virtual void BuildTransitionMatrixes( const ARM_PricingModel& model, const ARM_GP_VectorPtr& timeSteps, const ARM_PricingStatesPtr& states, double* timeStepStart = NULL, double* timeStepEnd = NULL);

	//// Constructors
	ARM_PDE1FExplicitNumericalScheme( size_t DiscretizationPointsNb, double StdDevNb, PDE_Truncator* Truncator ) : ARM_PDE1FNumericalScheme( DiscretizationPointsNb, StdDevNb, Truncator ) {}
	ARM_PDE1FExplicitNumericalScheme() : ARM_PDE1FNumericalScheme() {}
	ARM_PDE1FExplicitNumericalScheme( const ARM_PDE1FExplicitNumericalScheme& rhs ) : ARM_PDE1FNumericalScheme( rhs ) {}

	/// Standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new ARM_PDE1FExplicitNumericalScheme(*this); }

	virtual NumericalSchemeTypes getNumericalSchemeType() const { return ARM_PDENumericalScheme::Explicit1F; }
};

class ARM_PDE1FCrankNicholsonNumericalScheme : public ARM_PDE1FNumericalScheme
{
private:
	/// Matrixes are Tridiagonal (numerator and denominator)
	ARM_VectorPtrVector itsNumDiagonalElts;
	ARM_VectorPtrVector itsNumUpperDiagonalElts;
	ARM_VectorPtrVector itsNumLowerDiagonalElts;

	ARM_VectorPtrVector itsDenomDiagonalElts;
	ARM_VectorPtrVector itsDenomUpperDiagonalElts;
	ARM_VectorPtrVector itsDenomLowerDiagonalElts;

	/// New limit conditions for truncation./ Uncommment if you want to use it
	ARM_VectorPtrVector itsNumULLimitConditions;
	ARM_VectorPtrVector itsNumDRLimitConditions;
	ARM_VectorPtrVector itsDenomULLimitConditions;
	ARM_VectorPtrVector itsDenomDRLimitConditions;

	/// For faster computation, Part of the Thomas Algorithm computations are stored
	size_t itsLastMatrixIndex;
	ARM_VectorPtr itsNormalizationTerms;
	ARM_VectorPtr itsOtherCoeff;
	void UpdateDiagsForInversion( size_t MatrixIndex );
public:
	/// Reimplemented Especially CN Scheme
	virtual void Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, int toTimeIdx);
	/// With truncation
	// virtual void BuildTransitionMatrixes( const ARM_GP_VectorPtr& timeSteps,const ARM_GP_MatrixPtr& relativeDrifts,const ARM_GP_MatrixPtr& absoluteDrifts,const ARM_GP_MatrixPtr& vols,const ARM_GP_MatrixPtr& d1Vols,const ARM_GP_MatrixPtr& correls, const ARM_GP_T_Vector<size_t>& TruncationSizeVector, const ARM_PricingStatesPtr& states);
	/// Without truncation
	virtual void BuildTransitionMatrixes( const ARM_PricingModel& model, const ARM_GP_VectorPtr& timeSteps, const ARM_PricingStatesPtr& states, double* timeStepStart = NULL, double* timeStepEnd = NULL);

	void UpdateVectorWithMatrixByInverse( const ARM_GP_VectorPtr& DiagTerm, const ARM_GP_VectorPtr& UpperTerm, const ARM_GP_VectorPtr& LowerTerm, ARM_GP_Matrix& TransposedVector );
	void UpdateVectorWithMatrixByInverse( const ARM_GP_VectorPtr& DiagTerm, const ARM_GP_VectorPtr& UpperTerm, const ARM_GP_VectorPtr& LowerTerm, ARM_GP_Vector& Vector );
	void UpdateVectorWithMatrixByInverse( const ARM_GP_VectorPtr& DiagTerm, const ARM_GP_VectorPtr& UpperTerm, const ARM_GP_VectorPtr& LowerTerm, ARM_MemPool_Matrix::iterator vecBegin );

	//// Constructors
	ARM_PDE1FCrankNicholsonNumericalScheme( size_t DiscretizationPointsNb, double StdDevNb, PDE_Truncator* Truncator ) : ARM_PDE1FNumericalScheme( DiscretizationPointsNb, StdDevNb, Truncator ) {}
	ARM_PDE1FCrankNicholsonNumericalScheme() : ARM_PDE1FNumericalScheme() {}
	ARM_PDE1FCrankNicholsonNumericalScheme( const ARM_PDE1FCrankNicholsonNumericalScheme& rhs ) : ARM_PDE1FNumericalScheme( rhs ) {}

	/// Standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new ARM_PDE1FCrankNicholsonNumericalScheme(*this); }

	virtual NumericalSchemeTypes getNumericalSchemeType() const { return ARM_PDENumericalScheme::CN1F; }
};

CC_END_NAMESPACE()

#endif