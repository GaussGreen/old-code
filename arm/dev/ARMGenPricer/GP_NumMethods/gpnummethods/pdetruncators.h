/*!
 *
 * Copyright (c) CDC IXIS CM July 2005 Paris
 *
 *	\file pdetruncators.h
 *	\author  A Schauly
 *	\version 1.0
 *	\date July 2005
 */

#ifndef _INGPINFRA_PDETRUNCATORS_H
#define _INGPINFRA_PDETRUNCATORS_H

#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gpinfra/typedef.h"
#include "gpbase/typedef.h"

CC_BEGIN_NAMESPACE( ARM )

class PDE_Truncator : public ARM_RootObject
{
private: 
	ARM_GP_T_Vector<size_t> itsToTimeSizes;

public:

	/// Constructors/destructor
	PDE_Truncator() : itsToTimeSizes(0) {}
	PDE_Truncator( const PDE_Truncator& rhs ) : itsToTimeSizes(rhs.itsToTimeSizes) {}
	virtual ~PDE_Truncator() {}

	/// Standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const = NULL;
	virtual ARM_Object* Clone() const = NULL;

	/// Methods for truncation
	virtual void TruncateTranstionMatrixes( size_t toTimeIdx, const ARM_GP_VectorPtr& UpperTerms, const ARM_GP_VectorPtr& LowerTerms, const ARM_GP_VectorPtr& DiagTerms, 
		const ARM_GP_VectorPtr& newUpperLeftLimitConditions, const ARM_GP_VectorPtr& newLowerRightLimitConditions );
	virtual void TruncatePricingStates( size_t toTimeIdx, const ARM_PricingStatesPtr& states ) = NULL;

	/// Accessors
	inline const ARM_GP_T_Vector<size_t>& getToTimeSizes() const { return itsToTimeSizes; }
	inline const size_t getToTimeSize( size_t timeIdx ) const { return itsToTimeSizes[timeIdx]; }
	inline void setToTimeSizes( const ARM_GP_T_Vector<size_t>& toTimeSizes ) { itsToTimeSizes = ARM_GP_T_Vector<size_t>(toTimeSizes); }
};


class PDE_DummyTruncator : public PDE_Truncator
{
public:
	/// Constructors/destructor
	PDE_DummyTruncator() : PDE_Truncator() {}
	PDE_DummyTruncator( const PDE_DummyTruncator& rhs ) : PDE_Truncator(rhs){}
	virtual ~PDE_DummyTruncator() {};

	/// Standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new PDE_DummyTruncator(*this); }

	/// Methods for truncation
	virtual void TruncatePricingStates( size_t toTimeIdx, const ARM_PricingStatesPtr& states ) {}
	virtual void TruncateTranstionMatrixes( size_t toTimeIdx, const ARM_GP_VectorPtr& UpperTerms, const ARM_GP_VectorPtr& LowerTerms, const ARM_GP_VectorPtr& DiagTerms, 
		const ARM_GP_VectorPtr& newUpperLeftLimitConditions, const ARM_GP_VectorPtr& newLowerRightLimitConditions ) {}
};

class PDE1F_Truncator : public PDE_Truncator
{
private:
	void TruncateVector( size_t toTimeIdx, ARM_GP_Vector& vector );
	void TruncateMatrix( size_t toTimeIdx, ARM_GP_Matrix& matrix );

public:
	/// Constructors/destructor
	PDE1F_Truncator() : PDE_Truncator() {}
	PDE1F_Truncator( const PDE_DummyTruncator& rhs ) : PDE_Truncator(rhs) {}
	virtual ~PDE1F_Truncator() {}

	/// Standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new PDE1F_Truncator(*this); }

	/// Methods for truncation
	virtual void TruncatePricingStates( size_t toTimeIdx, const ARM_PricingStatesPtr& states );
};


class PDE2F_Truncator : public PDE_Truncator
{
private:
	void TruncateVector( size_t toTimeIdx, ARM_GP_Vector& vector );
	void TruncateMatrix( size_t toTimeIdx, ARM_GP_Matrix& matrix );

public:
	/// Constructors/destructor
	PDE2F_Truncator() : PDE_Truncator() {}
	PDE2F_Truncator( const PDE_DummyTruncator& rhs ) : PDE_Truncator(rhs) {}
	virtual ~PDE2F_Truncator() {}

	/// Standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new PDE2F_Truncator(*this); }

	/// Methods for truncation
	virtual void TruncatePricingStates( size_t toTimeIdx, const ARM_PricingStatesPtr& states );
};

CC_END_NAMESPACE()

#endif