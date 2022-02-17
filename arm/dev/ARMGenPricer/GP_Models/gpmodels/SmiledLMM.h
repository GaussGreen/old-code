/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
  */



#ifndef _INGPMODELS_SLMM_H
#define _INGPMODELS_SLMM_H

// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"

#include "gpmodels/SmiledMM.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_SmiledLMM : public ARM_SmiledMM
{
private:
	void			computeWeights();
	bool			FastCompute() const						{ return true;					};
	bool			IsOnSamePath(size_t i, size_t j) const	{ return ( itsWeight(i,j)==1 ); };
	inline void		setDelta( const ARM_GP_Vector& delta)	{ itsDelta = delta;				}

	virtual ARM_VectorPtr	DiscountFactorNoInterpol(	const string& curveName,
											double evalTime, 
											double maturityTime,
											const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr	ForwardDiscountFactor(const string& curveName,
											double evalTime, 
											double startTime,
											double endTime,
											const ARM_PricingStatesPtr& states) const;

	ARM_VectorPtr	ForwardDiscountFactorFromIdx(	const string& curveName,
											double evalTime,
											size_t IdxFrom,
											size_t IdxTo,
											size_t modelNb,
											const ARM_PricingStatesPtr& states) const;

public:
	ARM_SmiledLMM( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params=NULL, size_t timeStepsNb=500,size_t gridSize=501,double stdDevNb=6,bool skipPDE=false, bool allowInterpol=false, ARM_ModelParamsSmiled::CalibProxy calibProxy=ARM_ModelParamsSmiled::LocalVolatility);
	ARM_SmiledLMM( const ARM_SmiledLMM& rhs);
	virtual ~ARM_SmiledLMM();

	ARM_SmiledLMM& operator = (const ARM_SmiledLMM& rhs);
	virtual ARM_Object*		Clone() const { return new ARM_SmiledLMM(*this); };

	virtual void	setNumericalModelFitter(ARM_NumericalModelFitter*);
	virtual double	DZCDRate( double maturity, size_t i ) const;

	virtual void	MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	
	virtual ARM_VectorPtr	DiscountFactor(	const string& curveName,
											double evalTime, 
											double maturityTime,
											const ARM_PricingStatesPtr& states) const;

private:
	ARM_GP_Matrix		itsWeight;
	ARM_GP_Vector		itsDelta;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
