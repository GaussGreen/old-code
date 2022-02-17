/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
  */



#ifndef _INGPMODELS_SSMMCOL_H
#define _INGPMODELS_SSMMCOL_H

// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"

#include "gpmodels/SmiledMM.h"
#include "gpmodels/SmiledSMMcol.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_SmiledSMMcol : public ARM_SmiledMM
{
private:
	inline void		setTheta( const ARM_GP_Vector& theta)	{ itsTheta = theta;}
	inline void		setLinks(	const ARM_IntVector& prev,
								const ARM_IntVector& next,
								const ARM_IntVector& last,
								const ARM_IntVector& prev2,
								const ARM_IntVector& next2,
								const ARM_IntVector& last2,
								const ARM_IntVector& extra) { 
										itsPreviousRate = prev;
										itsNextRate = next;
										itsLastRate = last;
										itsPreviousRate2 = prev2;
										itsNextRate2 = next2;
										itsLastRate2 = last2;
										itsExtraLibor = extra;}

	inline	bool isDiag(int i)	const		{return (i==-1?false:(itsNextRate2[i]==itsExtraLibor[i]));}
	inline	bool isExtraLibor(int i) const	{return (i==-1?false:(itsLastRate[i]==-1));}

public:
	ARM_SmiledSMMcol( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params=NULL, size_t timeStepsNb=500,size_t gridSize=501,double stdDevNb=6,bool skipPDE=false, bool allowInterpol=false, ARM_ModelParamsSmiled::CalibProxy calibProxy=ARM_ModelParamsSmiled::LocalVolatility,bool cache=false);
	ARM_SmiledSMMcol( const ARM_SmiledSMMcol& rhs);
	virtual ~ARM_SmiledSMMcol();

	ARM_SmiledSMMcol& operator = (const ARM_SmiledSMMcol& rhs);
	virtual ARM_Object*		Clone() const { return new ARM_SmiledSMMcol(*this); };

	virtual void	setNumericalModelFitter(ARM_NumericalModelFitter*);
	void	checkNumericalModelFitter(ARM_NumericalModelFitter*);

	virtual void	computeDWeights();
	virtual double	DZCDRate( double maturity, size_t i ) const;

	virtual void	MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	void	MCModelStatesFromToNextTime_std(ARM_PricingStatesPtr& states,int timeIndex) const;
	void	MCModelStatesFromToNextTime_evalOnly(ARM_PricingStatesPtr& states,int timeIndex) const;

	virtual ARM_VectorPtr	DiscountFactor(	const string& curveName,
											double evalTime, 
											double maturityTime,
											const ARM_PricingStatesPtr& states) const;

	ARM_VectorPtr	ForwardDiscountFactorFromIdx(	const string& curveName,
													double evalTime,
													size_t IdxFrom,
													const ARM_PricingStatesPtr& states) const;
	ARM_VectorPtr	ForwardDiscountFactorFromIdx_std(	const string& curveName,
											double evalTime,
											size_t IdxFrom,
											const ARM_PricingStatesPtr& states) const;
	ARM_VectorPtr	ForwardDiscountFactorFromIdx_cache(	const string& curveName,
											double evalTime,
											size_t IdxFrom,
											const ARM_PricingStatesPtr& states) const;

	virtual	string			toString(const string& indent="",const string& nextIndent="") const;

private:
	ARM_GP_Vector	itsTheta;
	ARM_IntVector	itsPreviousRate;
	ARM_IntVector	itsNextRate;
	ARM_IntVector	itsLastRate;
	ARM_IntVector	itsPreviousRate2;
	ARM_IntVector	itsNextRate2;
	ARM_IntVector	itsLastRate2;
	ARM_IntVector	itsExtraLibor;
	ARM_GP_Matrix	itsWeight;
	ARM_GP_Matrix	itsWeightZC;

};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/