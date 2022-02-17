/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelparamsSFRMDiag.h
 *
 *  \brief Version diag of SFRM model params
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _INGPMODELS_MODELPARAMSSFRMDIAG_H
#define _INGPMODELS_MODELPARAMSSFRMDIAG_H

#include "gpbase/port.h"
#include "ModelParamsSFRM.h"

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsSFRMDiag
// \brief Interface class for model parameters of the SFRM model
//-----------------------------------------------------------------------------
class ARM_ModelParamsSFRMDiag : public ARM_ModelParamsSFRM
{
public:
	ARM_ModelParamsSFRMDiag( const ARM_ModelParamVector& params, ARM_IRIndex* index, size_t factorsNb = 1 );
	ARM_ModelParamsSFRMDiag( const ARM_ModelParamsSFRMDiag& rhs);
	ARM_ModelParamsSFRMDiag& operator=( const ARM_ModelParamsSFRMDiag& rhs );
	virtual ~ARM_ModelParamsSFRMDiag();
	
	/// pricing function
    virtual double VolatilityFunction(double t, double T) const;
	virtual double IntegratedLocalVariance(double s, double t) const;
	virtual double MaturityTerm(double T) const;
	virtual double MaturityTermSquared(double T) const;
	virtual double VolatilitySpotSquared(double s) const;

    virtual double VarianceToTime(double var,double minTime,double maxTime) const;

    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0);
	virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model ,int factorNb=0);

	/// Standard ARM object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

