/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: SABR_ModelParams.h,v $
 * Revision 1.1  2004/07/30 09:52:19  emezzine
 * Initial revision
 *
 *
 */

/*! \file SABR_ModelParams.h
 *
 *  \brief 
 *	\author  E.Ezzine
 *	\version 1.0
 *	\date February 2005
 */


#ifndef _INGPMODELS_SABR_MODELPARAMS_H
#define _INGPMODELS_SABR_MODELPARAMS_H

#include "gpbase/env.h"
#include "gpbase/assignop.h"
#include "gpbase/port.h"
#include "AnalyticModelParams.h"

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_SABR_ModelParams
// \brief Class for the model param of the SABR Model
//-----------------------------------------------------------------------------
class ARM_SABR_ModelParams : public ARM_AnalyticModelParams
{
private:
	void ValidateModelParams() const;
	bool itsIsSigmaParam;

public:
	ARM_SABR_ModelParams( const ARM_SABR_ModelParams& rhs );
	ARM_SABR_ModelParams( const ARM_ModelParamVector& params=ARM_ModelParamVector() );
	ASSIGN_OPERATOR(ARM_SABR_ModelParams)
	virtual ~ARM_SABR_ModelParams() {};

	/// How many factors?
    virtual size_t FactorCount() const { return 1; };

    double ImpliedVol(double undelying, 
        double strike, 
        double time,
        double tenor,
        int type,
        int IntegrationStep);
    double PartialDerivative(double underlying, 
       double strike, 
       double time,
       double tenor,
       int type,
       int IntegrationStep,
	   int modelParamType);

	/// calibration part
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model ,int factorNb=0);

	/// Standard ARM object support
	virtual ARM_Object* Clone() const;
	    /// Variance in [a,b] of the state variable
    virtual double StateLocalVariance(double a,double b) const;
	/// Variance in [a,b] of the state drift
	virtual double StateLocalDrift(double a,double b) const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
