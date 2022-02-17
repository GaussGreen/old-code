/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: SABR_SurfaceModelParams.h,v $
 * Revision 1.1  2004/07/30 09:52:19  ebenhamou
 * Initial revision
 *
 *
 */

/*! \file SABR_SurfaceModelParams.h
 *
 *  \brief 
 *	\author  O. Croissant
 *	\version 1.0
 *	\date October 2004
 */


#ifndef _INGPMODELS_SABR_SURFACEMODELPARAMS_H
#define _INGPMODELS_SABR_SURFACEMODELPARAMS_H

#include "gpbase/env.h"
#include "gpbase/assignop.h"
#include "gpbase/port.h"
#include "AnalyticModelParams.h"

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_SABR_SurfaceModelParams
// \brief Class for the model param of the SABR Model
//-----------------------------------------------------------------------------
class ARM_SABR_SurfaceModelParams : public ARM_AnalyticModelParams
{
private:
	void ValidateModelParams() const;
	bool itsIsSigmaParam;

public:
	ARM_SABR_SurfaceModelParams( const ARM_SABR_SurfaceModelParams& rhs );
	ARM_SABR_SurfaceModelParams( const ARM_ModelParamVector& params=ARM_ModelParamVector() );
	ASSIGN_OPERATOR(ARM_SABR_SurfaceModelParams)
	virtual ~ARM_SABR_SurfaceModelParams() {};

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
