/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelparamsEqHybrid.h
 *  \brief general class for model params in the equity hybrid model
 *	\author  A. Schauly
 *	\version 1.0
 *	\date March 2005
 */

#ifndef _INGPMODELS_MODELPARAMSEQHYBRID_H
#define _INGPMODELS_MODELPARAMSEQHYBRID_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/curve.h"

/// gpinfra
#include "gpinfra/modelparams.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/typedef.h"

/// gpmodels
#include "gpmodels/typedef.h"

/// stl
#include <vector>

CC_USING_NS(std,vector)


CC_BEGIN_NAMESPACE( ARM )


//-----------------------------------------------------------------------------
// \class ARM_ModelParamsInflationEquity
// \brief Interface class for model parameters of the Equity Hybrid model
//-----------------------------------------------------------------------------
class ARM_ModelParamsInflationEquity : public ARM_ModelParams 
{
private:
	/// validation of the model params (not virtual
	void ValidateModelParams() const;
	ARM_CurvePtr itsSquaredIntegratedVol;

public:
	ARM_ModelParamsInflationEquity( const ARM_ModelParamVector& params );
	ARM_ModelParamsInflationEquity( const ARM_ModelParamsInflationEquity& rhs);
	ARM_ModelParamsInflationEquity& operator=( const ARM_ModelParamsInflationEquity& rhs );
	virtual ~ARM_ModelParamsInflationEquity();

	/// Standard ARM object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const;

    inline ARM_Curve* GetVolCurve() const { return ((ARM_CurveModelParam& )GetModelParam(ARM_ModelParamType::Volatility)).GetCurve(); }

	/// How many factors?
	virtual size_t FactorCount() const { return 1; }

	/// pricing function
	/// exp(+Sum on [a,b] of sgm^2(s) ds)
	double IntegratedVolatility( double a, double b ) const;

	/// returns the multiplier (alpha)
	double getMultiplier() const;

    /// initialise the news parameters to calibrate (not pure virtual to allow the use of 
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) { PreComputeSquaredIntegratedVol(); };
	virtual iterator SetModelParamValue( int paramType, size_t i, double value, double time, double tenor = 0.0 );
	void PreComputeSquaredIntegratedVol();
	void InitSquaredIntegratedVol();

	double IntegratedLocalVarianceFromZero(double t) const;

    /// Variance in [a,b] of the state variable
    virtual double StateLocalVariance(double a,double b) const;
	/// Variance in [a,b] of the state drift
	virtual double StateLocalDrift(double a,double b) const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

