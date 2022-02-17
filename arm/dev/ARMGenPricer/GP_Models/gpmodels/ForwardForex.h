/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ForwardForex.h
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date July 2004
 */


#ifndef _INGPMODELS_FORWARDFOREX_H
#define _INGPMODELS_FORWARDFOREX_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
#include "gpbase/env.h"
#include "gpinfra/analyticpricingmodel.h"
#include "gpinfra/pricingfunctionequity.h"
#include "gpinfra/typedef.h"
#include <string>

#include <inst/forex.h>

CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE( ARM )

class ARM_PricingModel;
struct ARM_PricingContext;

//-----------------------------------------------------------------------------
// \class ARM_ForwardForex
// \brief
//  Model for a deterministic forex between two yield curves.
//  The model assumes that forex at time t equals it forward forex at time 0
//-----------------------------------------------------------------------------

class ARM_ForwardForex :	public ARM_AnalyticPricingModel,
							public ARM_PricingFunctionEquity
{
private :
    ARM_Forex           itsForex;
    /// Domestic curve is saved at ARM_Model level
    ARM_ZeroCurvePtr      itsForeignCurve;
    void CheckCurves();

public:
	/// 
	ARM_ForwardForex(const ARM_Forex& forex, 
		const ARM_ZeroCurvePtr& firstCurve,
		const ARM_ZeroCurvePtr& secondCurve);
	ARM_ForwardForex(const ARM_ForwardForex& rhs);
	virtual ~ARM_ForwardForex();
    ARM_ForwardForex& operator = (const ARM_ForwardForex& rhs);


    /// curve manangement
	void SetCurves( const ARM_ZeroCurvePtr& domesticCurve, const ARM_ZeroCurvePtr& foreignCurve);
    const ARM_Forex& GetForex() { return itsForex; }
    void SetForex(const ARM_Forex& forex) { itsForex = forex; CheckCurves(); }

    double Forward(double evalTime, double maturityTime) const;

	/// pricing part
    virtual ARM_VectorPtr DiscountFactor( 
		const string& curveName,
        double evalTime, 
		double maturityTime,
        const ARM_PricingStatesPtr& states) const;

	/// Forward function
	virtual ARM_VectorPtr Forward(
		const string& curveName, 
        double evalTime,
	    double expiryTime,
	    double settlementTime,
	    double payTime,
        const ARM_PricingStatesPtr& states) const;

	/// Call function (vectorial strike version)
	virtual ARM_VectorPtr CallVectorial(
		const string& curveName,
        double evalTime,
	    double expiryTime,
	    double settlementTime,
		const ARM_GP_Vector& strikePerState,
		int callPut,
	    double payTime,
        const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context) const;

    /// Always true because no model params
    virtual bool ValidateModelParams(const ARM_ModelParams& params) const {return true;}
	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod);

	virtual double ImpliedVol(const ARM_VanillaArg& arg) const ;

	/// convention support : calendar + gap support
	virtual string GetSettlementCalendar(const string& modelName="") const;
	virtual double GetSettlementGap(const string& modelName="") const;
	virtual int GetType() const;
    virtual size_t FactorCount() const{ return 0.0; } /// no model parameter in this non stochastic model!

    /// Standard ARM object support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LFXCL";}
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
