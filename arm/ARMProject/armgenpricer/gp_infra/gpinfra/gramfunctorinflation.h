/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gramfunctorir.h
 *
 *  \brief gram functor are functor for the grammar of the generic pricer
 *
 *	\author  A. Schauly
 *	\version 1.0
 *	\date October 2003
 */



#ifndef _INGPINFRA_GRAMFUNCTORINFLATION_H
#define _INGPINFRA_GRAMFUNCTORINFLATION_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gramfunctorbase.h"
#include "pricingfunctioninflation.h"
#include "gpbase/datestrip.h"

CC_BEGIN_NAMESPACE( ARM )

/// CPI Rate
struct ARM_GP_CPI : public ARM_GramFctor
{
	/// default constructor to allow its creation!
	ARM_GP_CPI() : ARM_GramFctor() {}

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_CPI( *this ) ); }

private:
	void GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );
	long ConvCPInterpMethod( const string& InfMethod );



	/// members to store
	ARM_PricingFuncInflation*	itsModelInflation;

	string						itsInfCurveName;
	string						itsDCFLag;
	long						itsDailyInterp;
	string						itsResetLag;
	double						itsCPITime;
	double						itsEvalTime;
	
	/// name for error writting
    static string itsFuncName;
};

struct ARM_GP_InfSwapCommon : public ARM_GramFctor
{
	/// default constructor to allow its creation!
	ARM_GP_InfSwapCommon() : ARM_GramFctor() {}

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes ) = 0;
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes ) = 0;
	virtual ARM_GramFctorPtr Clone() const = 0;

protected: 
	void ComputeDateStrips( ARM_PricingModel* mod );

	/// members to store
	ARM_PricingFuncInflation*	itsModelInflation;

	string						itsZcCurveName;
	string						itsInfCurveName;
	string						itsSwapType;
	double						itsEvalTime;
	double						itsFixingTime;
	double						itsSpread;
	double						itsStrike;
	double						itsCoupon;

	int							itsFloatResetGap;
	int							itsFloatPaymGap;
	int							itsFloatFreq;
	int							itsFloatDayCount;
	int							itsFixedFreq;
	int							itsFixedDayCount;

	ARM_Date					itsStartDate;
	ARM_Date					itsEndDate;

	ARM_DateStripPtr			itsInfNumDateStrip;
	ARM_DateStripPtr			itsInfDenomDateStrip;
	ARM_DateStripPtr			itsFixedLegDateStrip;

private:

	/// name for error writting
    static string itsFuncName;
};

/// Inflation Swap Rate
struct ARM_GP_InfSwapRate : public ARM_GP_InfSwapCommon
{
	/// default constructor to allow its creation!
	ARM_GP_InfSwapRate() : ARM_GP_InfSwapCommon() {}

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );

    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );

	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_InfSwapRate( *this ) ); }

private:
	void GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
		double evalDate, vector<ARM_ExpNodePtr>& nodes );

	/// name for error writting
    static string itsFuncName;
};


/// Inflation Swap
struct ARM_GP_InfSwap : public ARM_GP_InfSwapCommon
{
	/// default constructor to allow its creation!
	ARM_GP_InfSwap() : ARM_GP_InfSwapCommon() {}

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );

    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );

	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_InfSwap( *this ) ); }

private:
	void GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
		double evalDate, vector<ARM_ExpNodePtr>& nodes );

	/// name for error writting
    static string itsFuncName;
};


CC_END_NAMESPACE()

#endif