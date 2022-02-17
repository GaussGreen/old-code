/*
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *
 * $Log: gramfunctoreq.h,v $
 * Revision 1.1  2004/28/06 18:53:24  ebenhamou
 *  filtering of blank
 *
 *
 */


/*! \file gramfunctoreq.h
 *   
 *  \brief gram functor file for equity like models
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date June 2004
 */



#ifndef _INGPINFRA_GRAMFUNCTOREQ_H
#define _INGPINFRA_GRAMFUNCTOREQ_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gramfunctorbase.h"

#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_PricingFunctionEquity;
class ARM_EqFxBase;//for the call FX
class ARM_HybridIRFX;//for the Exotic Equity payoffs

//Spot operator
struct ARM_GP_SpotFctor : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );

    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );

	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_SpotFctor( *this ) ); }

private:
	void GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( double evalDate, ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	string		                itsModelName;
	double                      itsEvalTime;
	double                      itsSettlementTime;
	ARM_PricingFunctionEquity*  itsModelEquity;

	/// name for error writting
    static string itsFuncName;
};


//Fwd operator
struct ARM_GP_FwdFctor : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );

    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );

	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_FwdFctor( *this ) ); }

private:
	void GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	string                      itsModelName;
	double                      itsEvalTime;
	double						itsExpiryTime;
	double                      itsSettlementTime;
	double                      itsPayTime;
	ARM_PricingFunctionEquity*  itsModelEquity;

	/// name for error writting
    static string itsFuncName;
};

//Call operator
struct ARM_GP_CallFctor : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    
	virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );

	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_CallFctor( *this ) ); }

private:
	void GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	string		                itsFXname;
	double                      itsEvalTime;
	double						itsExpiryTime;
	double                      itsSettlementTime;
	ARM_PricingFunctionEquity*  itsModelEquity;
	double                      itsStrikeDouble;
	ARM_VectorPtr               itsStrikeVector;
	int                         itsCallPut;
	double                      itsPayTime;
	ARM_EqFxBase*				itsModelFX;
	ARM_HybridIRFX*				itsHybridModel;
	/// name for error writting
    static string itsFuncName;
};

//CallSpread operator
struct ARM_GP_CallSpreadFctor : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    
	virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );

	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_CallSpreadFctor( *this ) ); }

private:
	void GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	string		                itsFX1Name;
	string		                itsFX2Name;
	double                      itsEvalTime;
	double						itsExpiryTime;
	double                      itsSettlementTime1;
	double                      itsSettlementTime2;
	ARM_HybridIRFX*				itsHybridModel;
	ARM_EqFxBase*				itsModelFX1;
	ARM_EqFxBase*				itsModelFX2;
	ARM_VectorPtr               itsStrikeVector;
	double                      itsAlphaDouble;
	double                      itsBetaDouble;
	double                      itsPayTime;
	
	/// name for error writting
    static string itsFuncName;
};

//EqDigital operator
struct ARM_GP_EqDigitalFctor : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );

    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );

	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_EqDigitalFctor( *this ) ); }

private:
	void GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	string		                itsModelName;
	double                      itsEvalTime;
	double						itsExpiryTime;
	double                      itsSettlementTime;
	ARM_PricingFunctionEquity*  itsModelEquity;
	double                      itsStrikeDouble;
	ARM_VectorPtr               itsStrikeVector;
	double                      itsNotional;
	int                         itsCallPut;
	double                      itsPayTime;
	ARM_DigitType				itsDigitType;
	double						itsEpsilon;

	/// name for error writting
    static string itsFuncName;
};

//CallStrip operator
struct ARM_GP_CallStripFctor : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    
	virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );

	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_CallStripFctor( *this ) ); }

private:
	void GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	string		                itsModelName;
	double                      itsEvalTime;

	ARM_GP_Vector			    itsExpiryTimes;
	ARM_GP_Vector               itsSettlementTimes;
	ARM_GP_Vector               itsPayTimes;

	int                         itsCallPut;

	ARM_VectorPtr               itsStrikes;
	ARM_VectorPtr               itsNominals;

	ARM_PricingFunctionEquity*  itsModelEquity;

	/// name for error writting
    static string itsFuncName;
};

//RangeAccrual operator
struct ARM_GP_RangeAccrualFctor : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    
	virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );

	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_RangeAccrualFctor( *this ) ); }

private:
	void GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	string						itsCurveName;
	int							itsPAYindexType;
	double						itsPAYindexTerm;
	double						itsPAYindexSpread;
	string		                itsFXname;
	int							itsIRindexType;
	double                      itsEvalTime;
	double						itsStartTime;
	double                      itsEndTime;
	ARM_GP_Vector				itsFixingTimes;
	ARM_GP_Vector				itsIRindexResetTimes;
	ARM_GP_Vector				itsIRindexStartTimes;
	ARM_GP_Vector				itsIRindexEndTimes;
	ARM_GP_Vector				itsIRindexTerms;
	ARM_VectorPtr				itsFXDownBarrierVector;
	ARM_VectorPtr               itsFXUpBarrierVector;
	ARM_VectorPtr				itsIRDownBarrierVector;
	ARM_VectorPtr               itsIRUpBarrierVector;
	ARM_VectorPtr               itsNotionalVector;
	double                      itsPAYtime;
	ARM_PricingFunctionEquity*  itsModelEquity;

	/// name for error writting
    static string itsFuncName;
};


CC_END_NAMESPACE()

#endif


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

