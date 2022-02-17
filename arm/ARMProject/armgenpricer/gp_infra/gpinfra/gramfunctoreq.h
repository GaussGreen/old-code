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
/// namely for std::vector<double> and ARM_Matrix

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

// Equity Option and Greeks Base Class Functor
struct ARM_GP_EquityOptionAndGreeksFctor : public ARM_GramFctor
{
protected:
	explicit ARM_GP_EquityOptionAndGreeksFctor();

	virtual ~ARM_GP_EquityOptionAndGreeksFctor();

	ARM_GP_EquityOptionAndGreeksFctor(const ARM_GP_EquityOptionAndGreeksFctor& from);

	ARM_GP_EquityOptionAndGreeksFctor& operator=(const ARM_GP_EquityOptionAndGreeksFctor& from);


	virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );


	virtual void checkArgument( ARM_GramFctorArgVector& arg ) = 0;
	virtual void GrabMoreInputs( ARM_GramFctorArgVector& arg ) {;}

	void GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes, size_t statesSize = 1 );
	void SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	ARM_PricingFunctionEquity*  itsModelEquity;
	string		                itsFXname;
	double                      itsEvalTime;
	ARM_VectorPtr               itsExpiryTimeVector;
	ARM_VectorPtr               itsSettlementTimeVector;
	ARM_VectorPtr               itsStrikeVector;
	int                         itsCallPut;
	ARM_VectorPtr               itsPayTimeVector;

};

// Equity Call operator
struct ARM_GP_CallFctor : public ARM_GP_EquityOptionAndGreeksFctor
{
	explicit ARM_GP_CallFctor() : ARM_GP_EquityOptionAndGreeksFctor()
	{

	}

	virtual ~ARM_GP_CallFctor()
	{

	}

	ARM_GP_CallFctor(const ARM_GP_CallFctor& from)
		: ARM_GP_EquityOptionAndGreeksFctor(from)
	{

	}

	ARM_GP_CallFctor& operator=(const ARM_GP_CallFctor& from)
	{
		if ( &from != this) {
			ARM_GP_EquityOptionAndGreeksFctor::operator =(from);
		}
		return *this;
	}

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );

	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_CallFctor( *this ) ); }

private:
	void checkArgument( ARM_GramFctorArgVector& arg);

	/// name for error writting
    static string itsFuncName;
};


// Equity Call greeks
struct ARM_GP_GreekFctor : public ARM_GP_EquityOptionAndGreeksFctor
{
	explicit ARM_GP_GreekFctor() : ARM_GP_EquityOptionAndGreeksFctor(), greekType_("DELTA")
	{

	}

	virtual ~ARM_GP_GreekFctor()
	{

	}

	ARM_GP_GreekFctor(const ARM_GP_GreekFctor& from)
		: ARM_GP_EquityOptionAndGreeksFctor(from), greekType_(from.greekType_)
	{

	}

	ARM_GP_GreekFctor& operator=(const ARM_GP_GreekFctor& from)
	{
		if ( &from != this) {
			ARM_GP_EquityOptionAndGreeksFctor::operator =(from);
			greekType_ = from.greekType_;
		}
		return *this;
	}

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );

	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_GreekFctor( *this ) ); }

private:
	void checkArgument( ARM_GramFctorArgVector& arg );
	void GrabMoreInputs( ARM_GramFctorArgVector& arg );

	string greekType_;
	//int greekType_;

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

	std::vector<double>			    itsExpiryTimes;
	std::vector<double>               itsSettlementTimes;
	std::vector<double>               itsPayTimes;

	int                         itsCallPut;

	ARM_GP_VectorPtr               itsStrikes;
	ARM_GP_VectorPtr               itsNominals;

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
	std::vector<double>				itsFixingTimes;
	std::vector<double>				itsIRindexResetTimes;
	std::vector<double>				itsIRindexStartTimes;
	std::vector<double>				itsIRindexEndTimes;
	std::vector<double>				itsIRindexTerms;
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

