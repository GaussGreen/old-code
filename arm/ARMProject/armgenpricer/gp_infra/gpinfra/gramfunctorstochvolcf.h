/*
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *
 *	\file gramfunctorStochVolcf.h
 *   
 *  \brief gram functor file for closed form key words
 *
 *	\author  O.Croissant
 *	\version 1.0
 *	\date June 2005
 */



#ifndef _INGPINFRA_GRAMFUNCTORSTOCHVOLCF_H
#define _INGPINFRA_GRAMFUNCTORSTOCHVOLCF_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gramfunctorbase.h"


CC_BEGIN_NAMESPACE( ARM )

ARM_VectorPtr ARM_GramFctorArg_ConvertToVector( const ARM_GramFctorArg& arg, size_t vector );

//struct ARM_CF_StochVolDispatcher
//{
//	enum DistributionType
//	{
//		K_STOCHASTIC_BLACKSCHOLES_DIST,
//		K_STOCHASTIC_BLACKSCHOLES_GEOMETRIC_DIST,
//		K_STOCHASTIC_BLACKSCHOLES_ARITHMETIC_DIST,
//		K_SABR_DIST,
//		K_SABR_MC_DIST,
//	};
//	
//	enum GreekType
//	{
//		K_Delta,
//		K_Vega,
//		K_Gamma
//	};
//};

struct ARM_CF_StochVolCallFctor : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );

    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );

	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_CF_StochVolCallFctor( *this ) ); }

private:
	void GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( double evalDate, ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	double itsEvalTime;
	double itsMaturityTime;

	/// name for error writting
    static string itsFuncName;
};



struct ARM_CF_StochVolGreekFctor : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );

    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );

	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_CF_StochVolGreekFctor( *this ) ); }

private:
	void GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( double evalDate, ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	double itsEvalTime;
	double itsMaturityTime;

	/// name for error writting
    static string itsFuncName;
};



CC_END_NAMESPACE()

#endif


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

