/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pricingstatesprocess.h
 *
 *  \brief simple function to process gramnode
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#ifndef _INGPINFRA_PRICINGSTATESPROCESS_H
#define _INGPINFRA_PRICINGSTATESPROCESS_H

#include "gpbase/port.h"
#include "gramnode.h"
#include "pricingstates.h"
#include "pricingmodel.h"
#include "modelnamemap.h"

CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////////////////
///	Struct : PushBackAndSetFtor 
///	Action : Functor whose responsability is to take 
///				a statesPayoff push_back a double payoff
////			and set back to the currentState
////////////////////////////////////////////////////

struct PushBackAndSetFtor
{
	void ProcessStates( ARM_PricingStatesPtr& currentStates, size_t colIdx, ARM_VectorPtr& payoffs, ARM_PricingModel* model, double eventTime, bool isAZero )
	{
        size_t lastPayoff = currentStates->GetPayoffsSize();
		currentStates->resizePayoffs(lastPayoff+1);
        for(int i=0;i<currentStates->size();++i)
            currentStates->SetPayoff(i,lastPayoff,(*payoffs)[i]);
	}
};


////////////////////////////////////////////////////
///	Struct : AddToState0
///	Action : Functor whose responsability is to add to
///				state 0 the payoff 
////////////////////////////////////////////////////

struct AddToState0
{
    AddToState0(const string& payCurveName="") : itsPayCurveName(payCurveName) {}

	void ProcessStates( ARM_PricingStatesPtr& currentStates, size_t colIdx, ARM_VectorPtr& payoffs, ARM_PricingModel* model, double eventTime, bool isAZero )
	{
		bool isOtherPayoffs = currentStates->GetPayoffStates(colIdx).GetOtherPayoffsFlag();

		/// Save a snapshot of current payoff states before any numeraire treatment
		//if(isOtherPayoffs)
			//currentStates->GetPayoffStates(colIdx).push_backPayoffSnapshot( payoffs );

		// To prevent in place effect in the graph nodes
		ARM_VectorPtr newPayoffs( new std::vector<double>( *payoffs ) );

		if( !isAZero )
			model->ProcessPaidPayoffs(itsPayCurveName, newPayoffs, eventTime, currentStates );

		int intermediatePayoffPosition = currentStates->GetPayoffStates(colIdx).GetIntermediatePayoffsSize();
		if(isOtherPayoffs)
			currentStates->GetPayoffStates(colIdx).resizeIntermediatePayoffs(intermediatePayoffPosition+1);

		for(size_t i=0;i<currentStates->size();++i)
		{
			currentStates->GetPayoffStates(colIdx).SetPayoff(i,currentStates->GetPayoffStates(colIdx).GetPayoff(i) + (*newPayoffs)[i]); 
			if(isOtherPayoffs)
				currentStates->GetPayoffStates(colIdx).SetIntermediatePayoff(i,intermediatePayoffPosition,(*newPayoffs)[i]);
		}
	}

private:
    string itsPayCurveName;
};


////////////////////////////////////////////////////
///	Struct : CreateFirstPayoff
///	Action : Functor whose responsability is to create the
///				first payoff
////////////////////////////////////////////////////
struct CreateFirstPayoff
{
    CreateFirstPayoff(const string& payCurveName="") : itsPayCurveName(payCurveName) {}

	void ProcessStates( ARM_PricingStatesPtr& currentStates, size_t colIdx, ARM_VectorPtr& payoffs, ARM_PricingModel* model, double eventTime, bool isAZero )
	{
        bool isOtherPayoffs = currentStates->GetPayoffStates(colIdx).GetOtherPayoffsFlag();

        /// Save a snapshot of current payoff states before any numeraire treatment
        //if(isOtherPayoffs)
		  //  currentStates->GetPayoffStates(colIdx).push_backPayoffSnapshot(ARM_VectorPtr(new std::vector<double>(*payoffs)));

		// To prevent in place effect in the graph nodes
		ARM_VectorPtr newPayoffs( new std::vector<double>( *payoffs ) );

		if( !isAZero )
			model->ProcessPaidPayoffs(itsPayCurveName, newPayoffs, eventTime, currentStates );

        if(isOtherPayoffs)
		    currentStates->GetPayoffStates(colIdx).resizeIntermediatePayoffs(1,currentStates->size());

        for(size_t i=0;i<currentStates->size();++i)
		{
			currentStates->GetPayoffStates(colIdx).SetPayoff(i,(*newPayoffs)[i]); 
            if(isOtherPayoffs)
			    currentStates->GetPayoffStates(colIdx).SetIntermediatePayoff(i,0,(*newPayoffs)[i]);
		}
    }

private:
    string itsPayCurveName;
};


////////////////////////////////////////////////////
///	Struct : CreateFirstPayoffClosedForms
///	Action : Functor whose responsability is to create the
///				first payoff for closed forms
////////////////////////////////////////////////////
struct CreateFirstPayoffClosedForms
{
	void ProcessStates( ARM_PricingStatesPtr& currentStates, size_t colIdx, ARM_VectorPtr& payoffs, ARM_PricingModel* model, double eventTime, bool isAZero )
	{
        for(size_t i=0;i<currentStates->size();++i)
			currentStates->GetPayoffStates(colIdx).SetPayoff(i,(*payoffs)[i]);
    }
};


////////////////////////////////////////////////////
///	Routine: EvalChildNodeAndProcess
///	Returns: ARM_PricingStatesPtr 
///	Action : evaluate a child node and set it to the pricing states
///				( helper function for evaluation of child node )
////////////////////////////////////////////////////
template <typename ActionPolicy >
	void EvalChildNodeAndProcess( const ARM_ExpNodePtr& ChildNode, ARM_PricingModel* model, 
		ARM_PricingStatesPtr& states, size_t colIdx, ActionPolicy& Ftor, double fromTime, bool auxiliaryEval = false )
{
	// For the closed form num method;
	ARM_NumMethodPtr numMethod = model->GetNumMethod();

	bool fwdLoop = false;
		
	if (!numMethod.IsNull())
		fwdLoop = (numMethod->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING) && (numMethod->GetPricingDirCurrLoop() == ARM_NumMethod::GP_FWDLOOKING);

	ARM_GramFctorArg payoff;
	if (auxiliaryEval || fwdLoop)
	{
		ARM_ExpNodeRef* refNode = dynamic_cast<ARM_ExpNodeRef*>(&*ChildNode);
		if(refNode)
			payoff = refNode->EvalAux(model, states, true );
		else
			 throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": An auxiliary node should be a Ref Node. " );
	}
	else
	{
		payoff = ChildNode->Eval(model, states );
	}
	bool isAZero = false;    
    ARM_VectorPtr payoffVec;
    
    /// builds the payoff according to the return type
    switch( payoff.GetType() )
    {
    case GFAT_DOUBLE_TYPE:
		{
			if( fabs( payoff.GetDouble()  )< K_NEW_DOUBLE_TOL ) 
				isAZero = true;
			payoffVec = ARM_VectorPtr(new std::vector<double>(states->size(),payoff.GetDouble()));
		}
		break;
		
    case GFAT_VECTOR_TYPE:
        payoffVec = payoff.GetVector();
        break;
		
    default:
        throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "Trying to refer to something that is neither a double or a vector " );
    }

#ifdef __GP_STRICT_VALIDATION
    if( !states.IsNull() && (payoffVec->size() != states->size()) )
    {
        CC_Ostringstream os;
        os << "ModelStates size = "	<< states->size() 
            << " while payoff size = " << payoff.GetVector()->size();
        throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }
#endif

	double evalDate = fromTime;
	if (numMethod != ARM_NumMethodPtr(NULL))
		evalDate = numMethod->ConvertEvalDate(evalDate,0);

    Ftor.ProcessStates(states,colIdx,payoffVec,model,evalDate,isAZero);
};



CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

