/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: gramnode.cpp,v $
 * Revision 1.1  2003/13/08 16:41:37  ebenhamou
 * Initial revision
 *
 *
 */


/*! \file gramnode.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#include "gpinfra/gramnode.h"
#include "gpinfra/gramfunctor.h"
#include "gpinfra/gramfunctorsimple.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/functordef.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/exerciseboundary.h"

#include "gpbase/ostringstream.h"
#include "gpbase/env.h"

#include <glob/expt.h>
#include <functional>
#include <algorithm>


CC_USING_NS( std, string )

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_ExpNode 
///	Static Member Function: ModelDoesNotSupportInPlace
////////////////////////////////////////////////////
bool ARM_ExpNode::ModelDoesNotSupportInPlace( ARM_PricingModel* mod )
{
#if defined( __GP_STRICT_VALIDATION)
	if( !mod )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": model* is NULL!" ); 
#endif
	return 	mod->GetNumMethod() != ARM_NumMethodPtr(NULL)
		 && mod->GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING;
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNode 
///	Static Member: itsSeparator
////////////////////////////////////////////////////
const string ARM_ExpNode::itsSeparator = string("\t");


////////////////////////////////////////////////////
///	Class  : ARM_ExpNode 
///	Routine: ~ARM_ExpNode
///	Returns: nothing
///	Action : Destructor
////////////////////////////////////////////////////

ARM_ExpNode::~ARM_ExpNode()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeRef
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////

string ARM_ExpNodeRef::toString( const string& indent ) const 
{
	CC_Ostringstream os;
    os << indent << "<Reference-node ChildLocationString='" << itsChildLocationString;
	/// print if it is a PV node!
    if(itsIsAPVNode)
    	os << "' PV-Node, PayModelName='" << itsPayModelName << "' >\n";
    else
        os << "' >\n";

	/// if it is a shared node, display it with the number of references!
	if(itsCountor>1)
	{
		os  << indent << "=========> Shared Node (Nb=" << itsCountor << ")\n";

#if defined(__GP_SHOW_SHARED_NODE_COORDINATES)
		os	<< indent << "                       (Parents=";
		
		for( size_t i=0; i<itsParentCoords.size(); ++i )
		{
			os  << " [" << itsParentCoords[i].itsElem1 
				<< ","  << itsParentCoords[i].itsElem2
                << "]";
		}
		os << " )\n";
#endif

		os	<< itsChildNode->toString( indent + itsSeparator );
		os  << indent << "=========> End Shared Node\n";

	}
	else
		os	<< itsChildNode->toString( indent + itsSeparator ) ;

	os << indent << "</Reference-node>\n";
	
	return os.str();
}

////////////////////////////////////////////////////
///	Class  :ARM_ExpNodePtr
///	Routine: operator<<
///	Returns: void
///	Action : Display the contents of an ExpNode
////////////////////////////////////////////////////
ostream& operator<< (ostream& os, const ARM_ExpNode& expNode )
{ 
	return os; 
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeRef
///	Routine: Reset
///	Returns: void
///	Action : reset the object
////////////////////////////////////////////////////

bool ARM_ExpNodeRef::Reset(ResetLevel partialReset)
{ 
	bool retReset = true;

	// The IsReset flag prevent a "crazy" Reset which will take for ever
	// We reset the node in three cases:
	// _ ResetBucket
	// _ ResetPartial
	// _ ResetCountor
	// _ The child node has been reseted
	if ((!itsIsReseted && partialReset == ResetLoop) || (itsIsPreComputed && partialReset == ResetTotal) || (itsIsPreComputed && partialReset == ResetBucket) || ((itsEvalCountor!=-1) && partialReset == ResetCountor))
	{
		bool resetParent = itsChildNode->Reset(partialReset);
		if ( partialReset == ResetBucket || partialReset == ResetTotal || partialReset == ResetCountor  || resetParent)
		{
			itsIsPreComputed = false;
		}
		else
		{
			retReset = false;
		}

		itsIsReseted = true;
	}
	else
	{
		retReset = false;
	}

	if (partialReset == ResetCountor)
	{
		itsEvalCountor = -1;
	}
	else
	{
		itsEvalCountor = 0;
	}

	return retReset;
}



////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeRef
///	Routine: SetValue
///	Returns: void
///	Action : set the value taking into account if it is a PV node!
////////////////////////////////////////////////////
void ARM_ExpNodeRef::SetValue( const ARM_GramFctorArg& val, ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	itsValue		= val;
	itsIsPreComputed= true;
	itsIsReseted= false;
}



////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeRef
///	Routine: EvalAux
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type cell ref and 
/// handle the auxiliary node counting
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_ExpNodeRef::EvalAux( ARM_PricingModel* model, const ARM_PricingStatesPtr& states, bool auxiliaryEval )
{
	if ( (itsEvalCountor >= 0) && (itsEvalCountor >= itsCountor))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
						ARM_USERNAME + ": This refnode cannot be call more !" ); 
	}


	/// Check if this has been already precomputed!
	if(!itsIsPreComputed)
	{
		ARM_GramFctorArg val = itsChildNode->Eval(model,states);

		if( itsIsAPVNode )
		{
			double childTime = model->GetTimeFromDate( itsChildDate );
			ARM_VectorPtr payoff;

			switch( val.GetType() )
			{
			case GFAT_VECTOR_TYPE:
				{
					payoff = val.GetVector();
					break;
				}
			case GFAT_DOUBLE_TYPE:
				{
					payoff = ARM_VectorPtr( new std::vector<double>( states->size(), val.GetDouble() ) );
					break;
				}
			/// case GFAT_DATE_TYPE:
			/// case GFAT_STRING_TYPE:
			default:
				{
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
						ARM_USERNAME + ": Unsupported type... supported are double and vector" ); 
				}
			}
			model->ProcessPaidPayoffs(itsPayModelName, payoff, childTime, states );
			itsValue = ARM_GramFctorArg( payoff );
		}
		else
			itsValue=val;

		/// set that it is precomputed!
		itsIsPreComputed=true; 
	}

	itsIsReseted=false;

	if ((itsEvalCountor>=0) && !auxiliaryEval)
	{
		// Increase the evaluation countor for each evaluation of the ref node
		itsEvalCountor++;
	}

	ARM_GramFctorArg value;

	/// to avoid side-effect 
	/// we clone the value
	// if shared
	if(itsCountor>1)
		value = itsValue.Duplicate();
	else
		value = itsValue;

	// The value is destroyed at the last evaluation so the memory is released
	if (itsEvalCountor == itsCountor)
		itsValue.SetVector(ARM_GP_VectorPtr(NULL));
	
	return value;
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeRef
///	Routine: EvalAux
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type cell ref
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_ExpNodeRef::Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states)
{
	return EvalAux(model,states);
}

////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeRef
///	Routine: Eval
///	Returns: ARM_GramFctorArg
///	Action : Changes itsChildNode's Exercises to Trigger
////////////////////////////////////////////////////

void ARM_ExpNodeRef::ChangeExerciseIntoTrigger( ARM_DealDescription& dealDesc )
{
	if (!itsHasBeenChangeIntoTrigger)
	{
		itsChildNode->ChangeExerciseIntoTrigger( dealDesc );
		itsHasBeenChangeIntoTrigger = true;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeDateOrDouble
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////

string ARM_ExpNodeDateOrDouble::toString(const string& indent) const
{
	CC_Ostringstream os;
	os << indent << "<Date-or-double-node Double='" <<itsDouble << "' />\n";
	return os.str();
}



////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeDateOrDouble
///	Routine: Eval
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type Date Or Double
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_ExpNodeDateOrDouble::Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	return ARM_GramFctorArg(itsDouble);	
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeDoubleCst
///	Routine: Constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_ExpNodeDoubleCst::ARM_ExpNodeDoubleCst( const string& name, ARM_CstManagerPtr cstManager )
: itsCountor(1)
// FIXMEFRED: mig.vc8 (31/05/2007 10:42:04):iterator(null)
{
/// normally this is checked in the generic security, hence no use!
#ifdef __GP_STRICT_VALIDATION
	if( !cstManager->DoesCstNameExist( name ) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME +  ": " + name + "Cst does not exist in the cst manager!" ); 
#endif
	itsCstIterator = cstManager->find( name );
}

	

////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeDoubleCst
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////

string ARM_ExpNodeDoubleCst::toString(const string& indent) const
{
	CC_Ostringstream os;
	os << indent << "<Cst Manager linked double name='" << (*itsCstIterator).first
		<< "' value='" << (*itsCstIterator).second << "' countor='"
		<< itsCountor << "' />\n";
	return os.str();
}



////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeDoubleCst
///	Routine: Eval
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type a constant mananger linked double
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_ExpNodeDoubleCst::Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	return (*itsCstIterator->second);	
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeDate
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////

string ARM_ExpNodeDate::toString(const string& indent) const
{
	CC_Ostringstream os;
	os << indent << "<Date-node Date='"  << itsDate.toString() << "' />\n";
	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeDate
///	Routine: Eval
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type Date
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_ExpNodeDate::Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	return ARM_GramFctorArg(itsDate);	
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeFactor
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////

string ARM_ExpNodeFactor::toString(const string& indent) const
{
	CC_Ostringstream os;
	os << indent << "<Factor-node>\n"
		<< indent << itsSeparator << "(" << itspNode->toString( indent + itsSeparator ) << indent << ")\n"
		<< indent << "</Factor-node>\n";
	return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeFactor
///	Routine: ChangeExerciseIntoTrigger
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////

void ARM_ExpNodeFactor::ChangeExerciseIntoTrigger( ARM_DealDescription& dealDesc )
{
	itspNode->ChangeExerciseIntoTrigger( dealDesc );
}

////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeFactor
///	Routine: Eval
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type function
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_ExpNodeFactor::Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	return itspNode->Eval(model,states);
}



////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeString
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////

string ARM_ExpNodeString::toString(const string& indent) const
{
	CC_Ostringstream os;
	os << indent << "<String-node string='" << itsString << "' />\n";
	return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeString
///	Routine: Eval
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type string
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_ExpNodeString::Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	return ARM_GramFctorArg(itsString);	
}



////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeFunc
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////

string ARM_ExpNodeFunc::toString(const string& indent) const
{
	CC_Ostringstream os;

	os << indent << "<Function-node Function-Name='" << itsFunction.FuncName() << "' "
		<< "Eval-Date='" << ARM_Date( itsEvalDate ).toString();
    os << "', PayModelName='" << itsFunction.Fctor()->GetPayModelName() << "' >\n";

	size_t argSize = itsArgs.size();
	size_t i;
	for(i=0; i<argSize; ++i ) 
		os << itsArgs[i]->toString( indent + itsSeparator  );
	os << indent << "</Function-node>\n";

	return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeFunc
///	Routine: Eval
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type function
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_ExpNodeFunc::Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	/// we are forced to recomputed arguments 
	/// as when determining used time lags
	/// we use a dummmed pricing states!
	size_t argSize = itsArgs.size();
	size_t i;

	ARM_GramFctorArgVector args;
	args.reserve(argSize);
	for(i=0; i<argSize; ++i)
	{
		ARM_GramFctorArg arg = itsArgs[i]->Eval( model, states );

		/// should not be necessary 
		/// but is here to catch bug when populating additional function in the grammar!
#ifdef __GP_STRICT_VALIDATION
		GPAF_CheckArgReturnType( arg, itsFunction.Args()[i].Type(), itsArgs[i], itsFunction.FuncName() );
#endif
		args.push_back( arg );
	}

	// For the closed form num method;
	ARM_NumMethodPtr numMethod = model->GetNumMethod();

	double evalDate = itsEvalDate;

	if (numMethod != ARM_NumMethodPtr(NULL))
		evalDate = numMethod->ConvertEvalDate(evalDate,model->GetAsOfDate().GetJulian());


	return (*itsFunction.Fctor())( args, model, evalDate, states, itsArgs );
}



////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeFunc
///	Routine: GetUsedTimeLags
///	Returns: ARM_GramFctorArg
///	Action : finds the corresponding payTime of a function
///             given its argument, model and evaldate!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_ExpNodeFunc::GetUsedTimeLags( ARM_PricingModel* model)
{
	size_t argSize = itsArgs.size();
	size_t i;


	ARM_GramFctorArgVector args;
	args.reserve( argSize );
    ARM_PricingStatesPtr dumStates(new ARM_PricingStates());

	for(i=0; i<argSize; ++i )
	{
        /// should not need states, hence dummy states!
		ARM_GramFctorArg arg = itsArgs[i]->Eval( model, dumStates);
		
		/// should not be necessary 
		/// but is here to catch bug when populating additional function in the grammar!
#ifdef __GP_STRICT_VALIDATION
		GPAF_CheckArgReturnType( arg, itsFunction.Args()[i].Type(), itsArgs[i], itsFunction.FuncName() );
#endif
		args.push_back( arg );
	}

	// For the closed form num method;
	ARM_NumMethodPtr numMethod = model->GetNumMethod();

	double evalDate = itsEvalDate;
	double asOfDate = model->GetAsOfDate().GetJulian();
	if (numMethod != ARM_NumMethodPtr(NULL))
		evalDate = numMethod->ConvertEvalDate(evalDate,asOfDate);

	/// evaluate with the arg the used time lags
    return itsFunction.Fctor()->GetUsedTimeLags( args, model, evalDate, itsArgs );
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeFunc
///	Routine: Initialise
///	Returns: void
///	Action : Initialise all the argument nodes
////////////////////////////////////////////////////
bool ARM_ExpNodeFunc::Reset(ResetLevel partialReset)
{
	bool resetParent = false;
	bool result = false;
	/// reset nodes
	size_t i, argSize = itsArgs.size();

	for(i=0; i<argSize; ++i )
		resetParent |= itsArgs[i]->Reset(partialReset);

	if ( partialReset == ResetBucket || partialReset == ResetTotal || partialReset == ResetCountor || resetParent)
	{
		/// reset the gramFunction as well (contains pre-computed arguments!)
		if ( partialReset == ResetTotal )
			itsFunction.Reset( true );
		else 
			itsFunction.Reset( false );

		result = true;
	}

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeFunc
///	Routine: Eval
///	Returns: ARM_GramFctorArg
///	Action : Changes itsChildNode's Exercises to Trigger
////////////////////////////////////////////////////

void ARM_ExpNodeFunc::ChangeExerciseIntoTrigger( ARM_DealDescription& dealDesc )
{
	for( CC_STL_VECTOR( ARM_ExpNodePtr )::iterator iter = itsArgs.begin() ; 
			iter != itsArgs.end() ; iter++ )
		(*iter)->ChangeExerciseIntoTrigger( dealDesc );
}

////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeFunc
///	Routine: constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_ExpNodeFunc::ARM_ExpNodeFunc( const ARM_GramFunction& Function, const CC_STL_VECTOR( ARM_ExpNodePtr )& Args, double evalDate) 
	: itsFunction(Function), itsArgs(Args ), itsEvalDate(evalDate)
{
	CC_STL_VECTOR( ARM_ExpNodePtr )::const_iterator iter;

	for( iter = Args.begin() ; iter!=Args.end() ; iter++ )
	{
		(*iter)->AddParent( this );
		(*iter)->AddParentPointer( const_cast<ARM_ExpNodePtr*>(&(*iter)) );
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeExerciseFunc
///	Routine: Eval
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type function
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_ExpNodeExerciseFunc::Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	size_t argSize = itsArgs.size();
	ARM_GramFctorArgVector args;
	args.reserve( argSize );

	ARM_NumMethod::GP_PricingDirection pricingDir = model->GetNumMethod()->GetPricingDirCurrLoop();

	if (argSize < itsFunction.Args().size()-1 ) // because there is one "va arg" argument
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME +  ": " + "Invalid number of arguments in the exercise node!" );
	}

	ARM_GramFctorArg arg;
	if (pricingDir == ARM_NumMethod::GP_BCKWDLOOKING)
	{
		// Payoff 1
		arg = itsArgs[0]->Eval( model, states );
		args.push_back(arg);
		// Payoff 2
		arg = itsArgs[1]->Eval( model, states );
		args.push_back(arg);
	}
	else
	{
		// Payoff 1
		ARM_ExpNodeRef* nodeRef = dynamic_cast<ARM_ExpNodeRef*>(&*itsArgs[0]);
		if (nodeRef)
			arg = nodeRef->EvalAux(model,states,true);
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME +  ": " + "the payoff 1 argument of the exercise node is a not a referennce !" );
		args.push_back(arg);
		// Payoff 2
		nodeRef = dynamic_cast<ARM_ExpNodeRef*>(&*itsArgs[1]);
		if (nodeRef)
			arg = nodeRef->EvalAux(model,states,true);
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME +  ": " + "the payoff 2 argument of the exercise node is a not a referennce !" );
		args.push_back(arg);

	}
#ifdef __GP_STRICT_VALIDATION
	GPAF_CheckArgReturnType( arg, itsFunction.Args()[0].Type(), itsArgs[0], itsFunction.FuncName() );
	GPAF_CheckArgReturnType( arg, itsFunction.Args()[1].Type(), itsArgs[1], itsFunction.FuncName() );
#endif
	
	// Continuation value
	if (pricingDir == ARM_NumMethod::GP_BCKWDLOOKING)
	{
		arg = itsArgs[2]->Eval( model, states );
#ifdef __GP_STRICT_VALIDATION
		GPAF_CheckArgReturnType( arg, itsFunction.Args()[1].Type(), itsArgs[2], itsFunction.FuncName() );
#endif
		args.push_back(arg);
	}
	else
	{
		ARM_GramFctorArg emptyArg;
		args.push_back(emptyArg);
	}

	// Regression variables
	if ( model->GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING)
	{
		for ( size_t i=3 ; i<argSize ; i++ )
		{
			if (pricingDir == ARM_NumMethod::GP_BCKWDLOOKING)
			{
				arg = itsArgs[i]->Eval( model,states );
			}
			else
			{
				ARM_ExpNodeRef* nodeRef = dynamic_cast<ARM_ExpNodeRef*>(&*itsArgs[i]);
				if (nodeRef)
					arg = nodeRef->EvalAux(model,states,true);
				else
				{
					CC_Ostringstream os;

					os << ARM_USERNAME <<  ": " << "the " << i << "th argument of the exercise node is a not a referennce !";

					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str());
				}
			}
#ifdef __GP_STRICT_VALIDATION
			GPAF_CheckArgReturnType( arg, itsFunction.Args()[3].Type(), itsArgs[i], itsFunction.FuncName() );
#endif
			args.push_back(arg);
		}
	}

    return (*itsFunction.Fctor())( args, model, itsEvalDate, states, itsArgs );
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeExerciseFunc
///	Routine: Reset
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type function
////////////////////////////////////////////////////
bool ARM_ExpNodeExerciseFunc::Reset(ResetLevel partialReset)
{
	if ( partialReset == ResetLoop )
	{
		size_t argSize = itsArgs.size();
		if (argSize < itsFunction.Args().size()-1 )
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME +  ": " + "Invalid number of arguments in the exercise node!" );
		}

		// Reset of te continuation option argument.
		itsArgs[2]->Reset(partialReset);

		return true;
	}
	else
	{
		return ARM_ExpNodeFunc::Reset(partialReset);
	}
}

vector<string> ExtractArgs(string str)
{
	size_t leftparpos = str.find("(");
	size_t rightparpos = str.find(")");
	size_t curpos = leftparpos+1, lastcurpos = curpos;

	vector<string> args;

	while ((curpos = str.find(",",curpos)) != string::npos)
	{
		args.push_back(str.substr(lastcurpos,curpos-lastcurpos));
		lastcurpos = curpos+1;
		curpos = lastcurpos;
	}
	args.push_back(str.substr(lastcurpos,rightparpos-lastcurpos));

	return args;
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeExerciseFunc
///	Routine: ChangeExerciseIntoTrigger
///	Returns: void
///	Action : Changes itsChildNode's Exercises to Trigger
////////////////////////////////////////////////////

void ARM_ExpNodeExerciseFunc::ChangeExerciseIntoTrigger( ARM_DealDescription& dealDesc )
{
	// Gets Exercise Boundary
	ARM_GramFctorPtr fctor( itsFunction.Fctor() );
	ARM_GP_Exercise* exercfctor = dynamic_cast<ARM_GP_Exercise*>( &(*fctor) );
	ARM_ExerciseBoundary * exercbound = exercfctor->GetExerciseBoundary();

	if( ! exercbound ) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"exercise boundary missing. Did you run an amc on this security ? is resetflag=N ?" );

	ARM_GP_VectorPtr vec( exercbound->GetExerciseBoundary() );

	if( vec == ARM_GP_VectorPtr(NULL) )
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"exercise boundary params is null" );
	}

	// First we extract all the arguments of the Exercise keyword
	string str( dealDesc.GetElem( getRowInDealDesc(), getColInDealDesc() ) );
	vector<string> args = ExtractArgs(str);

    CC_Ostringstream os;

	ARM_ExpNodeFunc::ChangeExerciseIntoTrigger( dealDesc );

	os << "Trigger(";

	if(typeid(*exercbound) == typeid(ARM_AndersenExerciseBoundary ) )
	{
		if (vec->size() == 1)
		{
			double trigBarrier = (*vec)[0];
			string trigParam;
			// The trigger value has not been inputed by the user
			if (args.size() == 3)
			{
				trigParam = args[1];
			}
			// The trigger value has been inputed by the user
			else
			{
				// We use the 4th parameter as a trigger value
				trigParam = args[3];
			}

			bool  isGreaterThanValue = ((static_cast<ARM_AndersenExerciseBoundary*> (exercbound))->ExerciseIfPayoffIsGreaterThanValue());

			if (isGreaterThanValue)
			{
				os << trigParam << ",";
				os << std::fixed << (*vec)[0];
			}
			else
			{
				os << std::fixed << (*vec)[0] << ",";
				os <<  trigParam;
			}
			os << "," << args[0] << "," << args[1] << "," << args[2] << ")";

		}
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "the size of the andersen exercise boundary should be 1." );
		
	}
	else if(typeid(*exercbound) == typeid(ARM_LSExerciseBoundary ) )
	{
		os << args[1] << ",";

		// By default it's a 4 dimensions polynom of the payoff
		if (args.size() == 3)
		{
			os << (*vec)[0] << "+";
			os << (*vec)[1] << "*(" << args[1] << ")+";
			os << (*vec)[2] << "*POW(" << args[1] << ",2)+";
			os << (*vec)[3] << "*POW(" << args[1] << ",3)+";
			os << (*vec)[4] << "*POW(" << args[1] << ",4),";
		}
		// We use the parameter from the fourth
		else if ((args.size()-3) == vec->size())
		{
			for (size_t i = 0; i < vec->size(); ++i)
			{
				os << "(" << (*vec)[i] << ")" << "*" << args[3+i];

				if (i < vec->size()-1)
					os << "+";
			}
		}
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "the size of the longstaff exercise boundary is not compatible with the number of parameters of the keyword." );

		os << "," << args[0] << "," << args[1] << ","<< args[2] << ")";
	}
				

	//// updates deal description
	dealDesc.SetElem( getRowInDealDesc(), getColInDealDesc(), os.str(), 
		dealDesc.GetElemFormat( getRowInDealDesc(), getColInDealDesc() ) );
}



////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeLinearPVFunc
///	Routine: Eval
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type function
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_ExpNodeLinearPVFunc::Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	if ( model->GetNumMethod()->GetPricingDirection() != ARM_NumMethod::GP_FWDBCKWDLOOKING)
		return ARM_ExpNodeFunc::Eval( model, states );

	size_t argSize = itsArgs.size();
	ARM_GramFctorArgVector args;
	args.reserve( argSize );

	ARM_NumMethod::GP_PricingDirection pricingDir = model->GetNumMethod()->GetPricingDirCurrLoop();
	if (pricingDir != ARM_NumMethod::GP_BCKWDLOOKING )
		return ARM_GramFctorArg(ARM_VectorPtr(new std::vector<double>(states->size(),0.0)));

	ARM_GramFctorArg arg = itsArgs[0]->Eval( model, states );
#ifdef __GP_STRICT_VALIDATION
	GPAF_CheckArgReturnType( arg, itsFunction.Args()[0].Type(), itsArgs[0], itsFunction.FuncName() );
#endif
	args.push_back(arg);
	return (*itsFunction.Fctor())( args, model, itsEvalDate, states, itsArgs );
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeLinearPVFunc
///	Routine: Reset
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type function
////////////////////////////////////////////////////
bool ARM_ExpNodeLinearPVFunc::Reset(ResetLevel partialReset)
{
	if ( partialReset != ResetLoop )
		return ARM_ExpNodeFunc::Reset(partialReset);

	itsArgs[0]->Reset(partialReset);
	return true;
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeTriggerFunc
///	Routine: Eval
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type function
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_ExpNodeTriggerFunc::Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	size_t argSize = itsArgs.size();
	ARM_GramFctorArgVector args;
	args.reserve( argSize );
	ARM_GramFctorArg arg;

	ARM_NumMethod::GP_PricingDirection pricingDir = model->GetNumMethod()->GetPricingDirCurrLoop();

	if (argSize != itsFunction.Args().size() ) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME +  ": " + "Invalid number of arguments in the Trigger node!" );

	if( pricingDir == ARM_NumMethod::GP_BCKWDLOOKING )
	{
		// Var Test 1
		arg = itsArgs[0]->Eval( model, states );
#ifdef __GP_STRICT_VALIDATION
		GPAF_CheckArgReturnType( arg, itsFunction.Args()[0].Type(), itsArgs[0], itsFunction.FuncName() );
#endif
		args.push_back(arg);

		// Var Test 2
		arg = itsArgs[1]->Eval( model, states );
#ifdef __GP_STRICT_VALIDATION
		GPAF_CheckArgReturnType( arg, itsFunction.Args()[1].Type(), itsArgs[1], itsFunction.FuncName() );
#endif
		args.push_back(arg);

		// Payoff 1
		arg = itsArgs[2]->Eval( model, states );
#ifdef __GP_STRICT_VALIDATION
		GPAF_CheckArgReturnType( arg, itsFunction.Args()[2].Type(), itsArgs[2], itsFunction.FuncName() );
#endif
		args.push_back(arg);

		// Payoff 2
		arg = itsArgs[3]->Eval( model, states );
#ifdef __GP_STRICT_VALIDATION
		GPAF_CheckArgReturnType( arg, itsFunction.Args()[3].Type(), itsArgs[3], itsFunction.FuncName() );
#endif
		args.push_back(arg);
	}
	else
	{
		// Var Test 1
		ARM_ExpNodeRef* nodeRef = dynamic_cast<ARM_ExpNodeRef*>(&*itsArgs[0]);
		if (nodeRef)
			arg = nodeRef->EvalAux(model,states,true);
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME +  ": " + "the var test 1 argument of the exercise node is a not a referennce !" );
		args.push_back(arg);
		// Var Test 2
		nodeRef = dynamic_cast<ARM_ExpNodeRef*>(&*itsArgs[1]);
		if (nodeRef)
			arg = nodeRef->EvalAux(model,states,true);
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME +  ": " + "the var test 2 argument of the exercise node is a not a referennce !" );
		args.push_back(arg);
		// Payoff 1
		nodeRef = dynamic_cast<ARM_ExpNodeRef*>(&*itsArgs[2]);
		if (nodeRef)
			arg = nodeRef->EvalAux(model,states,true);
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME +  ": " + "the payoff 1 argument of the exercise node is a not a referennce !" );
		args.push_back(arg);
		// Payoff 2
		nodeRef = dynamic_cast<ARM_ExpNodeRef*>(&*itsArgs[3]);
		if (nodeRef)
			arg = nodeRef->EvalAux(model,states,true);
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME +  ": " + "the payoff 2 argument of the exercise node is a not a referennce !" );
		args.push_back(arg);
	}


	// Future Payoff
	if (pricingDir == ARM_NumMethod::GP_BCKWDLOOKING)
	{
		arg = itsArgs[4]->Eval( model, states );
#ifdef __GP_STRICT_VALIDATION
		GPAF_CheckArgReturnType( arg, itsFunction.Args()[4].Type(), itsArgs[4], itsFunction.FuncName() );
#endif
		args.push_back(arg);
	}

    return (*itsFunction.Fctor())( args, model, itsEvalDate, states, itsArgs );
}

////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeTriggerFunc
///	Routine: Reset
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type function
////////////////////////////////////////////////////
bool ARM_ExpNodeTriggerFunc::Reset(ResetLevel partialReset)
{
	if ( partialReset == ResetLoop )
	{
		size_t argSize = itsArgs.size();
		if (argSize != itsFunction.Args().size() )
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME +  ": " + "Invalid number of arguments in the Trigger node!" );
		}

		// Reset of the continuation option argument.
		itsArgs[4]->Reset(partialReset);

		return true;
	}
	else
	{
		return ARM_ExpNodeFunc::Reset(partialReset);
	}
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


