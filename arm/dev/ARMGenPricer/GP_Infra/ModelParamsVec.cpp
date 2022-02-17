/*!
 *
 * Copyright (c) IXIS CIB February 2004 Paris
 *
 *	\file modelparamsvec.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date February 2004
 */


/// gpinfra
#include "gpinfra/modelparamsvec.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrix.h"

/// gpcalib
//#include "gpcalib/modelfitter.h"

/// kernel headers
#include <inst/portfolio.h>

CC_BEGIN_NAMESPACE( ARM )


///////////////////////////////////////
///	Class  : ARM_ModelParamsVec 
///	Routine: CheckRange
///	Returns: void
///	Action : throw an exception if out of range!
///////////////////////////////////////
void ARM_ModelParamsVec::CheckRange( size_t factorNb ) const
{
	if( factorNb >= itsParamsVector.size() )
		ARM_THROW(ERR_INVALID_ARGUMENT,  ": out of range in ARM_ModelParamsVec::" );
}



///////////////////////////////////////
///	Class  : ARM_ModelParamsVec
///	Routine: GetModelParamName
///	Returns: Name of a parameter of the model
///	Action : 
///////////////////////////////////////

string ARM_ModelParamsVec::GetModelParamName( int paramType, size_t factorNb  ) const
{
#if defined(__GP_STRICT_VALIDATION)
	CheckRange( factorNb );
#endif
    return itsParamsVector[factorNb]->GetModelParamName( paramType, factorNb );
}


///////////////////////////////////////
///	Class  : ARM_ModelParamsVec 
///	Routine: DoesModelParamExist (const version)
///	Returns: a boolean
///	Action : Test if a model param with a 
/// specific type exists
///////////////////////////////////////
bool ARM_ModelParamsVec::DoesModelParamExist(int paramType, size_t factorNb  ) const
{
#if defined(__GP_STRICT_VALIDATION)
	CheckRange( factorNb );
#endif
    return itsParamsVector[factorNb]->DoesModelParamExist( paramType, factorNb );
}


///////////////////////////////////////
///	Class  : ARM_ModelParamsVec 
///	Routine: GetModelParam (const version)
///	Returns: A particular parameter
///	Action : 
///////////////////////////////////////

const ARM_ModelParam& ARM_ModelParamsVec::GetModelParam( int paramType, size_t factorNb  ) const
{ 
#if defined(__GP_STRICT_VALIDATION)
	CheckRange( factorNb );
#endif
    return itsParamsVector[factorNb]->GetModelParam( paramType, factorNb );
}


///////////////////////////////////////
///	Class  : ARM_ModelParamsVec 
///	Routine: GetModelParam (non const version)
///	Returns: A particular parameter
///	Action : 
///////////////////////////////////////

ARM_ModelParam& ARM_ModelParamsVec::GetModelParam( int paramType, size_t factorNb  )
{ 
#if defined(__GP_STRICT_VALIDATION)
	CheckRange( factorNb );
#endif
    return itsParamsVector[factorNb]->GetModelParam( paramType, factorNb );
}




///////////////////////////////////////
///	Class  : ARM_ModelParamsVec 
///	Routine: SetModelParam
///	Returns: 
///	Action : Remove the old param and Clone new param
///////////////////////////////////////
void ARM_ModelParamsVec::SetModelParam(ARM_ModelParam* param, size_t factorNb )
{
#if defined(__GP_STRICT_VALIDATION)
	CheckRange( factorNb );
#endif
    itsParamsVector[factorNb]->SetModelParam( param, factorNb );
}

///////////////////////////////////////
///	Class  : ARM_ModelParamsVec 
///	Routine: DeleteModelParam
///	Returns: 
///	Action : same as param type except 
///////////////////////////////////////
void ARM_ModelParamsVec::DeleteModelParam( int paramType, size_t factorNb  )
{
#if defined(__GP_STRICT_VALIDATION)
	CheckRange( factorNb );
#endif
    itsParamsVector[factorNb]->DeleteModelParam( paramType, factorNb );
}


///////////////////////////////////////
///	Class  : ARM_ModelParamsVec 
///	Routine: SetModelParamValue
///	Returns: ARM_ModelParam*
///	Action : update the value at break point i on the model param
///         nb ParamNum and return the new model param
///////////////////////////////////////
ARM_ModelParamVector::iterator ARM_ModelParamsVec::SetModelParamValue( 
	int paramType,
	size_t i,
	double value, 
	double time,
	double tenor,
	size_t factorNb )
{
#if defined(__GP_STRICT_VALIDATION)
	CheckRange( factorNb );
#endif
	return itsParamsVector[factorNb]->SetModelParamValue( paramType, i, value, 
		time, tenor, factorNb );
}


///////////////////////////////////////
///	Class  : ARM_ModelParamsVec 
///	Routine: MergeModelParam
///	Returns: 
///	Action : merge and replace a particular model parameter if not null 
///          element by elment and return by increasing strictly BreakPointTimes
///////////////////////////////////////
void ARM_ModelParamsVec::MergeModelParam(ARM_ModelParam* NewValue, size_t factorNb  )
{
#if defined(__GP_STRICT_VALIDATION)
	CheckRange( factorNb );
#endif
    itsParamsVector[factorNb]->MergeModelParam( NewValue, factorNb );
}

///////////////////////////////////////
///	Class  : ARM_ModelParamsVec 
///	Routine: PreProcessing
///	Returns: void
///	Action : PreProcessing
///////////////////////////////////////
void ARM_ModelParamsVec::PreProcessing(ARM_ModelFitter& modelFitter,int factorNb)
{
#if defined(__GP_STRICT_VALIDATION)
	CheckRange( factorNb );
#endif
	itsParamsVector[factorNb]->PreProcessing( modelFitter,factorNb );
}


///////////////////////////////////////
///	Class  : ARM_ModelParamsVec 
///	Routine: PreProcessing
///	Returns: void
///	Action : PreProcessing
///////////////////////////////////////
void ARM_ModelParamsVec::PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb )
{
#if defined(__GP_STRICT_VALIDATION)
	CheckRange( factorNb );
#endif
	itsParamsVector[factorNb]->PostProcessing( modelFitter, model, factorNb );
}



////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsVec
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_ModelParamsVec::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " Model Params Vector object\n";
	for( size_t i=0; i<itsParamsVector.size(); ++i )
		os << itsParamsVector[i]->toString( indent, nextIndent );
	return os.str();
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsVec
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
ARM_ModelParamVector ARM_ModelParamsVec::GetModelParams(size_t factorNb) const
{
#if defined(__GP_STRICT_VALIDATION)
	CheckRange(factorNb);
#endif
	return itsParamsVector[factorNb]->GetModelParams();
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

