/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file modelparams.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/// this headers has to come first
/// as firsttoinc defines pragma for warning to avoid


///gpbase
#include "gpbase/curve.h"
#include "gpbase/surface.h"
#include "gpbase/gpvector.h"

///gpinfra
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/surfacemodelparam.h"



CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////
///	Class  : ARM_ModelParams 
///	Routine: Constructor 
///	Returns:
///	Action : 
///////////////////////////////////////
ARM_ModelParams::ARM_ModelParams( const ARM_ModelParamVector& params )
:	itsParams( ARM_ModelParamType::Unknown )
{
    for( int i=0; i < ARM_ModelParamType::Unknown; ++i )
        itsParams[i] = NULL;

	size_t paramSize = params.size();
	ARM_GP_T_Vector<int> paramTypes;
	paramTypes.reserve(paramSize);

	/// validation to check that we do not have two parameters with the same type
    for( i=0; i<paramSize; ++i )
		paramTypes.push_back(params[i]->GetType());
	paramTypes.sort();
	paramTypes.unique();

	if( paramTypes.size() != paramSize )
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"There is two parameters with the same type: strictly forbidden!" );

    /// clone the params from params
    for( i=0; i<params.size(); ++i )
        itsParams[params[i]->GetType()] = ( (ARM_ModelParam*) params[i]->Clone() );
}


///////////////////////////////////////
///	Class  : ARM_ModelParams 
///	Routine: CopyNoCleanUp 
///	Returns:
///	Action : 
///////////////////////////////////////
void ARM_ModelParams::CopyNoCleanUp( const ARM_ModelParams& rhs )
{
	itsParams.resize( rhs.itsParams.size() );

	ARM_ModelParamVector::iterator
			NewParam = begin();

	ARM_ModelParams::const_iterator
			Param    = rhs.begin(),
			ParamEnd = rhs.end();

	for( ; Param != ParamEnd; ++NewParam, ++Param )
        *NewParam = *Param ? static_cast<ARM_ModelParam*>( (*Param)->Clone() ) : NULL;
}


///////////////////////////////////////
///	Class  : ARM_ModelParams 
///	Routine: Copy Constructor 
///	Returns:
///	Action : 
///////////////////////////////////////
ARM_ModelParams::ARM_ModelParams( const ARM_ModelParams &rhs  )
:	ARM_RootObject( rhs ), itsParams( rhs.size() )
{	
	CopyNoCleanUp( rhs );
}


///////////////////////////////////////
///	Class  : ARM_ModelParams 
///	Routine: CleanUp
///	Returns:
///	Action : 
///////////////////////////////////////
void ARM_ModelParams::CleanUp()
{
	ARM_ModelParamVector::iterator
			Param    = begin(),
			ParamEnd = end();

	for( ; Param != ParamEnd; ++Param )
		delete *Param;
}

///////////////////////////////////////
///	Class  : ARM_ModelParams 
///	Routine: Destructor 
///	Returns:
///	Action : 
///////////////////////////////////////
ARM_ModelParams::~ARM_ModelParams()
{
	CleanUp();
}



///////////////////////////////////////
///	Class  : ARM_ModelParams 
///	Routine: operator= 
///	Returns: object with new values
///	Action : 
///////////////////////////////////////

ARM_ModelParams& ARM_ModelParams::operator=( const ARM_ModelParams& rhs )
{
	if( this != &rhs )
	{
		ARM_RootObject::operator=( rhs );
		/// Delete the old ones.
		CleanUp();
		/// Copy in the new ones.
		CopyNoCleanUp( rhs );
	}

	return *this;
}

///////////////////////////////////////
///	Class  : ARM_ModelParams 
///	Routine: SetModelParamValue
///	Returns: ARM_ModelParam*
///	Action : update the value at break point i on the model param
///         nb ParamNum and return the new model param
///////////////////////////////////////
ARM_ModelParams::iterator ARM_ModelParams::SetModelParamValue( 
	int paramType,
	size_t i,
	double value, 
	double time,
	double tenor,
	size_t factorNb )
{
	ARM_ModelParam* foundParam = itsParams[paramType];
	
    /// find out whether we have already a similar model param according to its type!
    if(foundParam)
    {
        if( ARM_CurveModelParam* curveModelParam = dynamic_cast<ARM_CurveModelParam*>(foundParam) )
		{
			double lag = curveModelParam->GetCurve()->GetAbscisses()[i];
			if(fabs(lag - time)< K_FRM_TOL)
				curveModelParam->SetValueAtPoint(i,value);
			else
				curveModelParam->SetValue(time,value);
			
			return itsParams.begin()+paramType;
		}
		else if( ARM_SurfaceModelParam* surfaceModelParam = dynamic_cast<ARM_SurfaceModelParam*>(foundParam) )
		{
			ARM_SurfaceWithInterpol* surface = dynamic_cast<ARM_SurfaceWithInterpol*>(surfaceModelParam->GetSurface());
			if( !surface )
				throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					"Expected a discretised surface with an interpolator!" );
			
			size_t tenorSize	= surface->GetX2().size();
			size_t k			= i%tenorSize;
			size_t j			= i/tenorSize;
			double paramExpiry	= surface->GetX1()[j];
			double paramTenor	= surface->GetX2()[k];

			if(	fabs(paramExpiry - time)< K_FRM_TOL && 
				fabs(paramTenor  - tenor)< K_FRM_TOL )
				surfaceModelParam->SetValueAtPoint(j,k,value);
			else
				surfaceModelParam->SetValue(time,tenor,value);
			return itsParams.begin()+paramType;
		}
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "this function currently only supports ARM_CurveModelParam!" );
    }
    else
    {
        CC_Ostringstream msg;
        msg << ARM_USERNAME << ": could not find model param for type " 
            << ARM_ModelParamType::GetTypeString( (ARM_ModelParamType::ParamNb) paramType);
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg.str() );
    }
}

///////////////////////////////////////
///	Class  : ARM_ModelParams 
///	Routine: GetModelParamValidate
///	Returns: A particular parameter
///	Action : 
///////////////////////////////////////
void ARM_ModelParams::GetModelParamValidate( int paramType ) const
{ 
    if( paramType >= ARM_ModelParamType::Unknown )
	{
		CC_Ostringstream msg;
		msg << "ARM_ModelParams::ModelParam(): Parameter index " << paramType << " not in range 0.." << ARM_ModelParamType::Unknown;
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg.str() );
	}

	if( NULL == itsParams[ paramType ] )
	{
		CC_Ostringstream msg;
		msg << "ARM_ModelParams::ModelParam(): Parameter " << itsParams[ paramType ]->GetTypeString() << " has not been set.";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg.str() );
	}
}

///////////////////////////////////////
///	Class  : ARM_ModelParams 
///	Routine: GetModelParam (const version)
///	Returns: a boolean
///	Action : Test if a model param with a 
/// specific type exists
///////////////////////////////////////
// 
bool ARM_ModelParams::DoesModelParamExist(int paramType, size_t factorNb  ) const
{
    return itsParams[ paramType ] != 0;
}


///////////////////////////////////////
///	Class  : ARM_ModelParams 
///	Routine: GetModelParam (const version)
///	Returns: A particular parameter
///	Action : 
///////////////////////////////////////

const ARM_ModelParam& ARM_ModelParams::GetModelParam( int paramType, size_t factorNb  ) const
{ 
	GetModelParamValidate( paramType );
	return *( itsParams[ paramType ] );
}


///////////////////////////////////////
///	Class  : ARM_ModelParams 
///	Routine: GetModelParam (non const version)
///	Returns: A particular parameter
///	Action : 
///////////////////////////////////////

ARM_ModelParam& ARM_ModelParams::GetModelParam( int paramType, size_t factorNb  )
{ 
	GetModelParamValidate( paramType );
	return *( itsParams[ paramType ] );
}
///////////////////////////////////////
///	Class  : ARM_ModelParams 
///	Routine: SetModelParam
///	Returns: 
///	Action : Remove the old param and Clone new param
///////////////////////////////////////
void ARM_ModelParams::SetModelParam(ARM_ModelParam* param, size_t factorNb )
{
	ARM_ModelParam* found = itsParams[param->GetType()];

    /// find out whether we have already a similar model param according to its type!
    if(found)
        delete found;

    itsParams[param->GetType()] = static_cast<ARM_ModelParam*>(param->Clone());
}

///////////////////////////////////////
///	Class  : ARM_ModelParams 
///	Routine: DeleteModelParam
///	Returns: 
///	Action : same as param type except 
///////////////////////////////////////
void ARM_ModelParams::DeleteModelParam( int paramType, size_t factorNb  )
{
	delete itsParams[paramType];
	itsParams[paramType] = NULL;
}


///////////////////////////////////////
///	Class  : ARM_ModelParams 
///	Routine: MergeModelParam
///	Returns: 
///	Action : merge and replace a particular model parameter if not null 
///          element by elment and return by increasing strictly BreakPointTimes
///////////////////////////////////////

void ARM_ModelParams::MergeModelParam(ARM_ModelParam* NewValue, size_t factorNb  )
{
    ARM_ModelParam* found = itsParams[NewValue->GetType()];

    /// find out whether we have already a similar model param according to its type!
    if(found)
        found->MergeModelParam(NewValue);
    else
        itsParams[NewValue->GetType()] = static_cast<ARM_ModelParam*>(NewValue->Clone());
}

///////////////////////////////////////
///	Class  : ARM_ModelParams
///	Routine: GetModelParamName
///	Returns: Name of a parameter of the model
///	Action : 
///////////////////////////////////////

string ARM_ModelParams::GetModelParamName( int paramType, size_t factorNb  ) const
{
	GetModelParamValidate( paramType );
	return itsParams[ paramType ]->GetTypeString();
}



///////////////////////////////////////
///	Class  : ARM_ModelParams
///	Routine: GetModelParamNames
///	Returns: Names of all the parameters of the model
///	Action : 
///////////////////////////////////////

ARM_StringVector ARM_ModelParams::GetModelParamNames() const
{	
    ARM_StringVector Result;

    int  paramType = 0;
	
    for( ; paramType < ARM_ModelParamType::Unknown; ++paramType )
    {
        if (DoesModelParamExist(paramType))
		    Result.push_back(GetModelParamName( paramType ));
    }

	return Result;
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParams
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_ModelParams::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
	size_t j=0;

    os << "\n\n";
    os << indent << "ARM_ModelParams\n";
    os << indent << "---------------\n";

	for(int i=0;i<itsParams.size();++i)
	{
        if (itsParams[i] != NULL)
        {
		    os << indent << "Param #" << CC_NS(std,setw)(2) << ++j;
            os << itsParams[i]->toString(indent);
            os << "\n";
        }
	}
    
    return os.str();
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

