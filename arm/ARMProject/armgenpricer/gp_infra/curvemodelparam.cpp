/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file curvemodelparam.cpp
 *
 *  \brief file for a model param that is a curve
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"
#include "gpinfra/curvemodelparam.h"

/// gpbase
#include <glob/expt.h>			/// necessary for exception throwing
#include "gpbase/ostringstream.h"
#include "gpbase/curve.h"
#include "gpbase/stringmanip.h"

#include "gpinfra/calibparamcst.h"

#include <iomanip> /// for setprecision()

CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////
///	Class   : ARM_CurveModelParam
///	Routines: CopyNoCleanUp
///	Returns :
///	Action  : copy the member data
////////////////////////////////////////

void ARM_CurveModelParam::CopyNoCleanUp( const ARM_CurveModelParam& rhs )
{
    itsLowerBound	= rhs.itsLowerBound		? (ARM_GP_Vector*) rhs.itsLowerBound	: NULL;
    itsUpperBound	= rhs.itsUpperBound		? (ARM_GP_Vector*) rhs.itsUpperBound	: NULL;
    itsInitialCurve = rhs.itsInitialCurve	? (ARM_Curve*) rhs.itsInitialCurve		: NULL;
	itsCurve		= rhs.itsCurve			? (ARM_Curve*) rhs.itsCurve			: NULL;
	itsName			= rhs.itsName;
	itsCurrency		= rhs.itsCurrency;
}


///////////////////////////////////////
///	Class   : ARM_CurveModelParam
///	Routines: CleanUp
///	Returns :
///	Action  : Clean Up
////////////////////////////////////////

void ARM_CurveModelParam::CleanUp( )
{
    delete itsLowerBound;
    itsLowerBound = NULL;
    
	delete itsUpperBound;
    itsUpperBound = NULL;
    
	delete itsInitialCurve;
    itsInitialCurve = NULL;
	
	delete itsCurve;
	itsCurve = NULL;  
}

////////////////////////////////////////////////////
///	Class  : ARM_CurveModelParam
///	Routine: Default Constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_CurveModelParam::ARM_CurveModelParam( )
:
	ARM_ModelParam( ARM_ModelParamType::Unknown, true ),
    itsCurve(NULL),
	itsLowerBound(NULL),
    itsUpperBound(NULL),
    itsInitialCurve(NULL)
{}


////////////////////////////////////////////////////
///	Class  : ARM_CurveModelParam
///	Routine: Constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_CurveModelParam::ARM_CurveModelParam( 
	ParamType type, 
	ARM_GP_Vector* values, 
	ARM_GP_Vector* breakPointTimes,
	const string& paramName,			/// not used (here for legacy reason)
	const string& interpolatorName,
    ARM_GP_Vector* lowerBound,
    ARM_GP_Vector* upperBound,
	bool adviseBreakPointTimes,
	const string&	currency):	ARM_ModelParam(type,adviseBreakPointTimes),
	itsName(paramName),
    itsCurve(NULL),
	itsLowerBound(NULL),
    itsUpperBound(NULL),
    itsInitialCurve(NULL),
	itsCurrency(currency)
{
	ARM_Interpolator<double, double>* interpolator = NULL;
	
	/// to make it case insensitive!
	string interpolatorNameUpper = stringGetUpper(interpolatorName);
	if( interpolatorNameUpper == "LINEAR" )

		interpolator = new ARM_LinInterpCstExtrapolDble ;
	else if( interpolatorNameUpper == "STEPUPRIGHT" )

		interpolator = new ARM_StepUpRightOpenCstExtrapolDble ;
	else if( interpolatorNameUpper == "STEPUPLEFT" )

		interpolator = new ARM_StepUpLeftOpenCstExtrapolDble ;
	else
	{
		CC_Ostringstream msg;
		msg << " ARM ERR: unknown interpolator type: permitted are LINEAR, STEPUPRIGHT, STEPUPLEFT!" ;
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg.str() );
	}
	
    itsCurve = new ARM_Curve((breakPointTimes ? *breakPointTimes : ARM_GP_Vector(0)),
                             (values ? *values : ARM_GP_Vector(0)),interpolator);

	if(lowerBound && values)
        CC_NS(ARM_Check,CheckSameArgSize)(*lowerBound, *values, "lowerBound", "values",__LINE__,__FILE__);
    if(upperBound && values)
        CC_NS(ARM_Check,CheckSameArgSize)(*upperBound, *values, "upperBound", "values",__LINE__,__FILE__);
        
    itsLowerBound   = lowerBound ? new ARM_GP_Vector(*lowerBound) : new ARM_GP_Vector(size(), ARM_CalibParamCst::ModelParamStdLowerBound( GetType() ));
    itsUpperBound   = lowerBound ? new ARM_GP_Vector(*upperBound) : new ARM_GP_Vector(size(), ARM_CalibParamCst::ModelParamStdUpperBound( GetType() ));
    itsInitialCurve = itsCurve ? (ARM_Curve*)itsCurve->Clone() : NULL; 
};

	////////////////////////////////////////////////////
///	Class  : ARM_CurveModelParam
///	Routine: constructor
///	Returns: 
///	Action : To build a constructor mono-dim
////////////////////////////////////////////////////
ARM_CurveModelParam::ARM_CurveModelParam( ParamType type,
	double value,
	const string& paramName,	
	double lowerBound,
	double upperBound,
	bool adviseBreakPointTimes,
	const string &	currency):ARM_ModelParam(type, adviseBreakPointTimes),
	itsName(paramName),
    itsCurve(NULL),
	itsLowerBound(NULL),
    itsUpperBound(NULL),
    itsInitialCurve(NULL),
	itsCurrency(currency){
	itsCurve = new ARM_FlatCurve(value);
	
	itsLowerBound   = new ARM_GP_Vector(1,lowerBound); 
	itsUpperBound   = new ARM_GP_Vector(1,upperBound); 
	itsInitialCurve = itsCurve ? (ARM_Curve*)itsCurve->Clone() : NULL; 

}

////////////////////////////////////////////////////
///	Class  : ARM_CurveModelParam
///	Routine: Copy constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_CurveModelParam::ARM_CurveModelParam( const ARM_CurveModelParam& rhs )
:	ARM_ModelParam( rhs )
{
	CopyNoCleanUp( rhs );
}


////////////////////////////////////////////////////
///	Class  : ARM_CurveModelParam
///	Routine: assignment operator
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_CurveModelParam& ARM_CurveModelParam::operator =( const ARM_CurveModelParam& rhs )
{
	if( this !=	 &rhs )
	{
		CleanUp();
		ARM_ModelParam::operator=( rhs );
		CopyNoCleanUp( rhs );
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_CurveModelParam
///	Routine: ~ARM_CurveModelParam
///	Returns: 
///	Action : destructor
////////////////////////////////////////////////////
///////////////////////////////////////
/// destructor
////////////////////////////////////////
ARM_CurveModelParam::~ARM_CurveModelParam()
{
	CleanUp();
}


////////////////////////////////////////
///	Class  : ARM_CurveModelParam 
///	Routine: SetAndUpDate
///	Returns: void
///	Action : replace  and Update a particular itsValues from modelParam
///////////////////////////////////////
void ARM_CurveModelParam::SetAndUpDate(ARM_GP_Vector* values)
{
    if(values->size() != itsCurve->GetAbscisses().size())
    {            
        itsCurve->SetOrdinates(*values);
        ARM_Curve* curve = itsCurve->CptCurve(itsCurve->GetAbscisses()); 
        delete itsCurve;
        itsCurve = curve;
    }
    else
    {
        itsCurve->SetOrdinates(*values);
    }
}

////////////////////////////////////////
///	Class  : ARM_CurveModelParam 
///	Routine: UpdateValues
///	Returns: void
///	Action : 
///////////////////////////////////////
void ARM_CurveModelParam::UpdateValues(ARM_GP_Vector* breakPointTimes)
{
    if(breakPointTimes)
    {
		 /// to update itsLowerBound and itsUpperBound
		ARM_Curve tempcurve(itsCurve->GetAbscisses(), *itsLowerBound, itsInitialCurve->GetInterpolator()->Clone());
        ARM_Curve* curve = tempcurve.CptCurve(*breakPointTimes);
        delete itsLowerBound;
        itsLowerBound = new ARM_GP_Vector(curve->GetOrdinates());
        tempcurve.SetOrdinates(*itsUpperBound);
        delete curve;
        curve = tempcurve.CptCurve(*breakPointTimes);
        delete itsUpperBound;
        itsUpperBound = new ARM_GP_Vector(curve->GetOrdinates());
        delete curve;

        /// to update initial curve 
        curve = itsInitialCurve->CptCurve(*breakPointTimes);
        delete itsInitialCurve;
        itsInitialCurve = curve;

		 /// to update curve
        ARM_Curve* foritscurve = itsCurve->CptCurve(*breakPointTimes);
        delete itsCurve;
        itsCurve = foritscurve;
    }
}


////////////////////////////////////////
///	Class  : ARM_CurveModelParam 
///	Routine: SetValuesAndTimes
///	Returns: void
///	Action : changes the values and times for the given
///				breakpontTimes and values
///////////////////////////////////////

void ARM_CurveModelParam::SetValuesAndTimes(ARM_GP_Vector* breakPointTimes,ARM_GP_Vector* values)
{
    if(breakPointTimes && values)
    {
        itsCurve->SetAbscisses(*breakPointTimes);
        itsCurve->SetOrdinates(*values);

    }
}

///////////////////////////////////////
///	Class  : ARM_CurveModelParam
///	Routine: MergeModelParam
///	Returns: 
///	Action : merge a model param with another model param of same type
///         Warning! this assumes that the newValue has the same type
///////////////////////////////////////

void ARM_CurveModelParam::MergeModelParam(ARM_ModelParam* NewValue )
{
	ARM_CurveModelParam* curveModelParam = dynamic_cast<ARM_CurveModelParam*>(NewValue);

	if( !curveModelParam)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		    "trying to merge model param that is not a curve with a model param curve!" );

/// theoretically this should not be tested 
/// but in strict validation mode
/// we test and retest!
#if defined( __GP_STRICT_VALIDATION )
    if(curveModelParam->GetType() != GetType() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		    "trying to merge model param of a different type!" );
#endif

	ARM_FlatCurve* curve = dynamic_cast<ARM_FlatCurve*>(itsCurve);
	if(!curve)
	{
		/// these are the vector that will contain the 
		/// merged values!
		ARM_Curve Lowercurve(itsCurve->GetAbscisses(), *itsLowerBound);
		ARM_Curve Uppercurve(itsCurve->GetAbscisses(), *itsUpperBound);

		itsCurve->insert( curveModelParam->GetCurve() );

		ARM_Curve Lowertempcurve(curveModelParam->GetCurve()->GetAbscisses(), (*(curveModelParam->GetLowerBound()))); 
		Lowercurve.insert(&Lowertempcurve);
		delete itsLowerBound;
		itsLowerBound = new ARM_GP_Vector(Lowercurve.GetOrdinates());

		ARM_Curve Uppertempcurve(curveModelParam->GetCurve()->GetAbscisses(), (*(curveModelParam->GetUpperBound())));
		Uppercurve.insert(&Uppertempcurve);
		delete itsUpperBound;
		itsUpperBound = new ARM_GP_Vector(Uppercurve.GetOrdinates());
	}
	else
	{
		CleanUp();
		CopyNoCleanUp(*(ARM_CurveModelParam*)NewValue);
	}
}

///////////////////////////////////////
///	Class  : ARM_CurveModelParam
///	Routine: Clone
///	Returns: ARM_Object*
///	Action : Clone the object
///////////////////////////////////////

ARM_Object* ARM_CurveModelParam::Clone() const
{
	return new ARM_CurveModelParam(*this); 
}


////////////////////////////////////////////////////
///	Class   : ARM_CurveModelParam
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_CurveModelParam::toString(const string& indent, const string& nextIndent) const
{

    CC_Ostringstream os;
	os << "\n"<< indent << "Model Param Curve Type: " << GetTypeString() << " Interpolator Type: " << itsCurve->GetInterpolator()->toString() << "\n";
    
	if ( itsCurrency.size()!=0)
		os << "\n"<< indent << "Currency: " << itsCurrency << "\n";

	os << indent <<"\n";

    os << CC_NS(std,fixed);
	if( itsCurve )
	{
		if( itsCurve->GetAbscisses().size()>0)
		{
			os << indent << "Times \t Values(%)\t LowerBound \t UpperBound\n";
			for(int i=0;i<GetCurve()->GetAbscisses().size();++i)
			{
				os << indent << CC_NS(std,fixed)    << CC_NS(std,setw)(0) << CC_NS(std,setprecision)(0) << itsCurve->GetAbscisses()[i];
				os << "\t " << CC_NS(std,setw)(10)  << CC_NS(std,setprecision)(10) << itsCurve->GetOrdinates()[i]*100;
				os << "\t " << CC_NS(std,scientific)<<CC_NS(std,setw)(2) << CC_NS(std,setprecision)(2) << (*itsLowerBound)[i];
				os << "\t " << CC_NS(std,setw)(2)   << CC_NS(std,setprecision)(2) << (*itsUpperBound)[i];
				os << "\n";
			}
		}
        else
        {
	        os << indent << "No Times, no values set\n";
        }
	}
    return os.str();
}



////////////////////////////////////////////////////
///	Class  : ARM_CurveModelParam
///	Routine: SetCurve
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_CurveModelParam::SetCurve(ARM_Curve* curve )
{ 
	delete itsCurve; 
	itsCurve = curve;
}


////////////////////////////////////////////////////
///	Class  : ARM_CurveModelParam
///	Routine: GetData
///	Returns: ARM_GP_Vector*
///	Action : returns either the values or the break point times of a model param
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_CurveModelParam::GetData( ARM_DataType type ) const
{
	switch(type)
	{
		case ARM_ModelParamType::BreakPointTimes:
            return &(itsCurve->GetAbscisses());
        case ARM_ModelParamType::Values:
            return &(itsCurve->GetOrdinates());
        default:
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "Unknown type... permitted are Values and BreakPointTimes" );
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_CurveModelParam
///	Routine: size
///	Returns: size_t
///	Action : returns the size of the model param
////////////////////////////////////////////////////
size_t ARM_CurveModelParam::size() const
{
	return itsCurve ? itsCurve->size(): 0;  
}

////////////////////////////////////////////////////
///	Class  : ARM_CurveModelParam
///	Routine: rows
///	Returns: size_t
///	Action : returns the size of the model param
////////////////////////////////////////////////////
size_t ARM_CurveModelParam::rows() const
{
	return itsCurve ? itsCurve->size(): 0;  
}
////////////////////////////////////////////////////
///	Class  : ARM_CurveModelParam
///	Routine: cols
///	Returns: size_t
///	Action : returns the size of the model param
////////////////////////////////////////////////////
size_t ARM_CurveModelParam::cols() const
{
	return 1;  
}



////////////////////////////////////////////////////
///	Class  : ARM_CurveModelParam
///	Routine: GetValue
///	Returns: double
///	Action : Various getvalue functions
////////////////////////////////////////////////////

double ARM_CurveModelParam::GetValue( double x1 ) const
{	return itsCurve->Interpolate(x1); }

double ARM_CurveModelParam::GetValue( double x1, double x2 ) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetValue unimplemented for a curve!" ); }

double ARM_CurveModelParam::GetTimeAtPoint( size_t i ) const
{	return itsCurve->GetAbscisse(i); }

double ARM_CurveModelParam::GetTimeAtPoint( size_t i, size_t j ) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetTimeAtPoint(i,j) unimplemented for a curve!" ); }
double ARM_CurveModelParam::GetTenorAtPoint( size_t i ) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetTenorAtPoint(i) unimplemented for a curve!" ); }
double ARM_CurveModelParam::GetTenorAtPoint( size_t i, size_t j ) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetTenorAtPoint(i,j) unimplemented for a curve!" ); }

double ARM_CurveModelParam::GetValueAtPoint( size_t i ) const
{	return itsCurve->GetOrdinate(i); }

double ARM_CurveModelParam::GetValueAtPoint( size_t i, size_t j ) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetValue unimplemented for a curve!" ); }

double ARM_CurveModelParam::GetValue( double x1, double x2, double x3 ) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetValue with three arguments unimplemented for a surface!" ); }


////////////////////////////////////////////////////
///	Class  : ARM_CurveModelParam
///	Routine: SetValue
///	Returns: void
///	Action : Various SetValue functions
////////////////////////////////////////////////////

void ARM_CurveModelParam::SetValue( double x1, double value )
{	itsCurve->insert( x1, value ); }

void ARM_CurveModelParam::SetValue( double x1, double x2, double value  )
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetValue unimplemented for a curve!" ); }

void ARM_CurveModelParam::SetValueAtPoint( size_t i, double value )
{	itsCurve->GetOrdinate(i) = value; }

void ARM_CurveModelParam::SetValueAtPoint( size_t i, size_t j, double value )
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetValue unimplemented for a curve!" ); }

void ARM_CurveModelParam::SetValueAtPoint( size_t i, size_t j, size_t k,double value )
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetValue unimplemented for a curve!" ); }
 
void ARM_CurveModelParam::SetValue( double x1, double x2, double x3, double value )
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetValue with three arguments unimplemented for a curve!" ); }

double ARM_CurveModelParam::GetValueAtPoint( size_t i, size_t j, size_t k) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetValue with three arguments unimplemented for a curve!" ); }

///	Class  : ARM_CurveModelParam 
///	Routine: MultiplyOrdinates
///	Returns: void
///	Action : multiply the ordinates of the curve
///////////////////////////////////////


////////////////////////////////////////
///	Class  : ARM_CurveModelParam 
///	Routine: ToCurveModelParam
///	Returns: ARM_CurveModelParam
///	Action : downcast to an ARM_CurveModelParam
///////////////////////////////////////

ARM_CurveModelParam& ARM_CurveModelParam::ToCurveModelParam()
{ return *this; }

const ARM_CurveModelParam& ARM_CurveModelParam::ToCurveModelParam() const
{ return *this; }


///////////////////////////////////////
///	Class  : ARM_CurveModelParam 
///	Routine: GetSubCalibParm
///	Returns: ARM_CurveModelParam
///	Action : gets all elements of range [begin, end)
///////////////////////////////////////
ARM_CurveModelParam* ARM_CurveModelParam::GetCalibParam(size_t begin, size_t end)
{
   if(begin > end){size_t tmp = begin;begin = end;end = begin;}

    if((end > size()))
        throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
        " impossible to get any asset out mother range" );

    
    size_t i, k, size  = end - begin;
    ARM_GP_Vector values(size), times (size);
    ARM_GP_Vector lowerboud(size),upperbound(size);
    
    for( i=begin,k = 0; i<end; ++i, ++k)
    {
        times[k] = GetCurve()->GetAbscisses()[i];
        values[k]  = GetCurve()->GetOrdinates()[i];
        lowerboud[k]  = (*itsLowerBound)[i];
        upperbound[k] = (*itsUpperBound)[i];
    }

    ARM_CurveModelParam* calibparam = new ARM_CurveModelParam(GetType(),
        &times,
		&values,
        GetTypeString(),
		"STEPUPRIGHT",
        &lowerboud,
        &upperbound);

    return calibparam;
} 
////////////////////////////////////////
///	Class  : ARM_CurveModelParam 
///	Routine: BumpModelParam
///	Returns: void
///	Action : bumpa model param
///////////////////////////////////////
void ARM_CurveModelParam::BumpModelParam(size_t row, size_t column ,double shift, ARM_BumpParamType bumpType)
{
	switch(bumpType)
	{
		case ARM_ModelParamBump::isCumulative:
			{
				for(size_t i(0); i<row; ++i)
					itsCurve->GetOrdinate(i) += shift;
			}
			break;

        case ARM_ModelParamBump::isPerturbative:
			itsCurve->GetOrdinate(row-1) += shift;
			break;

        default:
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "Unknown type... permitted are Cumulative and Perturbative" );
	}
}
   

////////////////////////////////////////
///	Class  : ARM_CurveModelParam 
///	Routine: SetToInitialCurve
///	Returns: void
///	Action : puts back the curve to its initial value
///////////////////////////////////////

void ARM_CurveModelParam::SetToInitialCurve()
{
	SetCurve( itsInitialCurve );	
}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///	ARM_RangedCurveModelParam 
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_RangedCurveModelParam
///	Routine: Constructor
////////////////////////////////////////////////////

ARM_RangedCurveModelParam::ARM_RangedCurveModelParam( 
	ParamType type,
	ARM_GP_Vector* values,
	ARM_GP_Vector* breakPointTimes, 
	const string& paramName,
    const string& interpolatorName )
:	
	ARM_CurveModelParam( type, values, breakPointTimes, paramName, interpolatorName )
{
	for( size_t i=0; i<values->size(); ++i )
		SetValueAtPoint(i,(*values)[i]);
}


////////////////////////////////////////
///	Class  : ARM_RangedCurveModelParam 
///	Routine: SetValuesAndTimes
///	Returns: void
///	Action : changes the values and times for the given
///				breakpontTimes and values
///////////////////////////////////////
	void ARM_RangedCurveModelParam::SetValuesAndTimes(ARM_GP_Vector* breakPointTimes,ARM_GP_Vector* values)
{
	for( size_t i=0; i<values->size(); ++i )
		(*values)[i]=ModifiedValue((*values)[i]);
	ARM_CurveModelParam::SetValuesAndTimes( breakPointTimes,values);
}





CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

