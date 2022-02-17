/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file surfacemodelparam.cpp
 *
 *  \brief surface for model param
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2004
 */

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"
#include "gpinfra/surfacemodelparam.h"
#include <math.h>

/// gpbase
#include <glob/expt.h>			/// necessary for exception throwing
#include "gpbase/ostringstream.h"
#include "gpbase/surface.h"
#include "gpbase/surfacetypedef.h"
#include "gpbase/numericconstant.h"
#include "gpinfra/modelparamtype.h"

#include <iomanip> /// for setprecision()

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_SurfaceModelParam
///	Routine: Constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_SurfaceModelParam::ARM_SurfaceModelParam( 
	ParamType type, 
	ARM_Surface* surface, 
	const string& paramName,
	double lowerBound,
	double upperBound,
	bool adviseBreakPointTimes)
:
	ARM_ModelParam(type, adviseBreakPointTimes),
	itsSurface( surface? (ARM_Surface*) surface->Clone(): NULL ),
	itsLowerBound(lowerBound),
    itsUpperBound(upperBound),
	itsInitialSurface( surface ? (ARM_Surface*) surface->Clone() : NULL )
{
		if(		lowerBound == -ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER
			&&  upperBound ==  ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER )
		{
			itsLowerBound = ARM_CalibParamCst::ModelParamStdLowerBound( GetType() );
			itsUpperBound = ARM_CalibParamCst::ModelParamStdUpperBound( GetType() );
		}
}

	

////////////////////////////////////////////////////
///	Class  : ARM_SurfaceModelParam
///	Routine: Copy constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_SurfaceModelParam::ARM_SurfaceModelParam( const ARM_SurfaceModelParam& rhs )
:	ARM_ModelParam( rhs ), itsSurface( rhs.itsSurface ? (ARM_Surface*) rhs.itsSurface->Clone() : NULL ),
	itsLowerBound( rhs.itsLowerBound),
    itsUpperBound( rhs.itsUpperBound),
	itsInitialSurface( rhs.itsSurface ? (ARM_Surface*) rhs.itsSurface->Clone() : NULL )

{}

////////////////////////////////////////////////////
///	Class  : ARM_SurfaceModelParam
///	Routine: assignment operator
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_SurfaceModelParam& ARM_SurfaceModelParam::operator =( const ARM_SurfaceModelParam& rhs )
{
	if( this !=	 &rhs )
	{
		ARM_ModelParam::operator=( rhs );
		itsLowerBound		= rhs.itsLowerBound;
		itsUpperBound		= rhs.itsUpperBound;
		delete itsSurface;
		delete itsInitialSurface;
		itsSurface		  = rhs.itsSurface ? (ARM_Surface*) rhs.itsSurface->Clone() : NULL;
		itsInitialSurface = rhs.itsInitialSurface ? (ARM_Surface*) rhs.itsInitialSurface->Clone() : NULL;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_SurfaceModelParam
///	Routine: ~ARM_SurfaceModelParam
///	Returns: 
///	Action : destructor
////////////////////////////////////////////////////
///////////////////////////////////////
/// destructor
////////////////////////////////////////
ARM_SurfaceModelParam::~ARM_SurfaceModelParam()
{
	delete itsSurface;
	delete itsInitialSurface;
}


///////////////////////////////////////
///	Class  : ARM_SurfaceModelParam
///	Routine: Clone
///	Returns: ARM_Object*
///	Action : Clone the object
///////////////////////////////////////

ARM_Object* ARM_SurfaceModelParam::Clone() const
{
	return new ARM_SurfaceModelParam(*this); 
}


////////////////////////////////////////////////////
///	Class   : ARM_SurfaceModelParam
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_SurfaceModelParam::toString(const string& indent, const string& nextIndent) const
{

    CC_Ostringstream os;
	os << "\n"<< indent << "Model Param Surface Type: " << GetTypeString() << "\n";
	os << indent <<" lower bound = " << itsLowerBound << "\n";
	os << indent <<" upper bound = " << itsUpperBound << "\n";
	if( itsSurface )
	{
		if(itsSurface->GetX1().size() > 0 && itsSurface->GetX2().size() > 0)
		{
			double theVal = itsSurface->InterpolateAtPoint(0,0);
			if(*itsSurface == theVal)
				os << " constant values = " << theVal << "\n";
			else
				os << itsSurface->toString();
		}
		else
			os << " empty surface" << "\n";
	}
    return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_SurfaceModelParam
///	Routine: size
///	Returns: size_t
///	Action : returns the size of the model param
////////////////////////////////////////////////////
size_t ARM_SurfaceModelParam::size() const
{
	return itsSurface ? itsSurface->size(): 0;  
}

////////////////////////////////////////////////////
///	Class  : ARM_SurfaceModelParam
///	Routine: rows
///	Returns: size_t
///	Action : returns the size of the model param
////////////////////////////////////////////////////
size_t ARM_SurfaceModelParam::rows() const
{
	return GetSurface()->GetX1().size();  
}

////////////////////////////////////////////////////
///	Class  : ARM_SurfaceModelParam
///	Routine: cols
///	Returns: size_t
///	Action : returns the size of the model param
////////////////////////////////////////////////////
size_t ARM_SurfaceModelParam::cols() const
{
	return GetSurface()->GetX2().size();  
}


////////////////////////////////////////////////////
///	Class  : ARM_SurfaceModelParam
///	Routine: GetValue
///	Returns: double
///	Action : Various getvalue functions
////////////////////////////////////////////////////

double ARM_SurfaceModelParam::GetValue( double x1 ) const
{
    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetValue with only one argument unimplemented for a surface!" ); 
}

double ARM_SurfaceModelParam::GetValue( double x1, double x2 ) const
{
    return itsSurface->Interpolate(x1,x2); 
}

double ARM_SurfaceModelParam::GetValueAtPoint( size_t i ) const
{
    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetValue with only one argument unimplemented for a surface!" ); 
}

double ARM_SurfaceModelParam::GetTimeAtPoint( size_t i, size_t j ) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetTimeAtPoint(i ,j) unimplemented for a surface!" ); }
double ARM_SurfaceModelParam::GetTenorAtPoint( size_t i, size_t j ) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetTenorAtPoint(i,j) unimplemented for a curve!" ); }

double ARM_SurfaceModelParam::GetTimeAtPoint( size_t i) const
{   return itsSurface->GetX1()[i]; }
double ARM_SurfaceModelParam::GetTenorAtPoint( size_t i ) const
{	return itsSurface->GetX2()[i]; }

double ARM_SurfaceModelParam::GetValueAtPoint( size_t i, size_t j ) const
{
    return itsSurface->InterpolateAtPoint(i,j); 
}


double ARM_SurfaceModelParam::GetValue( double x1, double x2, double x3 ) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetValue with three arguments unimplemented for a surface!" ); }


////////////////////////////////////////////////////
///	Class  : ARM_SurfaceModelParam
///	Routine: SetValue
///	Returns: void
///	Action : Various SetValue functions
////////////////////////////////////////////////////

void ARM_SurfaceModelParam::SetValue( double x1, double value )
{	
    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetValue unimplemented for a surface!" ); 
}

void ARM_SurfaceModelParam::SetValue( double x1, double x2, double value  )
{	
    itsSurface->insert( x1, x2, value ); 
}

void ARM_SurfaceModelParam::SetValueAtPoint( size_t i, double value )
{
    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetValue unimplemented for a surface!" ); 
}

void ARM_SurfaceModelParam::SetValueAtPoint( size_t i, size_t j, double value )
{	
    itsSurface->insertAtPoint(i,j,value); 
}

void ARM_SurfaceModelParam::SetValueAtPoint( size_t i, size_t j, size_t k, double value )
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetValue unimplemented for a surface!" ); }


void ARM_SurfaceModelParam::SetValue( double x1, double x2, double x3, double value ) 
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetValue with three arguments unimplemented for a surface!" ); }

double ARM_SurfaceModelParam::GetValueAtPoint( size_t i, size_t j, size_t k) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetValue with three arguments unimplemented for a surface!" ); }


///////////////////////////////////////
///	Class  : ARM_SurfaceModelParam
///	Routine: SetAndUpDate
///	Returns: 
///	Action : Replace old values using interpolate if necessairy
///////////////////////////////////////
	
void ARM_SurfaceModelParam::SetAndUpDate(std::vector<double>* values)
{
	size_t k=0;
	for(size_t j1=0; j1<itsSurface->GetX3().rows(); ++j1 )
	{
		for(size_t j2=0; j2<itsSurface->GetX3().cols(); ++j2)
		{
			itsSurface->insert(itsSurface->GetX1()[j1],itsSurface->GetX2()[j2],(*values)[k]);
			++k;
		}
	}
}

///////////////////////////////////////
///	Class  : ARM_SurfaceModelParam
///	Routine: MergeModelParam
///	Returns: 
///	Action : merge a model param with another model param of same type
///         Warning! this assumes that the newValue has the same type
///////////////////////////////////////

void ARM_SurfaceModelParam::MergeModelParam(ARM_ModelParam* NewValue )
{
	ARM_SurfaceModelParam* surfaceModelParam = dynamic_cast<ARM_SurfaceModelParam*>(NewValue);

	if( !surfaceModelParam)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		    "trying to merge model param that is not a surface with a model param surface!" );

/// theoretically this should not be tested 
/// but in strict validation mode
/// we test and retest!
#if defined( __GP_STRICT_VALIDATION )
    if(surfaceModelParam->GetType() != GetType() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		    "trying to merge model param of a different type!" );
#endif
   ARM_FlatSurface* surface = dynamic_cast<ARM_FlatSurface*>(itsSurface);
   if(!surface)
	   itsSurface->insert( surfaceModelParam->GetSurface() );
   else
   {
	   delete itsSurface;
	   itsSurface = (ARM_Surface*)(surfaceModelParam->GetSurface())->Clone();
   }
}

///////////////////////////////////////
///	Class  : ARM_SurfaceModelParam
///	Routine: SetSurface
///	Returns: 
///	Action : Set the surface of a surface model param!
///////////////////////////////////////

void ARM_SurfaceModelParam::SetSurface( ARM_Surface* surface )
{ 
	delete itsSurface;
	itsSurface = surface; 
}

std::vector<double>* ARM_SurfaceModelParam::GetData( ARM_DataType type, long& rows, long& cols ) const
{
	return NULL;
	/*switch(type)
	{
		case ARM_ModelParamType::BreakPointTimes:
			{
				const std::vector<double>& result = itsSurface->GetX1();
				rows=result.size();
				cols=1;
				return new std::vector<double>(result);
			}
        case ARM_ModelParamType::Tenors:
            {
				const std::vector<double>& result = itsSurface->GetX2();
				rows=result.size();
				cols=1;
				return new std::vector<double>(result);
			}
        case ARM_ModelParamType::Values:
			{
				const ARM_GP_Matrix& values = itsSurface->GetX3();
				rows = values.rows();
				cols = values.cols();
				std::vector<double>& result = new std::vector<double>( values.size() );
				CC_NS(std,copy)(values.begin(), values.end(), result->begin());
				return result;
			}

        default:
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "Unknown type... permitted are values,BreakPointTimes or Tenors" );
	}*/
}


////////////////////////////////////////
///	Class  : ARM_SurfaceModelParam 
///	Routine: ToSurfaceModelParam
///	Returns: ARM_SurfaceModelParam
///	Action : downcast to an ARM_CurveModelParam
///////////////////////////////////////

ARM_SurfaceModelParam& ARM_SurfaceModelParam::ToSurfaceModelParam()
{ return *this; }

const ARM_SurfaceModelParam& ARM_SurfaceModelParam::ToSurfaceModelParam() const
{ return *this; }



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

