/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file correlmatparam.cpp
 *
 *  \brief a correlation matrix parameter
 *		based on curve model param but with matrix
 *		data
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpinfra/correlmatparam.h"

#include "gpbase/curve.h"
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/vectorarithmetic.h"
#include "gpbase/numericconstant.h"
#include "gpbase/gptrigomatrix.h"
#include "gpbase/datestrip.h"
 
CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_CorrelMatParam
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_CorrelMatParam::ARM_CorrelMatParam( 
	ParamType type,
	ARM_GP_Vector* breakPointTimes,
    ARM_GP_Vector* values,
	const string& interpolatorName,
	ARM_GP_Vector* lowerBound,
    ARM_GP_Vector* upperBound,
	bool adviseBreakPointTimes)
:	
	ARM_CurveModelParam( type, values, breakPointTimes, "Correlation", interpolatorName, lowerBound, upperBound, adviseBreakPointTimes), 
	itsCurves(NULL),
	itsRealizedCorrelMatrix(NULL)
{	
	CreateMultiCurves(interpolatorName); 
}


////////////////////////////////////////////////////
///	Class  : ARM_CorrelMatParam
///	Routine: Copy Constructor
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////
ARM_CorrelMatParam::ARM_CorrelMatParam( const ARM_CorrelMatParam& rhs )
:	
	ARM_CurveModelParam( rhs ), 
	itsCurves(rhs.itsCurves ? (ARM_MultiCurve*)rhs.itsCurves->Clone() : NULL),
	itsRealizedCorrelMatrix( rhs.itsRealizedCorrelMatrix? (ARM_GP_Matrix*) rhs.itsRealizedCorrelMatrix->Clone() : NULL )
{}



////////////////////////////////////////////////////
///	Class  : ARM_CorrelMatParam
///	Routine: operator=
///	Returns: 
///	Action : assignment operator
////////////////////////////////////////////////////

ARM_CorrelMatParam& ARM_CorrelMatParam::operator=( const ARM_CorrelMatParam& rhs )
{
	if( this != & rhs )
	{
		delete itsCurves;
		delete itsRealizedCorrelMatrix;
	    itsCurves				= rhs.itsCurves ? (ARM_MultiCurve*)rhs.itsCurves->Clone() : NULL;
		itsRealizedCorrelMatrix = rhs.itsRealizedCorrelMatrix? (ARM_GP_Matrix*) rhs.itsRealizedCorrelMatrix->Clone() : NULL;
	}
	return *this;
};


////////////////////////////////////////////////////
///	Class  : ARM_CorrelMatParam
///	Routine: ~ARM_CorrelMatParam
///	Returns: 
///	Action : destructor
////////////////////////////////////////////////////

ARM_CorrelMatParam::~ARM_CorrelMatParam()
{
	delete itsCurves;
	delete itsRealizedCorrelMatrix;
}

////////////////////////////////////////////////////
///	Class  : ARM_CorrelMatParam
///	Routine: SetMultiCurve
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_CorrelMatParam::SetMultiCurve(ARM_MultiCurve* mCurve)
{ 
	delete itsCurves;
	itsCurves=mCurve; 
}

////////////////////////////////////////////////////
///	Class  : ARM_CorrelMatParam
///	Routine: SetRealizedCorrel
///	Returns: 
///	Action : gives the realized correlation matrix
////////////////////////////////////////////////////
void ARM_CorrelMatParam::UpdateRealizedCorrel()
{
	delete itsRealizedCorrelMatrix;

	size_t size = itsCurves->size();
	itsRealizedCorrelMatrix = new ARM_GP_Matrix( size, size );
	size_t factorNb = (*itsCurves)[0].size();
	
	size_t i,j,k;
	ARM_Vector norm2(size,0.0);
	for( i=0; i<size; ++i )
	{
		for( k=0; k<factorNb; ++k )
			norm2[i] += (*itsCurves)[i][k] * (*itsCurves)[i][k];
	}

	for( i=0; i<size; ++i )
	{
		for( j=0; j<=i; ++j )
		{
			double sum = 0.0;
			for( k=0; k<factorNb; ++k )
				sum += (*itsCurves)[i][k] * (*itsCurves)[j][k];
			(*itsRealizedCorrelMatrix)(i,j) = (*itsRealizedCorrelMatrix)(j,i) = sum/norm2[i];
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_CorrelMatParam
///	Routine: CreateMultiCurves
///	Returns: 
///	Action : create the multi curve with the correct interpolator
////////////////////////////////////////////////////

void ARM_CorrelMatParam::CreateMultiCurves( const string& interpolatorName )
{
	itsCurves = new ARM_MultiCurve;
	ARM_MultiCurveInterpolator* interpolator;

	if( interpolatorName != "" )
	{
		if( interpolatorName == "LINEAR" )
			interpolator = new ARM_LinInterpCstExtrapolVec ;
		else if( interpolatorName == "STEPUPRIGHT" )
			interpolator = new ARM_StepUpRightOpenCstExtrapolVec ;
		else if( interpolatorName == "STEPUPLEFT" )
			interpolator = new ARM_StepUpLeftOpenCstExtrapolVec ;
		else
		{
			CC_Ostringstream msg;
			msg << " ARM ERR: unknown interpolator type: permitted are LINEAR, STEPUPRIGHT, STEPUPLEFT!" ;
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg.str() );
		}
	}
	else
	{
		interpolator = new ARM_StepUpRightOpenCstExtrapolVec ;
	}
	itsCurves->SetInterpolator(interpolator);
}


////////////////////////////////////////////////////
///	Class  : ARM_CorrelMatParam
///	Routine: GetData
///	Returns: 
///	Action : get multicurve datas
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_CorrelMatParam::GetData( ARM_DataType type, long& rows, long& cols ) const
{
	switch(type)
	{
		case ARM_ModelParamType::BreakPointTimes:
			{
				const ARM_GP_Vector& result = GetCurve()->GetAbscisses();
				rows=result.size();
				cols=1;
				return new ARM_GP_Vector(result);
			}
        case ARM_ModelParamType::Values:
			{
				if((rows=itsCurves->size())<1)
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "No ordinates found" );
				cols = (itsCurves->GetOrdinates())[0].size();
				ARM_GP_Vector* result = new ARM_GP_Vector(rows*cols,0.0);
				for(size_t i=0,k=0;i<rows;++i)
				{
					const ARM_GP_Vector& val = (itsCurves->GetOrdinates())[i];
					for(size_t j=0;j<cols;++j,++k)
						(*result)[k] = val[j];
				}

				return result;
			}

        default:
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
				"Unknown type... permitted are Values or BreakPointTimes" );
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_CorrelMatParam
///	Routine: toString
///	Returns: 
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_CorrelMatParam::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;

	os << ARM_CurveModelParam::toString(indent, nextIndent);

	if( itsCurves && itsCurves->size() )
	{
		os << "\n";
		os << indent << "Multi Curves\n";
		os << indent << "------------\n";
		os << itsCurves->toString( indent, nextIndent);
	}

	if( itsRealizedCorrelMatrix )
	{
		os << "\n";
		os << indent << "Realized correlation matrix\n";
		os << indent << "---------------------------\n";
		os << itsRealizedCorrelMatrix->toString(indent, nextIndent );
	}

	return os.str();
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	Class  : ARM_CorrelACPMatParam	////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_CorrelACPMatParam
///	Routine: ARM_CorrelACPMatParam
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_CorrelACPMatParam::ARM_CorrelACPMatParam( 
	ARM_GP_Vector* breakPointTimes,
	ARM_GP_Vector* values,
	const string& interpolatorName,
	ARM_GP_Vector* lowerBound,
	ARM_GP_Vector* upperBound,
	bool adviseBreakPointTimes)
:	
	ARM_CorrelMatParam( ARM_ModelParamType::Correlation, breakPointTimes, values, interpolatorName, lowerBound, upperBound, adviseBreakPointTimes ),
	itsInputCorrelMatrix(NULL),
	itsACPMatrix( NULL ),
	itsFactorCount(-1)
{}


////////////////////////////////////////////////////
///	Class  : ARM_CorrelACPMatParam
///	Routine: ARM_CorrelACPMatParam
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////
ARM_CorrelACPMatParam::ARM_CorrelACPMatParam( const ARM_CorrelACPMatParam& rhs )
:	
	ARM_CorrelMatParam(rhs), 
	itsInputCorrelMatrix( rhs.itsInputCorrelMatrix? (ARM_GP_Matrix*) rhs.itsInputCorrelMatrix->Clone() : NULL ),
	itsACPMatrix( rhs.itsACPMatrix? (ARM_GP_Matrix*) rhs.itsACPMatrix->Clone() : NULL ),
	itsFactorCount( rhs.itsFactorCount )
{}



////////////////////////////////////////////////////
///	Class  : ARM_CorrelACPMatParam
///	Routine: operator=
///	Returns: 
///	Action : assignment operator
////////////////////////////////////////////////////
ARM_CorrelACPMatParam& ARM_CorrelACPMatParam::operator=(const ARM_CorrelACPMatParam& rhs )
{
	if( this != & rhs )
	{
		ARM_CorrelMatParam::operator =(rhs);
		delete itsInputCorrelMatrix;
		delete itsACPMatrix;
		itsInputCorrelMatrix	= rhs.itsInputCorrelMatrix? (ARM_GP_Matrix*) rhs.itsInputCorrelMatrix->Clone() : NULL;
		itsACPMatrix			= rhs.itsACPMatrix? (ARM_GP_Matrix*) rhs.itsACPMatrix->Clone() : NULL;
		itsFactorCount			= rhs.itsFactorCount;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_CorrelACPMatParam
///	Routine: Destructor
///	Returns: 
////////////////////////////////////////////////////
ARM_CorrelACPMatParam::~ARM_CorrelACPMatParam()
{
	delete itsInputCorrelMatrix;
	delete itsACPMatrix;
}


////////////////////////////////////////////////////
///	Class  : ARM_CorrelACPMatParam
///	Routine: UpdateMultiCurves
///	Returns: ARM_MultiCurve*
///	Action : computes the output matrix
////////////////////////////////////////////////////
void ARM_CorrelACPMatParam::SetInputCorrelMatrix( ARM_GP_Matrix* inputCorrelMatrix )
{
	delete itsInputCorrelMatrix;
	itsInputCorrelMatrix = inputCorrelMatrix;
}


////////////////////////////////////////////////////
///	Class  : ARM_CorrelACPMatParam
///	Routine: UpdateMultiCurves
///	Returns: ARM_MultiCurve*
///	Action : computes the output matrix
////////////////////////////////////////////////////
void ARM_CorrelACPMatParam::UpdateMultiCurves( ARM_GP_Vector* breakPointTimes )
{
	size_t size = itsInputCorrelMatrix->rows();

#if defined __GP_STRICT_VALIDATION
	if( itsFactorCount == -1 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "unitialised factor count!" );
	if( !breakPointTimes )
		ARM_THROW( ERR_INVALID_ARGUMENT, "!breakPointTimes" );
	if( size != breakPointTimes->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, " size != breakPointTimes->size()" );
#endif

	size_t i,j;
	ARM_GP_Vector eigenValues(size);
	ARM_GP_Matrix eigenVectors(itsFactorCount,size,0.0);
	itsACPMatrix = ACPTransformation( itsInputCorrelMatrix, eigenValues);

	/// get all the curveFactors and renormalize them!
	for(i=0; i<itsFactorCount; ++i)
	{
		if( fabs(eigenValues[i])> K_NEW_DOUBLE_TOL )
		{
			if(eigenValues[i] < -K_NEW_DOUBLE_TOL)	
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					" there is a nonpositive Eignvalue ...!" );

			double eigenSqrt = sqrt(eigenValues[i]);
			for(j=0; j<size; ++j) 
				eigenVectors(i,j) = (*itsACPMatrix)(j,i)*eigenSqrt;
		}
	}

	/// Create the multi Curve
#if defined __GP_STRICT_VALIDATION
	if( !GetMultiCurve()->GetInterpolator() )
		ARM_THROW( ERR_INVALID_ARGUMENT, " no interpolator set!" );
#endif

	ARM_MultiCurveInterpolator* interpolator = (ARM_MultiCurveInterpolator*) GetMultiCurve()->GetInterpolator()->Clone();

	ARM_GP_T_Vector< ARM_GP_Vector > ordinates(breakPointTimes->size() );
	ARM_MultiCurve* mCurve = new ARM_MultiCurve( *breakPointTimes, ordinates, interpolator );
	for(i=0;i<size; ++i)
	{
        ARM_GP_Vector factor(itsFactorCount);
		for(j=0;j<itsFactorCount;++j)        
			factor[j] = eigenVectors(j,i);
        ordinates[i] = factor;
    } 
	mCurve->SetOrdinates(ordinates);
	SetMultiCurve(mCurve);
}


////////////////////////////////////////////////////
///	Class  : ARM_CorrelACPMatParam
///	Routine: PrintMatrix
///	Returns: print the matrix in the ostringstream object
///	Action : void
////////////////////////////////////////////////////
void ARM_CorrelACPMatParam::PrintMatrix( CC_Ostringstream& os, const string& matrixName, ARM_GP_Matrix* matrix, const string& indent, const string& nextIndent ) const
{
	os << matrixName;
	if( matrix )
		os << matrix->toString( indent, nextIndent ) << "\n";
	else
		os << "\n";
}


////////////////////////////////////////////////////
///	Class  : ARM_CorrelACPMatParam
///	Routine: toString
///	Returns: 
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_CorrelACPMatParam::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
    os << "\n\n=========> ACP CorrelMatrix <========: \n"; 
	PrintMatrix( os, " Input  Matrix ", itsInputCorrelMatrix, indent, nextIndent );
	PrintMatrix( os, " ACP    Matrix ", itsACPMatrix, indent, nextIndent );
	os << " Factor Count " << itsFactorCount << "\n\n";
	os << ARM_CorrelMatParam::toString( indent, nextIndent );
	return os.str();
}



////////////////////////////////////////////////////
///	Class  : ARM_CorrelACPMatParam
///	Routine: SetFactorCount
///	Returns: void
///	Action : set factor count
////////////////////////////////////////////////////
void ARM_CorrelACPMatParam::SetFactorCount( size_t factorCount )
{
	itsFactorCount=factorCount;
}



////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	Class  : ARM_CorrelTrigoMatParam	////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_CorrelTrigoMatParam
///	Routine: ARM_CorrelTrigoMatParam
///	Returns: 
///	Action : constructor
////////////////////////////////////////////////////

ARM_CorrelTrigoMatParam::ARM_CorrelTrigoMatParam( 
	double theta,
	const ARM_DateStrip& datestrip,
	double julianAsOfDate,
	const string& interpolatorName,
	ARM_GP_Vector* lowerBound,
	ARM_GP_Vector* upperBound,
	bool adviseBreakPointTimes )
:	
	ARM_CorrelACPMatParam( &ARM_GP_Vector(1,0.0), &ARM_GP_Vector(1,theta), interpolatorName, lowerBound, upperBound ), 
	itsTheta(theta),
	itsCorrelBreakPointTimes( NULL )
{
	itsCorrelBreakPointTimes = (ARM_GP_Vector*) datestrip.GetResetDates()->Clone();
	*itsCorrelBreakPointTimes -= julianAsOfDate;

	size_t n = itsCorrelBreakPointTimes->size();
	ARM_GP_Matrix* inputMatrix = TrigoMatrix(n,theta);
	SetInputCorrelMatrix( inputMatrix );
}


ARM_CorrelTrigoMatParam::ARM_CorrelTrigoMatParam( 
	double theta,
	const ARM_GP_Vector& resetDates,
	double julianAsOfDate,
	const string& interpolatorName,
	ARM_GP_Vector* lowerBound,
	ARM_GP_Vector* upperBound,
	bool adviseBreakPointTimes )
:	
	ARM_CorrelACPMatParam( &ARM_GP_Vector(1,0.0), &ARM_GP_Vector(1,theta), interpolatorName, lowerBound, upperBound ), 
	itsTheta(theta),
	itsCorrelBreakPointTimes( NULL )
{
	itsCorrelBreakPointTimes = (ARM_GP_Vector*) resetDates.Clone();
	*itsCorrelBreakPointTimes -= julianAsOfDate;

	size_t n = itsCorrelBreakPointTimes->size();
	ARM_GP_Matrix* inputMatrix = TrigoMatrix(n,theta);
	SetInputCorrelMatrix( inputMatrix );
}

////////////////////////////////////////////////////
///	Class  : ARM_CorrelTrigoMatParam
///	Routine: ARM_CorrelTrigoMatParam
///	Returns: 
///	Action : copy constructor
////////////////////////////////////////////////////
ARM_CorrelTrigoMatParam::ARM_CorrelTrigoMatParam( const ARM_CorrelTrigoMatParam& rhs )
:	ARM_CorrelACPMatParam(rhs),
	itsTheta(rhs.itsTheta),
	itsCorrelBreakPointTimes( rhs.itsCorrelBreakPointTimes? (ARM_GP_Vector*) rhs.itsCorrelBreakPointTimes->Clone(): NULL )
{}



////////////////////////////////////////////////////
///	Class  : ARM_CorrelTrigoMatParam
///	Routine: operator=
///	Returns: ARM_CorrelTrigoMatParam&
///	Action : operator=
////////////////////////////////////////////////////

ARM_CorrelTrigoMatParam& ARM_CorrelTrigoMatParam::operator=(const ARM_CorrelTrigoMatParam& rhs )
{
	if( this != & rhs )
	{
		ARM_CorrelACPMatParam::operator =(rhs);
		itsTheta = rhs.itsTheta;
		delete itsCorrelBreakPointTimes;
		itsCorrelBreakPointTimes = rhs.itsCorrelBreakPointTimes? (ARM_GP_Vector*) rhs.itsCorrelBreakPointTimes->Clone(): NULL;
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_CorrelTrigoMatParam
///	Routine: ~ARM_CorrelTrigoMatParam
///	Returns: 
///	Action : destructor
////////////////////////////////////////////////////

ARM_CorrelTrigoMatParam::~ARM_CorrelTrigoMatParam()
{
	delete itsCorrelBreakPointTimes;
}


////////////////////////////////////////////////////
///	Class  : ARM_CorrelTrigoMatParam
///	Routine: SetFactorCount
///	Returns: 
///	Action : set the nb of factors
////////////////////////////////////////////////////
void ARM_CorrelTrigoMatParam::SetFactorCount( size_t factorCount )
{
	ARM_CorrelACPMatParam::SetFactorCount(factorCount);
	UpdateMultiCurves( itsCorrelBreakPointTimes );
}



////////////////////////////////////////////////////
///	Class  : ARM_CorrelTrigoMatParam
///	Routine: MergeModelParam
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_CorrelTrigoMatParam::MergeModelParam( ARM_ModelParam* NewValue )
{
	if( !dynamic_cast<ARM_CorrelTrigoMatParam*>(NewValue) )
		ARM_THROW( ERR_INVALID_ARGUMENT, " cannot merge a CorrelTrigoMatParam with something different than a CorrelTrigoMatParam!" );

	ARM_CurveModelParam::MergeModelParam( NewValue );
	ARM_GP_Vector& values = GetCurve()->GetOrdinates();

#if defined __GP_STRICT_VALIDATION
	if( values.size() != 1 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "values->size() != 1" );
	if( !itsCorrelBreakPointTimes )
		ARM_THROW( ERR_INVALID_ARGUMENT, " !itsCorrelBreakPointTimes" );
#endif

	itsTheta = values[0]; 
	size_t n = itsCorrelBreakPointTimes->size();
	ARM_GP_Matrix* inputMatrix = TrigoMatrix(n,itsTheta);
	SetInputCorrelMatrix( inputMatrix );
	UpdateMultiCurves( itsCorrelBreakPointTimes );
}


///////////////////////////////////////////////////
///	Class  : ARM_CorrelTrigoMatParam
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////

string ARM_CorrelTrigoMatParam::toString( const string& indent, const string& nextIndent ) const
{
#if defined __GP_STRICT_VALIDATION
	if( !itsCorrelBreakPointTimes )
		ARM_THROW( ERR_INVALID_ARGUMENT, "!itsCorrelBreakPointTimes" );
#endif

	CC_Ostringstream os;
	os << "\n Correl Trigo Matrix \n\n";
	os << " CorrelBreakPointTimes : " << itsCorrelBreakPointTimes->toString(indent,nextIndent) << "\n";
	os << " Theta                 : " << itsTheta << "\n";
	os << ARM_CorrelACPMatParam::toString( indent, nextIndent );
	return os.str();
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	Class  : ARM_CorrelMatNoCalibParam ////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_CorrelMatNoCalibParam
///	Routine: ARM_CorrelMatNoCalibParam::
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_CorrelMatNoCalibParam::ARM_CorrelMatNoCalibParam( 
	ARM_GP_Vector* breakPointTimes,
	ARM_GP_Vector* values,
	const string& interpolatorName,
	ARM_GP_Vector* lowerBound,
	ARM_GP_Vector* upperBound,
	bool adviseBreakPointTimes)
:	
	ARM_CorrelACPMatParam( breakPointTimes, breakPointTimes, interpolatorName, lowerBound, upperBound, adviseBreakPointTimes )
{
	if( !values )
		ARM_THROW( ERR_INVALID_ARGUMENT, "!values " );
	if( !breakPointTimes )
		ARM_THROW( ERR_INVALID_ARGUMENT, "!breakPointTimes " );
	if( values->size() != breakPointTimes->size() * breakPointTimes->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "values->size() != breakPointTimes->size() * breakPointTimes->size()" );

	ARM_GP_Matrix* correlMatrix = new ARM_GP_Matrix( breakPointTimes->size(), breakPointTimes->size(), *values );
	SetInputCorrelMatrix( correlMatrix );
}


////////////////////////////////////////////////////
///	Class  : ARM_CorrelMatNoCalibParam
///	Routine: ARM_CorrelMatNoCalibParam::
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////

ARM_CorrelMatNoCalibParam::ARM_CorrelMatNoCalibParam( const ARM_CorrelMatNoCalibParam& rhs )
:	ARM_CorrelACPMatParam(rhs){}


////////////////////////////////////////////////////
///	Class  : ARM_CorrelMatNoCalibParam
///	Routine: operator=
///	Returns: ARM_CorrelMatNoCalibParam&
///	Action : assignment operator
////////////////////////////////////////////////////
ARM_CorrelMatNoCalibParam& ARM_CorrelMatNoCalibParam::operator=(const ARM_CorrelMatNoCalibParam& rhs )
{
	if( this!= &rhs )
		ARM_CorrelACPMatParam::operator =(rhs);
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_CorrelMatNoCalibParam
///	Routine: SetFactorCount
///	Returns: 
///	Action : set the nb of factors
////////////////////////////////////////////////////
void ARM_CorrelMatNoCalibParam::SetFactorCount( size_t factorCount )
{
	ARM_CorrelACPMatParam::SetFactorCount(factorCount);
	UpdateMultiCurves( &GetCurve()->GetAbscisses() );
}


////////////////////////////////////////////////////
///	Class  : ARM_CorrelMatNoCalibParam
///	Routine: ComputeMultiCurves
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_CorrelMatNoCalibParam::MergeModelParam( ARM_ModelParam* NewValue )
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ": cannot this correlation matrix with another model param!" );
}



////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	Class  : ARM_BrownianCorrelationParam ////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_BrownianCorrelationParam
///	Routine: ARM_BrownianCorrelationParam
///	Returns: 
///	Action : constructor
////////////////////////////////////////////////////

ARM_BrownianCorrelationParam::ARM_BrownianCorrelationParam( 
	ARM_GP_Vector* breakPointTimes,
	ARM_GP_Vector* values,
	const string& interpolatorName,
	ARM_GP_Vector* lowerBound,
	ARM_GP_Vector* upperBound,
	bool adviseBreakPointTimes )
:	ARM_CorrelMatParam( ARM_ModelParamType::BrownianCorrelation, breakPointTimes, values, interpolatorName, lowerBound, upperBound, adviseBreakPointTimes )
{
	ComputeMultiCurves( values );
}



////////////////////////////////////////////////////
///	Class  : ARM_BrownianCorrelationParam
///	Routine: ARM_BrownianCorrelationParam
///	Returns: 
///	Action : copy constructor
////////////////////////////////////////////////////
ARM_BrownianCorrelationParam::ARM_BrownianCorrelationParam( const ARM_BrownianCorrelationParam& rhs)
:	ARM_CorrelMatParam(rhs) {}



////////////////////////////////////////////////////
///	Class  : ARM_BrownianCorrelationParam
///	Routine: operator=
///	Returns: 
////////////////////////////////////////////////////

ARM_BrownianCorrelationParam& ARM_BrownianCorrelationParam::operator=( const ARM_BrownianCorrelationParam& rhs )
{
	if( this != &rhs )
		ARM_CorrelMatParam::operator =( rhs);
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_BrownianCorrelationParam
///	Routine: ARM_BrownianCorrelationParam
///	Returns: 
///	Action : destructor
////////////////////////////////////////////////////
ARM_BrownianCorrelationParam::~ARM_BrownianCorrelationParam()
{}


////////////////////////////////////////////////////
///	Class  : ARM_BrownianCorrelationParam
///	Routine: SetFactorCount
///	Returns: 
///	Action : set the nb of factors
////////////////////////////////////////////////////
void ARM_BrownianCorrelationParam::SetFactorCount( size_t factorCount )
{
	if( factorCount != 2 )
		ARM_THROW( ERR_INVALID_ARGUMENT, " factorCount != 2" );
	ARM_GP_Vector& values = GetCurve()->GetOrdinates();
	ComputeMultiCurves(&values);
}



////////////////////////////////////////////////////
///	Class  : ARM_BrownianCorrelationParam
///	Routine: MergeModelParam
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_BrownianCorrelationParam::MergeModelParam( ARM_ModelParam* NewValue )
{
	if( !dynamic_cast<ARM_BrownianCorrelationParam*>(NewValue) )
		ARM_THROW( ERR_INVALID_ARGUMENT, " cannot merge a CorrelVecTrigoMatParam with something different than a CorrelVecTrigoMatParam !" );
	ARM_CurveModelParam::MergeModelParam( NewValue );

	ComputeMultiCurves( &GetCurve()->GetOrdinates() );
}

////////////////////////////////////////////////////
///	Class  : ARM_BrownianCorrelationParam
///	Routine: ComputeMultiCurves
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_BrownianCorrelationParam::ComputeMultiCurves( ARM_GP_Vector* values )
{
	const ARM_GP_Vector& times	= GetCurve()->GetAbscisses();
	size_t size	= values->size();

#if defined __GP_STRICT_VALIDATION
	if( times.size() != values->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, " newValues->size() != values->size()" );
#endif

	ARM_MultiCurveInterpolator* interpolator = (ARM_MultiCurveInterpolator*) GetMultiCurve()->GetInterpolator()->Clone();
	ARM_GP_T_Vector< ARM_GP_Vector > ordinates(times.size() );
	ARM_MultiCurve* mCurve = new ARM_MultiCurve( times, ordinates, interpolator );
	
	for( size_t i=0; i<size; ++i )
	{
        ARM_GP_Vector factors(2);
		factors[0] = cos(ARM_NumericConstants::ARM_PI*(*values)[i]);
		factors[1] = sin(ARM_NumericConstants::ARM_PI*(*values)[i]);
        ordinates[i] = factors;
    }
	mCurve->SetOrdinates( ordinates );
	SetMultiCurve(mCurve);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

