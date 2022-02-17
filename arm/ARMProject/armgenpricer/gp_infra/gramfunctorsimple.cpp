/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gramfunctor.cpp
 *
 *  \brief gramfunctor are functor corresponding to grammar functions
 *	the functor design allows to have a context and a unified command like
 *	interface
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#include <iostream>
#include <iomanip>
#include <fstream>



#include "gpbase/env.h"

#include "gpinfra/gramfunctorsimple.h"

/// gpbase
#include "gpbase/eventviewerfwd.h"
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/checkinputs.h"
#include "gpbase/ostringstream.h"
#include "gpbase/gpvector.h"

/// gpinfra
#include "gpinfra/gramfunctorarg.h"
#include "gpinfra/gramfunctorargcheck.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/functorop.h"
#include "gpinfra/argconvdefault.h"
#include "gpinfra/gramfunctorconv.h"
#include "gpinfra/gramnode.h"
#include "gpinfra/exerciseboundary.h"

/// for easy debugging of shared nodes!
#if defined(__GP_SHOW_SHARED_NODE_COORDINATES)
#include "gpbase/pair.h"
#endif


CC_BEGIN_NAMESPACE( ARM )

/// the corresponding static string
string ARM_GP_PowVector::itsFuncName= "Pow function";
string ARM_GP_IfVector::itsFuncName	= "If function";
string ARM_GP_DCF::itsFuncName		= "DCF (day count fraction)";
string ARM_GP_PV::itsFuncName		= "PV";
string ARM_GP_UnPay::itsFuncName	= "UnPay";
string ARM_GP_Exercise::itsFuncName = "Exercise";
string ARM_GP_Frontier::itsFuncName = "Frontier";
string ARM_GP_Trigger::itsFuncName  = "Trigger";
string ARM_GP_SumSerieVector::itsFuncName = "Sum Serie 1/(1+x)";
/// Stat
string								ARM_GP_StatAverageVector::itsFuncName				= "Statistics Average";
std::vector< std::vector<double> >	ARM_GP_StatAverageVector::initialValuesPerStates_	= std::vector< std::vector<double> >();

string								ARM_GP_StatStdDevVector::itsFuncName				= "Statistics Standard Deviation";
std::vector< std::vector<double> >	ARM_GP_StatStdDevVector::initialValuesPerStates_	= std::vector< std::vector<double> >();

////////////////////////////////////////////////////
/// \function: If in place: 
/// \param v1 first vector
/// \param v2 second vector
/// \param v3 third vector
///
/// \brief
//		- if the first vector has a size of 1, 
///			affect to v1 v2 if true, v3 otherwise
///		- else, v1,v2,v3 should have the same size
///			and foreach elem do v1[i]=if(v1[i]) v2[i]; else v3[i];
////////////////////////////////////////////////////

void IfInPlace( ARM_GP_VectorPtr& v1, const ARM_VectorPtr& v2, const ARM_VectorPtr& v3 )
{
	std::vector<double>::iterator v1Elem		= (*v1).begin();
	std::vector<double>::const_iterator v2Elem	= (*v2).begin();
	std::vector<double>::const_iterator v3Elem	= (*v3).begin();
	std::vector<double>::iterator end			= (*v1).end();
	
	for( ; v1Elem != end; ++v1Elem, ++v2Elem, ++v3Elem )
		*v1Elem = *v1Elem? *v2Elem : *v3Elem;
}



////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: operator()
///	Returns: ARM_GramFctorArg 
///	Action : uset itsOp to tansfor the initial inputs
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_UnaryOp::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
											double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking if we are in a fwdbckwd nummetho or not. 
	if( ARM_ExpNode::ModelDoesNotSupportInPlace( mod ) )
		return OperatorWithClone(arg, mod, evalDate, states, nodes );
	else
		return OperatorInPlace(arg, mod, evalDate, states, nodes );
}



////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: OperatorInPlace
///	Returns: ARM_GramFctorArg 
///	Action : uset itsOp to tansfor the initial inputs
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_UnaryOp::OperatorInPlace( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
											double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 1, itsFuncName );
	GPAF_CheckArgType(  arg[0], GFAT_VECTOR_TYPE, itsFuncName );

	/// use function InPlace for efficiency reason
	if( GFAT_DOUBLE_TYPE == arg[0].GetType() )
		arg[0].SetDouble( itsOp( arg[0].GetDouble() ) );
	else
		FuncUnaryInPlace( arg[0].GetVector(), itsOp );
	return arg[0];
}


////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: OperatorWithClone
///	Returns: ARM_GramFctorArg 
///	Action : uset itsOp to tansfor the initial inputs
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_UnaryOp::OperatorWithClone( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
											double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 1, itsFuncName );
	GPAF_CheckArgType(  arg[0], GFAT_VECTOR_TYPE, itsFuncName );

	if( GFAT_DOUBLE_TYPE == arg[0].GetType() )
		return ARM_GramFctorArg( itsOp( arg[0].GetDouble() ) );
	else
	{
		ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );//( new std::vector<double>( *(arg[0].GetVector()->Clone()) );
		FuncUnaryInPlace( newvec, itsOp );
		return ARM_GramFctorArg( newvec );
	}
	return arg[0];
}

////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: operator()
///	Returns: ARM_GramFctorArg 
///	Action : uset itsOp to tansfor the initial inputs
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_PowVector::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
											  double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	if( ARM_ExpNode::ModelDoesNotSupportInPlace( mod ) )
		return OperatorWithClone(arg, mod, evalDate, states, nodes );
	else
		return OperatorInPlace(arg, mod, evalDate, states, nodes );

}

////////////////////////////////////////////////////
///	Class  : ARM_GP_PowVector
///	Routine: operatorInPlace
///	Returns: ARM_GramFctorArg 
///	Action : uset itsOp to tansfor the initial inputs
////////////////////////////////////////////////////

/// Power operator the second argument has to be a double!
ARM_GramFctorArg ARM_GP_PowVector::OperatorInPlace( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
											  double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 2, ARM_GP_PowVector::itsFuncName );
	GPAF_CheckArgType(  arg[0], GFAT_VECTOR_TYPE, ARM_GP_PowVector::itsFuncName );
	
	if( GFAT_DOUBLE_TYPE != arg[1].GetType() )
	{
		CC_Ostringstream os;
		os << "Pow function takes only exponent as simple deterministic double! " << ARM_USERNAME 
			<< " please advise!";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}


	if( GFAT_DOUBLE_TYPE == arg[0].GetType() )
		arg[0].SetDouble( pDbleBinaryFunc( pow )( arg[0].GetDouble(), arg[1].GetDouble() ) );
	else
		FuncUnaryInPlace( arg[0].GetVector(), CC_NS( std, bind2nd )( pDbleBinaryFunc( pow ), ( arg[1].GetDouble() ) ) );

	return arg[0];
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_PowVector
///	Routine: operatorWithClone
///	Returns: ARM_GramFctorArg 
///	Action : uset itsOp to tansfor the initial inputs
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_PowVector::OperatorWithClone( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
											  double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 2, ARM_GP_PowVector::itsFuncName );
	GPAF_CheckArgType(  arg[0], GFAT_VECTOR_TYPE, ARM_GP_PowVector::itsFuncName );

	if( GFAT_DOUBLE_TYPE != arg[1].GetType() )
	{
		CC_Ostringstream os;
		os << "Pow function takes only exponent as simple deterministic double! " << ARM_USERNAME 
			<< " please advise!";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}


	if( GFAT_DOUBLE_TYPE == arg[0].GetType() )
		return ARM_GramFctorArg( pDbleBinaryFunc( pow )( arg[0].GetDouble(), arg[1].GetDouble() ) );
	else
	{
		ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );//( new std::vector<double>( *(arg[0].GetVector() ) ) );
		FuncUnaryInPlace( newvec, CC_NS( std, bind2nd )( pDbleBinaryFunc( pow ), ( arg[1].GetDouble() ) ) );
		return ARM_GramFctorArg( newvec );
	}
}



/// If operator
ARM_GramFctorArg ARM_GP_IfVector::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
											 double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 3, ARM_GP_IfVector::itsFuncName );
	GPAF_CheckArgType(  arg[0], GFAT_VECTOR_TYPE, ARM_GP_IfVector::itsFuncName );
	GPAF_CheckArgType(  arg[1], GFAT_DATE_VECTOR_OR_DOUBLE_TYPE, ARM_GP_IfVector::itsFuncName );
	GPAF_CheckArgType(  arg[2], GFAT_DATE_VECTOR_OR_DOUBLE_TYPE, ARM_GP_IfVector::itsFuncName );

	if( ARM_ExpNode::ModelDoesNotSupportInPlace( mod ) )
		return OperatorWithClone(arg, mod, evalDate, states, nodes );
	else
		return OperatorInPlace(arg, mod, evalDate, states, nodes );

}

ARM_GramFctorArg ARM_GP_IfVector::OperatorWithClone( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
											 double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// handle all the cases
	if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
		&& arg[1].GetType() == GFAT_VECTOR_TYPE
		&& arg[2].GetType() == GFAT_VECTOR_TYPE )
	{	
		ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );//( ARM_GP_VectorPtr( new std::vector<double>( *(arg[0].GetVector()->Clone()) ) );
		GPAF_CheckArgVecSameSize( arg, ARM_GP_IfVector::itsFuncName );
		IfInPlace( newvec, arg[1].GetVector(), arg[2].GetVector() );
		return ARM_GramFctorArg( newvec );
	}
	else
		if( arg[0].GetType() == GFAT_DOUBLE_TYPE )
		{
			GPAF_CheckArgVecSameSize( arg, ARM_GP_IfVector::itsFuncName );
			return ARM_GramFctorArg( arg[0].GetDouble()? arg[1] : arg[2] );
		}
		else
			if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
				&& arg[1].GetType() == GFAT_DOUBLE_TYPE
				&& arg[2].GetType() == GFAT_VECTOR_TYPE )
			{
				GPAF_CheckArgVecSameSize( arg, ARM_GP_IfVector::itsFuncName );
				ARM_VectorPtr vec( new std::vector<double>( arg[0].GetVector()->size(), arg[1].GetDouble() ) );
				arg[1].SetVector( vec );
				ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );//( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );
				IfInPlace( newvec, arg[1].GetVector(), arg[2].GetVector() );
				return ARM_GramFctorArg( newvec );
			}
			else
				if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
					&& arg[1].GetType() == GFAT_VECTOR_TYPE
					&& arg[2].GetType() == GFAT_DOUBLE_TYPE )
				{
					GPAF_CheckArgVecSameSize( arg, ARM_GP_IfVector::itsFuncName );
					ARM_VectorPtr vec( new std::vector<double>( arg[0].GetVector()->size(), arg[2].GetDouble() ) );
					arg[2].SetVector( vec );
					ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );//( new std::vector<double>( *(arg[0].GetVector()->Clone()) );
					IfInPlace( newvec, arg[1].GetVector(), arg[2].GetVector() );
					return ARM_GramFctorArg( newvec );
				}
				else
					if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
						&& arg[1].GetType() == GFAT_DOUBLE_TYPE
						&& arg[2].GetType() == GFAT_DOUBLE_TYPE )
					{
						ARM_VectorPtr vec1( new std::vector<double>( arg[0].GetVector()->size(), arg[1].GetDouble() ) );
						ARM_VectorPtr vec2( new std::vector<double>( arg[0].GetVector()->size(), arg[2].GetDouble() ) );
						arg[1].SetVector( vec1 );
						arg[2].SetVector( vec2 );
						ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );//( new std::vector<double>( *(arg[0].GetVector()->Clone()) );
						IfInPlace( newvec, arg[1].GetVector(), arg[2].GetVector() );
						return ARM_GramFctorArg( newvec );
					}
					////////////////////
					//// date section
					////////////////////
					else
						if(	arg[0].GetType() == GFAT_VECTOR_TYPE 
							&& arg[1].GetType() == GFAT_DATE_TYPE
							&& arg[2].GetType() == GFAT_VECTOR_TYPE )
						{
							GPAF_CheckArgVecSameSize( arg, ARM_GP_IfVector::itsFuncName );
							ARM_VectorPtr vec( new std::vector<double>( arg[0].GetVector()->size(), arg[1].GetDate().GetJulian() ) );
							arg[1].SetVector( vec );
							ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );//( new std::vector<double>( *(arg[0].GetVector()->Clone()) );
							IfInPlace( newvec, arg[1].GetVector(), arg[2].GetVector() );
							return ARM_GramFctorArg( newvec );
						}
						else
							if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
								&& arg[1].GetType() == GFAT_VECTOR_TYPE
								&& arg[2].GetType() == GFAT_DATE_TYPE )
							{
								GPAF_CheckArgVecSameSize( arg, ARM_GP_IfVector::itsFuncName );
								ARM_VectorPtr vec( new std::vector<double>( arg[0].GetVector()->size(), arg[2].GetDate().GetJulian() ) );
								arg[2].SetVector( vec );
								ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );//( new std::vector<double>( *(arg[0].GetVector()->Clone()) );
								IfInPlace( newvec, arg[1].GetVector(), arg[2].GetVector() );
								return ARM_GramFctorArg( newvec );
							}
							else
								if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
									&& arg[1].GetType() == GFAT_DATE_TYPE
									&& arg[2].GetType() == GFAT_DATE_TYPE )
								{
									ARM_VectorPtr vec1( new std::vector<double>( arg[0].GetVector()->size(), arg[1].GetDate().GetJulian() ) );
									ARM_VectorPtr vec2( new std::vector<double>( arg[0].GetVector()->size(), arg[2].GetDate().GetJulian() ) );
									arg[1].SetVector( vec1 );
									arg[2].SetVector( vec2);
									ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );//( new std::vector<double>( *(arg[0].GetVector()->Clone()) );
									IfInPlace( newvec, arg[1].GetVector(), arg[2].GetVector() );
									return ARM_GramFctorArg( newvec );
								}
								else 
									throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": the condition in the if cannot be a date!" );
}

ARM_GramFctorArg ARM_GP_IfVector::OperatorInPlace( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
											 double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// handle all the cases
	if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
		&& arg[1].GetType() == GFAT_VECTOR_TYPE
		&& arg[2].GetType() == GFAT_VECTOR_TYPE )
	{	
		IfInPlace( arg[0].GetVector(), arg[1].GetVector(), arg[2].GetVector() );
	}
	else
		if( arg[0].GetType() == GFAT_DOUBLE_TYPE )
		{
			GPAF_CheckArgVecSameSize( arg, ARM_GP_IfVector::itsFuncName );
			arg[0] = arg[0].GetDouble()? arg[1] : arg[2];
		}
		else
			if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
				&& arg[1].GetType() == GFAT_DOUBLE_TYPE
				&& arg[2].GetType() == GFAT_VECTOR_TYPE )
			{
				GPAF_CheckArgVecSameSize( arg, ARM_GP_IfVector::itsFuncName );
				ARM_VectorPtr vec( new std::vector<double>( arg[0].GetVector()->size(), arg[1].GetDouble() ) );
				arg[1].SetVector( vec );
				IfInPlace( arg[0].GetVector(), arg[1].GetVector(), arg[2].GetVector() );
			}
			else
				if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
					&& arg[1].GetType() == GFAT_VECTOR_TYPE
					&& arg[2].GetType() == GFAT_DOUBLE_TYPE )
				{
					GPAF_CheckArgVecSameSize( arg, ARM_GP_IfVector::itsFuncName );
					ARM_VectorPtr vec( new std::vector<double>( arg[0].GetVector()->size(), arg[2].GetDouble() ) );
					arg[2].SetVector( vec );
					IfInPlace( arg[0].GetVector(), arg[1].GetVector(), arg[2].GetVector() );
				}
				else
					if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
						&& arg[1].GetType() == GFAT_DOUBLE_TYPE
						&& arg[2].GetType() == GFAT_DOUBLE_TYPE )
					{
						ARM_VectorPtr vec1( new std::vector<double>( arg[0].GetVector()->size(), arg[1].GetDouble() ) );
						ARM_VectorPtr vec2( new std::vector<double>( arg[0].GetVector()->size(), arg[2].GetDouble() ) );
						arg[1].SetVector( vec1 );
						arg[2].SetVector( vec2 );
						IfInPlace( arg[0].GetVector(), arg[1].GetVector(), arg[2].GetVector() );
					}
					////////////////////
					//// date section
					////////////////////
					else
						if(	arg[0].GetType() == GFAT_VECTOR_TYPE 
							&& arg[1].GetType() == GFAT_DATE_TYPE
							&& arg[2].GetType() == GFAT_VECTOR_TYPE )
						{
							GPAF_CheckArgVecSameSize( arg, ARM_GP_IfVector::itsFuncName );
							ARM_VectorPtr vec( new std::vector<double>( arg[0].GetVector()->size(), arg[1].GetDate().GetJulian() ) );
							arg[1].SetVector( vec );
							IfInPlace( arg[0].GetVector(), arg[1].GetVector(), arg[2].GetVector() );
						}
						else
							if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
								&& arg[1].GetType() == GFAT_VECTOR_TYPE
								&& arg[2].GetType() == GFAT_DATE_TYPE )
							{
								GPAF_CheckArgVecSameSize( arg, ARM_GP_IfVector::itsFuncName );
								ARM_VectorPtr vec( new std::vector<double>( arg[0].GetVector()->size(), arg[2].GetDate().GetJulian() ) );
								arg[2].SetVector( vec );
								IfInPlace( arg[0].GetVector(), arg[1].GetVector(), arg[2].GetVector() );
							}
							else
								if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
									&& arg[1].GetType() == GFAT_DATE_TYPE
									&& arg[2].GetType() == GFAT_DATE_TYPE )
								{
									ARM_VectorPtr vec1( new std::vector<double>( arg[0].GetVector()->size(), arg[1].GetDate().GetJulian() ) );
									ARM_VectorPtr vec2( new std::vector<double>( arg[0].GetVector()->size(), arg[2].GetDate().GetJulian() ) );
									arg[1].SetVector( vec1 );
									arg[2].SetVector( vec2);
									IfInPlace( arg[0].GetVector(), arg[1].GetVector(), arg[2].GetVector() );
								}
								else 
									throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": the condition in the if cannot be a date!" );

	return arg[0];
}






/// binary operator like minimum and maximum!

ARM_GramFctorArg ARM_GP_BinaryOpWithDates::OperatorInPlace( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
													  double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{

	/// handle all the cases
	if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
		&& arg[1].GetType() == GFAT_VECTOR_TYPE )
	{
		GPAF_CheckArgVecSameSize( arg, itsFuncName );
		FuncBinaryInPlace(arg[0].GetVector(),arg[1].GetVector(), itsOp );	
	}
	else
		if(    arg[0].GetType() == GFAT_DOUBLE_TYPE 
			&& arg[1].GetType() == GFAT_VECTOR_TYPE )
		{
			FuncUnaryInPlace( arg[1].GetVector(), CC_NS( std, bind1st )( itsOp, arg[0].GetDouble() ) );
			arg[0]=arg[1];
		}
		else
			if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
				&& arg[1].GetType() == GFAT_DOUBLE_TYPE )
			{
				FuncUnaryInPlace( arg[0].GetVector(), CC_NS( std, bind2nd )( itsOp, arg[1].GetDouble() ) );
			}
			else
				if(    arg[0].GetType() == GFAT_DOUBLE_TYPE 
					&& arg[1].GetType() == GFAT_DOUBLE_TYPE )
				{
					arg[0].SetDouble( ( itsOp( arg[0].GetDouble(), arg[1].GetDouble() ) ) );
				}
				////////////////////
				//// date section
				////////////////////
				else
					if(    arg[0].GetType() == GFAT_DATE_TYPE 
						&& arg[1].GetType() == GFAT_VECTOR_TYPE )
					{
						FuncUnaryInPlace( arg[1].GetVector(), CC_NS( std, bind1st )( itsOp, arg[0].GetDate().GetJulian() ) );
						arg[0]=arg[1];
					}
					else
						if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
							&& arg[1].GetType() == GFAT_DATE_TYPE )
							FuncUnaryInPlace( arg[0].GetVector(), CC_NS( std, bind2nd )( itsOp, arg[1].GetDate().GetJulian() ) );
						else
							if(    arg[0].GetType() == GFAT_DATE_TYPE 
								&& arg[1].GetType() == GFAT_DATE_TYPE )
								arg[0].SetDate( ARM_Date( itsOp( arg[0].GetDate().GetJulian(), arg[1].GetDate().GetJulian() ) ) );
							else 
								throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": permitted types are dates, double, vectors!" );
							
	return arg[0];
}

ARM_GramFctorArg ARM_GP_BinaryOpWithDates::OperatorWithClone( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
													  double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// handle all the cases
	if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
		&& arg[1].GetType() == GFAT_VECTOR_TYPE )
	{
		GPAF_CheckArgVecSameSize( arg, itsFuncName );
		ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );//( new std::vector<double>( *(arg[0].GetVector()->Clone()) );
		FuncBinaryInPlace(newvec,arg[1].GetVector(), itsOp );
		return ARM_GramFctorArg( newvec );
	}
	else
		if(    arg[0].GetType() == GFAT_DOUBLE_TYPE 
			&& arg[1].GetType() == GFAT_VECTOR_TYPE )
		{
			ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );//( new std::vector<double>( *(arg[1].GetVector()->Clone()) );
			FuncUnaryInPlace( newvec, CC_NS( std, bind1st )( itsOp, arg[0].GetDouble() ) );
			return ARM_GramFctorArg( newvec );
		}
		else
			if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
				&& arg[1].GetType() == GFAT_DOUBLE_TYPE )
			{
				ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );//( new std::vector<double>( *(arg[0].GetVector()->Clone()) );
				FuncUnaryInPlace( newvec, CC_NS( std, bind2nd )( itsOp, arg[1].GetDouble() ) );
				return ARM_GramFctorArg( newvec );
			}
			else
				if(    arg[0].GetType() == GFAT_DOUBLE_TYPE 
					&& arg[1].GetType() == GFAT_DOUBLE_TYPE )
					return ARM_GramFctorArg( ( itsOp( arg[0].GetDouble(), arg[1].GetDouble() ) ) );

				////////////////////
				//// date section
				////////////////////
				else
					if(    arg[0].GetType() == GFAT_DATE_TYPE 
						&& arg[1].GetType() == GFAT_VECTOR_TYPE )
					{
						ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );//( new std::vector<double>( *(arg[1].GetVector()->Clone()) );
						FuncUnaryInPlace( newvec, CC_NS( std, bind1st )( itsOp, arg[0].GetDate().GetJulian() ) );
						return ARM_GramFctorArg( newvec );
					}
					else
						if(	   arg[0].GetType() == GFAT_VECTOR_TYPE 
							&& arg[1].GetType() == GFAT_DATE_TYPE )
						{
							ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );//( new std::vector<double>( *(arg[0].GetVector()->Clone()) );
							FuncUnaryInPlace( newvec, CC_NS( std, bind2nd )( itsOp, arg[1].GetDate().GetJulian() ) );
							return ARM_GramFctorArg( newvec );
						}
						else
							if(    arg[0].GetType() == GFAT_DATE_TYPE 
								&& arg[1].GetType() == GFAT_DATE_TYPE )
								return ARM_GramFctorArg( ARM_Date( itsOp( arg[0].GetDate().GetJulian(), arg[1].GetDate().GetJulian() ) ) );
							else 
								throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": permitted types are dates, double, vectors!" );
							
}


ARM_GramFctorArg ARM_GP_BinaryOpWithDates::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
													  double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 2, itsFuncName );
	GPAF_CheckArgType(  arg[0], GFAT_DATE_VECTOR_OR_DOUBLE_TYPE, itsFuncName );
	GPAF_CheckArgType(  arg[1], GFAT_DATE_VECTOR_OR_DOUBLE_TYPE, itsFuncName );

	if( ARM_ExpNode::ModelDoesNotSupportInPlace( mod ) )
		return OperatorWithClone(arg, mod, evalDate, states, nodes );
	else
		return OperatorInPlace(arg, mod, evalDate, states, nodes );

}


/// Day Count Fraction
ARM_GramFctorArg ARM_GP_DCF::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
										double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 3, ARM_GP_DCF::itsFuncName );
	GPAF_CheckArgType(  arg[0], GFAT_DATE_OR_VECTOR_TYPE,	ARM_GP_DCF::itsFuncName );
	GPAF_CheckArgType(  arg[1], GFAT_DATE_OR_VECTOR_TYPE,	ARM_GP_DCF::itsFuncName );
	GPAF_CheckArgType(  arg[2], GFAT_STRING_TYPE,			ARM_GP_DCF::itsFuncName );

	/// checking if we are in a fwdbckwd nummethod or not. 
	bool doesSupportInPlace = ARM_ExpNode::ModelDoesNotSupportInPlace( mod );
;
	
	int dayCount = ARM_ArgConv_DayCount.GetNumber( arg[2].GetString() );
	if(	   arg[0].GetType() == GFAT_DATE_TYPE 
		&& arg[1].GetType() == GFAT_DATE_TYPE )
	{	
		double dcf = CountYearsWithoutException( dayCount, arg[0].GetDate(), arg[1].GetDate() );
		return ARM_GramFctorArg( dcf );
	}
	else
		if(	   arg[0].GetType() == GFAT_DATE_TYPE 
			&& arg[1].GetType() == GFAT_VECTOR_TYPE )
		{
			double dateJulian	= arg[0].GetDate().GetJulian();
			ARM_GP_VectorPtr v1;
			if( doesSupportInPlace )
				v1 = ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(arg[1].GetVector()->Clone()) );
			else
				v1 = arg[1].GetVector();

			for( size_t i=0; i<v1->size(); ++i )
				(*v1)[i] = CountYearsWithoutException( dayCount, dateJulian, (*v1)[i] );
			return ARM_GramFctorArg( v1 );
		}
		else
			if(    arg[0].GetType() == GFAT_VECTOR_TYPE
				&& arg[1].GetType() == GFAT_DATE_TYPE )
			{
				ARM_VectorPtr v0;
				if( doesSupportInPlace )
					v0 = ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );
				else
					v0 = arg[0].GetVector();

				double dateJulian	= arg[1].GetDate().GetJulian();
				for( size_t i=0; i<v0->size(); ++i )
					(*v0)[i] = CountYearsWithoutException( dayCount, (*v0)[i], dateJulian );
				return ARM_GramFctorArg( v0 );
			}
			else 
				if(    arg[0].GetType() == GFAT_VECTOR_TYPE
					&& arg[1].GetType() == GFAT_VECTOR_TYPE ) 
					
				{
					GPAF_CheckArgVecSameSize( arg, ARM_GP_DCF::itsFuncName );
					ARM_VectorPtr v0;
					if( doesSupportInPlace )
						v0 = ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );
					else
						v0 = arg[0].GetVector();
					ARM_VectorPtr v1 = arg[1].GetVector();
					for( size_t i=0; i<v0->size(); ++i )
						(*v0)[i] = CountYearsWithoutException( dayCount, (*v0)[i], (*v1)[i]);
					return ARM_GramFctorArg( v0 );
				}
				else 
					if(	   arg[0].GetType() == GFAT_DOUBLE_TYPE 
						&& arg[1].GetType() == GFAT_DOUBLE_TYPE )
					{	
						double dcf = CountYearsWithoutException( dayCount, arg[0].GetDouble(), arg[1].GetDouble() );
						return ARM_GramFctorArg( dcf );
					}
					else
						if(	   arg[0].GetType() == GFAT_DOUBLE_TYPE 
							&& arg[1].GetType() == GFAT_VECTOR_TYPE )
						{
							double dateJulian	= arg[0].GetDouble();
							ARM_GP_VectorPtr v1;
							if( doesSupportInPlace )
								v1 = ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(arg[1].GetVector()->Clone()) );
							else
								v1 = arg[1].GetVector();
							for( size_t i=0; i<v1->size(); ++i )
								(*v1)[i] = CountYearsWithoutException( dayCount, dateJulian, (*v1)[i] );
							return ARM_GramFctorArg( v1 );
						}
						else
							if(    arg[0].GetType() == GFAT_VECTOR_TYPE
								&& arg[1].GetType() == GFAT_DOUBLE_TYPE )
							{
								ARM_GP_VectorPtr v0;
								if( doesSupportInPlace )
									v0 = ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );
								else
									v0 = arg[0].GetVector();
								double dateJulian	= arg[1].GetDouble();
								for( size_t i=0; i<v0->size(); ++i )
									(*v0)[i] = CountYearsWithoutException( dayCount, (*v0)[i], dateJulian );
								return ARM_GramFctorArg( v0 );
							}
							else
								if(	   arg[0].GetType() == GFAT_DOUBLE_TYPE 
									&& arg[1].GetType() == GFAT_DATE_TYPE )
								{
									double dcf = CountYearsWithoutException( dayCount, arg[0].GetDouble(), arg[1].GetDate().GetJulian() );
									return ARM_GramFctorArg( dcf );
								}
								else
									if(    arg[0].GetType() == GFAT_DATE_TYPE
										&& arg[1].GetType() == GFAT_DOUBLE_TYPE )
									{
										double dcf = CountYearsWithoutException( dayCount, arg[0].GetDate().GetJulian(), arg[1].GetDouble() );
										return ARM_GramFctorArg( dcf );
									}
									else
										throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": dates can be dates|vector|double!" );
}


/// PV operator
ARM_GramFctorArg ARM_GP_PV::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
									   double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// PV operator is very special as it is handled payoff paid at one date and
	/// evaluated at a different date.
	GPAF_CheckArgSize( arg, 1, "PV" );
	GPAF_CheckArgType(  arg[0], GFAT_VECTOR_TYPE, "PV" );
	
	double evalTime = GetAndCheckEvalTime( evalDate, mod, ARM_GP_PV::itsFuncName );
	
#if defined( __GP_STRICT_VALIDATION)
	if( dynamic_cast<ARM_ExpNodeRef*>( &*nodes[0] ) == NULL )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": this should be a node reference!" );
	
	ARM_ExpNodeRef& nodeRef = (ARM_ExpNodeRef&) *nodes[0];
	double parentTime		= mod->GetTimeFromDate( nodeRef.GetParentDate() );
	
	if( evalTime != parentTime )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": a pv should only be executed at parentTime!" );
#endif
	
	ARM_GP_VectorPtr payoff;
	
	if( ARM_ExpNode::ModelDoesNotSupportInPlace( mod ) )
		payoff = ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );
	else
		payoff = arg[0].GetVector();
	
	/// unpaid the payoff
	mod->ProcessUnPaidPayoffs( GetPayModelName(), payoff, evalTime, states );
	
	return ARM_GramFctorArg(payoff);
}




/// PV operator
ARM_GramFctorArg ARM_GP_UnPay::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
										  double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// PV operator is very special as it is handled payoff paid at one date and
	/// evaluated at a different date.
	GPAF_CheckArgSize( arg, 1, "PV" );
	GPAF_CheckArgType(  arg[0], GFAT_VECTOR_TYPE, "PV" );
	
	double evalTime = GetAndCheckEvalTime( evalDate, mod, ARM_GP_UnPay::itsFuncName );
	ARM_VectorPtr payoff;
	
	switch( arg[0].GetType() )
	{
	case GFAT_VECTOR_TYPE:
		{
			if( ARM_ExpNode::ModelDoesNotSupportInPlace( mod ) )
				payoff = ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );
			else 
				payoff = arg[0].GetVector();
			/// unpaid the payoff
			mod->ProcessUnPaidPayoffs( GetPayModelName(), payoff, evalTime, states );
			break;
		}
	case GFAT_DOUBLE_TYPE:
		{
			payoff = ARM_VectorPtr( new std::vector<double>( states->size(), arg[0].GetDouble() ) );
			/// unpaid the payoff
			mod->ProcessUnPaidPayoffs( GetPayModelName(), payoff, evalTime, states );
			/// SetTheValue
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
	/// other case do nothing!
	return ARM_GramFctorArg(payoff);
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_Exercise
///	Routine: ARM_GP_Exercise
///	Returns: 
///	Action : default constructor
///////////////////////////////////////////////
ARM_GP_Exercise::ARM_GP_Exercise() : ARM_GramFctor()
{
	itsExerciseBoundary = NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_Exercise
///	Routine: ARM_GP_Exercise
///	Returns: 
///	Action : copy constructor
///////////////////////////////////////////////
ARM_GP_Exercise::ARM_GP_Exercise(const ARM_GP_Exercise& rhs) : ARM_GramFctor(rhs) 
{
	if ( rhs.itsExerciseBoundary )
		itsExerciseBoundary = static_cast<ARM_ExerciseBoundary *>( itsExerciseBoundary->Clone() );
	else
		itsExerciseBoundary = NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_Exercise
///	Routine: ARM_GP_Exercise
///	Returns: 
///	Action : default destructor
///////////////////////////////////////////////
ARM_GP_Exercise::~ARM_GP_Exercise()
{
	if ( itsExerciseBoundary ) 
		delete itsExerciseBoundary;
	itsExerciseBoundary = NULL;
}


////////////////////////////////////////////////////
///	Routine: ComputeCombin
///	Returns: 
///	Action : Compute all the monome of n variables
/// with the degree at least equal to m
///////////////////////////////////////////////

vector<vector<double> > ComputeCombin(int n, int m)
{
	vector<double> init(n,0);
	vector<double> cur;
	vector<double> dest;
	vector<vector<double> > v;
	vector<vector<double> > prevv;
	vector<vector<double> > newv;

	v.push_back(init);
	prevv.push_back(init);

	int i, j, k, start;

	for (i = 0; i < m; ++i)
	{
		newv.resize(0);
		for (j = 0; j < prevv.size(); ++j)
		{
			cur = prevv[j];
			start = n-1;
			while ((cur[start] == 0) && (start>0))
				start--;

			for (k = start; k < n; ++k)
			{
				dest = cur;
				dest[k]++;
				newv.push_back(dest);
			}
		}
		for (j = 0; j < newv.size(); ++j)
			v.push_back(newv[j]);
		prevv = newv;
	}

	return v;
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_Exercise
///	Routine: operator()
///	Returns: ARM_GramFctorArg 
///	Action : Execise operator
///////////////////////////////////////////////
ARM_GramFctorArg ARM_GP_Exercise::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
									   double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// PV operator is very special as it is handled payoff paid at one date and
	/// evaluated at a different date.
	
	ARM_NumMethod::GP_PricingDirection pricingDir = mod->GetNumMethod()->GetPricingDirCurrLoop();

	if ( mod->GetNumMethod()->GetPricingDirection() != ARM_NumMethod::GP_FWDBCKWDLOOKING )
		if (pricingDir == ARM_NumMethod::GP_FWDLOOKING)
		{
			GPAF_CheckArgSize( arg, 2, ARM_GP_Exercise::itsFuncName );
		}
		else if (pricingDir == ARM_NumMethod::GP_BCKWDLOOKING)
		{
			GPAF_CheckArgSize( arg, 3, ARM_GP_Exercise::itsFuncName );
		}

	GPAF_CheckArgType(  arg[0], GFAT_VECTOR_TYPE, ARM_GP_Exercise::itsFuncName );
	GPAF_CheckArgType(  arg[1], GFAT_VECTOR_TYPE, ARM_GP_Exercise::itsFuncName );
	if (pricingDir == ARM_NumMethod::GP_BCKWDLOOKING)
	{
		GPAF_CheckArgType(  arg[2], GFAT_VECTOR_TYPE, ARM_GP_Exercise::itsFuncName );
	}
	
	double evalTime = GetAndCheckEvalTime( evalDate, mod, ARM_GP_Exercise::itsFuncName );

#if defined( __GP_STRICT_VALIDATION)
/**** no more valid because 2nd argument may be any node and not only a reference one
	if( dynamic_cast<ARM_ExpNodeRef*>( &*nodes[1] ) == NULL )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": this should be a node reference!" );
****/
	
	ARM_ExpNodeRef* nodeRef=NULL;
    if( (nodeRef = dynamic_cast<ARM_ExpNodeRef*>(&*nodes[2])) )
    {
       double parentTime = mod->GetTimeFromDate( nodeRef->GetParentDate() );
	    
	    if( evalTime != parentTime )
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": a pv should only be executed at parentTime!" );
    }

#endif

	ARM_GP_VectorPtr payoff1,payoff2;

	if (arg[0].GetType() == GFAT_VECTOR_TYPE)
	{
		payoff1 = ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );
	}
	else if (arg[0].GetType() == GFAT_DOUBLE_TYPE)
	{
		payoff1 = ARM_GP_VectorPtr( new ARM_GP_Vector( states->size(), arg[0].GetDouble() ) );
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": execise payoff 1 option should be double or vector!" );
	}


	if (arg[1].GetType() == GFAT_VECTOR_TYPE)
	{
		payoff2 = ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(arg[1].GetVector()->Clone()) );
	}
	else if (arg[1].GetType() == GFAT_DOUBLE_TYPE)
	{
		payoff2 = ARM_GP_VectorPtr( new ARM_GP_Vector( states->size(), arg[1].GetDouble() ) );
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": execise payoff 2 option should be double or vector!" );
	}

	// We compute the continuation option and evaluate the node just in the backward loop
	if (pricingDir == ARM_NumMethod::GP_BCKWDLOOKING)
	{
		ARM_GP_VectorPtr contOpt;

		if (arg[2].GetType() == GFAT_VECTOR_TYPE)
		{
			contOpt = arg[2].GetVector();
		}
		else if (arg[2].GetType() == GFAT_DOUBLE_TYPE)
		{
			contOpt = ARM_GP_VectorPtr( new ARM_GP_Vector( states->size(), arg[2].GetDouble() ) );
		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": execise continuation option should be double or vector!" );
		}

		/// unpaid the contination option
		mod->ProcessUnPaidPayoffs( GetPayModelName(), contOpt, evalTime, states );


		///////////////////////////////////////////////////////////////////
		/// begin exercise boundary calculation 
		///////////////////////////////////////////////////////////////////

		if ( itsExerciseBoundary == NULL ) /// A new exercise boundary is computed
			itsExerciseBoundary = mod->GetNumMethod()->ComputeExercise( payoff2, contOpt, itsMat );


		///////////////////////////////////////////////////////////////////
		/// end exercise boundary calculation 
		///////////////////////////////////////////////////////////////////
		if ( itsExerciseBoundary != NULL ) // Which is the case for AMC only
		{
			itsExerciseBoundary->EvalInPlace( payoff2, contOpt, itsMat );
		}
		else  // Which is the case for bckwdlooking nummethod such as trees. 
			//FIXMEFRED: mig.vc8 : CC_Max -> CC_Max_c
			FuncBinaryInPlace( payoff2, contOpt, pDbleBinaryFunc( CC_Max_c<double> ) );

		CC_NS( ARM_GP, plus )<double> plusOp;
		FuncBinaryInPlace( payoff2, payoff1, plusOp );

		return ARM_GramFctorArg(payoff2);
	}
	else
	{
		ARM_ExerciseBoundary* exerciseBoundary = mod->GetNumMethod()->ComputeExercise(ARM_GP_VectorPtr(NULL),ARM_GP_VectorPtr(NULL),ARM_GP_MatrixPtr(NULL));
		if (exerciseBoundary && (exerciseBoundary->IsAutomatic()))
		{
			ARM_GP_MatrixPtr variables;
			ARM_GP_VectorPtr vec;
			if (arg.size() > 3)
			{
				variables = ARM_GP_MatrixPtr(new ARM_GP_Matrix(arg.size()-3, payoff2->size()));
				for ( size_t j=0 ; j<arg.size()-3 ; j++ )
				{
					if (arg[j+3].GetType() == GFAT_VECTOR_TYPE)
					{
						vec = arg[j+3].GetVector(); // FIXME
						for ( size_t i= 0 ; i<payoff2->size() ; i++ )
							(*variables)(j,i) = vec->Elt(i);
					}
					else if (arg[j+3].GetType() == GFAT_DOUBLE_TYPE)
					{
						for ( size_t i= 0 ; i<payoff1->size() ; i++ )
							(*variables)(j,i) = arg[j+3].GetDouble();
					}
					else
					{
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": execise continuation option should be double or vector!" );
					}
				}
			}
			else
			{
				variables = states->GetModelStates();
			}

			int nbStates = variables->cols();
			int n = variables->rows();
			int m = exerciseBoundary->Degree();

			vector<vector<double> > combin = ComputeCombin(n,m);

			int nbVars = combin.size();

			int i, j, k, l;

			itsMat = ARM_GP_MatrixPtr(new ARM_GP_Matrix(nbStates,nbVars));

			for (i = 0; i < nbStates; ++i)
				for (j = 0; j < nbVars; ++j)
					(*itsMat)(i,j) = 1.0;

			for (i = 0; i < nbVars; ++i)
			{
				for (j = 0; j < n; ++j)
				{
					for (k = 0; k < combin[i][j]; ++k)
					{
						for (l = 0; l < nbStates; ++l)
							(*itsMat)(l,i) *= (*variables)(j,l);
					}
				}
			}

			delete exerciseBoundary;
		}
		else
		{
			/// are we in andersen method or with user defined parametrization of the exercise boundary?
			if(	arg.size() > 3 || !mod->GetNumMethod()->NeedToCreateDefaultArgument() )
			{
				/// in case regression functions are specified in the gensec
				ARM_GP_Matrix* mat = new ARM_GP_Matrix( payoff2->size(), arg.size()-3 );
				ARM_GP_VectorPtr vec( NULL );

				for ( size_t j=0 ; j<arg.size()-3 ; j++ )
				{
					if (arg[j+3].GetType() == GFAT_VECTOR_TYPE)
					{
						vec = arg[j+3].GetVector(); // FIXME
						for ( size_t i= 0 ; i<payoff2->size() ; i++ )
							(*mat)(i,j) = vec->Elt(i);
					}
					else if (arg[j+3].GetType() == GFAT_DOUBLE_TYPE)
					{
						for ( size_t i= 0 ; i<payoff1->size() ; i++ )
							(*mat)(i,j) = arg[j+3].GetDouble();
					}
					else
					{
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": execise continuation option should be double or vector!" );
					}
				}

				itsMat = ARM_GP_MatrixPtr( mat );

			}
			else
			{
				/// part for the longstaff schwartz algorithm
				ARM_GP_Matrix * mat = new ARM_GP_Matrix( payoff2->size(), 5, 1.0 );

				for (size_t i = 0 ; i < payoff2->size() ; i++ )
					for (size_t j = 1 ; j<5 ; j++ )
						(*mat)(i,j) = (*mat)(i,j-1) * (*payoff2)[i];

				itsMat = ARM_GP_MatrixPtr( mat );
			}
		}

		return ARM_GramFctorArg(ARM_VectorPtr(new std::vector<double>(states->size(),0.0)));
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_Frontier
///	Routine: ARM_GP_Frontier
///	Returns: 
///	Action : default constructor
///////////////////////////////////////////////
ARM_GP_Frontier::ARM_GP_Frontier() : ARM_GramFctor()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_Frontier
///	Routine: ARM_GP_Frontier
///	Returns: 
///	Action : copy constructor
///////////////////////////////////////////////
ARM_GP_Frontier::ARM_GP_Frontier(const ARM_GP_Frontier& rhs) : ARM_GramFctor(rhs) 
{
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_Frontier
///	Routine: ARM_GP_Frontier
///	Returns: 
///	Action : default destructor
///////////////////////////////////////////////
ARM_GP_Frontier::~ARM_GP_Frontier()
{
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_Frontier
///	Routine: operator()
///	Returns: ARM_GramFctorArg 
///	Action : Frontier operator
///////////////////////////////////////////////
ARM_GramFctorArg ARM_GP_Frontier::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
									   double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	if (mod->GetNumMethod()->SupportFrontierComputation())
	{
		ARM_NumMethod::GP_PricingDirection pricingDir = mod->GetNumMethod()->GetPricingDirCurrLoop();

		GPAF_CheckArgType(  arg[0], GFAT_VECTOR_TYPE, ARM_GP_Frontier::itsFuncName );
		GPAF_CheckArgType(  arg[1], GFAT_VECTOR_TYPE, ARM_GP_Frontier::itsFuncName );
		GPAF_CheckArgType(  arg[2], GFAT_VECTOR_TYPE, ARM_GP_Frontier::itsFuncName );

		double evalTime = GetAndCheckEvalTime( evalDate, mod, ARM_GP_Frontier::itsFuncName );

		ARM_GP_VectorPtr index;

		if (arg[2].GetType() == GFAT_VECTOR_TYPE)
		{
			index = ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(arg[2].GetVector()->Clone()) );
		}
		else if (arg[2].GetType() == GFAT_DOUBLE_TYPE)
		{
			index = ARM_GP_VectorPtr( new ARM_GP_Vector( states->size(), arg[2].GetDouble() ) );
		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": execise payoff 2 option should be double or vector!" );
		}


		// We compute the continuation option and evaluate the node just in the backward loop
		if (pricingDir == ARM_NumMethod::GP_BCKWDLOOKING)
		{
			ARM_GP_VectorPtr intrinsic;

			if (arg[0].GetType() == GFAT_VECTOR_TYPE)
			{
				intrinsic = arg[0].GetVector();
			}
			else if (arg[0].GetType() == GFAT_DOUBLE_TYPE)
			{
				intrinsic = ARM_GP_VectorPtr( new ARM_GP_Vector( states->size(), arg[0].GetDouble() ) );
			}
			else
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": execise continuation option should be double or vector!" );
			}

			ARM_VectorPtr continuation;

			if (arg[1].GetType() == GFAT_VECTOR_TYPE)
			{
				continuation = arg[1].GetVector();
			}
			else if (arg[1].GetType() == GFAT_DOUBLE_TYPE)
			{
				continuation = ARM_GP_VectorPtr( new ARM_GP_Vector( states->size(), arg[1].GetDouble() ) );
			}
			else
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": exercise continuation option should be double or vector!" );
			}

			/// unpaid the contination option
			mod->ProcessUnPaidPayoffs( GetPayModelName(), continuation, evalTime, states );

			int size = intrinsic->size();
			std::vector<double>::iterator iterIdx = index->begin();
			std::vector<double>::iterator iterIV = intrinsic->begin();
			std::vector<double>::iterator iterCV = continuation->begin();
			double sgn = (*iterIV)-(*iterCV);
			double value;
			
			std::vector<double> idxFrontier;
			int i = 0;
			for ( ; iterIV != intrinsic->end(); ++iterIV , ++iterCV, ++iterIdx, ++i)
			{
				if (((*iterIV)-(*iterCV))*sgn < 0 )
				{
					idxFrontier.push_back(i);
					sgn = (*iterIV)-(*iterCV);
				}
			}
			int nbFrontier = idxFrontier.size();
			int minDist = size;
			int idx = size/2;
			if (nbFrontier>0)
			{
				if (nbFrontier>1)
				{
					CC_Ostringstream os;
					os << ARM_USERNAME << " : warning: multi-valued frontier (" << nbFrontier << ")\n";
					ARM_TheEventViewer.Instance()->AddToMessage(os.str());
				}
				for (int j=0;j<nbFrontier;j++)
				{
					int dist = (idxFrontier[j]-size/2)>0?idxFrontier[j]-size/2:size/2-idxFrontier[j];
					if (dist<minDist)
					{
						idx = idxFrontier[j];
						minDist = dist;
					}
				}
				double f1 = (*intrinsic)[idx-1]-(*continuation)[idx-1];
				double f2 = (*intrinsic)[idx]-(*continuation)[idx];
				double df = f2 - f1;
				if (fabs(df)<K_NEW_DOUBLE_TOL)
					value = (*index)[idx];
				else
					value = ( (*index)[idx]*f1 - (*index)[idx-1]*f2 ) / ( f1 - f2 );
			}
			else
				value = index->Elt(size/2);

			for ( iterIdx = index->begin() ; iterIdx != index->end() ; ++iterIdx )
				(*iterIdx) = value;

			/// unpaid the contination option
			mod->ProcessUnPaidPayoffs( GetPayModelName(), index, evalTime, states );

			return ARM_GramFctorArg(index);
		}
		else
		{
			return ARM_GramFctorArg(ARM_VectorPtr(new std::vector<double>(states->size(),0.0)));
		}
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": num method does not support frontier!" );
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_Trigger
///	Routine: ARM_GP_Trigger
///	Returns: 
///	Action : default constructor
///////////////////////////////////////////////
ARM_GP_Trigger::ARM_GP_Trigger() : ARM_GramFctor()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_Trigger
///	Routine: ARM_GP_Trigger
///	Returns: 
///	Action : copy constructor
///////////////////////////////////////////////
ARM_GP_Trigger::ARM_GP_Trigger(const ARM_GP_Trigger& rhs) : ARM_GramFctor(rhs) 
{
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_Trigger
///	Routine: ARM_GP_Trigger
///	Returns: 
///	Action : default destructor
///////////////////////////////////////////////
ARM_GP_Trigger::~ARM_GP_Trigger()
{
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_Trigger
///	Routine: operator()
///	Returns: ARM_GramFctorArg 
///	Action : Execise operator
///////////////////////////////////////////////
ARM_GramFctorArg ARM_GP_Trigger::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
									   double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	size_t i;
	
	ARM_NumMethod::GP_PricingDirection pricingDir = mod->GetNumMethod()->GetPricingDirCurrLoop();


	if (pricingDir == ARM_NumMethod::GP_FWDLOOKING)
		GPAF_CheckArgSize( arg, 4, ARM_GP_Trigger::itsFuncName );

	else if (pricingDir == ARM_NumMethod::GP_BCKWDLOOKING)
		GPAF_CheckArgSize( arg, 5, ARM_GP_Trigger::itsFuncName );


	for ( i=3 ; i<4 ; i++ )
		GPAF_CheckArgType(  arg[i], GFAT_VECTOR_TYPE, ARM_GP_Trigger::itsFuncName );

	if (pricingDir == ARM_NumMethod::GP_BCKWDLOOKING)
		GPAF_CheckArgType(  arg[4], GFAT_VECTOR_TYPE, ARM_GP_Trigger::itsFuncName );

	
	double evalTime = GetAndCheckEvalTime( evalDate, mod, ARM_GP_Trigger::itsFuncName );
	

#if defined( __GP_STRICT_VALIDATION)
	if( dynamic_cast<ARM_ExpNodeRef*>( &*nodes[4] ) == NULL )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": this should be a node reference!" );
	
	ARM_ExpNodeRef& nodeRef = (ARM_ExpNodeRef&) *nodes[4];
	double parentTime		= mod->GetTimeFromDate( nodeRef.GetParentDate() );
	
	if( evalTime != parentTime )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": a pv should only be executed at parentTime!" );
#endif

	if (pricingDir == ARM_NumMethod::GP_BCKWDLOOKING)
	{
		ARM_GP_VectorPtr contOpt;
		ARM_GP_VectorPtr testVar1, payoff1, payoff2;

		//// gets all the vector needed for evaluation

		// testVar1
		if (arg[0].GetType() == GFAT_VECTOR_TYPE)
			testVar1 = arg[0].GetVector();
		else if (arg[0].GetType() == GFAT_DOUBLE_TYPE)
			testVar1 = ARM_GP_VectorPtr( new ARM_GP_Vector( states->size(), arg[0].GetDouble() ) );
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": execise testVar1 option should be double or vector!" );

		/// testVar2
		if (arg[1].GetType() == GFAT_VECTOR_TYPE)
			testVar2 = arg[1].GetVector();
		else if (arg[1].GetType() == GFAT_DOUBLE_TYPE)
			testVar2 = ARM_GP_VectorPtr( new ARM_GP_Vector( states->size(), arg[1].GetDouble() ) );
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": execise testVar2 option should be double or vector!" );

		// payoff1
		if (arg[2].GetType() == GFAT_VECTOR_TYPE)
			payoff1 = arg[2].GetVector();
		else if (arg[2].GetType() == GFAT_DOUBLE_TYPE)
			payoff1 = ARM_GP_VectorPtr( new ARM_GP_Vector( states->size(), arg[2].GetDouble() ) );
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": execise payoff1 option should be double or vector!" );

		// payoff2
		if (arg[3].GetType() == GFAT_VECTOR_TYPE)
			payoff2 = arg[3].GetVector();
		else if (arg[3].GetType() == GFAT_DOUBLE_TYPE)
			payoff2 = ARM_GP_VectorPtr( new ARM_GP_Vector( states->size(), arg[3].GetDouble() ) );
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": execise payoff2 option should be double or vector!" );

		// contopt
		if (arg[4].GetType() == GFAT_VECTOR_TYPE)
			contOpt = arg[4].GetVector();
		else if (arg[4].GetType() == GFAT_DOUBLE_TYPE)
			contOpt = ARM_GP_VectorPtr( new ARM_GP_Vector( states->size(), arg[4].GetDouble() ) );
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": execise continuation option should be double or vector!" );

		ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(payoff2->Clone() ) ) ;

		/// unpaid the contination option
		mod->ProcessUnPaidPayoffs( GetPayModelName(), contOpt, evalTime, states );

		for ( i=0 ; i< contOpt->size() ; i++ )
			payoff2->Elt(i) = payoff1->Elt(i) + (testVar1->Elt(i) >= testVar2->Elt(i) ? payoff2->Elt(i) :  contOpt->Elt(i));

		return payoff2;
	}
	else
	{
		return ARM_GramFctorArg(ARM_VectorPtr(new std::vector<double>(states->size(),0.0)));
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_MFactorFctor
///	Routine: ARM_GP_MFactorFctor
///	Returns: 
///	Action : default constructor
///////////////////////////////////////////////
ARM_GP_MFactorFctor::ARM_GP_MFactorFctor() : ARM_GramFctor()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_MFactorFctor
///	Routine: ARM_GP_MFactorFctor
///	Returns: 
///	Action : copy constructor
///////////////////////////////////////////////
ARM_GP_MFactorFctor::ARM_GP_MFactorFctor(const ARM_GP_MFactorFctor& rhs) : ARM_GramFctor(rhs) 
{
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_MFactorFctor
///	Routine: ~ARM_GP_MFactorFctor
///	Returns: 
///	Action : default destructor
///////////////////////////////////////////////
ARM_GP_MFactorFctor::~ARM_GP_MFactorFctor()
{
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_MFactorFctor
///	Routine: operator()
///	Returns: ARM_GramFctorArg 
///	Action : gets i^th model factor
///////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_MFactorFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
									double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	if (arg[0].GetType() == GFAT_DOUBLE_TYPE)
	{
		int paramIdx = (size_t) (arg[0].GetDouble());
		ARM_GP_MatrixPtr modelStates = states->GetModelStates();
		if ( modelStates->rows() > paramIdx )
			return ARM_GramFctorArg( ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>( modelStates->GetRow( paramIdx )->Clone() ) ) );
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "Bad param for ModelFactor" );

	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "parameters of ModelFactor are supposed to be int!" );
}

ARM_GramFctorArg ARM_GP_SumSerieVector::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
											  double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	if( ARM_ExpNode::ModelDoesNotSupportInPlace( mod ) )
		return OperatorWithClone(arg, mod, evalDate, states, nodes );
	else
		return OperatorInPlace(arg, mod, evalDate, states, nodes );

}

////////////////////////////////////////////////////
///	Class  : ARM_GP_PowVector
///	Routine: operatorInPlace
///	Returns: ARM_GramFctorArg 
///	Action : uset itsOp to tansfor the initial inputs
////////////////////////////////////////////////////

/// Power operator the second argument has to be a double!
ARM_GramFctorArg ARM_GP_SumSerieVector::OperatorInPlace( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
											  double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 2, ARM_GP_SumSerieVector::itsFuncName );
	GPAF_CheckArgType(  arg[0], GFAT_VECTOR_TYPE, ARM_GP_SumSerieVector::itsFuncName );
	
	if( GFAT_DOUBLE_TYPE != arg[1].GetType() )
	{
		CC_Ostringstream os;
		os << "Pow function takes only exponent as simple deterministic double! " << ARM_USERNAME 
			<< " please advise!";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}


	if( GFAT_DOUBLE_TYPE == arg[0].GetType() )
	{
		int k, size = (int)(arg[1].GetDouble());
		double res = 0., x = arg[0].GetDouble();
		for(k = 1; k < size; k++)
		{
			res += 1./pow(1.+x,(double)k);
		}

		arg[0].SetDouble(res);
	}
	else
	{
		int k, size = (int)(arg[1].GetDouble());
		int n, nbStates = arg[0].GetVector()->size();
		for(n = 0; n < nbStates; n++)
		{
			double res = 0., x = (*(arg[0].GetVector()))[n];
			for(k = 1; k < size; k++)
			{
				res += 1./pow(1.+x,(double)k);
			}
			(*(arg[0].GetVector()))[n] = res;
		}
	}

	return arg[0];
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_PowVector
///	Routine: operatorWithClone
///	Returns: ARM_GramFctorArg 
///	Action : uset itsOp to tansfor the initial inputs
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_SumSerieVector::OperatorWithClone( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
											  double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 2, ARM_GP_SumSerieVector::itsFuncName );
	GPAF_CheckArgType(  arg[0], GFAT_VECTOR_TYPE, ARM_GP_SumSerieVector::itsFuncName );

	if( GFAT_DOUBLE_TYPE != arg[1].GetType() )
	{
		CC_Ostringstream os;
		os << "Pow function takes only exponent as simple deterministic double! " << ARM_USERNAME 
			<< " please advise!";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}

	if( GFAT_DOUBLE_TYPE == arg[0].GetType() )
	{
		int k, size = (int)(arg[1].GetDouble());
		double res = 0., x = arg[0].GetDouble();
		for(k = 1; k < size; k++)
		{
			res += 1./pow(1.+x,(double)k);
		}

		return ARM_GramFctorArg(res);
	}
	else
	{
		int k, size = (int)(arg[1].GetDouble());
		int n, nbStates = arg[0].GetVector()->size();
		ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(arg[0].GetVector()->Clone()) );//(new std::vector<double>(nbStates,0.));

		for(n = 0; n < nbStates; n++)
		{
			double x = (*(arg[0].GetVector()))[n];
			for(k = 1; k < size; k++)
			{
				(*newvec)[n] += 1./pow(1.+x,(double)k);
			}
		}

		return ARM_GramFctorArg( newvec );
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_StatAverageVector
///	Routine: operator()
///	Returns: ARM_GramFctorArg 
///	Action : Compute average of a sequence
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_GP_StatAverageVector::operator()(ARM_GramFctorArgVector& arg, 
													  ARM_PricingModel* mod,
													  double evalDate, 
													  const ARM_PricingStatesPtr& states, 
													  vector< ARM_ExpNodePtr >& nodes )
{
	/*
	std::string fileName = "c:/Debug/StatAverageVector_operator.txt";

	std::ofstream os(fileName.c_str());

	if (!os) {
		std::string msg = "Failed to open file : " + fileName;
		os.close();
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg.c_str() );
	}
	*/


	/// checking of the size and type
	GPAF_CheckArgSize( arg, 4, ARM_GP_StatAverageVector::itsFuncName );
	// New Value
	GPAF_CheckArgType(  arg[0], GFAT_VECTOR_TYPE, ARM_GP_StatAverageVector::itsFuncName );
	// Initial Values
	GPAF_CheckArgType(  arg[1], GFAT_VECTOR_TYPE, ARM_GP_StatAverageVector::itsFuncName );
	// Average size
	GPAF_CheckArgType(  arg[2], GFAT_DOUBLE_TYPE, ARM_GP_StatAverageVector::itsFuncName );
	// Reset Initial values
	GPAF_CheckArgType(  arg[3], GFAT_DOUBLE_TYPE, ARM_GP_StatAverageVector::itsFuncName );


	ARM_GP_VectorPtr primePerStates = arg[0].GetVector();
	size_t histoSize = static_cast<size_t>(arg[2].GetDouble());
	isInitialized_   = ( static_cast<int>(arg[3].GetDouble()) == 1? false : true);

	/*
	os << "ARM_GP_StatAverageVector::operator()(...)" << std::endl;
	os << "histo size : " << histoSize << std::endl;
	*/

	if (!isInitialized_)
	{
		ARM_GP_VectorPtr histoPrimes = arg[1].GetVector();
		if ( histoPrimes->size() >= histoSize )
		{
			size_t statesSize = primePerStates->size();

			//os << "states size : " << statesSize << std::endl;

			initialValuesPerStates_ = std::vector< std::vector<double> >(histoSize, std::vector<double>());
			//os << "[";
			for (size_t i = 0; i < histoSize; ++i )
			{
				double value = (*histoPrimes)[(histoSize-1)-i];
				/*
				os.precision(4);
				os.width(10);
				os << value << std::endl;
				*/
				initialValuesPerStates_[i] = std::vector<double>(statesSize, value);
			}
			//os << "]" << std::endl;
		}
		else {
			CC_Ostringstream oss;
			oss << "StatAverage function : vector of historic primes must be greater then " << histoSize << ", " << ARM_USERNAME 
				<< " please advise!";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, oss.str() );
		}
		isInitialized_ = true;
	}
	else {
		if (initialValuesPerStates_.size() == 0) {
			CC_Ostringstream oss;
			oss << "StatAverage function : vector of initial values is not initialized" << ARM_USERNAME 
				<< " please advise!";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, oss.str() );
		}
	}
	
	//os.close();

	bool cloneResult = ARM_ExpNode::ModelDoesNotSupportInPlace( mod );
	return pushNewValuesAndComputeAverage(arg, histoSize, cloneResult);
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_StatAverageVector
///	Routine: pushNewValuesAndComputeAverage
///	Returns: ARM_GramFctorArg 
///	Action : Add new value to sequence
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_GP_StatAverageVector::pushNewValuesAndComputeAverage(ARM_GramFctorArgVector& arg, size_t histoSize, bool cloneResult)
{
	ARM_GP_VectorPtr result(NULL);
	ARM_GP_VectorPtr primePerStates = arg[0].GetVector();
	size_t nbStates					= primePerStates->size();

	if (cloneResult)
		result = ARM_GP_VectorPtr(new ARM_GP_Vector(nbStates,0.));
	else
		result = primePerStates;

	// push new primes per states
	std::vector<double> pps(nbStates,0.);
	size_t i, j;
	for(i = 0; i < pps.size(); ++i)
	{
		pps[i] = (*primePerStates)[i];
	}
	initialValuesPerStates_.erase(initialValuesPerStates_.begin());
	initialValuesPerStates_.push_back(pps);

	for(i = 0; i < pps.size(); ++i)
	{
		double average = 0.;
		for(j = 0; j < histoSize; ++j) 
		{
			average += initialValuesPerStates_[j][i];
		}
		average /= static_cast<double>(histoSize);
		(*result)[i] = average;
	}

	if (cloneResult)
		return ARM_GramFctorArg(result);
	else
		return arg[0];
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_StatStdDevVector
///	Routine: operator()
///	Returns: ARM_GramFctorArg 
///	Action : Compute Standart Deviation of a sequence
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_GP_StatStdDevVector::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
											  double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 5, ARM_GP_StatStdDevVector::itsFuncName );
	// New Value
	GPAF_CheckArgType(  arg[0], GFAT_VECTOR_TYPE, ARM_GP_StatStdDevVector::itsFuncName );
	// Average
	GPAF_CheckArgType(  arg[1], GFAT_VECTOR_TYPE, ARM_GP_StatStdDevVector::itsFuncName );
	// Initial Values
	GPAF_CheckArgType(  arg[2], GFAT_VECTOR_TYPE, ARM_GP_StatStdDevVector::itsFuncName );
	// Average size
	GPAF_CheckArgType(  arg[3], GFAT_DOUBLE_TYPE, ARM_GP_StatStdDevVector::itsFuncName );
		// Reset Initial values
	GPAF_CheckArgType(  arg[4], GFAT_DOUBLE_TYPE, ARM_GP_StatStdDevVector::itsFuncName );
	

	ARM_GP_VectorPtr primePerStates		= arg[0].GetVector();
	size_t histoSize					= static_cast<size_t>(arg[3].GetDouble());
	isInitialized_						= ( static_cast<int>(arg[4].GetDouble()) == 1? false : true);

	if (!isInitialized_)
	{
		ARM_GP_VectorPtr histoPrimes = arg[2].GetVector();
		if ( histoPrimes->size() >= histoSize )
		{
			size_t statesSize = primePerStates->size();

			initialValuesPerStates_ = std::vector< std::vector<double> >(histoSize, std::vector<double>());
			for (size_t i = 0; i < histoSize; ++i )
			{
				double value = (*histoPrimes)[(histoSize-1)-i];
				initialValuesPerStates_[i] = std::vector<double>(statesSize, value);
			}
		}
		else {
			CC_Ostringstream os;
			os << "StatStdDev function : vector of historic primes must be greater then " << histoSize << ", " << ARM_USERNAME 
				<< " please advise!";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
		isInitialized_ = true;
	}
	else {
		if (initialValuesPerStates_.size() == 0) {
			CC_Ostringstream oss;
			oss << "StatAverage function : vector of initial values is not initialized" << ARM_USERNAME 
				<< " please advise!";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, oss.str() );
		}
	}
	
	bool cloneResult = ARM_ExpNode::ModelDoesNotSupportInPlace( mod );
	return pushNewValuesAndComputeStdDev(arg, histoSize, cloneResult);
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_StatStdDevVector
///	Routine: pushNewValuesAndComputeStdDev
///	Returns: ARM_GramFctorArg 
///	Action : Add new value to sequence
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_GP_StatStdDevVector::pushNewValuesAndComputeStdDev(ARM_GramFctorArgVector& arg, size_t histoSize, bool cloneResult)
{
	ARM_GP_VectorPtr result(NULL);

	ARM_GP_VectorPtr primePerStates		= arg[0].GetVector();
	ARM_GP_VectorPtr averagePerStates	= arg[1].GetVector();
	size_t nbStates						= primePerStates->size();

	if (cloneResult)
		result = ARM_GP_VectorPtr(new ARM_GP_Vector(nbStates,0.));
	else
		result = primePerStates;

	// push new primes per states
	std::vector<double> pps(primePerStates->size(),0.);
	size_t i, j;
	for(i = 0; i < pps.size(); ++i)
	{
		pps[i] = (*primePerStates)[i];
	}
	initialValuesPerStates_.erase(initialValuesPerStates_.begin());
	initialValuesPerStates_.push_back(pps);

	for(i = 0; i < pps.size(); ++i)
	{
		double average = (*averagePerStates)[i],
			stdDev = 0.;
		for(j = 0; j < histoSize; ++j) 
		{
			double value = initialValuesPerStates_[j][i];
			stdDev += value*value;
		}
		stdDev /= static_cast<double>(histoSize);
		stdDev = sqrt(stdDev-average*average);
		(*result)[i] = stdDev;
	}

	if (cloneResult)
		return ARM_GramFctorArg(result);
	else
		return arg[0];
}





CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

