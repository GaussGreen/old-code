/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParamsQ1F.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */


/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"

#include "gpmodels/ModelParamsQNF.h"

/// gpbase
#include "gpbase/numericconstant.h"
#include "gpbase/curve.h"
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/cloneutilityfunc.h"

/// gpinfra
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"

/// gpnumlib
#include "gpnumlib/normalinvcum.h"
#include "gpnumlib/gaussiananalytics.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQNF
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsQNF::ARM_ModelParamsQNF( const ARM_CurveModelParam& QParam, const vector<ARM_ModelParamsHW1F*>& factorParams, const ARM_GP_Matrix& correlMatrix )
:	ARM_ModelParams(), itsFactorParams(factorParams.size(), NULL ), itsQParam(QParam), itsCorrelMatrix(correlMatrix)
{	
	DuplicateCloneablePointorAndNullVectorInPlace<ARM_ModelParamsHW1F>( factorParams, itsFactorParams );
	ValidateModelParams();
}



////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsQNF
///	Routines: Validate
///	Returns :
///	Action  : validate the model params to check that this is compatible with the Q1F model
////////////////////////////////////////////////////
void ARM_ModelParamsQNF::ValidateModelParams() const
{	
    /// checks that the q model is of size 1 since the current model is with cst q!
	if( ARM_ModelParamType::QParameter != itsQParam.GetType() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": ARM_ModelParamsQNF should have a Q model parameter!");

	/// check only a single Q
	if( itsQParam.size() != 1)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": ARM_ModelParamsQNF should have a single Q model parameter!");

	/// check the matrix
	if( itsCorrelMatrix.rows() != itsFactorParams.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": ARM_ModelParamsQNF: CorrelMatParam.rows() != Nb Of factors!");
	
	if( itsCorrelMatrix.cols() != itsFactorParams.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": ARM_ModelParamsQNF: CorrelMatParam.cols() != Nb Of factors!");

	/// Check correlation number
	size_t i,j;
	for( i=0; i<itsCorrelMatrix.cols(); ++i )
		if( fabs(itsCorrelMatrix(i,i)-1.0) > K_NEW_DOUBLE_TOL )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": diagonal != 1.0 !");
	
	for( i=0; i<itsCorrelMatrix.cols(); ++i )
	{
		for( j=0; j<i; ++j )
		{
			if( fabs(itsCorrelMatrix(i,j)-itsCorrelMatrix(j,i)) > K_NEW_DOUBLE_TOL )
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not symmetric!");

// FIXMEFRED: mig.vc8 (25/05/2007 15:30:01):take care
			if( abs(itsCorrelMatrix(i,j) > 1.0 ) )
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": correlation term >= 1.0!");
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQNF
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsQNF::ARM_ModelParamsQNF( const ARM_ModelParamsQNF& rhs )
:	ARM_ModelParams(rhs), itsFactorParams( rhs.itsFactorParams.size(), NULL ), itsQParam( rhs.itsQParam ), itsCorrelMatrix(rhs.itsCorrelMatrix)
{
	DuplicateCloneablePointorAndNullVectorInPlace<ARM_ModelParamsHW1F>( rhs.itsFactorParams, itsFactorParams );
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQNF
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsQNF::~ARM_ModelParamsQNF()
{
    DeletePointorVector<ARM_ModelParamsHW1F>( itsFactorParams );
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQNF
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelParamsQNF& ARM_ModelParamsQNF::operator=(const ARM_ModelParamsQNF& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParams::operator=(rhs);
	    DeletePointorVector<ARM_ModelParamsHW1F>( itsFactorParams );
		DuplicateCloneablePointorAndNullVectorInPlace<ARM_ModelParamsHW1F>( rhs.itsFactorParams, itsFactorParams );
		itsQParam			= rhs.itsQParam;
		itsCorrelMatrix		= rhs.itsCorrelMatrix;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsQNF
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ModelParamsQNF::Clone() const
{
	return new ARM_ModelParamsQNF(*this);
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsQNF
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_ModelParamsQNF::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "ARM_ModelParamsQNF\n";
    os << "----------------------\n";
	os << " Dimension : " << itsFactorParams.size() << "\n";
	os << " Q Param   :\n" << itsQParam.toString() << "\n";
	for( size_t i=0; i<itsFactorParams.size(); ++i )
		os << i+1 << ") Factor:\n " << itsFactorParams[i]->toString() << "\n";
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsQNF
///	Routines: GetModelParam
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
const ARM_ModelParamsHW1F* ARM_ModelParamsQNF::GetModelParamsPerDim( size_t dim ) const
{
#ifdef __GP_STRICT_VALIDATION
	if( dim >= itsFactorParams.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": ARM_ModelParamsQNF out of bound!");
#endif
	return itsFactorParams[dim];
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

