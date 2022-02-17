/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file compositegen.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */


#include "gpnumlib/compositegen.h"
#include "expt.h"

/// gpbase
#include "gpbase/ostringstream.h"
/// gpnumlib
#include "gpnumlib/boxmuller.h"
#include "gpnumlib/invcum.h"
#include "gpnumlib/condinvcum.h"
#include "gpnumlib/transposer.h"
#include "gpnumlib/skipper.h"
#include "gpnumlib/mixtegen.h"
#include "gpnumlib/uniform.h"
#include "gpnumlib/stddevfunc.h"

/// standard library
#include <cmath>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///			 ARM_NormalCompositeGen			 ///////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_NormalCompositeGen
///	Routine: constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_NormalCompositeGen::ARM_NormalCompositeGen( 
const ARM_RandomGeneratorPtr& randomGen1,
const ARM_RandomGeneratorPtr& randomGen2,
ARM_DrawMethod method, 
double nbStdDevs,
int firstNbTimes,
int firstNbDims,
int order,
int firstSimulations )
{
	switch(method)
	{
	case BoxMuller_Method:
		itsFunction = ARM_ManipulatorMethodPtr ( new ARM_BoxMuller(randomGen1) );
		break;
	case InvCumDistribution:
		itsFunction = ARM_ManipulatorMethodPtr ( new ARM_InvCumFunction(randomGen1, ARM_InvCumFunction::INV_ERF_MORRO ) );
		break;
	case InvCumDistributionFast:
		itsFunction = ARM_ManipulatorMethodPtr ( new ARM_InvCumFunction(randomGen1, ARM_InvCumFunction::INV_ERF ) );
		break;
	case CondInvCumDistribution:
		itsFunction = ARM_ManipulatorMethodPtr ( new ARM_CondInvCumFunction(randomGen1, ARM_InvCumFunction::INV_ERF_MORRO, nbStdDevs ) );
		break;
	case Transposer:
		itsFunction = ARM_ManipulatorMethodPtr ( new ARM_Transposer(randomGen1, (ARM_Transposer::TransposeOrder) order, firstSimulations) );
		break;
	case Skipper:
		itsFunction = ARM_ManipulatorMethodPtr ( new ARM_SkipperFunction(randomGen1,nbStdDevs) );
		break;
	case MixteGen:
		itsFunction = ARM_ManipulatorMethodPtr ( new ARM_MixteGen(randomGen1, randomGen2, firstNbTimes, firstNbDims) );
		break;
	default:
	    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Unknown method type!");
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_NormalCompositeGen  
///	Routine: reset
///	Returns: void
///	Action : reset the random number generator
////////////////////////////////////////////////////*
void ARM_NormalCompositeGen::reset(size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb)
{
	itsFunction->reset(dim, nbOfPointsList, factorNb);
}

void ARM_NormalCompositeGen::reset( const ARM_GP_T_Vector<size_t>& nbOfPointsList, const ARM_GP_T_Vector<size_t>& factorNb )
{
	itsFunction->reset(nbOfPointsList, factorNb );
}

////////////////////////////////////////////////////
///	Class  : ARM_NormalCompositeGen
///	Routine: constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
double ARM_NormalCompositeGen::DrawOne()
{
	return (*itsFunction)();
}



////////////////////////////////////////////////////
///	Class  : ARM_NormalCompositeGen
///	Routine: toString
///	Returns: 
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_NormalCompositeGen::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " Normal Composite generator " << CC_NS(std,endl);
	os << itsFunction->toString();
	return os.str();
}




/////////////////////////////////////////////////////////////////
///	Class  : ARM_NormalCompositeGen
///	Routine: copy constructor, assignment operator
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_NormalCompositeGen::ARM_NormalCompositeGen( const ARM_NormalCompositeGen& rhs )
:	ARM_RandomGenerator( rhs ), 
	itsFunction( rhs.itsFunction->Clone() )
{
}


ARM_NormalCompositeGen& ARM_NormalCompositeGen::operator =(const ARM_NormalCompositeGen& rhs )
{
	if( this != & rhs )
	{
		ARM_RandomGenerator::operator =( rhs );
		itsFunction		= ARM_ManipulatorMethodPtr( rhs.itsFunction->Clone() );
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_NormalCompositeGen
///	Routine: desctructor (necessary even for pure virtual case)
///	Action : draw a vector of uniform random numbers
////////////////////////////////////////////////////
ARM_NormalCompositeGen::~ARM_NormalCompositeGen()
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_NormalCompositeGen
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_NormalCompositeGen::Clone() const
{
	return new ARM_NormalCompositeGen(*this);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_NormalCompositeGen
///	Routine: ARM_MomentFuncPtr
///	Returns: StdDevComputationFunc
///	Action : gives the functor to compute the std dev
/////////////////////////////////////////////////////////////////
ARM_MomentFuncPtr ARM_NormalCompositeGen::StdDevComputationFunc() const
{
	return itsFunction->GetBaseGen()->StdDevComputationFunc(); 
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

