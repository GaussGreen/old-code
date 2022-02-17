/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file invcum.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */


#include "gpnumlib/invcum.h"
#include "gpbase/ostringstream.h"
#include "gpbase/cloneutilityfunc.h"
#include "expt.h"
#include "gpnumlib/normalinvcum.h"
#include "gpnumlib/uniform.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///				ARM_InvCumFunction			////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_InvCumFunction
///	Routine: constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_InvCumFunction::ARM_InvCumFunction(const ARM_RandomGeneratorPtr& randomGen, invCumAlgoType algoType )
: ARM_ManipulatorMethod(randomGen)
{
#if defined(__GP_STRICT_VALIDATION)
	if(randomGen == ARM_RandomGeneratorPtr(NULL) )
	    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "pointor on random number generator is NULL!");

	if( !dynamic_cast<ARM_UniformGenerator*>(&*randomGen) 
		&& randomGen->GetDistributionType() != ARM_RandomGenerator::ARM_Uniform )
	    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "random number base generator is not uniform!!");

#endif
	switch(algoType)
	{
	case INV_ERF_MORRO:
		itsFunc = &ARM_NormalInvCum::Inverse_erf_Moro;
		break;
	case INV_ERF:
		itsFunc = &ARM_NormalInvCum::Inverse_erf;
		break;
	default:
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": unknown inverse cumulative normal algorithm!" );

	}
}

////////////////////////////////////////////////////
///	Class  : ARM_InvCumFunction  
///	Routine: reset
///	Returns: void
///	Action : reset the random number generator
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_InvCumFunction
///	Routine: operator()() const
///	Returns: 
///	Action : creates a normal from a uniform
////////////////////////////////////////////////////
double ARM_InvCumFunction::operator()() const
{
	return (*itsFunc)( (*itsRandomGen)() );
}


string ARM_InvCumFunction::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << " Normal Generator  : " << "Inverse Cumulative function" << CC_NS(std,endl);
	os << " Base Generator    : " << itsRandomGen->toString();
	return os.str();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_InvCumFunction
///	Routine: copy constructor, assignment operator
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_InvCumFunction::ARM_InvCumFunction( const ARM_InvCumFunction& rhs )
:	ARM_ManipulatorMethod( CreateClonedPtr(&*rhs.itsRandomGen) ), /// clone it for independence!
	itsFunc(rhs.itsFunc)
{
	
}


////////////////////////////////////////////////////
///	Class  : ARM_InvCumFunction
///	Routine: desctructor (necessary even for pure virtual case)
///	Action : draw a vector of uniform random numbers
////////////////////////////////////////////////////
ARM_InvCumFunction::~ARM_InvCumFunction()
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_InvCumFunction
///	Routine: Clone
///	Returns: ARM_ManipulatorMethod*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_ManipulatorMethod* ARM_InvCumFunction::Clone() const
{
	return new ARM_InvCumFunction(*this);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

