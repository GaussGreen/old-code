/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file manipulator.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#include "gpnumlib/manipulator.h"
#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpnumlib/random.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ManipulatorMethod
///	Routine: desctructor (necessary even for pure virtual case)
///	Action : method to transform uniform random gen into other
///			 distribution
////////////////////////////////////////////////////
ARM_ManipulatorMethod::ARM_ManipulatorMethod(const ARM_RandomGeneratorPtr& randomGen) : itsRandomGen(randomGen)
{
}

ARM_ManipulatorMethod::ARM_ManipulatorMethod(const ARM_ManipulatorMethod& rhs)
{
	itsRandomGen = ARM_RandomGeneratorPtr( (ARM_RandomGenerator*) rhs.itsRandomGen->Clone() );
}

ARM_ManipulatorMethod::~ARM_ManipulatorMethod()
{}

void ARM_ManipulatorMethod::reset(size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb)
{
	ARM_GP_T_Vector<size_t> factorsNb;
	if(factorNb > 0) factorsNb.resize(dim/factorNb,factorNb);
	reset(nbOfPointsList, factorsNb);
}

void ARM_ManipulatorMethod::reset(const ARM_GP_T_Vector<size_t>& nbOfPointsList, const ARM_GP_T_Vector<size_t>& factorNb)
{
	if(itsRandomGen != ARM_RandomGeneratorPtr(NULL)) itsRandomGen->reset(nbOfPointsList, factorNb);
	resetwork(nbOfPointsList, factorNb);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

