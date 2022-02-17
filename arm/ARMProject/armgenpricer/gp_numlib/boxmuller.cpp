/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file boxmuller.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#include "gpnumlib/boxmuller.h"
#include "gpnumlib/uniform.h"
#include "gpnumlib/quasirandom.h"
#include "gpbase/ostringstream.h"
#include "expt.h"
#include <cmath>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///				ARM_BoxMuller				////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_BoxMuller
///	Routine: constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_BoxMuller::ARM_BoxMuller(const ARM_RandomGeneratorPtr& randomGen )
: ARM_ManipulatorMethod(randomGen), 	itsIset( false ), itsGset( 0 )
{
	if(randomGen == ARM_RandomGeneratorPtr(NULL) )
	    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "pointor on random number generator is NULL!");
	if( dynamic_cast<ARM_QuasiRandom*>(&*randomGen) )
	    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "box muller incompatible with quasi random.!");
	if( !dynamic_cast<ARM_UniformGenerator*>(&*randomGen) 
		|| randomGen->GetDistributionType() != ARM_RandomGenerator::ARM_Uniform )
	    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "random number base generator is not random uniform.!");
}


////////////////////////////////////////////////////
///	Class  : ARM_BoxMuller  
///	Routine: reset
///	Returns: void
///	Action : reset the random number generator
////////////////////////////////////////////////////
void ARM_BoxMuller::resetwork( const ARM_GP_T_Vector<size_t>& nbOfPointsList, const ARM_GP_T_Vector<size_t>& factorsNb ) 
{
	itsIset = false;
	itsGset = 0;
}

////////////////////////////////////////////////////
///	Class  : ARM_BoxMuller
///	Routine: operator()() const
///	Returns: 
///	Action : uses the box muller method to compute a normal from uniform
////////////////////////////////////////////////////
double ARM_BoxMuller::operator()() const
{
	double fac,rsq,v1,v2;

	if ( !itsIset ) 
	{
		do 
		{
			v1  = 2.0 * (*itsRandomGen)() - 1.0;
			v2  = 2.0 * (*itsRandomGen)() - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} 
		while ( ( rsq >= 1.0 ) || ( rsq == 0.0 ) );

		fac		= sqrt( -2.0 * log( rsq ) / rsq );
		CC_MUTABLE( ARM_BoxMuller, itsGset ) = v1 * fac;
		CC_MUTABLE( ARM_BoxMuller, itsIset ) = true;
		return v2 * fac;
	} 
	else 
	{
		CC_MUTABLE( ARM_BoxMuller, itsIset ) = false;
		return itsGset;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_NormalCompositeGen
///	Routine: BoxMuller
///	Returns: 
///	Action : uses the box muller method to compute a normal
////////////////////////////////////////////////////
string ARM_BoxMuller::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << " Normal Generator  : " << "Box Muller method" << CC_NS(std,endl);
	os << " Base Generator    : " << itsRandomGen->toString();
	return os.str();
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_BoxMuller
///	Routine: copy constructor, assignment operator
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_BoxMuller::ARM_BoxMuller( const ARM_BoxMuller& rhs )
:	ARM_ManipulatorMethod(rhs), /// clone it for independence!
	itsIset( rhs.itsIset ),
	itsGset( rhs.itsGset )
{}


ARM_BoxMuller& ARM_BoxMuller::operator =(const ARM_BoxMuller& rhs )
{
	if( this != & rhs )
	{
		ARM_ManipulatorMethod::operator =( rhs );
		itsRandomGen	= ARM_RandomGeneratorPtr( (ARM_RandomGenerator*) rhs.itsRandomGen->Clone() ); /// clone it for independence!
		itsIset			= rhs.itsIset;
		itsGset			= rhs.itsGset;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_BoxMuller
///	Routine: desctructor (necessary even for pure virtual case)
///	Action : draw a vector of uniform random numbers
////////////////////////////////////////////////////
ARM_BoxMuller::~ARM_BoxMuller()
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_BoxMuller
///	Routine: Clone
///	Returns: ARM_ManipulatorMethod*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_ManipulatorMethod* ARM_BoxMuller::Clone() const
{
	return new ARM_BoxMuller(*this);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

