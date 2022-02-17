/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file skipper.cpp
 *  \brief 
 *	\author  R. GUILLEMOT
 *	\version 1.0
 *	\date Novvember 2005
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpnumlib/skipper.h"
#include "gpnumlib/random.h"
#include "gpbase/ostringstream.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_SkipperFunction
///	Routine: ARM_SkipperFunction
///	Returns: void
///	Action : Constructor
////////////////////////////////////////////////////

ARM_SkipperFunction::ARM_SkipperFunction(const ARM_RandomGeneratorPtr& randomGen, double nbStdDev) :
ARM_ManipulatorMethod(randomGen),
itsNbStdDev(nbStdDev)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_SkipperFunction
///	Routine: ARM_SkipperFunction
///	Returns: void
///	Action : Copy Constructor
////////////////////////////////////////////////////

ARM_SkipperFunction::ARM_SkipperFunction(const ARM_SkipperFunction& rhs)
: ARM_ManipulatorMethod(rhs),
itsNbStdDev(rhs.itsNbStdDev)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_SkipperFunction
///	Routine: Clone()
///	Returns: void
///	Action : Create a copy of the object
////////////////////////////////////////////////////
ARM_ManipulatorMethod* ARM_SkipperFunction::Clone() const
{
	return new ARM_SkipperFunction(*this);
}

////////////////////////////////////////////////////
///	Class  : ARM_SkipperFunction
///	Routine: ~ARM_SkipperFunction
///	Returns: void
///	Action : Destructor
////////////////////////////////////////////////////

ARM_SkipperFunction::~ARM_SkipperFunction()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_SkipperFunction  
///	Routine: reset
///	Returns: void
///	Action : reset the random number generator
////////////////////////////////////////////////////
//void ARM_SkipperFunction::reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorsNb ) 
//{
//	itsRandomGen->reset(dim, nbOfPointsList, factorsNb);
//}

////////////////////////////////////////////////////
///	Class  : ARM_SkipperFunction
///	Routine: operator()
///	Returns: void
///	Action : Reset the transoser
////////////////////////////////////////////////////

double ARM_SkipperFunction::operator()() const
{
	double x;

	while (((x = (*itsRandomGen)()) >itsNbStdDev) || (x < -itsNbStdDev))
		x = 1.0;
	return x;
}

////////////////////////////////////////////////////
///	Class  : ARM_SkipperFunction
///	Routine: toString()
///	Returns: void
///	Action : Display the contents
////////////////////////////////////////////////////

string ARM_SkipperFunction::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << " Skipper Function:" << CC_NS(std,endl);
	os << indent << " Nb Std Dev:" << itsNbStdDev;

	return os.str();
}

CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----