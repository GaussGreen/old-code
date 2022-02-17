/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file condinvcum.cpp
 *  \brief 
 *	\author  R. GUILLEMOT
 *	\version 1.0
 *	\date Novvember 2005
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpnumlib/condinvcum.h"

#include "gpnumlib/random.h"
#include "gpnumlib/gaussiananalytics.h"
#include <math.h>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_CondInvCumFunction
///	Routine: ARM_CondInvCumFunction
///	Returns: void
///	Action : Constructor
////////////////////////////////////////////////////

ARM_CondInvCumFunction::ARM_CondInvCumFunction(const ARM_RandomGeneratorPtr& randomGen, invCumAlgoType algoType, double nbStdDev) :
ARM_InvCumFunction(randomGen,algoType),
itsNbStdDev(nbStdDev)
{
	itsMinProba = ARM_GaussianAnalytics::cdfNormal2(-itsNbStdDev);
	itsMaxProba = ARM_GaussianAnalytics::cdfNormal2(itsNbStdDev);
	size_t dim = randomGen->dim();

// FIXMEFRED: mig.vc8 (22/05/2007 18:12:55):cast
	itsProbaCorrect = pow(itsMaxProba-itsMinProba,static_cast<int>(dim));
}

////////////////////////////////////////////////////
///	Class  : ARM_CondInvCumFunction
///	Routine: ARM_CondInvCumFunction
///	Returns: void
///	Action : Copy Constructor
////////////////////////////////////////////////////

ARM_CondInvCumFunction::ARM_CondInvCumFunction(const ARM_CondInvCumFunction& rhs)
: ARM_InvCumFunction(rhs),
itsNbStdDev(rhs.itsNbStdDev),
itsMinProba(rhs.itsMinProba),
itsMaxProba(rhs.itsMaxProba),
itsProbaCorrect(rhs.itsProbaCorrect)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_CondInvCumFunction
///	Routine: Clone()
///	Returns: void
///	Action : Create a copy of the object
////////////////////////////////////////////////////

ARM_ManipulatorMethod* ARM_CondInvCumFunction::Clone() const
{
	return new ARM_CondInvCumFunction(*this);
}

////////////////////////////////////////////////////
///	Class  : ARM_CondInvCumFunction
///	Routine: ~ARM_CondInvCumFunction
///	Returns: void
///	Action : Destructor
////////////////////////////////////////////////////

ARM_CondInvCumFunction::~ARM_CondInvCumFunction()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_CondInvCumFunction
///	Routine: operator()
///	Returns: void
///	Action : Reset the transoser
////////////////////////////////////////////////////

double ARM_CondInvCumFunction::operator()() const
{
	double x = (*itsRandomGen)();

	return itsFunc(x*itsMinProba+(1.0-x)*itsMaxProba);
}

////////////////////////////////////////////////////
///	Class  : ARM_CondInvCumFunction
///	Routine: toString()
///	Returns: void
///	Action : Display the contents
////////////////////////////////////////////////////

string ARM_CondInvCumFunction::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << " Cond Inv Cum Function:" << CC_NS(std,endl);
	os << indent << " Min Proba:" << itsMinProba << CC_NS(std,endl);
	os << indent << " Max Proba:" << itsMaxProba << CC_NS(std,endl);
	os << indent << " Proba Correct:" << itsProbaCorrect;

	return os.str();
}

CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----