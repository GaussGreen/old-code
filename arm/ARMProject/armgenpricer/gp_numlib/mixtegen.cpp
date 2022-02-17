/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file mixtegen.cpp
 *  \brief 
 *	\author  R. GUILLEMOT
 *	\version 1.0
 *	\date December 2005
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/ostringstream.h"
#include "gpbase/cloneutilityfunc.h"
#include "gpnumlib/mixtegen.h"
#include "gpnumlib/random.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_MixteGen
///	Routine: ARM_MixteGen
///	Returns: void
///	Action : Constructor
////////////////////////////////////////////////////

ARM_MixteGen::ARM_MixteGen(
const ARM_RandomGeneratorPtr& randomGen1,
const ARM_RandomGeneratorPtr& randomGen2,
int firstNbTimes,
int firstNbDims) :
ARM_ManipulatorMethod(randomGen1),
itsRandomGen2(randomGen2),
itsFirstNbTimes(firstNbTimes),
itsFirstNbDims(firstNbDims),
itsNbTimes(0),
itsNbOfPointsList(0),
itsCountor(0),
itsDimCountor(0),
itsTimeCountor(0),
itsBucketCountor(0)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_MixteGen
///	Routine: ARM_MixteGen
///	Returns: void
///	Action : Copy Constructor
////////////////////////////////////////////////////

ARM_MixteGen::ARM_MixteGen(const ARM_MixteGen& rhs)
: 
ARM_ManipulatorMethod(rhs),
itsRandomGen2(rhs.itsRandomGen2),
itsFirstNbTimes(rhs.itsFirstNbTimes),
itsFirstNbDims(rhs.itsFirstNbDims),
itsNbTimes(rhs.itsNbTimes),
itsNbOfPointsList(rhs.itsNbOfPointsList),
itsCountor(rhs.itsCountor),
itsDimCountor(rhs.itsDimCountor),
itsBucketCountor(rhs.itsBucketCountor)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_Transposer
///	Routine: Clone()
///	Returns: void
///	Action : Create a copy of the object
////////////////////////////////////////////////////

ARM_ManipulatorMethod* ARM_MixteGen::Clone() const
{
	return new ARM_MixteGen(*this);
}

////////////////////////////////////////////////////
///	Class  : ARM_MixteGen
///	Routine: ~ARM_MixteGen
///	Returns: void
///	Action : Destructor
////////////////////////////////////////////////////

ARM_MixteGen::~ARM_MixteGen()
{
}


////////////////////////////////////////////////////
///	Class  : ARM_MixteGen
///	Routine: reset
///	Returns: void
///	Action : Reset the mixte generator
////////////////////////////////////////////////////

void ARM_MixteGen::reset(const ARM_GP_T_Vector<size_t>& nbOfPointsList, const ARM_GP_T_Vector<size_t>& factorNb)
{
	resetwork(nbOfPointsList, factorNb);
}

////////////////////////////////////////////////////
///	Class  : ARM_MixteGen
///	Routine: reset
///	Returns: void
///	Action : Reset the mixte generator
////////////////////////////////////////////////////

void ARM_MixteGen::resetwork(const ARM_GP_T_Vector<size_t>& nbOfPointsList, const ARM_GP_T_Vector<size_t>& factorNb )
{
	itsNbTimes = factorNb.size();
	itsStepNbFactor = factorNb;
	itsNbOfPointsList = nbOfPointsList;

	itsCountor = 0;
	itsDimCountor = 0;
	itsTimeCountor = 0;
	itsBucketCountor = 0;


	ARM_GP_T_Vector<size_t> factorNb1;
	ARM_GP_T_Vector<size_t> factorNb2;
	
	int firstTotalNbDims = 0;
	int totalNbDims = 0;

	for(int k = 0; (k < itsFirstNbTimes) && (k < itsNbTimes); k++)
	{
		totalNbDims += factorNb[k];
		if (itsFirstNbDims > 0)
		{
			factorNb1.push_back(factorNb[k] > itsFirstNbDims ? itsFirstNbDims : factorNb[k]);
			firstTotalNbDims += (factorNb[k] > itsFirstNbDims ? itsFirstNbDims : factorNb[k]);
		}

		if (factorNb[k] > itsFirstNbDims)
			factorNb2.push_back(factorNb[k] - itsFirstNbDims);
	}

	for(int k = itsFirstNbTimes; k < itsNbTimes; k++)
	{
		totalNbDims += factorNb[k];
		factorNb2.push_back(factorNb[k]);
	}

	itsRandomGen->reset(nbOfPointsList, factorNb1);
	if (factorNb2.size())
		itsRandomGen2->reset(totalNbDims-firstTotalNbDims,nbOfPointsList, factorNb[0]);
}

////////////////////////////////////////////////////
///	Class  : ARM_MixteGen
///	Routine: operator()
///	Returns: void
///	Action : Reset the transoser
////////////////////////////////////////////////////0

double ARM_MixteGen::operator()() const
{
	double value;

	if ((itsTimeCountor < itsFirstNbTimes) && (itsDimCountor < itsFirstNbDims))
	{
		value = (*itsRandomGen)();
	}
	else
		value = (*itsRandomGen2)();
	
	itsDimCountor = (itsDimCountor+1)%itsStepNbFactor[itsTimeCountor];
	if (!itsDimCountor)
		itsCountor = (itsCountor+1)%itsNbOfPointsList[itsBucketCountor];
	if (!itsDimCountor&&!itsCountor)
		itsTimeCountor = (itsTimeCountor+1)%(itsNbTimes);
	if (!itsCountor&&!itsDimCountor&&!itsTimeCountor)
		itsBucketCountor++;

	return value;
}

////////////////////////////////////////////////////
///	Class  : ARM_MixteGen
///	Routine: toString()
///	Returns: void
///	Action : Display the contents
////////////////////////////////////////////////////

string ARM_MixteGen::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;

	os << " MixteGen:\n";
	os << " " << itsFirstNbTimes << " first times " << itsFirstNbDims << " first factors:\n";
	os << itsRandomGen->toString() << "\n";

	size_t remainNbFactor = 0;

	for(int k = 0; (k < itsFirstNbTimes) && (k < itsNbTimes); k++)
	{
		if (itsStepNbFactor[k] < itsFirstNbDims)
			remainNbFactor += itsStepNbFactor[k] - itsFirstNbDims;
	}

	for(int k = itsFirstNbTimes; k < itsNbTimes; k++)
	{
		remainNbFactor += itsStepNbFactor[k];
	}
	

	if (remainNbFactor > 0)
	{
		os << " then "<< remainNbFactor << " : " << "\n";
		os << itsRandomGen2->toString() << "\n";
	}

	return os.str();
}

CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----