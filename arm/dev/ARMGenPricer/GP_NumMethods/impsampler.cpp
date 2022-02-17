/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file impsampler.cpp
 *
 *  \brief sampler object are used to discretise in time according to 
 *      various hypotheses
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date October 2005
 */

#include "gpbase/removeidentifiedwarning.h"

#include "gpnummethods/impsampler.h"

// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/stringmanip.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/gpvector.h"

// gpnummethods
#include "gpnummethods/sampler.h"

CC_USING_NS( std, pair )
CC_USING_NS( std, make_pair )
CC_USING_NS( std, endl )


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// The Dummy importance sampler class
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_DummyImpSampler
///	Routine: Constructor
///	Returns:
///	Action : 
////////////////////////////////////////////////////

ARM_DummyImpSampler::ARM_DummyImpSampler() :
ARM_ImpSampler()
{
}


////////////////////////////////////////////////////
///	Class  : ARM_DummyImpSampler
///	Routine: Copy Constructor
///	Returns:
///	Action : 
////////////////////////////////////////////////////

ARM_DummyImpSampler::ARM_DummyImpSampler(const ARM_DummyImpSampler& rhs) :
ARM_ImpSampler(rhs)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_DummyImpSamper
///	Routine: operator=
///	Returns:
///	Action : 
////////////////////////////////////////////////////

ARM_DummyImpSampler& ARM_DummyImpSampler::operator=(const ARM_DummyImpSampler& rhs)
{
	if( this != & rhs )
	{
		this->~ARM_DummyImpSampler();
        new (this) ARM_DummyImpSampler(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_DummyImpSampler
///	Routine: operator=
///	Returns:
///	Action : Compute the drift to add at the process
/// It does nothing
////////////////////////////////////////////////////

ARM_VectorPtrVector ARM_DummyImpSampler::ComputeDrifts(
	int factorNb,
	const ARM_GP_Vector& timeSteps,
	ARM_SamplerBase* sampler)
{
	ARM_VectorPtrVector dummyDrifts;

	return dummyDrifts;
}

////////////////////////////////////////////////////
///	Class  : ARM_DummyImpSampler
///	Routine: toString
///	Returns:
///	Action : Display the content of the object
////////////////////////////////////////////////////
string ARM_DummyImpSampler::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;

	os << "Dummy (it does nothing)" << CC_NS(std,endl);

	return os.str();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// The Proportional importance sampler class
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_PropImpSampler
///	Routine: Constructor
///	Returns:
///	Action : 
////////////////////////////////////////////////////

ARM_PropImpSampler::ARM_PropImpSampler(const ARM_MultiCurve* alpha) :
ARM_ImpSampler(),
itsAlpha(static_cast<ARM_MultiCurve*>(alpha->Clone()))
{
}


////////////////////////////////////////////////////
///	Class  : ARM_PropImpSampler
///	Routine: Copy Constructor
///	Returns:
///	Action : 
////////////////////////////////////////////////////

ARM_PropImpSampler::ARM_PropImpSampler(const ARM_PropImpSampler& rhs) :
ARM_ImpSampler(rhs),
itsAlpha(static_cast<ARM_MultiCurve*>(rhs.itsAlpha->Clone()))
{
}

////////////////////////////////////////////////////
///	Class  : ARM_PropImpSampler
///	Routine: operator=
///	Returns:
///	Action : 
////////////////////////////////////////////////////

ARM_PropImpSampler& ARM_PropImpSampler::operator=(const ARM_PropImpSampler& rhs)
{
	if( this != & rhs )
	{
		this->~ARM_PropImpSampler();
        new (this) ARM_PropImpSampler(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_PropImpSampler
///	Routine: operator=
///	Returns:
///	Action : Compute the drift to add at the process
////////////////////////////////////////////////////

ARM_VectorPtrVector ARM_PropImpSampler::ComputeDrifts(
	int factorNb,
	const ARM_GP_Vector& timeSteps,
	ARM_SamplerBase* sampler)
{
	ARM_SamplerNDBase* samplerND = sampler->ToSamplerNDBase();

	int nbTimeSteps = timeSteps.size();

	ARM_VectorPtrVector drifts(nbTimeSteps-1);

	ARM_GP_Vector meanRevCoeff(factorNb,1.0);
	
	size_t timeIdx, factorIdx;

	ARM_GP_Vector alpha(factorNb,0.0);

	for(timeIdx = 0; timeIdx < nbTimeSteps-1; ++timeIdx)
	{
		alpha = itsAlpha->Interpolate(timeSteps[timeIdx+1]);

		// if the alpha curve don't have the good size we resize it
		if (alpha.size() < factorNb)
			alpha = ARM_GP_Vector(factorNb,0.0);

		drifts[timeIdx] = ARM_GP_VectorPtr(new ARM_GP_Vector(factorNb));

		for (factorIdx = 0; factorIdx < factorNb; ++factorIdx)
		{
			meanRevCoeff[factorIdx] /= samplerND->GetRelativeDrifts()(timeIdx,factorIdx);
			(*drifts[timeIdx])[factorIdx] = samplerND->GetLocalVar(timeIdx)[factorIdx]*alpha[factorIdx]*meanRevCoeff[factorIdx];
		}
	}

	return drifts;
}

void ARM_PropImpSampler::BootstrapAlpha(
	int factorNb,
	const ARM_GP_Vector& timeSteps,
	const ARM_GP_Vector& rowSteps,
	const ARM_GP_T_Vector<ARM_GP_Vector>& globalAlpha,
	ARM_SamplerBase* sampler)
{
	ARM_SamplerNDBase* samplerND = sampler->ToSamplerNDBase();

	int nbTimeSteps = timeSteps.size();
	int nbRowSteps = rowSteps.size();

	ARM_GP_Vector lastDrift(factorNb,0.0);
	ARM_GP_Vector lastDriftWithAlpha(factorNb,0.0);
	ARM_GP_Vector incDrift(factorNb,0.0);

	ARM_GP_Vector meanRevCoeff(factorNb,1.0);
	
	size_t rowIdx, timeIdx, factorIdx;

	ARM_GP_Vector alpha(factorNb,0.0);

	for(rowIdx = 0, timeIdx = 0; rowIdx < nbRowSteps; ++rowIdx)
	{
		for (factorIdx = 0; factorIdx < factorNb; ++factorIdx)
			incDrift[factorIdx] = 0.0;


		while((timeIdx < nbTimeSteps-1) && timeSteps[timeIdx]<=rowSteps[rowIdx])
		{
			for (factorIdx = 0; factorIdx < factorNb; ++factorIdx)
			{
				meanRevCoeff[factorIdx] /= samplerND->GetRelativeDrifts()(timeIdx,factorIdx);
				incDrift[factorIdx] += samplerND->GetLocalVar(timeIdx)[factorIdx]*meanRevCoeff[factorIdx];
			}
			++timeIdx;
		}

		timeIdx--;

		for (factorIdx = 0; factorIdx < factorNb; ++factorIdx)
		{
			alpha[factorIdx] = ((globalAlpha[rowIdx])[factorIdx]*(lastDrift[factorIdx]+incDrift[factorIdx])-lastDriftWithAlpha[factorIdx])/incDrift[factorIdx];
			lastDrift[factorIdx] += incDrift[factorIdx];
			lastDriftWithAlpha[factorIdx] += alpha[factorIdx]*incDrift[factorIdx];
		}

		itsAlpha->insert(rowSteps[rowIdx], alpha);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_PropImpSampler
///	Routine: operator=
///	Returns:
///	Action : Compute the exponential to multiply with the payoff
////////////////////////////////////////////////////
 
void ARM_PropImpSampler::ComputeExps(
	int nbStates,
	int factorNb,
	const ARM_GP_Vector& timeSteps,
	const ARM_MatrixPtrVector& processStates,
	ARM_SamplerBase* sampler)
{
	// Remove all the element of the map
	// Very important for multi bucket
	itsExpProcessStates.clear();

	ARM_SamplerNDBase* samplerND = sampler->ToSamplerNDBase();

	int nbTimeSteps = timeSteps.size();

	size_t timeIdx, stateIdx, factorIdx;

	ARM_GP_Vector alpha(factorNb,0.0);
	ARM_GP_Vector meanRevCoeff(factorNb,1.0);
	ARM_GP_Vector expMat(nbStates,1.0);

	double var;

	for(timeIdx = 0; timeIdx < nbTimeSteps-1; ++timeIdx)
	{
		alpha = itsAlpha->Interpolate(timeSteps[timeIdx+1]);

		// if the alpha curve don't have the good size we resize it
		if (alpha.size() < factorNb)
			alpha = ARM_GP_Vector(factorNb,0.0);
			
		for (factorIdx = 0; factorIdx < factorNb; ++factorIdx)
		{
			meanRevCoeff[factorIdx] /= samplerND->GetRelativeDrifts()(timeIdx,factorIdx);
			for (stateIdx = 0; stateIdx < nbStates; ++stateIdx)
			{
				var = samplerND->GetLocalVar(timeIdx)[factorIdx]*alpha[factorIdx]*alpha[factorIdx]*meanRevCoeff[factorIdx]*meanRevCoeff[factorIdx];
				expMat(stateIdx) *= exp(-0.5*var-alpha[factorIdx]*meanRevCoeff[factorIdx]*(*processStates[timeIdx])(factorIdx,stateIdx));
			}
		}

		itsExpProcessStates.insert(make_pair(timeSteps[timeIdx+1], ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(expMat.Clone()))));
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_PropImpSampler
///	Routine: ProcessPaidPayoffs 
///	Returns: 
///	Action : This functions applies the importance 
/// sampling exponential
////////////////////////////////////////////////////

void ARM_PropImpSampler::ProcessPaidPayoffs(ARM_VectorPtr& payoffs, double evalTime) const
{	
	ARM_ImpSampler::ExpMap::const_iterator iter;

	iter = itsExpProcessStates.find( evalTime );

	if (iter == itsExpProcessStates.end())
	{
		CC_Ostringstream os;
		os << "Importance Sampling exponential cannot be found in the map for the date << evalTime" << endl;
		ARM_THROW( ERR_INVALID_ARGUMENT, os.str());
	}

	for (size_t stateIdx = 0; stateIdx < payoffs->size(); ++stateIdx)
		(*payoffs)[stateIdx] *= (*((*iter).second))[stateIdx];
}

////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: ProcessUnPaidPayoffs
///	Returns: 
///	Action : This function removes the importance 
/// sampling exponential
////////////////////////////////////////////////////

void ARM_PropImpSampler::ProcessUnPaidPayoffs(ARM_VectorPtr& payoffs, double evalTime) const
{
	ARM_ImpSampler::ExpMap::const_iterator iter;
	
	iter = itsExpProcessStates.find( evalTime );

	if (iter == itsExpProcessStates.end())
	{
		CC_Ostringstream os;
		os << "Importance Sampling exponential cannot be found in the map for the date << evalTime" << endl;
		ARM_THROW( ERR_INVALID_ARGUMENT, os.str());
	}

	for (size_t stateIdx = 0; stateIdx < payoffs->size(); ++stateIdx)
		(*payoffs)[stateIdx] /= (*((*iter).second))[stateIdx];
}

////////////////////////////////////////////////////
///	Class  : ARM_PropImpSampler
///	Routine: toString
///	Returns:
///	Action : Display the content of the object
////////////////////////////////////////////////////
string ARM_PropImpSampler::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;

	os << "Proportional" << CC_NS(std,endl);
	os << indent << "Alpha:" << CC_NS(std,endl);
	os << indent << itsAlpha->toString();

	return os.str();
}

CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

