/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file impsampleropt.cpp
 *  \brief
 *
 *	\author  R Guillemot
 *	\version 1.0
 *	\date December 2005
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/singleton.h"
#include "gpbase/eventviewer.h"
#include "gpbase/ostringstream.h"


#include "gpinfra/pricingmodel.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/gramfunctorarg.h"

#include "gpnumlib/random.h"
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/optimizernd.h"

#include "gpnummethods/typedef.h"
#include "gpnummethods/impsampleropt.h"
#include "gpnummethods/mcmethod.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ImpSamplerOptFunc
/// Desc: This functor is used to optimze the drift
/// of the importance sampling
////////////////////////////////////////////////////
class ImpSamplerOptFunc : public NagNDOptimizer::NDFunc
{
public:
	ImpSamplerOptFunc(
		ARM_GenSecurity* genSec,
		ARM_PricingModel* pricingModel,
		ARM_GP_MultiCurvePtr alpha,
		bool withMC,
		const ARM_GP_VectorPtr& rowTimes=ARM_GP_VectorPtr(new ARM_GP_Vector(NULL)),
		int timeStep=-1)
		:
	itsGenSec(genSec),
	itsPricingModel(pricingModel),
	itsAlpha(alpha),
	itsWithMC(withMC),
	itsRowTimes(rowTimes),
	itsTimeStep(timeStep)
	{
	}

	virtual double operator()(std::vector<double> x) const;

private:
	ARM_GenSecurity* itsGenSec;
	ARM_PricingModel* itsPricingModel;
	ARM_GP_MultiCurvePtr itsAlpha;
	bool itsWithMC;
	ARM_GP_VectorPtr itsRowTimes;
	int itsTimeStep;
};

double ImpSamplerOptFunc::operator()(std::vector<double> x) const
{
	if (itsTimeStep == -1)
	{
		// Insert the new alpha curve point
		itsAlpha->insert(0, x);
	}
	else
	{
		itsAlpha->insert((*itsRowTimes)[itsTimeStep], x);
	}

	// Evaluate the current row of the gensec
	ARM_GenPricer genPricer(
		itsGenSec,
		itsPricingModel);

	double ret;

	genPricer.Price();
	if (itsTimeStep == -1)
	{
		if (itsWithMC)
			// We want the minimim of variance
			ret = genPricer.GetPricerInfo()->GetContents().GetData("StdDev").GetDouble();
		else
			// We want the maximum of value
			ret = -genPricer.GetPricerInfo()->GetContents().GetData("Price").GetDouble();
	}
	else
	{
		if (itsWithMC)
			// We want the minimim of variance
			ret = (*genPricer.GetPricerInfo()->GetContents().GetData("IntermediateStdDevs").GetVector())[itsTimeStep];
		else
			// We want the maximum of value
			ret = -(*genPricer.GetPricerInfo()->GetContents().GetData("IntermediatePrices").GetVector())[itsTimeStep];
	}

	return ret;
}


////////////////////////////////////////////////////
///	Class  : ARM_ImpSamplerOpt
///	Routine: ComputeModel
///	Returns:
///	Action : Optimize the model for the
/// importance sampling
////////////////////////////////////////////////////
ARM_PricingModel* ARM_ImpSamplerOpt::ComputeModel(
		ARM_GenSecurity* genSec,
		ARM_PricingModel* model,
		const ARM_MultiCurve& initGuess,
		const ARM_MultiCurve& lowerBounds,
		const ARM_MultiCurve& upperBounds,
		bool withMC,
		long nbSteps,
		bool bootstrap)
{
	// Create a copy of the model
	ARM_PricingModel* newModel = static_cast<ARM_PricingModel*>(model->Clone());

	// Get the MC method of the model
	ARM_MCMethod* MCMethod = dynamic_cast<ARM_MCMethod*>(&*(model->GetNumMethod()));

	if (!MCMethod)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": the model doesn't contain a MC num method.");

	ARM_MCMethod* optMCMethod = static_cast<ARM_MCMethod*>(MCMethod->Clone());

	if (withMC)
	{
		optMCMethod->Initialize(nbSteps,optMCMethod->GetRandGenVector(),nbSteps);
	}
	else
	{
		// Create a fake mc method with one path and a "Null" generator for the optimization
		ARM_RandomGeneratorPtr  unifRandomGen( ARM_RandGenFactory.Instance()->CreateRandGen( 
			ARM_RandGenFactoryImp::Null,
			ARM_RandGenFactoryImp::UnknownTransformAlgo,
			ARM_RandomGeneratorPtr(NULL),
			ARM_RandomGeneratorPtr(NULL),
			-1,
			1,
			1,
			ARM_GP_T_Vector<size_t>(1,10),
			0.5) );

		/// Normal
		ARM_RandomGeneratorPtr normRandGen( ARM_RandGenFactory.Instance()->CreateRandGen( 
			ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
			ARM_RandGenFactoryImp::InvNormCum,
			unifRandomGen ) );

		vector<ARM_RandomGeneratorPtr> randGenList(1,normRandGen);
		optMCMethod->Initialize(nbSteps,randGenList,nbSteps);
	}

	ARM_NumMethodPtr optNumMethod(optMCMethod);
	newModel->SetNumMethod(optNumMethod);

	// Retreive the alpha curve of the importance sampler

	ARM_PropImpSampler* propImpSampler = dynamic_cast<ARM_PropImpSampler*>(&*optMCMethod->GetImpSampler());

	if (!propImpSampler)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": The MC method doesn't contain a proportional importance sampler.");

	ARM_GP_MultiCurvePtr alpha = propImpSampler->GetAlpha();
	size_t nbFactor = model->FactorCount();

	double asOfDate = model->GetAsOfDate().GetJulian();

	if (bootstrap)
	{
		ARM_GP_VectorPtr rowSteps = ARM_GP_VectorPtr(new ARM_GP_Vector(*genSec->GetRowDates()));

		*(rowSteps) -= newModel->GetAsOfDate().GetJulian();

		// For each time steps of the gensec
		for (size_t i = 0; i < rowSteps->size(); ++i)
		{
			ImpSamplerOptFunc samplerOptFunc(
			genSec,
			newModel,
			alpha,
			withMC,
			rowSteps,
			i);	

			// The optimizer
			NagNDOptimizer optim(&samplerOptFunc);

			// Optimize the alpha curve for each row of the gensec
			std::vector<double> guess = initGuess.Interpolate((*rowSteps)[i]);
			std::vector<double> lowerBound = lowerBounds.Interpolate((*rowSteps)[i]);
			std::vector<double> upperBound = upperBounds.Interpolate((*rowSteps)[i]);

			// fill missing with the last value

			if (guess.size() < nbFactor)
			{
				size_t size = guess.size();
				for (size_t i = 0; i < nbFactor-size; ++i)
					guess.push_back(guess[size-1]);
			}
			if (lowerBound.size() < nbFactor)
			{
				size_t size = lowerBound.size();
				for (size_t i = 0; i < nbFactor-size; ++i)
					lowerBound.push_back(lowerBound[size-1]);
			}
			if (upperBound.size() < nbFactor)
			{
				size_t size = upperBound.size();
				for (size_t i = 0; i < nbFactor-size; ++i)
					upperBound.push_back(upperBound[size-1]);
			}

			std::vector<double> res;

			res = optim.Optimize(
							guess,
							lowerBound,
							upperBound);

			alpha->insert((*rowSteps)[i], res);
		}
	}
	else
	{
		// Buid the multi dim optimize
		// The function
		ImpSamplerOptFunc samplerOptFunc(
			genSec,
			newModel,
			alpha,
			withMC);
		// The optimizer
		NagNDOptimizer optim(&samplerOptFunc);

		// Optimize the alpha curve for each row of the gensec
		std::vector<double> guess = initGuess.Interpolate(0);
		std::vector<double> lowerBound = lowerBounds.Interpolate(0);
		std::vector<double> upperBound = upperBounds.Interpolate(0);

		// fill missing with the last value

		if (guess.size() < nbFactor)
		{
			size_t size = guess.size();
			for (size_t i = 0; i < nbFactor-size; ++i)
				guess.push_back(guess[size-1]);
		}
		if (lowerBound.size() < nbFactor)
		{
			size_t size = lowerBound.size();
			for (size_t i = 0; i < nbFactor-size; ++i)
				lowerBound.push_back(lowerBound[size-1]);
		}
		if (upperBound.size() < nbFactor)
		{
			size_t size = upperBound.size();
			for (size_t i = 0; i < nbFactor-size; ++i)
				upperBound.push_back(upperBound[size-1]);
		}

		std::vector<double> res;

		res = optim.Optimize(
						guess,
						lowerBound,
						upperBound);

		alpha->insert(0, res);
	}

	// Create a new MC method which is a copy of the previous one with the optimized
	// alpha curve
	ARM_MCMethod* newMCMethod = static_cast<ARM_MCMethod*>(MCMethod->Clone());
	ARM_NumMethodPtr newNumMethod(newMCMethod);

	propImpSampler = dynamic_cast<ARM_PropImpSampler*>(&*newMCMethod->GetImpSampler());
	if (!propImpSampler)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": The MC method doesn't contain a proportional importance sampler.");

	propImpSampler->SetAlpha(ARM_GP_MultiCurvePtr(static_cast<ARM_MultiCurve*>(alpha->Clone())));

	newModel->SetNumMethod(newNumMethod);

	return newModel;
}


CC_END_NAMESPACE()