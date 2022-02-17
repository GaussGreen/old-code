//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : CCMCopulaModel.hpp
//
//   Description : Implementation of class CCMCopulaModel, which implements the 
//				   CopulaModel interface as a mixture of an RFL copula, an
//                 independence copula and a dependece copula
//
//   Date        : Feb 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CCMCopulaModel.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/CcmOnlyParameters.hpp"
#include "edginc/CDSHelper.hpp"
#include "edginc/MaturityPeriod.hpp"

DRLIB_BEGIN_NAMESPACE

#define CCM_MC_TINY  1e-10


void CCMCopulaModel::Update()
{
	depenCopula.valueDate = valueDate;
	rflCopula.valueDate = valueDate;
	indepCopula.valueDate = valueDate;

	depenCopula.names = names;
	rflCopula.names = names;
	indepCopula.names = names;
};

void
CCMCopulaModel::addNames(const CreditAssetWrapperArray& xnames)
{
	const string method = "CCMCopulaModel::addNames";
	try
	{
		CopulaModelBase::addNames(xnames);
		floorCurves.clear();
		int numNames = xnames.size();

		for (int i = 0; i < numNames; i++)
		{
			SingleCreditAssetConstSP namei = 
				SingleCreditAssetConstSP::dynamicCast(names[i]->getSP());

			CreditEngineParametersConstSP 
				cep(namei->getEngineParams(CcmOnlyParameters::TYPE));

			const CcmOnlyParameters *ccmp = 
				dynamic_cast<const CcmOnlyParameters *>(cep.get());

			if (!ccmp)
			{
				throw ModelException("Failed to obtain CcmOnlyParameters");
			}
			else
			{
				if (ccmp->seniorCurve.get())
					floorCurves.push_back(ccmp->seniorCurve);
			}
		}
	}
	catch (exception &e)
	{
		throw ModelException(e, method);
	}
}


CCMCopulaModel::Simulation::Simulation(
	const DateTime& valueDate,
	CopulaModelBase &model,
	DateTimeArray timeLine,
	int nSamples) 
	: BaseSimulation(valueDate, model, timeLine, nSamples)
{
	CCMCopulaModel *model1 = dynamic_cast<CCMCopulaModel*>(&model);

	if (!model1)
	{
		throw ModelException("Expected CCMCopula");
	}

	model1->populateInnerCopulas();

	// each simulation constructor computes surv. probabilities to each 
	// of the dates on the timeline, and resets its random num generator

	model1->Update();

	model1->ComputeCurves(timeLine);

	dcSim = model1->depenCopula.simulation(
		timeLine,
		nSamples);

	indSim = model1->indepCopula.simulation(
		timeLine, 
		nSamples);

	rflSim = model1->rflCopula.simulation(
		timeLine, 
		nSamples);

	randsPerPath =  rflSim->getNumberRands() + 
					indSim->getNumberRands() + 
					dcSim->getNumberRands();

};

void CCMCopulaModel::Simulation::UpdateSim()
{
	dcSim->UpdateSim(); 
    indSim->UpdateSim(); 
    rflSim->UpdateSim(); 
}

void CCMCopulaModel::Simulation::setRandoms(IMCRandomSP rand)
{
	randoms = rand->getRandomNumbers()[0];

	// also, set randoms of inner copulas
/*
	rflSim->setRandoms(
				randoms);

	indSim->setRandoms(
				randoms + 
				rflSim->getNumberRands());

	dcSim->setRandoms(
				randoms + 
				rflSim->getNumberRands() + 
				indSim->getNumberRands());
*/
	dcSim->setRandoms(
				randoms);

	indSim->setRandoms(
				randoms + 
				dcSim->getNumberRands());

	rflSim->setRandoms(
				randoms + 
				dcSim->getNumberRands() + 
				indSim->getNumberRands());

}

void
CCMCopulaModel::updateBetas(const vector<double>& betas)
{
	const string method = "CCMCopulaModel::updateBetas";

	DoubleArray localBetas;

	try
	{
		int i;
		for (i = 0; i < this->getNumAssets(); i++)
		{
			SingleCreditAssetConstSP namei = 
				SingleCreditAssetConstSP::dynamicCast(names[i]->getSP());

			CreditEngineParametersConstSP 
				cep(namei->getEngineParams(CcmOnlyParameters::TYPE));

			const CcmOnlyParameters *ccmp 
				= dynamic_cast<const CcmOnlyParameters *>(cep.get());

			double beta = betas[i];

			beta = beta + (1-beta)*ccmp->getBetaTweak();

			localBetas.push_back(beta);
		}

		rflCopula.updateBetas(localBetas);
	}
	catch (exception &e)
	{
		throw ModelException(e, method);
	}
}

void
CCMCopulaModel::populateInnerCopulas()
{
	const string method = "CCMCopulaModel::getMarket";

	try{
//first clear the existing containers - check for memory leaks
		
		defaultRatesArray->clear();

		depenCopula.hasDefaulted.clear();
		depenCopula.defaultDates.clear();
		depenCopula.recoveries.clear();

		depenCopula.defaultRatesArray->clear();

		indepCopula.hasDefaulted.clear();
		indepCopula.defaultDates.clear();
		indepCopula.recoveries.clear();
		indepCopula.defaultRatesArray->clear();

		rflCopula.hasDefaulted.clear();
		rflCopula.defaultDates.clear();
		rflCopula.recoveries.clear();
		rflCopula.defaultRatesArray->clear();

		ratios.clear();
		cataRecFactor.clear();
		betaTweak.clear();

		//populate the recoveries of the CCMCopula with the market recovery 
		// rate of the names
		int numAssets = getNumAssets();
		recoveries.resize(numAssets);
		depenCopula.recoveries.resize(numAssets);

		//then populate with new ones
		for (int i = 0; i < numAssets; i++)
		{
			SingleCreditAssetConstSP namei;
			DefaultRatesSP nameDefaultRates;
			try
			{
				namei = 
					SingleCreditAssetConstSP::dynamicCast(names[i]->getSP());

				nameDefaultRates = namei->cdsParSpreads->defaultRates();
			}
			catch(exception& e)
			{
				stringstream message;
				message << "Unable to compute the defaultRates";
				message << " object for the name " << namei->getName();

				throw ModelException(e, method, message.str());
			}
			catch(...)
			{
				stringstream message;
				message << "Unable to compute the defaultRates object";
				message << " for the name " << namei->getName();
				throw ModelException(method, message.str());
			}

			defaultRatesArray->push_back(
				CleanSpreadCurve::convertFromDefaultRates(nameDefaultRates));
			
			recoveries[i].clear();
			recoveries[i].push_back(namei->recoveryRate());
			depenCopula.hasDefaulted.push_back(namei->hasDefaulted());
			depenCopula.defaultDates.push_back(namei->defaultDate());

			CreditEngineParametersConstSP 
				cep(namei->getEngineParams(CcmOnlyParameters::TYPE));

			const CcmOnlyParameters *ccmp = 
				dynamic_cast<const CcmOnlyParameters *>(cep.get());

			if (!ccmp)
			{
				throw ModelException("Failed to obtain CcmOnlyParameters");
			};

			ratios.push_back(ccmp->getIndependenceFactor());

			cataRecFactor.push_back(ccmp->getCatastrophicRecoveryFactor());

			betaTweak.push_back(ccmp->getBetaTweak());

			depenCopula.recoveries[i].
				push_back(namei->recoveryRate()
				*ccmp->getCatastrophicRecoveryFactor());

			if (ccmp->seniorCurve.get())
			{
				DefaultRatesSP defaultRates = floorCurves[i]->defaultRates();

				CleanSpreadCurveSP cleanCurve = 
					CleanSpreadCurve::convertFromDefaultRates(defaultRates);

				depenCopula.defaultRatesArray->push_back(cleanCurve);
			}
			else
			{
				depenCopula.defaultRatesArray->
					push_back(CleanSpreadCurveSP(0));
			}
		};

		indepCopula.hasDefaulted = 
			rflCopula.hasDefaulted = depenCopula.hasDefaulted;

		indepCopula.defaultDates = 
			rflCopula.defaultDates = depenCopula.defaultDates;
	}
	catch (exception &e)
	{
		throw ModelException(e, method);
	}
}


// split default probability among the 3 copulas for each name
// once we know the dates we want to get right (in terms of default prob),
// we create curves with exactly those dates for the individual copulas
void CCMCopulaModel::ComputeCurves(const DateTimeArray& dates)
{
	const string method = "CCMCopulaModel::ComputeCurves";

	try
	{
		int numAssets = this->getNumAssets();
		rflCopula.recoveries.resize(numAssets);
		indepCopula.recoveries.resize(numAssets);

		for (int i=0; i< numAssets; ++i)
		{
			//rflCopula.recoveries[i].resize(dates.size(), 1.); 
			//by default loss = 0 //revisit
			rflCopula.recoveries[i].resize(dates.size(), recoveries[i][0]);
		}

		for (int nameIdx = 0; nameIdx < numAssets; ++nameIdx)
		{
			DoubleArray indepProbArray, rflProbArray;

			int size = defaultRatesArray->size();

			DateTimeArray usedDates;

			//survival probs
			double sProb, dependentProb, indepenProb, rflProb; 

			sProb = dependentProb = indepenProb = rflProb = 1;

			//survival probs
			double lastSProb, lastDependentProb, lastIndepenProb, lastRflProb; 

			lastSProb  = lastDependentProb = lastIndepenProb = lastRflProb = 1;

			for (int point = 0; point<dates.size(); ++point)
			{
				DateTime  date = dates[point];

				if (date > valueDate)
				{
					usedDates.push_back(date);


					try
					{
						sProb = (*defaultRatesArray)[nameIdx]->
									getDefaultPV(valueDate, date);
					}
					catch(exception& e)
					{
						stringstream message;
						message << 
							"Unable to compute the survival probability for " <<
							"the name " << nameIdx << " for the date " << 
							date.toString();

						throw ModelException(e, method, message.str());
					}
					catch(...)
					{
						stringstream message;
						message << 
							"Unable to compute survival probability for name "
							<< nameIdx << " for the date " << date.toString();
						throw ModelException(method, message.str());
					}

					if ((*depenCopula.defaultRatesArray)[nameIdx].get())
					{
						try
						{
							dependentProb = 
								(*depenCopula.defaultRatesArray)[nameIdx]->
								getDefaultPV(valueDate, date);
						}
						catch(exception& e)
						{
							stringstream message;
							message << 
								"Unable to compute survival probability for" <<
								" dependent copula curve " << nameIdx << 
								" for date " << date.toString();

							throw ModelException(e, method, message.str());
						}
						catch(...)
						{
							stringstream message;
							message << "Unable to compute survival probability"<<
								" for dependent copula curve " << nameIdx << 
								" for date " << date.toString();

							throw ModelException(method, message.str());
						}
					}
					else
						dependentProb = 1;

					indepenProb = exp( ratios[nameIdx] * 
									log(sProb/dependentProb));

					if (indepenProb > 1.)
						throw ModelException(method, 
							string("CCM Calibration Error : ")+
							string("Independent copula survival probs > 1 "));

					if (indepProbArray.size() > 0)
					{
						if (indepenProb > indepProbArray.back())
							throw ModelException(method, 
								string("CCM Calibration Error: Independent ") +
							    string("copula survival probs are increasing"));
					}
					indepProbArray.push_back(indepenProb);

					rflProb = exp( (1. - ratios[nameIdx]) * 
						log(sProb/dependentProb) );

					if (rflProb > 1.)
						throw ModelException(method, 
							string("CCM Calibration Error : RFL probs > 1 "));

					if (rflProbArray.size() > 0)
					{
						if (rflProb > rflProbArray.back())
							throw ModelException(method, 
								string("CCM Calibration Error: RFL copula ")+
								string("survival probs are increasing"));
					}

					rflProbArray.push_back(rflProb);

					//compute non-catastrophic recovery and check for 
					// legitimate bounds

					double lastNonCatasProb = lastSProb/lastDependentProb;

					double nonCatasProb =sProb/dependentProb;

					//because this appears in the denominator
					if ((lastNonCatasProb - nonCatasProb)< CCM_MC_TINY)
					{
						stringstream message;
						message << "At the date " << date.toString() << 
							" for the name idx " <<  nameIdx << 
							" CCM recovery calibration fails " ;

						throw ModelException(method, message.str());
					}

					//semi-analytical case - not
					//double nonCatastrophicRecovery = 1. + 
					// ((1.-dependentProb)*
					// (1.-depenCopula.recoveries[nameIdx][0]) - 
					// (1.-sProb)*(1. - recoveries[nameIdx][0]))/
					// (dependentProb - sProb);

					double nonCatastrophicRecovery = 
						(1 - recoveries[nameIdx][0])*(lastSProb - sProb)
						- (1- depenCopula.recoveries[nameIdx][0]) * 
						(lastDependentProb - dependentProb)*lastNonCatasProb;

					nonCatastrophicRecovery = 1 - 
						nonCatastrophicRecovery /
							(dependentProb*(lastNonCatasProb - nonCatasProb));

					if (fabs(nonCatastrophicRecovery) < CCM_MC_TINY) 
						nonCatastrophicRecovery = 0;

					if ((nonCatastrophicRecovery > 1) && 
						( nonCatastrophicRecovery < (1 + CCM_MC_TINY))) 
						
						nonCatastrophicRecovery = 1;

					if ((nonCatastrophicRecovery > 1.) || 
						 (nonCatastrophicRecovery < 0))
					{
						stringstream message;
						message << "At date " << date.toString() << 
							    " for name idx " <<  nameIdx << 
								" nonCatastrophicRecovery is " 
								<< nonCatastrophicRecovery;

						message << 
							" Therefore, CCM recovery calibration fails " ;

						throw ModelException(method, message.str());
					}
					else
					{
						rflCopula.recoveries[nameIdx][point] = 
							nonCatastrophicRecovery;  
					}

					lastSProb = sProb;
					lastDependentProb = dependentProb;
					lastIndepenProb = indepenProb;
					lastRflProb = rflProb;
				}
			}

			indepCopula.recoveries = rflCopula.recoveries;

			//convert from spots to forwards
			indepCopula.defaultRatesArray->push_back(
				CleanSpreadCurve::cleanSpreadCurveFromDF(
							valueDate, 
							usedDates,
							indepProbArray
							));

			rflCopula.defaultRatesArray->push_back(
				CleanSpreadCurve::cleanSpreadCurveFromDF(
							valueDate, 
							usedDates,
							rflProbArray
							));
		}
	}
	catch(exception& e)
	{
		throw ModelException(e, method);
	}
};


static int min(int x, int y, int z)
{
	return Maths::min(x, Maths::min(y,z));
}


void
CCMCopulaModel::Simulation::computeSample(long idx)
{

	dcSim->computeSample(idx);
	indSim->computeSample(idx);
	rflSim->computeSample(idx);

	// compute member default info using min of all default times
	for (int i = 0; i < getNumAssets(); i++)
	{
		simState.setRecovery(i, 1); //ie no loss

		simState.setDefaultIndex(i, 
			min(dcSim->simulationState().getDefaultIndices()[i],
			    indSim->simulationState().getDefaultIndices()[i],
				rflSim->simulationState().getDefaultIndices()[i])
		);

		if (simState.getDefaultIndices()[i] < timeLine.size()-1)
		{
			if (simState.getDefaultIndices()[i] == 
				 dcSim->simulationState().getDefaultIndices()[i])   
				 //default due to catastrophic default
				simState.setRecovery(i,
					dcSim->simulationState().getRecoveries()[i]);
			else
				 //rfl and ind recovery rate would be the same
				simState.setRecovery(i, 
					rflSim->simulationState().getRecoveries()[i]);
		}
	}
}


DRLIB_END_NAMESPACE
