//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : IndependentCopulaModel.hpp
//
//   Description : Implementation of for class IndependentCopulaModel, a class 
//                 for the independence copula in the context of MC CCM
//
//   Date        : Feb 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/IndependentCopulaModel.hpp"
#include "edginc/SingleCreditAsset.hpp"

DRLIB_BEGIN_NAMESPACE



COPULA_MODEL
IndependentCopulaModel::type() const
{
	return DEPENDENT;
}


void 
IndependentCopulaModel::Simulation::UpdateSim()
{
	BaseSimulation::UpdateSim();

	for (int i = 0; i < model->getNumAssets(); i++)
	{
		for (int j = 0; j < timeLine.size()-1; j++)
		{
			zScoresInv[i][j] = (probabilities[i][j] > TINY) ? 
				imsl_d_normal_inverse_cdf(probabilities[i][j]) : ZERO_EVENT;
		};
	};
};

void
IndependentCopulaModel::Simulation::computeSample(long idx)
{
	int i, j;

	const IndependentCopulaModel* mod = 
					dynamic_cast<const IndependentCopulaModel*>(model);

	// set all names to defaulted at numDates()-1
	// number of dates includes value date or a date in the past

	for (i = 0; i < mod->getNumAssets(); i++)
	{
		simState.setDefaultIndex(i, timeLine.size()-1); 
		//simState.setRecoveries()[i] = 1; //corresponding to no loss; //revisit
		simState.setRecovery(i, mod->recoveries[i][0]);
	}

	for (i = 0; i < mod->getNumAssets(); i++)
	{
		if (mod->hasDefaulted[i])
		{
			for (j = 0; j < timeLine.size(); j++)
			{
				if ( timeLine[j] > mod->defaultDates[i] )
				{	
					simState.setDefaultIndex(i, j-1);
					simState.setRecovery(i, mod->recoveries[i][j]);
					break;
				}
			}
		} else
		{
			for (j = 1; j < timeLine.size(); j++)
			{
				if ( (simState.getDefaultIndices()[i] >= j) && (randoms[i] > zScoresInv[i][j-1]) )
				{	
					simState.setDefaultIndex(i, j-1);
					simState.setRecovery(i, mod->recoveries[i][j]);
					break;
				}
			}
		}
	}
}

IndependentCopulaModel::Simulation::Simulation(
	const DateTime& valueDate,
	CopulaModelBase &model, 
	DateTimeArray timeLine,
	long nSamples) :	
		BaseSimulation(valueDate, model, timeLine, nSamples)
{

	// the baseSimulation constructor computes default probs to timeline points


	// here, we would do any initialization of the random number generation

	randsPerPath = model.getNumAssets();
	zScoresInv = CDoubleMatrix(model.getNumAssets(), timeLine.size());
}


DRLIB_END_NAMESPACE





