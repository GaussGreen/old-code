//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : IndependentCopulaModel.hpp
//
//   Description : Implementation for class RFLCopulaModel, a class 
//                 implementing the rfl copula in the context of 
//				   MC CCM
//
//   Date        : Feb 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/RFLCopulaModel.hpp"

DRLIB_BEGIN_NAMESPACE

void 
RFLCopulaModel::Simulation::computeSample(long idx)
{
	int i, j;
	const RFLCopulaModel* mod = dynamic_cast<const RFLCopulaModel*>(model);

	// set all names to defaulted at numDates()
	for (i = 0; i < model->getNumAssets(); i++)
	{
		simState.setDefaultIndex(i, timeLine.size()-1);
		//simState.setRecoveries()[i] = 1; //corresponding to no loss; //revisit
		simState.setRecovery(i, mod->recoveries[i][0]);
	}

	for (i = 0; i < model->getNumAssets(); i++)
	{
		double beta = mod->betas[i];
		double betaInv = sqrt(1-beta*beta);
		
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
		} 
		else
		{
			double rand = randoms[i+1]*betaInv + randoms[0]*beta;
			
			// first date is either value date or in the past
			for (j = 1; j < timeLine.size(); j++) 
			{
				if ( (simState.getDefaultIndices()[i] >= j) 
					&& (rand > zScoresInv[i][j-1]) )
				{	
					simState.setDefaultIndex(i, j-1);
					simState.setRecovery(i, mod->recoveries[i][j]);
					break;
				}
			}
		}
	}
}

void
RFLCopulaModel::getMarket(const CModel* model, const CMarketDataSP market)
{
	int i;
	for (i=0; i < getNumAssets(); ++i)
	{
		ICDSParSpreadsWrapper temp;
		ICDSParSpreads::getMarketData(model, market.get(), temp );
	}
}

void RFLCopulaModel::Simulation::UpdateSim()
{
	BaseSimulation::UpdateSim();
	
	for (int i = 0; i < model->getNumAssets(); i++)
	{
		for (int j = 0; j < timeLine.size()-1; j++)
		{
			if ( ( probabilities[i][j] < 1. - TINY) 
				&& ( probabilities[i][j] > TINY) )

				zScoresInv[i][j] = 
					imsl_d_normal_inverse_cdf(probabilities[i][j]);

			else 
				zScoresInv[i][j] = 
					(probabilities[i][j] < 1. - TINY) ? 
						- ZERO_EVENT : ZERO_EVENT;
		}
	}
}

RFLCopulaModel::Simulation::Simulation(
	const DateTime& valueDate,
	CopulaModelBase &mod, 
	DateTimeArray timeLine,
	int nSamples) 
	:
	BaseSimulation(valueDate, mod, timeLine, nSamples)
{

	const RFLCopulaModel* mod1 = dynamic_cast<const RFLCopulaModel*>(model);

	randsPerPath = mod1->getNumAssets() + 1;

	zScoresInv = CDoubleMatrix(mod1->getNumAssets(), timeLine.size());
}

DRLIB_END_NAMESPACE
