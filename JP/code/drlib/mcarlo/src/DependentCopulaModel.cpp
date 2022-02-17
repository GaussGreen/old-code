//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : IndependentCopulaModel.hpp
//
//   Description : Implementation for class DependentCopulaModel, a class 
//                 implementing the dependence copula in the context of 
//				   MC CCM
//
//   Date        : Feb 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/DependentCopulaModel.hpp"
#include "edginc/SingleCreditAsset.hpp"

DRLIB_BEGIN_NAMESPACE

COPULA_MODEL
DependentCopulaModel::type() const
{
	return DEPENDENT;
}

void 
DependentCopulaModel::Simulation::UpdateSim()
{
	BaseSimulation::UpdateSim();

	for (int i = 0; i < model->getNumAssets(); i++)
	{
		for (int j = 0; j < timeLine.size()-1; j++)
		{
			zScoresInv[i][j] 
				= (probabilities[i][j] > TINY) ? 
				imsl_d_normal_inverse_cdf(probabilities[i][j]) : ZERO_EVENT;
		};
	};
};

void
DependentCopulaModel::Simulation::computeSample(long idx)
{
	int i, j;
	const DependentCopulaModel* mod = dynamic_cast<const DependentCopulaModel*>(model);

	// set all names to defaulted at numDates()

	for (i = 0; i < mod->getNumAssets(); i++)
	{
		simState.setDefaultIndex(i, timeLine.size()-1);
		simState.setRecovery(i, mod->recoveries[i][0]); 
	};

	for (i = 0; i < mod->getNumAssets(); i++)
	{
		if (mod->hasDefaulted[i])
		{
			for (j = 0; j < timeLine.size(); j++)
			{
				if ( timeLine[j] > mod->defaultDates[i] )
				{	
					simState.setDefaultIndex(i,j-1); // default date is known, but in the future
					break;
				}
			}
		} else
		{
			for (j = 1; j <  timeLine.size(); j++)         // first point in timeLine is either valDate or in the past
			{
				if ((simState.getDefaultIndices()[i] >= j) &&  
					(randoms[0] > zScoresInv[i][j-1]) )
				{	
					simState.setDefaultIndex(i, j-1);
					break;
				}
			}
		}
	}
}

DependentCopulaModel::Simulation::Simulation(
	const DateTime& valueDate,
	CopulaModelBase &mod, 
	DateTimeArray timeLine,
	long nSamples) 
	:	BaseSimulation(valueDate, mod, timeLine, nSamples)
{
	// the baseSimulation constructor computes default probs to timeline points

	zScoresInv = CDoubleMatrix(mod.getNumAssets(), timeLine.size()-1);
	
	randsPerPath = 1;

}

DRLIB_END_NAMESPACE





