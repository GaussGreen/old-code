//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : BaseSimulation.cpp
//
//   Description : Declares the CopulaModelBase, data and functions common to
//                 all copulas in CCM
//
//   Date        : Feb 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/DefaultRates.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/imsl.h"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/BaseSimulation.hpp"

DRLIB_BEGIN_NAMESPACE


const DateTimeArray&
BaseSimulation::timeline() const
{
	return timeLine;
}

void BaseSimulation::SimulationState::setDefaultIndex(int i, int value)
{
	defIdx[i] = value;	
};


void BaseSimulation::SimulationState::setRecovery(int i, double value)
{
	recoveries[i] = value;	
};


BaseSimulation::SimulationState::~SimulationState()
{
};

BaseSimulation::BaseSimulation(
	const DateTime& valueDate,
	const ICopulaModelBase &mod,
	DateTimeArray  dates,
	long  nSamples) : 
	valueDate(valueDate),
	model(&mod),
	numSamples(nSamples),
	simState(mod.getNumAssets()),
	timeLine(dates)
{
	probabilities = CDoubleMatrix(mod.getNumAssets(), dates.size()-1);
}

void BaseSimulation::UpdateSim()
{
	int timelineSize = timeLine.size();
	for (int i = 0; i < model->getNumAssets(); i++)
	{
		if (model->nameHasDefaulted(i))
		{
			for (int j = 1; j < timelineSize; j++)
				probabilities[i][j-1] =	0;
		}
		else
		{
			if (model->getCurves()[i].get()) // JCP what is this handling?
			{
				for (int j = 1; j < timelineSize; j++)
				{
					if ( valueDate <= timeLine[j] )
					{
						probabilities[i][j-1] =	(model->getCurves()[i])->getDefaultPV(
															valueDate, timeLine[j]);
					} else {
						probabilities[i][j-1] = 1;
					}
				};
			}
			else
			{
				for (int j = 1; j < timelineSize; j++)
					probabilities[i][j-1] = 1;
			}
		}
	}
}

DRLIB_END_NAMESPACE

