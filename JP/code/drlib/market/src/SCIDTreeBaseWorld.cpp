#include "edginc/config.hpp"
#include "edginc/SCIDTreeBaseWorld.hpp"

DRLIB_BEGIN_NAMESPACE


void TreeBaseWorld::InitializeRNG(long seed) 
{
	rngUniformFast.initializeObject(seed,1); 
	rngNormalFast.initializeObject(seed+1,1); 
	IpathFast = 0; 
	rngUniformFull.initializeObject(seed,1); 
	rngNormalFull.initializeObject(seed+1,1); 
	IpathFull = 0; 
};

double TreeBaseWorld::getRandomNb(bool fast, bool uniform)
{
	Sequence *seq;
	if (fast)
	{
		if (uniform) seq = &rngUniformFast;
			else seq = &rngNormalFast;
		seq->populateVector(IpathFast);
		IpathFast++;
		return *(seq->getVector());
	}
	else
	{
		if (uniform) seq = &rngUniformFull;
		else seq = &rngNormalFull;
		seq->populateVector(IpathFast);
		IpathFast++;
		return *(seq->getVector());
	}

}
UniformRandomSequence rngUniformFast, rngUniformFull;
GaussianRandomSequence rngNormalFast, rngNormalFull;
long IpathFast, IpathFull;




void TreeBaseWorld::setDiscretizationTimesCF(DoubleArray &_times) 
{ 
	cfCIR->setDiscretizationTimes(_times);
}
void TreeBaseWorld::setDiscretizationTimesCF(double lastTime, double delta) 
{ 
	cfCIR->setDiscretizationTimes(lastTime,delta);
}

void TreeBaseWorld::setDiscretizationTimesIDIO(DoubleArray &_times) 
{ 
	for (size_t i=0; i<idioCIR.size(); i++) 
		idioCIR[i]->setDiscretizationTimes(_times); 
}

void TreeBaseWorld::setDiscretizationTimesIDIO(double lastTime, double delta) 
{ 
	for (size_t i=0; i<idioCIR.size(); i++) 
		idioCIR[i]->setDiscretizationTimes(lastTime,delta);
}

void TreeBaseWorld::setUnderlyingProcesses(const int nbNames, 
										   const double volDiff, 
										   const double decayDiff,
										   const DoubleArray &jumpDecay,
										   const DoubleArray &jumpImpact, 
										   const DoubleArray &jumpFreq,
										   const DoubleArray &jumpRatios,
										   const DoubleArrayArrayArray &jumpMeans,
										   const DoubleArrayArrayArray &jumpWeights)
{
	static const string method = "TreeBaseWorld::setUnderlyingProcesses";
	try
	{
		QLIB_VERIFY(jumpDecay.size()==jumpImpact.size() 
				&& jumpImpact.size()==jumpFreq.size() 
				&& jumpImpact.size()==jumpRatios.size(), "wrong jump input");
		if (nbNames>0)
		{
			IntArray trueJumpMarkets;
			for (int k=0; k<jumpDecay.size(); ++k)
			{
				if (jumpDecay[k]*jumpFreq[k]*jumpImpact[k]*jumpRatios[k]>1e-15)
					trueJumpMarkets.push_back(k);
			}
			nbMarkets = trueJumpMarkets.size();
			DoubleArray freq(nbMarkets);
			for (int k=0; k<nbMarkets; ++k)
			{
				int index = trueJumpMarkets[k];
				freq[k] = jumpFreq[index];
			}
			poissonDiffuseFast.Initialize(freq);
			poissonDiffuseFull.Initialize(freq);

			hasIdVol = (volDiff>0);
			cfCIR = CIRAffineSP(new CIRAffine(1.0,volDiff,decayDiff,1.0));
			idioCIR.resize(nbNames);
			jumpMarket.resize(nbNames);
			
			bool parametricDist = (jumpMeans.size()==jumpFreq.size()) && (jumpWeights.size()==jumpFreq.size());
			for (int m=0; m<nbNames; m++)
			{
				idioCIR[m] = CIRAffineSP(new CIRAffine(1.0,volDiff,decayDiff,1.0));

				jumpMarket[m].resize(nbMarkets);
				for (int k=0; k<nbMarkets; ++k)
				{
					int index = trueJumpMarkets[k];
					try
					{
						if (parametricDist)
						{
							DoubleArraySP pieceT = DoubleArraySP (new DoubleArray(jumpMeans[index].size()));
							DoubleArrayArray jm = jumpMeans[index], jw = jumpWeights[index];

							NonParametricDistributionSP jumpDist = NonParametricDistributionSP(
								new NonParametricDistribution(pieceT ,jm,jw));
							jumpMarket[m][k] = JumpAffineSP(new JumpAffine(0.0,jumpDecay[index], jumpFreq[index], jumpDist));
						}
						else
						{
							ExponentialJumpDistributionSP jumpDist = ExponentialJumpDistributionSP(new ExponentialJumpDistribution(1.0));
							jumpMarket[m][k] = JumpAffineSP(new JumpAffine(0.0,jumpDecay[index], jumpFreq[index], jumpDist));
						}
					}
					catch (exception& e)
					{
						throw ModelException(e, "error in definition of the jump distribution");
					}
					// NO IMPACTS FOR NOW!!!

				}
			}
		}
	}
	catch (exception& e)
	{
		throw ModelException(e, method);
	}
}

void TreeBaseWorld::CalibrateSingleNames(DateTime today,
										DateTimeArray &singleNameDates,
										DoubleMatrix &survProbas,  // size (nbNames,singleNameDates.size())
										double idioRatio,
										double cfRatio,
										DoubleArray &jumpRatios)
{
	static const string method = "TreeBaseWorld::CalibrateSingleNames";
	try
	{
		QLIB_VERIFY(survProbas.numCols()==idioCIR.size() && survProbas.numRows()==singleNameDates.size(), "wrong number of survival probabilities");
		//	QLIB_VERIFY(ratios.size()== "" NB JUMPS "", "not enough ratios provided");
		QLIB_VERIFY(idioRatio>=0, "negative idiosyncratic ratio!");
		QLIB_VERIFY(cfRatio>=0,   "negative common factor ratio!");
		QLIB_VERIFY(jumpRatios.size()>=nbMarkets, "jumpRatios of the wrong size!");
		for (int i=0; i<jumpRatios.size(); i++) QLIB_VERIFY(jumpRatios[i] >=0, "negative Jump ratio!");

		double sumRatios = idioRatio + cfRatio;
		for (int i=0; i<nbMarkets; i++) sumRatios += jumpRatios[i];
		QLIB_VERIFY(abs(sumRatios-1.0)<1e-10, "negative common factor ratio!");
		hasCF = (cfRatio>0);

		DoubleArray singleNameTimes(singleNameDates.size());
		for (int k=0; k<singleNameDates.size(); k++) 
			singleNameTimes[k] = today.yearFrac(singleNameDates[k]);
		int nbNames = survProbas.numCols();
		DoubleArray survProbaForIdio(singleNameTimes);
		DoubleArray y;
		DoubleArray meanOneTime;
		DoubleArrayArrayArray means;

		for (int m=0; m<nbNames; m++)
		{
			fill(survProbaForIdio.begin(),survProbaForIdio.end(),0.0);
			for (int k=0; k<nbMarkets; ++k)
			{
				jumpMarket[m][k]->calibrate(singleNameTimes, survProbas[m], jumpRatios[k]); 
				for (int i=0; i<singleNameTimes.size(); ++i)
					survProbaForIdio[i] += jumpMarket[m][k]->minusLogSurvProba(singleNameTimes[i]);
			}
			if (hasCF) 
			{
				cfCIR->setWhichg(m,false);
				cfCIR->calibrate(singleNameTimes, survProbas[m], cfRatio);
				for (int i=0; i<singleNameTimes.size(); ++i)
					survProbaForIdio[i] += cfCIR->minusLogSurvProba(singleNameTimes[i]);
			}
			for (int i=0; i<singleNameTimes.size(); ++i)
				survProbaForIdio[i] = survProbas[m][i]/exp(-survProbaForIdio[i]);
			idioCIR[m]->calibrate(singleNameTimes, &survProbaForIdio[0], 1.0);
		}
	}
	catch (exception& e)
	{
		throw ModelException(e, method);
	}
}

void TreeBaseWorld::SimulateCFFast(double futureTime, double lambdaFutureTime)
{
	if (hasCF)
		cfCIR->diffuse(&rngNormalFast,IpathFast,futureTime,lambdaFutureTime);
}

double TreeBaseWorld::SimulateJumpsFast(double lastTime,
									  bool conditional, 
									  double unif,
									  double futureTime)
{
	if (nbMarkets>0)
	{
		if ((unif<=0)||(unif>=1))
		{
			rngUniformFast.populateVector(IpathFast);
			IpathFast++;
			unif = *(rngUniformFast.getVector());
		}
		poissonDiffuseFast.GetPoissonNbJumps(lastTime, conditional, unif, &rngUniformFast, IpathFast, futureTime);
		for (int k=0; k<nbMarkets; ++k)
			poissonDiffuseFast.simulateJumpTimeCondNbJumps(lastTime, k, &rngUniformFast, IpathFast, futureTime);
		for (size_t m=0; m<idioCIR.size(); ++m)
			for (int k=0; k<nbMarkets; ++k)
				jumpMarket[m][k]->setSimulatedJumpTimes(poissonDiffuseFast.jumpTimes[k],
														poissonDiffuseFast.nbJumps[k]);
		return poissonDiffuseFast.getFirstJump();
	}
	return 0.0;
}

void TreeBaseWorld::SimulateJumpsFull(double lastTime,
									  bool conditional)
{
	if (nbMarkets>0)
	{
		rngUniformFull.populateVector(IpathFull);
		IpathFull++;
		double unif = *(rngUniformFull.getVector());
		poissonDiffuseFull.GetPoissonNbJumps(lastTime,conditional, unif, &rngUniformFull, IpathFull);
		for (int i=0; i<nbMarkets; ++i)
			poissonDiffuseFull.simulateJumpTimeCondNbJumps(lastTime,i,&rngUniformFull, IpathFull);
		for (size_t m=0; m<idioCIR.size(); ++m)
			for (int k=0; k<nbMarkets; ++k)
				jumpMarket[m][k]->setSimulatedJumpTimes(poissonDiffuseFull.jumpTimes[k],
														poissonDiffuseFull.nbJumps[k]);
	}
}

void TreeBaseWorld::SimulateSpreadFull(baseWorldSpread &bws)
{
	if (hasIdVol)
		for (size_t m=0; m<idioCIRFull.size(); ++m)
			idioCIRFull[m]->diffuse(&rngNormalFull, IpathFull);
	for (size_t m=0; m<idioCIRFull.size(); ++m)
		idioCIRFull[m]->readLambdaT(bws.intensity[m]);

	if (hasCF)
	{
		DoubleArray cfspread(bws.nbTimes);
		cfCIRFull->diffuse(&rngNormalFull,IpathFull);
		for (size_t m=0; m<idioCIRFull.size(); ++m)
		{
			cfCIRFull->setWhichg(m);
			cfCIRFull->readLambdaT(cfspread);
			for (int i=0; i<bws.nbTimes; ++i)
				bws.intensity[m][i] += cfspread[i];
		}
	}
	if (nbMarkets>0)
	{
		double firstJump = poissonDiffuseFull.getFirstJump();
		int firstJumpIndex = lower_bound(bws.discTimesSpread.begin(),bws.discTimesSpread.end(),firstJump) 
			- bws.discTimesSpread.begin();
		for (int k=0; k<nbMarkets; ++k)
			for (size_t m=0; m<idioCIRFull.size(); ++m)
			{
				jumpMarketFull[m][k]->simulateJumpSizes(&rngUniformFull,IpathFull);
				for (int i=firstJumpIndex; i<bws.nbTimes; ++i)
					bws.intensity[m][i] += jumpMarketFull[m][k]->readLambdaT(bws.discTimesSpread[i]);
			}
	}
}


void TreeBaseWorld::setWorlds(vector<LinearInterpolantSP> _f)
{
	f=_f;
}



/* retrieve unconditional survival probabilities **/
void TreeBaseWorld::minusLogSurvProba( 
			   int world,  // x axis of world assumed to be the same as cfF and idioF x axis
			   int name,
			   const DoubleArray &bigT,
			   double *minusLogSurvProb,  // (O) must be of size bigT
			   bool add,
			   double futureTime,
			   double lambdaFutureTimeIdio,
			   double lambdaFutureTimeCf,
			   double *lambdaFutureTimeJump)
{
	QLIB_VERIFY( futureTime == 0 || lambdaFutureTimeJump!=0, "lambdaFutureTimeJump should not be NULL");
	minusLogSurvProbaIdio( world, name, bigT, minusLogSurvProb, add, futureTime, lambdaFutureTimeIdio);
	if (hasCF) 
	{
		cfCIR->setWhichg(name);
		cfCIR->setFs(f[world]);
		for (int i=0; i<bigT.size(); ++i)
		{
			double logsp = cfCIR->minusLogSurvProba(bigT[i],futureTime,lambdaFutureTimeIdio);
			minusLogSurvProb[i]+=logsp;
		}
	}
	for (int market=0; market<nbMarkets; ++market)
	{
		jumpMarket[name][market]->setFs(f[world]);
		double lambdaFutureTimeJumpMarket;
		if (futureTime!=0) lambdaFutureTimeJumpMarket = lambdaFutureTimeJump[market];
		for (int i=0; i<bigT.size(); ++i)
		{
			double logsp = jumpMarket[name][market]->minusLogSurvProba(bigT[i],futureTime,lambdaFutureTimeJumpMarket);
			minusLogSurvProb[i]+=logsp;
		}
	}
};


/* retrieve unconditional idiosyncratic survival probabilities **/
void TreeBaseWorld::minusLogSurvProbaIdio( 
				   int world,
				   int name,
				   const DoubleArray &bigT,
				   double *minusLogSurvProb,  // (O) must be of size bigT
				   bool add,
				   double futureTime,
				   double lambdaFutureTimeIdio)
{
	idioCIR[name]->setFs(f[world]);
	for (int i=0; i<bigT.size(); ++i)
	{
		double logsp = idioCIR[name]->minusLogSurvProba(bigT[i],futureTime,lambdaFutureTimeIdio);
		if (add) minusLogSurvProb[i]+=logsp;
			else minusLogSurvProb[i]=logsp;
	}
};

/* retrieve stochastic survival probability for common factor **/
void TreeBaseWorld::minusLogCondSurvProbaCF(
					 int world,
					 int name,
					 const DoubleArray &bigT,
					 double *minusLogSurvProb,  // (O) must be of size bigT
					 bool add,
					 double futureTime)
{
	if (hasCF)
	{
		cfCIR->setFs(f[world]);
		cfCIR->setWhichg(name);
		cfCIR->read_F_MinusLogSurvProba(futureTime, bigT, minusLogSurvProb, add);
	}
};


void TreeBaseWorld::minusLogCondSurvProbaJump(
							int world,
							int name,
							const DoubleArray &bigT,
							double *minusLogSurvProb,  // (O) must be of size bigT
							bool fastMC,
							bool add,
							double futureTime,
							double lambdaFutureTimeJumpMarket)
{
	if (!add)
		for (int i=0; i<bigT.size(); ++i)minusLogSurvProb[i]=0;
	for (int market=0; market<nbMarkets; ++market)
	{
		jumpMarket[name][market]->setFs(f[world]);
		for (int i=0; i<bigT.size(); ++i)
			minusLogSurvProb[i]+=jumpMarket[name][market]->read_F_MinusLogSurvProba(!fastMC, bigT[i],
																					futureTime,
																					lambdaFutureTimeJumpMarket);
	}
};

DRLIB_END_NAMESPACE

