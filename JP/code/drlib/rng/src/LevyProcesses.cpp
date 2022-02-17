#include "edginc/LevyProcesses.h"

CORE_BEGIN_NAMESPACE


BrownianMotion::BrownianMotion(double *timesBM, int timesBMsize, long seed)
{
	setParameters(timesBM,timesBMsize,seed);
}	

BrownianMotion::BrownianMotion(BrownianMotion &bm)
{
	std::vector<double> time;
	bm.getDiscretization(time);
	setParameters(&time[0],time.size(),bm.getSeed());
}


void BrownianMotion::setParameters(double *timesBM, int timesBMsize, long seed)
{
	m_t.resize(timesBMsize-1);
	m_sqrtt.resize(timesBMsize-1);
	for (unsigned int k=0; k<m_t.size(); k++)
	{
		m_t[k] = timesBM[k+1] - timesBM[k];
		m_sqrtt[k] = sqrt(m_t[k]);
	}
	m_BM.initializeObject(seed,m_t.size());
}

void BrownianMotion::GetBMincrements(double *bmIncrements)
{
	double *BMI = m_BM.getVector();
	for (unsigned int k=0; k<m_t.size(); k++) bmIncrements[k] = m_sqrtt[k]*BMI[k];
}


//////////////////////////////////////////////////////////////////////
// Poisson Process
//////////////////////////////////////////////////////////////////////


PoissonProcess::PoissonProcess(double *freq, int freqSize, double maturity, long seed)
{
	setParameters(freq,freqSize,maturity,seed);
}

PoissonProcess::PoissonProcess(PoissonProcess &pp)
{
	std::vector<double> freq;
	pp.getFrequency(freq);
	setParameters( &freq[0], freq.size(), pp.getMaturity(),pp.getSeed());
}


void PoissonProcess::setParameters(double *freq, int freqSize, double maturity, long seed)
{
	m_maturity = maturity;
	m_sumFreq = 0;
	m_freq.resize(freqSize); 
	for (int i=0; i<freqSize; i++)
	{
		m_freq[i]=freq[i];
		m_sumFreq += freq[i];
	}
	m_poisson.initializeObject(seed,1);
}



int PoissonProcess::GetPoissonNbJumps(bool cond,int *nbJump)
{
	return GetPoissonNbJumps(cond, nbJump, *m_poisson.getVector());
}

int PoissonProcess::GetPoissonNbJumps(bool cond,int *nbJump, double unif) // if cond false, nb Jumps conditionned on at least one jump...
																		// unif is the uniform rv which will give the total nb of jumps
																		// (does not have to be random)
{
	double proba, Sproba;
	unsigned int i;
	int nbJumpSum;

	if (cond==true)
	{
		nbJumpSum = 1;
		proba = Sproba = m_sumFreq * m_maturity / (exp(m_sumFreq*m_maturity) - 1.0);
	}
	else
	{
		nbJumpSum = 0;
		proba = Sproba = exp(-m_sumFreq*m_maturity);
	}
	
	while ((unif>Sproba)&&(Sproba<1-1e-20)) // second test to make sure we do not spend the rest of our life in this loop
	{
		nbJumpSum++;
		proba *= m_sumFreq * m_maturity / nbJumpSum;
		Sproba += proba;
	}; 
    int answer = nbJumpSum;
	// we have simulated the nb Jumps of the sum of the Poisson processes...
	// now we need to decide which one has it
	for (i=0; i<m_freq.size(); i++) nbJump[i]=0;

	double remainingFreq = m_sumFreq, alpha;
	for (i=0; (i<m_freq.size())&&(nbJumpSum>0); i++)
	{
		alpha = m_freq[i] / remainingFreq;
		remainingFreq -= m_freq[i];
		if (alpha>=1-1e-9) 
		{
			nbJump[i] = nbJumpSum;
			nbJumpSum = 0;
		}
		else
		{
			unif = *m_poisson.getVector();
			proba = Sproba = pow(1-alpha,nbJumpSum);
			while ((unif>Sproba)&&(nbJump[i]<nbJumpSum))
			{
				proba *= alpha*(nbJumpSum-nbJump[i])/((1-alpha)*(nbJump[i]+1));
				Sproba += proba;
				nbJump[i]++;
			}
			nbJumpSum -= nbJump[i];
		}
	}
	return answer;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////



CORE_END_NAMESPACE
