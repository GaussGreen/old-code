/*************************************************************************
 LevyProcesses.h
 Author: N.Victoir, Credit QR
 
 Simulates incremement of Brownian Motion and nb Jumps of Poisson processes, possibly conditional on the sum of them
 having at least one jump. Used in the sCID quick fix project, should not be used elsewhere
 *************************************************************************/
#ifndef _LEVYPROCESSES__H_
#define _LEVYPROCESSES__H_


#include "edginc/UniformRandomSequence.h"
#include "edginc/GaussianRandomSequence.h"

CORE_BEGIN_NAMESPACE



class RNG_DLL Uniform
{
public:
	Uniform(){};
	Uniform(Uniform &u) { initializeObject(u.getSeed(), u.getDimension()); };
	Uniform(long seed, int dimension):m_urs(seed, dimension), m_iPath(0) {};
	void initializeObject(long seed, int dimension) { m_urs.initializeObject(seed, dimension); m_iPath=0; }
	double * getVector()
	{
		m_urs.populateVector(m_iPath);
		m_iPath++;
		return m_urs.getVector();
	}
	long getIpath(){ return m_iPath; };
	void setIpath(long ip) { m_iPath = ip; };
	long getSeed() { return m_urs.getSeed();}
	int getDimension() {return m_urs.getDimension(); };

private:
	UniformRandomSequence m_urs;
	long m_iPath;
};

class RNG_DLL Uniform1
{
public:
	Uniform1(){};
	void initializeObject(long seed) { int dim = 1; m_urs.initializeObject(seed, dim); m_iPath=0; }
	double getVector()
	{
		m_urs.populateVector(m_iPath);
		m_iPath++;
		return *m_urs.getVector();
	}
	long getIpath(){ return m_iPath; };
	void setIpath(long ip) { m_iPath = ip; };
	long getSeed() { return m_urs.getSeed();}
private:
	UniformRandomSequence m_urs;
	long m_iPath;
};

// just Gaussian random number, but deals with iPath on its own, the way I like it
class RNG_DLL Gaussian
{
public:
	Gaussian(){};
	Gaussian(long seed, int dimension):m_gauss(seed, dimension), m_iPath(0) {};
	void initializeObject(long seed, int dimension) { m_gauss.initializeObject(seed, dimension); m_iPath=0; }
	double * getVector()
	{
		m_gauss.populateVector(m_iPath);
		m_iPath++;
		return m_gauss.getVector();
	}
	long getIpath(){return m_iPath;};
	void setIpath(long ip) { m_iPath = ip; };
	long getSeed() { return m_gauss.getSeed();}
private:
	GaussianRandomSequence m_gauss;
	long m_iPath;
};

// as above, in dim 1
class RNG_DLL Normal
{
public:
	Normal(){};
	void initializeObject(long seed) { int dim=1; m_gauss.initializeObject(seed, dim); m_iPath=0; }
	double getRandomNumber()
	{
		m_gauss.populateVector(m_iPath);
		m_iPath++;
		return *m_gauss.getVector();
	}
	long getIpath(){return m_iPath;};
	void setIpath(long ip) { m_iPath = ip; };
	long getSeed() { return m_gauss.getSeed();}
private:
	GaussianRandomSequence m_gauss;
	long m_iPath;
};


// simulate the increments of Brownian Motion
class RNG_DLL BrownianMotion
{
public:
	BrownianMotion(){};
	BrownianMotion(BrownianMotion &bm);
	BrownianMotion(double *timesBM, int timesBMsize, long seed);

	void setParameters(double *timesBM, int timesBMsize, long seed);
	void GetBMincrements(double *bmIncrements);

	void setIpath(long ip) { m_BM.setIpath(ip); };
	long getIpath() { return m_BM.getIpath(); };
	long getSeed() {return m_BM.getSeed(); };
	void getDiscretization(std::vector<double> &time) { time = m_t; };
private:
	std::vector<double> m_t, m_sqrtt;
	Gaussian m_BM;
};


// simulate the jumps of Poisson process
// possibly conditionally on the sum of the poisson process being >0 at maturity

class RNG_DLL PoissonProcess // just number of jumps
{
public:
	PoissonProcess(){};
	PoissonProcess(PoissonProcess &pp);
	PoissonProcess(double *freq, int freqSize, double maturity, long seed);

	void setParameters(double *freq, int freqSize, double maturity, long seed);
	void changeMaturity(double newMaturity) { m_maturity = newMaturity; }
	void setIpath(long ip) { m_poisson.setIpath(ip); };
	double probaNoJump() { return exp(-m_sumFreq*m_maturity); };
	int  GetPoissonNbJumps(bool cond,int *nbJump); // if cond false, nb Jumps conditionned on at least one jump...
	int  GetPoissonNbJumps(bool cond,int *nbJump, double unif); // if cond false, nb Jumps conditionned on at least one jump...
														// unif is used instead of a random number

	long getSeed() { return m_poisson.getSeed(); };
	long getIpath() { return m_poisson.getIpath(); };
	double getMaturity() { return m_maturity; };
	void getFrequency(std::vector<double> &freq) { freq = m_freq; };

private:
	double m_maturity, m_sumFreq;
	std::vector<double> m_freq;
	Uniform m_poisson;
};

CORE_END_NAMESPACE

#endif