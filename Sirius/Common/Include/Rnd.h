#ifndef _RND_H_
#define _RND_H_

#include "smart.h"
#include "INTIntegrator.h"

class Cholesky;
 class randomGenerator : public RCObject
 {

 protected:

	 int m_icounter;
	 long m_seed;
	 CMatrix m_filledRandoms;//[iscenario][ivec]

 public:

	 int m_numberDimensions;
	 int m_nPaths;


	 virtual void generateRandoms(CVector& randoms,int idim=0);
	 virtual void generateUniformRandoms(CVector& randoms,int idim=0);

	 virtual void generateRandoms(GVector < CVector >& randoms,GVector < int > numberOfFactors, Cholesky & cholesky,int ipath,CVector& rndTemp, bool genGauss = true);
	 virtual void generateUniformRandoms(GVector < CVector >& randoms,GVector < int > numberOfFactors, Cholesky & cholesky,int ipath,CVector& rndTemp );


	 virtual double generateRandom();
	 virtual double generateUniformRandom();
	 void throwNext(CVector& randoms,int fillflag=0,int useRandoms=0);
	 double throwNext(int fillflag=0,int useRandoms=0);

	 void initialize(long idum=-123,int numberScenariosToBeStored=0,int numberRandomsPerScenario=0);
	 randomGenerator(long idum=-123,int numberScenariosToBeStored=0,int numberRandomsPerScenario=0);

	 virtual void average(CVector& results,int numberPaths);
	 virtual void completeValue(double& payoff,int ipath){};
	 virtual int offset(){return 0;};
	 virtual ~randomGenerator();
 };


class Diophantine: public randomGenerator
{


protected:

	RCPtr < CAlpha > m_pAlpha ;	
	RCPtr < CWeights > m_pWeights ;
	RCPtr < CDomain > m_pDomain ;
	RCPtr < CPeriodization > m_pPeriodization ;
	RCPtr < CNTintegrator > m_pIntegrator ;
	RCPtr < CNTparameters > m_pParameters ;


	int m_AmIsConst;
	int m_JacobiIsConst;
	double m_jacobi;
	vector<double>  m_line ;//[idimesion]
	pair<long, long> m_Nrange;

	long m_Neff;


public:


	void generateRandoms(CVector& randoms,int ipath=0);
	void generateUniformRandoms(CVector& randoms,int ipath=0);

	void average(CVector& results,int numberPaths);
	void completeValue(double& payoff,int ipath);
	int  offset();

	void initialize(	
						int ndimension,
						int npaths,
						CAlpha * pAlpha ,	
						CWeights * pWeights,
						CDomain * pDomain,
						CPeriodization * pPeriodization,
						CNTintegrator* pIntegrator,
						CNTparameters* pParameters);


	Diophantine(int ndimension,
				int npaths,
				CAlpha * pAlpha ,	
				CWeights * pWeights,
				CDomain * pDomain,
				CPeriodization * pPeriodization,
				CNTintegrator* pIntegrator,
				CNTparameters* pParameters);

	Diophantine(){};

	double generateRandom(){assert(false);return -1e99;};

};


#endif