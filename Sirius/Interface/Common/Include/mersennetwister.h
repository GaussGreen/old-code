#ifndef MERSENNETWISTER_H
#define MERSENNETWISTER_H


#pragma once

#include "smart.h"
#include "rnd.h"

class CMersenneTwister	:	public randomGenerator 
{
public:
	CMersenneTwister(){}
	CMersenneTwister( int dim );
	~CMersenneTwister(){};

	virtual void generateRandoms( CVector& nRandVect, int idim = 0 );
	virtual void generateRandoms( double* nptr );
	virtual void generateUniformRandoms( CVector& uRandVect, int idim = 0 );

	virtual double generateRandom();
	virtual double generateUniformRandom();

	void initialize(int dimension);
	void initialize(int dimension, int numberScenariosToBeStored,int numberRandomsPerScenario);

//	void generate_numbers();

//	static double *setup_random_numbers();	

//	double* next();

protected:

	
	/* initializes mt[N] with a seed */
	void init_genrand(unsigned long s);

	/* initialize by an array with array-length */
	/* init_key is the array for initializing keys */
	/* key_length is its length */
	/* slight change for C++, 2004/2/26 */
	void init_by_array(unsigned long init_key[], int key_length);

	/* generates a random number on [0,0xffffffff]-interval */
	unsigned long genrand_int32(void);

	/* generates a random number on [0,0x7fffffff]-interval */
	long genrand_int31(void);

	/* generates a random number on [0,1]-real-interval */
	double genrand_real1(void);

	/* generates a random number on [0,1)-real-interval */
	double genrand_real2(void);

	/* generates a random number on (0,1)-real-interval */
	double genrand_real3(void);

	/* generates a random number on [0,1) with 53-bit resolution*/
	double genrand_res53(void) ;

protected:

	void	generateRandoms(GVector < CVector >& randoms,GVector < int > numberOfFactors, Cholesky & cholesky,int ipath,CVector& rndTemp,bool generateGauss);
	void	initCorrelationParameters(GVector < CVector >& randoms,GVector < int > numberOfFactors, Cholesky & cholesky,int ipath,CVector& rndTemp,bool generateGauss);

	int		m_nAssets ;
	int		m_nFactors ;
	int		m_nSteps ;

	CVector m_correl;
};


//extern bool mt_initialized;
//extern unsigned long mt_nbCalls;
//extern double *pRandomNumber;
//extern double random_numbers[];

#endif
