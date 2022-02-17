/* file: SciSobol.h */
/* Sobol sequences up to d = 350 */

#ifndef SCISOBOL_H
#define SCISOBOL_H

#include "MlEqMaths.h"
#include "Rnd.h"

// commenting out following line removes an important warning

 class Sobol : public randomGenerator
 {

	friend double DISTfunc(double x,void* vp);

	class SciSobolSeqStruct 
	{
	 public:
		int dimension;
		int counter;
		int leap;
		double fac;
		unsigned long *iv;
		unsigned long *ix;
	}m_SciSobolSeqStruct;

	void SciSobolInit(
		int dimension,
		int skip,
		int leap);

	void SciSobolNext(
		double *x,
		void* vp,
		double (*DISTfunc)(double x,void* vp)
		);

	void SciSobolFree();
	CVector m_rand;//[m_dimesion]

public:

	int m_dimension;

	Sobol(int dimension,int skip,int leap,int numberScenariosToBeStored,int numberRandomsPerScenario);
	void initialize(int dimension,int skip=1010,int leap=1,int numberScenariosToBeStored=0,int numberRandomsPerScenario=0);
	void virtual generateRandoms(CVector& random,int ipath);
//	void virtual generateRandoms(GVector < CVector > & randoms,int ipath);
	double generateRandom();

	virtual  ~Sobol(); 
};

#endif
