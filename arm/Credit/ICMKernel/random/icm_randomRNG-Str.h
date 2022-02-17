#ifndef _ICM_RANDOM_RNG_STR_H_
#define _ICM_RANDOM_RNG_STR_H_

#include "ARMKernel\glob\linalg.h"
#include "ICMKernel\random\icm_RandomRan1.h"

/*
From RNG STR used in CorrelIntertranches 
initia seed and CurrentSeed
*/


class ICM_RandomRNG_Str : public ICM_RandomRan1
{
protected :
	void Init(long InitialSeed);

public : 

	ICM_RandomRNG_Str(); 
	ICM_RandomRNG_Str(long InitialSeed);
	ICM_RandomRNG_Str(const ICM_RandomRNG_Str& ref); 

	~ICM_RandomRNG_Str();

	ICM_RandomRNG_Str& operator=(const ICM_RandomRNG_Str& ref);
	
	ARM_Object * Clone(void);
	// Pure virtual functions to be defined in derived classes : 

	// function to generate a number
	void GenerateRandomRNG_Str(ARM_Vector& RandomVector) ;
	/// function to generate a vector
	virtual void GenerateRandoms(ARM_Vector& RandomVector) ;
	// reset enables to reset the random number generator!
	// gives the dimension and indicates the nb of points to generate!
	virtual void reset();
};

inline void ICM_RandomRNG_Str::GenerateRandoms(ARM_Vector& RandomVector) 
{	
	GenerateRandomRNG_Str(RandomVector);
}

inline void ICM_RandomRNG_Str::GenerateRandomRNG_Str(ARM_Vector& RandomVector) 
{
	//double itsCurrentSeed=itsInitialSeed; 
	double * dvect = RandomVector.GetElt();
	for( int i=0;i<RandomVector.size();i++) 
	{
		// double resultat = 0.;
		// double x,u,a, test = 0.;
		// a = seed;
		double test = itsCurrentSeed * 16807./2147483647.;
		if (test > 0)   // Fix = Floor
		{	itsCurrentSeed = (itsCurrentSeed  * 16807.) - floor(test) * 2147483647;}
		else			// Fix = Floor + 1 
		{	itsCurrentSeed = (itsCurrentSeed  * 16807.) - (floor(test) + 1) * 2147483647;}
		dvect[i] = itsCurrentSeed/2147483647.;
	}
}


#endif 