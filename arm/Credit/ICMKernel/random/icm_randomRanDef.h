#ifndef _ICM_RANDOM_RANDEF_H_
#define _ICM_RANDOM_RANDEF_H_

#include "ARMKernel\glob\linalg.h"
#include "ICMKernel\random\icm_RandomGenerator.h"

class ICM_RandomRanDef : public ICM_RandomGenerator
{

private :
	int			itsInitialSeed; 
	bool		itsIsRandomSeed;
	//double		itsRandom;

protected :

	void Init(int InitialSeed, bool IsRandomSeed);

public : 

	ICM_RandomRanDef(); 
	ICM_RandomRanDef(int InitialSeed,  bool IsRandomSeed);
	ICM_RandomRanDef(const ICM_RandomRanDef& ref); 

	~ICM_RandomRanDef();

	ICM_RandomRanDef& operator=(const ICM_RandomRanDef& ref);
	
	// Set - Get Functions 
	int getInitialSeed() const { return itsInitialSeed;}
	bool getIsRandomSeed() const { return itsIsRandomSeed;}
	//double getRandom() const { return itsRandom;}

	/*
	void setInitialSeed(const int i) ;
	void setIsRandomSeed(const bool b) ;
	*/
	// ----------------------------
	//	Copy of members data
	// ----------------------------
/*	void BitwiseCopy(const ARM_Object * src)
	{
	    ICM_RandomRanDef * Generator = (ICM_RandomRanDef *) src ;

		itsInitialSeed      = Generator->itsInitialSeed ;
		itsCurrentSeed      = Generator->itsCurrentSeed ;
	
	}


	// -------------
	//	Copy Method 
	// -------------
	
	void Copy(const ARM_Object* src)
	{
		ICM_RandomGenerator::Copy(src) ;
 
		BitwiseCopy(src) ;
	}
*/
	// --------------
	//	Clone Method
	// --------------

	ARM_Object * Clone(void);

	virtual void View(char* id, FILE* ficOut);
	
	// Pure virtual functions to be defined in derived classes : 

	/// function to generate a vector
	virtual void GenerateRandoms(ARM_Vector& RandomVector);
	
	// reset enables to reset the random number generator!
	// gives the dimension and indicates the nb of points to generate!
	virtual void reset();
	virtual void setParameters(const std::string& sParamName, double dParamValue);
};

inline void ICM_RandomRanDef::GenerateRandoms(ARM_Vector& RandomVector) 
{
	long Dim = RandomVector.size(); 
	double * dVect = RandomVector.GetElt();
	for (int i = 0 ; i<Dim; i++ ) 
		dVect[i] = (double)rand()/(RAND_MAX);
}

#endif 