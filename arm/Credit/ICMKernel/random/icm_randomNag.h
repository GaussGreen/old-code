#ifndef _ICM_RANDOM_NAG_H_
#define _ICM_RANDOM_NAG_H_

#include "ARMKernel\glob\linalg.h"
#include "ICMKernel\random\icm_RandomGenerator.h"
#include <nag.h>
#include <nagg05.h>

class ICM_RandomNag : public ICM_RandomGenerator
{

private :
	long		itsInitialSeed; 
	bool		itsIsRandomSeed;
	//double		itsRandom;

protected :

	void Init(long InitialSeed, bool IsRandomSeed);

public : 

	ICM_RandomNag();
	ICM_RandomNag(long InitialSeed, bool IsRandomSeed);
	ICM_RandomNag(const ICM_RandomNag& ref);

	~ICM_RandomNag();

	ICM_RandomNag& operator=(const ICM_RandomNag& ref);
	
	// Set - Get Functions 
	long getInitialSeed() const{return itsInitialSeed;}
	bool getIsRandomSeed() const{return itsIsRandomSeed;}
	//double getRandom() const{return itsRandom;}



	// ----------------------------
	//	Copy of members data
	// ----------------------------
/*	void BitwiseCopy(const ARM_Object * src)
	{
	    ICM_RandomNag * Generator = (ICM_RandomNag *) src ;

		itsInitialSeed      = Generator->itsInitialSeed ;
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

	// function to generate a number
	virtual void GenerateRandoms(ARM_Vector& RandomVector) ;
	/// function to generate a vector
	
	// reset enables to reset the random number generator!
	// gives the dimension and indicates the nb of points to generate!
	virtual void reset();

	virtual void setParameters(const std::string& sParamName, double dParamValue);


};

inline void ICM_RandomNag::GenerateRandoms(ARM_Vector& RandomVector) 
{	
	long Dim = RandomVector.size(); 
	double * dVect = RandomVector.GetElt();
	for (int i = 0 ; i<Dim; i++ ) 
		dVect[i] = g05cac();
}

#endif 