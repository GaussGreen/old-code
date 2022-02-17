#ifndef _ICM_RANDOM_RANMAR_H_
#define _ICM_RANDOM_RANMAR_H_

#include "ARMKernel\glob\linalg.h"
#include "ICMKernel\random\icm_RandomGenerator.h"

class ICM_RandomRanmar : public ICM_RandomGenerator
{

private :
	//double		itsRandom;
	int			its_ij, its_kl;

// non settable
	float its_u[98];
	float its_c1;
	float its_cd; 
	float its_cm;
	int its_i97;
	int its_j97;
	bool its_testvar;
protected :

	void Init(int _ij, int _kl);

public : 

	ICM_RandomRanmar(); 
	ICM_RandomRanmar(int _ij, int _kl);
	ICM_RandomRanmar(const ICM_RandomRanmar& ref); 

	~ICM_RandomRanmar();

	ICM_RandomRanmar& operator=(const ICM_RandomRanmar& ref);
	
	// Set - Get Functions 
	int getij() const{return its_ij;}
	int getkl() const{return its_kl;}
	//double getRandom() const{return itsRandom;}



	// ----------------------------
	//	Copy of members data
	// ----------------------------
/*	void BitwiseCopy(const ARM_Object * src)
	{
	    ICM_RandomRanmar * Generator = (ICM_RandomRanmar *) src ;

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

	void rmarin(int ij, int kl);

	ARM_Object * Clone(void);

	virtual void View(char* id, FILE* ficOut);
	
	// Pure virtual functions to be defined in derived classes : 

	// function to generate a number
	void GenerateRannmar(ARM_Vector& RandomVector)  ;
	/// function to generate a vector
	virtual void GenerateRandoms(ARM_Vector& RandomVector) ;
	// reset enables to reset the random number generator!
	// gives the dimension and indicates the nb of points to generate!
	virtual void reset();
	virtual void setParameters(const std::string& sParamName, double dParamValue);
};

inline void ICM_RandomRanmar::GenerateRandoms(ARM_Vector& RandomVector) 
{
	GenerateRannmar(RandomVector);
}

inline void ICM_RandomRanmar::GenerateRannmar(ARM_Vector& RandomVector) 
{
	float uni;
	double* dVect = RandomVector.GetElt();
	for (int i = 0 ; i<RandomVector.size(); i++ ) {		
		uni =its_u[its_i97] -its_u[its_j97];
		if (uni < 0.0) uni += 1.0;
		its_u[its_i97] = uni;
		its_i97--;
		if (its_i97==0) its_i97 = 97;
		its_j97--;
		if (its_j97==0) its_j97 = 97;
		its_c1 -=its_cd;
		if (its_c1<0.0)its_c1 +=its_cm;
		uni -=its_c1;
		if (uni<0.0) uni += 1.0;
		dVect[i] = uni;
	}
}

#endif 