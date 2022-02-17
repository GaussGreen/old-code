#ifndef _ICM_RANDOM_GENERATOR_H_
#define _ICM_RANDOM_GENERATOR_H_

#include "ARMKernel\glob\linalg.h"


class ICM_RandomGenerator : public ARM_Object
{

protected :

	virtual void Init();

public : 

	ICM_RandomGenerator();

	ICM_RandomGenerator(const ICM_RandomGenerator& ref);
	
	virtual ~ICM_RandomGenerator();

	// ----------------------------
	//	Copy of members data
	// ----------------------------
/*	void BitwiseCopy(const ARM_Object * src)
	{
	    ICM_RandomGenerator * Generator = (ICM_RandomGenerator *) src ;

		itsDim      = Generator->getDim() ;
	}

	// -------------
	//	Copy Method 
	// -------------
	
	void Copy(const ARM_Object* src)
	{
		ARM_Object::Copy(src) ;
 
		BitwiseCopy(src) ;
	}

*/
	virtual void View(char* id, FILE* ficOut);

	// Pure virtual functions to be defined in derived classes : 

	/// function to generate a vector
	double GenerateOneRandom() ;

	virtual	void GenerateRandoms(ARM_Vector& RandomVector) = 0;
	
	// reset enables to reset the random number generator!
	// gives the dimension and indicates the nb of points to generate!
	virtual void reset()  = 0;

	virtual void setParameters(const std::string& sParamName, double dParamValue) =0;
};

inline double ICM_RandomGenerator::GenerateOneRandom(){
	ARM_Vector v(1);
	GenerateRandoms(v);
	return v.Elt(0);
}


#endif 