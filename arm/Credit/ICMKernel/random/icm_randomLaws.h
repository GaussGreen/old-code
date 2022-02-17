#ifndef _ICM_RANDOM_LAWS_H_
#define _ICM_RANDOM_LAWS_H_

#include "ARMKernel\glob\linalg.h"
#include "ICMKernel\random\icm_RandomGenerator.h"
#include "ICMKernel\util\icm_macro.h"

/* 
Abstract class for random variables constructed from uniform variables.
like Normal or Beta...
/!\ uniform vector is const. it's NEVER CLONED.
*/



class ICM_RandomLaws : public ICM_RandomGenerator
{	
protected :
	ICM_RandomGenerator* its_pRandomUnif; // cloned
	void Init();

public : 

	ICM_RandomLaws(); 
	ICM_RandomLaws(const ICM_RandomGenerator& pRandomUnif);
	ICM_RandomLaws(const ICM_RandomLaws& ref); 
	~ICM_RandomLaws();

	virtual void View(char* id, FILE* ficOut);
	
	void SetUniformGenerator(const ICM_RandomGenerator& pRandomUnif);
	const ICM_RandomGenerator* getUniformGenerator() const ;
	ARM_Object * Clone(void);
	/// function to generate a vector
	void GenerateRandoms(ARM_Vector& RandomVectorRes);
	virtual void setParameters(const std::string& sParamName, double dParamValue) {};

	virtual void reset() {
	 if (its_pRandomUnif) 
		 its_pRandomUnif->reset();
	}
};

inline void ICM_RandomLaws::GenerateRandoms(ARM_Vector& RandomVectorRes)
{
	if(its_pRandomUnif)
		its_pRandomUnif->GenerateRandoms(RandomVectorRes);
	else 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_RandomLaws::GenerateRandoms its_pRandomUnif NULL");
}
#endif 