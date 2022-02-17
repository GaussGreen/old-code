#ifndef _ICM_RANDOM_INVNORMMORO_H_
#define _ICM_RANDOM_INVNORMMORO_H_

#include "ARMKernel\glob\linalg.h"
#include "ICMKernel\random\icm_RandomLaws.h"

/* from 
 ep::MathSrv::cumNorm*/


class ICM_RandomInvNormMoro : public ICM_RandomLaws
{

private :
	// unsetted
	double its_a1;
	double its_a2;
	double its_a3;
	double its_a4;

	double its_b1;
	double its_b2;
	double its_b3;
	double its_b4;

	double its_c1;
	double its_c2;
	double its_c3;
	double its_c4;
	double its_c5;
	double its_c6;
	double its_c7;
	double its_c8;
	double its_c9;

protected :
	void Init();

public : 

	ICM_RandomInvNormMoro(); 
	ICM_RandomInvNormMoro(const ICM_RandomGenerator& pRandomUnif);
	ICM_RandomInvNormMoro(const ICM_RandomInvNormMoro& ref); 
	~ICM_RandomInvNormMoro();
	ICM_RandomInvNormMoro& operator=(const ICM_RandomInvNormMoro& ref);
	ARM_Object * Clone(void);
	// Pure virtual functions to be defined in derived classes : 

	// function to generate  a vector
	void GenerateRandomInvNormMoro(ARM_Vector& RandomVector) const ;
	/// function to generate a vector
	virtual void GenerateRandoms(ARM_Vector& RandomVector) ;
};

inline void ICM_RandomInvNormMoro::GenerateRandoms(ARM_Vector& RandomVector) 
{	
	ICM_RandomLaws::GenerateRandoms(RandomVector);
	GenerateRandomInvNormMoro(RandomVector);
}

inline void ICM_RandomInvNormMoro::GenerateRandomInvNormMoro (ARM_Vector& RandomVector) const {

		//  Replaces NormSInv for quasi-random sequences (eg Faure)
	//  See Moro (1995)	
	double* dRandom = RandomVector.GetElt();
	double resultat = 0.;
	double y=0. ,p = 0.;
	  
	for (int i=0; i< RandomVector.size(); i++) {
		y = dRandom[i] - 0.5;
		
		if (fabs(y)<0.42)
		{	
			p = y*y;
			p = y * (((its_a4 * p + its_a3) * p + its_a2) * p + its_a1) / ((((its_b4 * p + its_b3) * p + its_b2) * p + its_b1) * p + 1);
		}
		else
		{
			if (y > 0)	{p = log(-log(1 - dRandom[i]));}
			if (y <= 0) {p = log(-log(dRandom[i]));}
			p = its_c1 + p * (its_c2 + p * (its_c3 + p * (its_c4 + p * (its_c5 + p * (its_c6 + p * (its_c7 + p * (its_c8 + p * its_c9)))))));
			if (y <= 0) {p = -p;}
		}

		dRandom[i] = p;
	}
}

#endif 