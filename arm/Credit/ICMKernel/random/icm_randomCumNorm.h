#ifndef _ICM_RANDOM_CUMNORMSVR_H_
#define _ICM_RANDOM_CUMNORMSVR_H_

#include "ARMKernel\glob\linalg.h"
#include "ICMKernel\random\icm_RandomLaws.h"

/* from 
 ep::MathSrv::cumNorm
 fonction de repartition d'une variable normale centrée reduite.

*/

class ICM_RandomCumNorm : public ICM_RandomLaws
{
public : 

	ICM_RandomCumNorm(); 
	ICM_RandomCumNorm(const ARM_Vector* pvUnifRand);
	ICM_RandomCumNorm(const ICM_RandomCumNorm& ref); 

	~ICM_RandomCumNorm();
	ICM_RandomCumNorm& operator=(const ICM_RandomCumNorm& ref);
	ARM_Object * Clone(void);
	
	// Pure virtual functions to be defined in derived classes : 

	// function to generate  a vector
	void GenerateRandomCumNorm(ARM_Vector& RandomVector) const ;
	/// function to generate a vector
	virtual void GenerateRandoms(ARM_Vector& RandomVector) ;
};

inline void ICM_RandomCumNorm::GenerateRandoms(ARM_Vector& RandomVector) 
{	
	ICM_RandomLaws::GenerateRandoms(RandomVector);
	GenerateRandomCumNorm(RandomVector);
}


inline void ICM_RandomCumNorm::GenerateRandomCumNorm(ARM_Vector& RandomVector) const 
{
	
	double	t=0;
	double	z=0;	
	double 	ans=0;
	double poly=0;
	const double * pUnifVect = its_pUniformRandom->GetElt();
	double * pRandVect = RandomVector.GetElt();
	for( int i=0; i< RandomVector.size(); i++) {
		z = fabs(pUnifVect[i]/sqrt(2));
		t = 1.0/(1.0+0.5*z);
		poly = -1.26551223 +
			t * ( 1.00002368 +
			t * ( 0.37409196 +
			t * ( 0.09678418 +
    		t * (-0.18628806 +
			t * ( 0.27886807 +
			t * (-1.13520398 +
			t * ( 1.48851587 +
    		t * (-0.82215223 +
 			t *   0.17087277 ) ) ) ) ) ) ) );
		ans = t * exp( -z*z + poly);

		if (pUnifVect[i] <=0.0)
			pRandVect[i] = 0.5 * ans;
		else
			pRandVect[i] = 1.0 - 0.5 * ans;
	}
}
#endif 