#ifndef _ICM_RANDOM_RAN2_H_
#define _ICM_RANDOM_RAN2_H_

#include "ARMKernel\glob\linalg.h"
#include "ICMKernel\random\icm_RandomGenerator.h"

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

class ICM_RandomRan2 : public ICM_RandomGenerator
{

private :
	long		itsInitialSeed; 
	long		itsCurrentSeed; 
	//double		itsRandom;

	// unsetted
	long its_iy;
	long its_iv[NTAB];
	long its_idum2;

protected :

	void Init(long InitialSeed);

public : 

	ICM_RandomRan2(); 
	ICM_RandomRan2(long InitialSeed);
	ICM_RandomRan2(const ICM_RandomRan2& ref); 

	~ICM_RandomRan2();

	ICM_RandomRan2& operator=(const ICM_RandomRan2& ref);
	
	// Set - Get Functions 
	long getInitialSeed() const{return itsInitialSeed;}
	long getCurrentSeed() const{return itsCurrentSeed;}
	//double getRandom() const{return itsRandom;}

	ARM_Object * Clone(void);

	virtual void View(char* id, FILE* ficOut);
	
	// Pure virtual functions to be defined in derived classes : 

	// function to generate a numbe
	void GenerateRandomRan2(ARM_Vector& RandomVector);
	/// function to generate a vector
	virtual void GenerateRandoms(ARM_Vector& RandomVector) ;

	
	// reset enables to reset the random number generator!
	// gives the dimension and indicates the nb of points to generate!
	virtual void reset();

	virtual void setParameters(const std::string& sParamName, double dParamValue);
};

inline void ICM_RandomRan2::GenerateRandoms(ARM_Vector& RandomVector) 
{
	GenerateRandomRan2(RandomVector);
}

inline void ICM_RandomRan2::GenerateRandomRan2(ARM_Vector& RandomVector) 
{
	
	
	double * dVect = RandomVector.GetElt();
	for(int i =0; i< RandomVector.size(); i++) {
		int j;
		long k;
		double temp;
		
		if (itsCurrentSeed <= 0) 
		{
			if (-(itsCurrentSeed) < 1) itsCurrentSeed=1;
			else itsCurrentSeed = -(itsCurrentSeed);
			its_idum2=(itsCurrentSeed);
			for (j=NTAB+7;j>=0;j--) 
			{
				k=(itsCurrentSeed)/IQ1;
				itsCurrentSeed=IA1*(itsCurrentSeed-k*IQ1)-k*IR1;
				if (itsCurrentSeed < 0) itsCurrentSeed += IM1;
				if (j < NTAB) its_iv[j] = itsCurrentSeed;
			}
			its_iy=its_iv[0];
		}

		k=(itsCurrentSeed)/IQ1;
		itsCurrentSeed=IA1*(itsCurrentSeed-k*IQ1)-k*IR1;
		if (itsCurrentSeed < 0) itsCurrentSeed += IM1;
		k=its_idum2/IQ2;
		its_idum2=IA2*(its_idum2-k*IQ2)-k*IR2;
		if (its_idum2 < 0) its_idum2 += IM2;
		j=its_iy/NDIV;
		its_iy=its_iv[j]-its_idum2;
		its_iv[j] = itsCurrentSeed;
		if (its_iy < 1) its_iy += IMM1;
		if ((temp=AM*its_iy) > RNMX) 
		{	
			dVect[i] = RNMX;}
		else 
		{	
			dVect[i] = temp;}
	}
	

}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

#endif 