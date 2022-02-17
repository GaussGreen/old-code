#ifndef _ICM_RANDOM_RAN1_H_
#define _ICM_RANDOM_RAN1_H_

#include "ARMKernel\glob\linalg.h"
#include "ICMKernel\random\icm_RandomGenerator.h"
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPSr 1.2e-7
#define RNMX (1.0-EPSr)

class ICM_RandomRan1 : public ICM_RandomGenerator
{
private :
	long its_iy;
	long its_iv[NTAB];

protected :
	long		itsInitialSeed; 
	long		itsCurrentSeed; 
	//double		itsRandom;

protected :

	void Init(long InitialSeed);

public : 

	ICM_RandomRan1(); 
	ICM_RandomRan1(long InitialSeed);
	ICM_RandomRan1(const ICM_RandomRan1& ref); 

	~ICM_RandomRan1();

	ICM_RandomRan1& operator=(const ICM_RandomRan1& ref);
	
	// Set - Get Functions 
	long getInitialSeed() const{return itsInitialSeed;}
	long getCurrentSeed() const{return itsCurrentSeed;}
	//double getRandom() const{return itsRandom;}



	// ----------------------------
	//	Copy of members data
	// ----------------------------
/*	void BitwiseCopy(const ARM_Object * src)
	{
	    ICM_RandomRan1 * Generator = (ICM_RandomRan1 *) src ;

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

	// function to generate a number
	void GenerateRandomRan1(ARM_Vector& RandomVector) ;
	/// function to generate a vector
	virtual void GenerateRandoms(ARM_Vector& RandomVector) ;

	
	// reset enables to reset the random number generator!
	// gives the dimension and indicates the nb of points to generate!
	virtual void reset();

	virtual void setParameters(const std::string& sParamName, double dParamValue);


};

inline void ICM_RandomRan1::GenerateRandoms(ARM_Vector& RandomVector) 
{	
	GenerateRandomRan1(RandomVector);
}

inline void ICM_RandomRan1::GenerateRandomRan1(ARM_Vector& RandomVector) 
{

	double * dVect = RandomVector.GetElt();
	for(int i =0; i< RandomVector.size(); i++) {
			int j;
			long k;	
			double temp;
		/// Initialize.
		if (itsCurrentSeed <= 0 || !its_iy)
		{
			/// Be sure to prevent itsCurrentSeed = 0.
			if (-(itsCurrentSeed) < 1) 
				itsCurrentSeed=1;
			else 
				itsCurrentSeed = -(itsCurrentSeed);
			for (j=NTAB+7;j>=0;j--)
			{
				/// Load the shuffle table (after 8 warm-ups).
				k=(itsCurrentSeed)/IQ;
				itsCurrentSeed=IA*(itsCurrentSeed-k*IQ)-IR*k;
				if (itsCurrentSeed < 0) 
					itsCurrentSeed += IM;
				if (j < NTAB) 
					its_iv[j] = itsCurrentSeed;
			}
			its_iy=its_iv[0];
		}
		
		/// Start here when not initializing.
		k=(itsCurrentSeed)/IQ;
		
		/// Compute itsCurrentSeed=(IAitsCurrentSeed) % IM without overows by Schrage's method. 
		itsCurrentSeed=IA*(itsCurrentSeed-k*IQ)-IR*k;
		if (itsCurrentSeed < 0) 
			itsCurrentSeed += IM;
		j=its_iy/NDIV;
		
		/// Will be in the range 0..NTAB-1.
		its_iy=its_iv[j];
		/// Output previously stored value and reffill the shuffle table
		its_iv[j] = itsCurrentSeed;
		
		if ((temp=AM*its_iy) > RNMX) 
			/// Because users don't expect endpoint values.
		{
			dVect[i] =  RNMX;
		}else {	
			dVect[i] =  temp;
		}
	}
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB 
#undef NDIV
#undef EPSr
#undef RNMX
#endif 