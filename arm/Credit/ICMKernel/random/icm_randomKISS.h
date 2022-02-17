#ifndef _ICM_RANDOM_KISS_H_
#define _ICM_RANDOM_KISS_H_

#include "ARMKernel\glob\linalg.h"
#include "ICMKernel\random\icm_RandomGenerator.h"

/*
Subject: Re: KISS random number generator
Author: George Marsaglia <geo@stat.fsu.edu> 
Date Posted: Jan 17 2003 12:00:35:000PM
(...)
Four random seeds are required,
0<=x<2^32, 0<y<2^32, 0<=z<2^32, 0<=c<698769069.
If the static x,y,z,c is placed outside the KISS proc, a seed-set
routine may be added to allow the calling program to set the seeds.
In gaming machines these would presumably be set by clocks or
registers based on cumulative previous results.
The period is 2^32*(2^32-1)*(698769069*2^31-1)> 2^124 or 10^37.

For variations: In the congruential part, 69069 may be replaced
by any multiplier that is 3 or 5 mod 8; for 12345, any odd constant.
For the 3-shifts, there are 161 triples besides (13,17,5) that
provide maximal period (if interested, send for article), and
for the multiply-with-carry sequence, choose any multiplier 'a'
such as 698769069 for which both a*2^32-1 and a*2^31-1 are primes.
*/

/*
LONG_MAX 
 Maximum value for a variable of type long.
 2147483647

ULONG_MAX 
 Maximum value for a variable of type unsigned long.
 4294967295 (0xffffffff)
 
  */

#define ARM_ULONG_MAX 4294967295;

class ICM_RandomKISS : public ICM_RandomGenerator
{

private :
	
	//double		itsRandom;
	unsigned long			its_x;
	unsigned long			its_y;
	unsigned long			its_z;
	unsigned long			its_c;

	// unsettable "unsigned long long"
	unsigned __int64		its_a; 



protected :

	void Init(unsigned long x, unsigned long y, unsigned long z, unsigned long c);

public : 

	ICM_RandomKISS(); 
	ICM_RandomKISS(unsigned long x, unsigned long y, unsigned long z, unsigned long c);
	ICM_RandomKISS(const ICM_RandomKISS& ref); 

	~ICM_RandomKISS();

	ICM_RandomKISS& operator=(const ICM_RandomKISS& ref);
	
	// Set - Get Functions 
	int getX() const{return its_x;}
	int getY() const{return its_y;}
	int getZ() const{return its_z;}
	int getC() const{return its_c;}
	//double getRandom() const{return itsRandom;}



	// ----------------------------
	//	Copy of members data
	// ----------------------------
/*	void BitwiseCopy(const ARM_Object * src)
	{
	    ICM_RandomKISS * Generator = (ICM_RandomKISS *) src ;

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
	void GenerateRandomKISS(ARM_Vector& RandomVector) ;
	/// function to generate a vector
	virtual void GenerateRandoms(ARM_Vector& RandomVector) ;
	// reset enables to reset the random number generator!
	// gives the dimension and indicates the nb of points to generate!
	virtual void reset();

	virtual void setParameters(const std::string& sParamName, double dParamValue);


};

inline void ICM_RandomKISS::GenerateRandoms(ARM_Vector& RandomVector) 
{	
	GenerateRandomKISS(RandomVector);
}

inline void ICM_RandomKISS::GenerateRandomKISS(ARM_Vector& RandomVector) 
{
	unsigned __int64 t=0;
	double * dVect = RandomVector.GetElt();
	int Dim = RandomVector.size();
	for (int i = 0 ; i<Dim; i++ ) {		
		its_x=69069*its_x+12345;
		its_y= pow(static_cast<double>(its_y),static_cast<double>((its_y<<13))); 
		its_y= pow(static_cast<double>(its_y),static_cast<double>((its_y>>17)));
		its_y= pow(static_cast<double>(its_y),static_cast<double>((its_y<<5)));
		t=its_a*its_z+its_c;
		its_c=(t>>32);
		dVect[i]= (double)(its_x+its_y+(its_z=t))/ (double)ARM_ULONG_MAX;
	}
}
#undef ARM_ULONG_MAX 

#endif 