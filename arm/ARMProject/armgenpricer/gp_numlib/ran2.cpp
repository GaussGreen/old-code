/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file ran2.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */


#include "gpnumlib/ran2.h"
#include "gpbase/ostringstream.h"
#include "expt.h"
#include <cmath>	/// for floor

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
////////////////////////////////////////////////////
//////////	ARM_RandUniform_NRRan2  ////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_NRRan2
///	Routine: constructor
///	Returns: built object
///	Action : 
////////////////////////////////////////////////////
ARM_RandUniform_NRRan2::ARM_RandUniform_NRRan2(long seed)
:	ARM_UniformGenerator(), itsCurrentSeed(seed), itsOriginalSeed(seed)
{
	if(seed >= 0 )
	    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,  
			"Initial Seed has to be negative!" );
};


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_NRRan2  
///	Routine: reset
///	Returns: void
///	Action : reset the random number generator
////////////////////////////////////////////////////
void ARM_RandUniform_NRRan2::reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb )
{
	itsCurrentSeed = itsOriginalSeed;
}


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_NRRan2
///	Routine: DrawOne
///	Returns: double
///	Action : do a next seed and get next number!
///		based on Numerical Recipess ran2
///		see http://www.library.cornell.edu/nr/bookcpdf.html
///		(C) Copr. 1986-92 Numerical Recipes Software #$#]2
///		Long period (> 2 × 1018) random number generator of L’Ecuyer with Bays-Durham shu.e
///		and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
///		the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
///		itsSeed between successive deviates in a sequence. RNMX should approximate the largest .oating
///	value that is less than 1.
////////////////////////////////////////////////////
double ARM_RandUniform_NRRan2::DrawOne()
{
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
	
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;
	
	if (itsCurrentSeed <= 0) 
	{
		if (-(itsCurrentSeed) < 1) itsCurrentSeed=1;
		else itsCurrentSeed = -(itsCurrentSeed);
		idum2=(itsCurrentSeed);
		for (j=NTAB+7;j>=0;j--) 
		{
			k=(itsCurrentSeed)/IQ1;
			itsCurrentSeed=IA1*(itsCurrentSeed-k*IQ1)-k*IR1;
			if (itsCurrentSeed < 0) itsCurrentSeed += IM1;
			if (j < NTAB) iv[j] = itsCurrentSeed;
		}
		iy=iv[0];
	}

	k=(itsCurrentSeed)/IQ1;
	itsCurrentSeed=IA1*(itsCurrentSeed-k*IQ1)-k*IR1;
	if (itsCurrentSeed < 0) itsCurrentSeed += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = itsCurrentSeed;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
	
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
};


/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_NRRan2
///	Routine: copy constructor, assignment operator
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_RandUniform_NRRan2::ARM_RandUniform_NRRan2( const ARM_RandUniform_NRRan2& rhs )
:	ARM_UniformGenerator( rhs ), 
	itsCurrentSeed( rhs.itsCurrentSeed ), 
	itsOriginalSeed( rhs.itsOriginalSeed )
{}


ARM_RandUniform_NRRan2& ARM_RandUniform_NRRan2::operator =(const ARM_RandUniform_NRRan2& rhs )
{
	if( this != & rhs )
	{
		ARM_UniformGenerator::operator =( rhs );
		itsCurrentSeed	= rhs.itsCurrentSeed;
		itsOriginalSeed	= rhs.itsOriginalSeed;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_NRRan2
///	Routine: desctructor (necessary even for pure virtual case)
///	Action : draw a vector of uniform random numbers
////////////////////////////////////////////////////
ARM_RandUniform_NRRan2::~ARM_RandUniform_NRRan2()
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_NRRan2
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_RandUniform_NRRan2::Clone() const
{
	return new ARM_RandUniform_NRRan2(*this);
}


CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----