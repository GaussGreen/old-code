/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file ran1.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#include "gpnumlib/ran1.h"
#include "gpbase/ostringstream.h"
#include "expt.h"
#include <cmath>	/// for floor

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
////////////////////////////////////////////////////
//////////	ARM_RandUniform_NRRan1  ////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_NRRan1
///	Routine: constructor
///	Returns: built object
///	Action : 
////////////////////////////////////////////////////
ARM_RandUniform_NRRan1::ARM_RandUniform_NRRan1(long seed)
:	ARM_UniformGenerator(), itsCurrentSeed(seed), itsOriginalSeed(seed)
{
	if(seed >= 0 )
	    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,  
			"Initial Seed has to be negative!" );
};


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_NRRan1  
///	Routine: reset
///	Returns: void
///	Action : reset the random number generator
////////////////////////////////////////////////////
void ARM_RandUniform_NRRan1::reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb )
{
	itsCurrentSeed = itsOriginalSeed;
}


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_NRRan1
///	Routine: DrawOne
///	Returns: double
///	Action : do a next seed and get next number!
///		based on Numerical Recipess ran1
///		see http://www.library.cornell.edu/nr/bookcpdf.html
///		(C) Copr. 1986-92 Numerical Recipes Software #$#]2
///
///
/// Minimal random number generator of Park and Miller with Bays-Durham shuffle and added
/// safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
/// values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
/// successive deviates in a sequence. RNMX should approximate the largest floating value that is
/// less than 1.
////////////////////////////////////////////////////

double ARM_RandUniform_NRRan1::DrawOne()
{
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	/// Initialize.
	if (itsCurrentSeed <= 0 || !iy)
	{
		/// Be sure to prevent idum = 0.
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
				iv[j] = itsCurrentSeed;
		}
		iy=iv[0];
	}
	
	/// Start here when not initializing.
	k=(itsCurrentSeed)/IQ;
	
	/// Compute idum=(IAitsCurrentSeed) % IM without overows by Schrage's method. 
	itsCurrentSeed=IA*(itsCurrentSeed-k*IQ)-IR*k;
	if (itsCurrentSeed < 0) 
		itsCurrentSeed += IM;
	j=iy/NDIV;
	
	/// Will be in the range 0..NTAB-1.
	iy=iv[j];
	/// Output previously stored value and reffill the shuffle table
	iv[j] = itsCurrentSeed;
	
	if ((temp=AM*iy) > RNMX) 
		/// Because users don't expect endpoint values.
		return RNMX;
	else
		return temp;

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
};


/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_NRRan1
///	Routine: copy constructor, assignment operator
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_RandUniform_NRRan1::ARM_RandUniform_NRRan1( const ARM_RandUniform_NRRan1& rhs )
:	ARM_UniformGenerator( rhs ), 
	itsCurrentSeed( rhs.itsCurrentSeed ), 
	itsOriginalSeed( rhs.itsOriginalSeed )
{}


ARM_RandUniform_NRRan1& ARM_RandUniform_NRRan1::operator =(const ARM_RandUniform_NRRan1& rhs )
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
///	Class  : ARM_RandUniform_NRRan1
///	Routine: desctructor (necessary even for pure virtual case)
///	Action : draw a vector of uniform random numbers
////////////////////////////////////////////////////
ARM_RandUniform_NRRan1::~ARM_RandUniform_NRRan1()
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_NRRan1
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_RandUniform_NRRan1::Clone() const
{
	return new ARM_RandUniform_NRRan1(*this);
}


CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----
