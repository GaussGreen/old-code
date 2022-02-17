/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file Ran3.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#include "gpnumlib/ran3.h"
#include "gpbase/ostringstream.h"
#include "expt.h"
#include <cmath>	/// for floor

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
////////////////////////////////////////////////////
//////////	ARM_RandUniform_NRRan3  ////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_NRRan3
///	Routine: constructor
///	Returns: built object
///	Action : 
////////////////////////////////////////////////////
ARM_RandUniform_NRRan3::ARM_RandUniform_NRRan3(long seed)
:	ARM_UniformGenerator(), itsCurrentSeed(seed), itsOriginalSeed(seed)
{
	if(seed >= 0 )
	    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,  
			"Initial Seed has to be negative!" );
};


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_NRRan3  
///	Routine: reset
///	Returns: void
///	Action : reset the random number generator
////////////////////////////////////////////////////
void ARM_RandUniform_NRRan3::reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb ) 
{
	itsCurrentSeed = itsOriginalSeed;
}


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_NRRan3
///	Routine: DrawOne
///	Returns: double
///	Action : do a next seed and get next number!
///		based on Numerical Recipess Ran3
///		see http://www.library.cornell.edu/nr/bookcpdf.html
///		(C) Copr. 1986-92 Numerical Recipes Software #$#]2
///
///	Returns a uniform random deviate between 0.0 and 1.0. Set itsSeed
///	to any negative value to initialize or reinitialize the sequence.
////////////////////////////////////////////////////

double ARM_RandUniform_NRRan3::DrawOne()
{
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
	
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;
	
	if (itsCurrentSeed < 0 || iff == 0)
	{
		iff=1;
		mj=MSEED-(itsCurrentSeed < 0 ? -itsCurrentSeed : itsCurrentSeed);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++)
		{
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++)
			{
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
			inext=0;
			inextp=31;
			itsCurrentSeed=1;
	}
	if (++inext == 56) 
		inext=1;
	if (++inextp == 56) 
		inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) 
		mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
};


/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_NRRan3
///	Routine: copy constructor, assignment operator
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_RandUniform_NRRan3::ARM_RandUniform_NRRan3( const ARM_RandUniform_NRRan3& rhs )
:	ARM_UniformGenerator( rhs ), 
	itsCurrentSeed( rhs.itsCurrentSeed ), 
	itsOriginalSeed( rhs.itsOriginalSeed )
{}


ARM_RandUniform_NRRan3& ARM_RandUniform_NRRan3::operator =(const ARM_RandUniform_NRRan3& rhs )
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
///	Class  : ARM_RandUniform_NRRan3
///	Routine: desctructor (necessary even for pure virtual case)
///	Action : draw a vector of uniform random numbers
////////////////////////////////////////////////////
ARM_RandUniform_NRRan3::~ARM_RandUniform_NRRan3()
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_NRRan3
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_RandUniform_NRRan3::Clone() const
{
	return new ARM_RandUniform_NRRan3(*this);
}


CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----
