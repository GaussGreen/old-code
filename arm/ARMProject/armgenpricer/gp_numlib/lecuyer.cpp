/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file lecuyer.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#include "gpnumlib/lecuyer.h"
#include <cmath>		/// for floor


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
////////////////////////////////////////////////////
//////////	ARM_RandUniform_Lecuyer   ////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Lecuyer
///	Routine: DrawOne
///	Returns: double
///	Action : do a next seed and get next number!
///				generator of L'Ecuyer with Bayes & Durham shuffling
///				Combination of two short periods LCG to obtain a longer 
///				period generator. The period is the least common multiple 
///				of the 2 others.
////////////////////////////////////////////////////

double ARM_RandUniform_Lecuyer::DrawOne()
{
	static long A1		= 40014;			/// multiplier of the 1st generator
	static long A2		= 40692;			/// multiplier of the 2nd generator
	static long M1		= 2147483647;		/// 2**31 - 1  
	static double InvM1	= 1./2147483647.;	/// inverse of M1
	static long M2		= 2147483399;		/// 2**31 - 249
	static long Q1		= 53668;			/// m1 div a1  
	static long Q2		= 52774;			/// m2 div a2  
	static long R1		= 12221;			/// m1 mod a1  
	static long R2		= 3791;				/// m2 mod a2  
	long N1;
	
	static long x;
	static long y= 978543162;
	int j;
	long hi;               /// high order bit
	static long z= 0;
	static long t[32];     //// 32 is the size of a computer word
	
	N1= (M1/32);
	
	/// First call to the sequence 
	if (itsFirstUse) 
    {
		x= 437651926;
		y= x;
		
		///  After 8 "warm-ups", initialisation of the shuffle table
		for (j= 39; j>= 0; j--)
		{
			/// Park & Miller's generator
			hi= x/Q1;  
			x= A1*(x-hi*Q1) - R1*hi;
			if (x < 0) 
				x+= M1;
			if (j < 32)
				t[j]= x;
		}
		z= t[0];
		itsFirstUse = false;
    }
	
	/// For each call to the sequence, computation of a new point
	/// First generator */
	hi= x/Q1;
	x= A1*(x-hi*Q1) - R1*hi;
	if (x < 0) 
		x+= M1;
	
	/// Second generator */
	hi= y/Q2;
	y= A2*(y-hi*Q2) - R2*hi;
	if (y < 0) 
		y+= M2;
	
	/// Shuffling procedure of Bayes & Durham
	/// Index j dependent on the last point
	j= z/N1;
	
	/// Next point dependent on j
	z= t[j]- y; 
	t[j]= x;
	
	
	/// To avoid 0 value 
	if(z<1) 
		z+= M1-1;

	return z*InvM1;
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Lecuyer
///	Routine: copy constructor, assignment operator
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_RandUniform_Lecuyer::ARM_RandUniform_Lecuyer( const ARM_RandUniform_Lecuyer& rhs )
:	ARM_UniformGenerator( rhs ), 
	itsFirstUse( rhs.itsFirstUse)
{}


ARM_RandUniform_Lecuyer& ARM_RandUniform_Lecuyer::operator =(const ARM_RandUniform_Lecuyer& rhs )
{
	if( this != & rhs )
	{
		ARM_UniformGenerator::operator =( rhs );
		itsFirstUse		= rhs.itsFirstUse;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Lecuyer
///	Routine: desctructor (necessary even for pure virtual case)
///	Action : draw a vector of uniform random numbers
////////////////////////////////////////////////////
ARM_RandUniform_Lecuyer::~ARM_RandUniform_Lecuyer()
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Lecuyer
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_RandUniform_Lecuyer::Clone() const
{
	return new ARM_RandUniform_Lecuyer(*this);
}




CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----