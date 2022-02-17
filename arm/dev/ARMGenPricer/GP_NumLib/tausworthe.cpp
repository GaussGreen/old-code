/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file tausworthe.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#include "gpnumlib/tausworthe.h"

CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////
////////////////////////////////////////////////////
//////////	ARM_RandUniform_Tausworthe   ///////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Tausworthe
///	Routine: Bit_random
///	Returns: double
///	Action : Generation of a random bit
///				Algorithm based on a prime polynomial : 
///				here we choose x^18 + x^5 + x^2 + x + 1
///				It is described in 'Numerical Recipes in C' page 296.
////////////////////////////////////////////////////

int ARM_RandUniform_Tausworthe::Bit_random()
{
	int degre = 18;
	static unsigned long a;
	unsigned long new_bit;
	
	/// Initialisation for the n first values
	/// random number over [1, 2^18] ; 2^18= 262144
	if(itsFirstUse)
    {
		a = 176355;
		itsFirstUse = false;
    }
    
	/// Next bit calculation by the recurrence relation
	new_bit= (a & (1<<17)) >> 17
		^ (a & (1<<4)) >> 4
		^ (a & (1<<1)) >> 1
		^ (a & 1);
	a <<= 1;
	/// The new bit is shift on the right
	a ^= new_bit;
	/// The most left bit is put to 0
	a ^= (1 << degre);
	
	return((int)new_bit);
}


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Tausworthe
///	Routine: DrawOne
///	Returns: double
///	Action : Generation of a word of k random bits
////////////////////////////////////////////////////

unsigned long ARM_RandUniform_Tausworthe::Random_word(int k)
{
	int i, bit;
	unsigned long mot = 0;
	for(i= 0; i< k; i++)
    {
		bit= Bit_random();
		mot= (mot<<1) ^ bit;
    }
	return(mot);
}


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Tausworthe
///	Routine: DrawOne
///	Returns: double
///	Action : do a next seed and get next number!
///			based on MRGK3 random number generator
///			see premia
////////////////////////////////////////////////////
double ARM_RandUniform_Tausworthe::DrawOne()
{
	/// maximal number of combined generators
#define TAUS_MAX 10
	int L= 32;
	int J= 3;
	
	int i;
	static unsigned long u[TAUS_MAX];
	static unsigned long c[TAUS_MAX];
	unsigned long b;
	static int k[TAUS_MAX], q[TAUS_MAX];
	static int s[TAUS_MAX], r[TAUS_MAX], t[TAUS_MAX];
	unsigned long v= 0;
	
	/// First call to the sequence. Initialisation
	if (itsFirstUse) 
    {
		/// Choice of the parameters to have ME-CF generators (cf L'Ecuyer)
		k[0]= 31; q[0]= 13; s[0]= 12;
		r[0]= k[0]- q[0]; t[0]= k[0]- s[0];
		
		k[1]= 29; q[1]= 2; s[1]= 4;
		r[1]= k[1]- q[1]; t[1]= k[1]- s[1];
		
		k[2]= 28; q[2]= 3; s[2]= 17;
		r[2]= k[2]- q[2]; t[2]= k[2]-s[2];
		
		/// constant c : k bits to one and (L-k) bits to zero
		/// c[j]= 2^32 - 2^(L-k[j])
		c[0]= 4294967294ul;
		c[1]= 4294967288ul;
		c[2]= 4294967280ul;
		
		/// initialisation of each generator
		u[0]= 0; u[1]= 0; u[2]= 0;
		for(i= 0; i< J; i++)
		{
			/// The first k bits are chosen randomly
			u[i]= Random_word(k[i]);
			/// The next L-k bits are initially fixed to zero
			u[i] <<= (L- k[i]);  
			/// Now they are computed with the recurrence on u
			b= u[i] << q[i];
			b ^= u[i];
			b >>= k[i];
			u[i] ^= b;
		}
    }
	
	/// For each call to the sequence, computation of a new point
	for(i= 0; i< J; i++)
    { 
		/// Calculus of the next point for the J generators
		/// 6 steps explained by L'Ecuyer
		b=(u[i]<<q[i]) ^ u[i];         /// Steps 1 and 2
		b >>= t[i];                    /// Step 3
		u[i]= (u[i] & c[i]) << s[i];   /// Steps 4 et 5
		u[i] ^= b;                     /// Step 6
		/// Combination : XOR between the J generators
		v^= u[i];
    }
	
	/// Normalization by 1/2^32
	return v * 2.3283064365e-10;
#undef TAUS_MAX
}




/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Tausworthe
///	Routine: copy constructor, assignment operator
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_RandUniform_Tausworthe::ARM_RandUniform_Tausworthe( const ARM_RandUniform_Tausworthe& rhs )
:	ARM_UniformGenerator( rhs ), 
	itsFirstUse( rhs.itsFirstUse)
{}


ARM_RandUniform_Tausworthe& ARM_RandUniform_Tausworthe::operator =(const ARM_RandUniform_Tausworthe& rhs )
{
	if( this != & rhs )
	{
		ARM_UniformGenerator::operator =( rhs );
		itsFirstUse		= rhs.itsFirstUse;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Tausworthe
///	Routine: desctructor (necessary even for pure virtual case)
///	Action : draw a vector of uniform random numbers
////////////////////////////////////////////////////
ARM_RandUniform_Tausworthe::~ARM_RandUniform_Tausworthe()
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Tausworthe
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_RandUniform_Tausworthe::Clone() const
{
	return new ARM_RandUniform_Tausworthe(*this);
}




CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----