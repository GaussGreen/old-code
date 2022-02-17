/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file knuth.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#include "gpnumlib/knuth.h"
#include <cmath>		/// for floor


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
////////////////////////////////////////////////////
//////////	ARM_RandUniform_Knuth   ////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Knuth
///	Routine: DrawOne
///	Returns: double
///	Action : do a next seed and get next number!
///			based on Knuth random number generator 
///			(see premia software). It is  based on MRG 
///			and it uses a substractive method
////////////////////////////////////////////////////
double ARM_RandUniform_Knuth::DrawOne()
{
	static long M		= 1000000000;
	static double InvM  = 0.000000001;
	static long SEED	= 161803398;
	
	/// Initialize the sequence with a positive seed
	static long alea= 1; 
	static int inc1, inc2;
	static long t_alea[56];
	long X_n, y_k;
	int i, ii, l;
		
	/// First call to the sequence
	if(itsFirstUse)
    {
		X_n= SEED- alea;
		X_n%= M;
		t_alea[55]= X_n;
		y_k= 1;
		
		/// Initialization of the table
		for(i= 1; i<= 54; i++)
		{
		ii= (21*i)%55; 
		
		/// 21 was chosen to alleviate initial
		/// nonrandomness problems
		t_alea[ii]= y_k;
		y_k= X_n - y_k;
		if(y_k < 0) 
			y_k+= M;
		X_n= t_alea[ii];
		}
		
		/// Randomization of the elements of the table
		for(l=  1; l<= 4; l++)
		{
			for(i= 1; i<= 55; i++)
			{
				t_alea[i]-= t_alea[1+(i+30)%55];
				if(t_alea[i] < 0) 
					t_alea[i]+= M;
			}
		}
		inc1= 0;
		inc2= 31;  /// 31 is a special value of Knuth : 31= 55-24
		alea= 1;
		itsFirstUse= false;
    }
	
	/// For each call to the sequence, computation of a new point
	if(++inc1 == 56) 
		inc1= 1;
	if(++inc2 == 56) 
		inc2= 1;
	/// Substractive method
	X_n= t_alea[inc1] - t_alea[inc2];
	
	if(X_n < 0) 
		X_n+= M;
	t_alea[inc1]= X_n;
	
	/// Normalized value
	return X_n*InvM;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Knuth
///	Routine: copy constructor, assignment operator
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_RandUniform_Knuth::ARM_RandUniform_Knuth( const ARM_RandUniform_Knuth& rhs )
:	ARM_UniformGenerator( rhs ), 
	itsFirstUse( rhs.itsFirstUse)
{}


ARM_RandUniform_Knuth& ARM_RandUniform_Knuth::operator =(const ARM_RandUniform_Knuth& rhs )
{
	if( this != & rhs )
	{
		ARM_UniformGenerator::operator =( rhs );
		itsFirstUse		= rhs.itsFirstUse;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Knuth
///	Routine: desctructor (necessary even for pure virtual case)
///	Action : draw a vector of uniform random numbers
////////////////////////////////////////////////////
ARM_RandUniform_Knuth::~ARM_RandUniform_Knuth()
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Knuth
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_RandUniform_Knuth::Clone() const
{
	return new ARM_RandUniform_Knuth(*this);
}


CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----