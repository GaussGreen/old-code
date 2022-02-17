/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file Niederreiter.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpnumlib/niederreiter.h"
#include "gpnumlib/niederreiterinit.h"
#include "gpbase/ostringstream.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////	ARM_Niederreiter   ///////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_QuasiRandom
///	Routine: constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_Niederreiter::ARM_Niederreiter(int firstSimulations)
:	ARM_QuasiRandom(firstSimulations), 
	itsFirstUse(true)
{}




////////////////////////////////////////////////////
///	Class  : ARM_Niederreiter
///	Routine: DrawAll
///	Returns: 
///	Action : get the next random number
////////////////////////////////////////////////////
void ARM_Niederreiter::DrawAll()
{
	int i, j;
	unsigned long saut, gray, dim;
	static double facteur;
	static unsigned long initial_d, initialX_n[CC_NS(ARM_NIEDERREITER,DIM_MAX_NIED)+1];
	
	/// First call to the sequence 
	if (itsFirstUse) 
    {
		
		/// Initialization of initX_n[]
		for (i=1; i<=CC_NS(ARM_NIEDERREITER,DIM_MAX_NIED); i++) 
			initialX_n[i]= 0;
		
		facteur= 1.0/(double)(1UL << (CC_NS(ARM_NIEDERREITER,BIT_MAX_NIED)+1));
		saut = 1L << 12;
		initial_d = saut;
		
		/// Gray code of saut
		gray = saut^(saut >> 1); 
		for (i=0; i<=CC_NS(ARM_NIEDERREITER,BIT_MAX_NIED); i++) 
		{
			if (gray == 0) 
				break;
			for (j= 1; j<= CC_NS(ARM_NIEDERREITER,DIM_MAX_NIED); j++) 
			{
				if ((gray & 1) == 0) 
					break;
				/// XOR sum
				initialX_n[j] ^= CC_NS(ARM_NIEDERREITER,C)[i][j];
				
			}
			gray >>= 1;
		}
		itsFirstUse = false;
	}
	
	/// Calculation of a new quasi-random vector on each call
	dim= initial_d++;
	
	/// Research of the rightmost 0 bit
	for (i=0; i<=CC_NS(ARM_NIEDERREITER,BIT_MAX_NIED); i++) 
	{
		if ((dim & 1) == 0) 
			break;
		dim >>= 1;
	}
	/// Computation of the term n from the term n-1
	for (j=1; j<= itsDim; j++) 
	{
		itsCurrentValues[itsCurrentBucketPos][itsCurrentPos][j-1]= (double)initialX_n[j]*facteur;
		initialX_n[j] ^= CC_NS(ARM_NIEDERREITER,C)[i][j];
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_Niederreiter
///	Routine: toString
///	Returns: 
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_Niederreiter::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << "Niederreiter Sequence with dim " << itsDim;
	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_Niederreiter
///	Routine: SetDim
///	Returns: 
///	Action : reset the dimension!
////////////////////////////////////////////////////
void ARM_Niederreiter::SetDim( size_t dim )
{
	itsDim		= dim;
	itsFirstUse = true;
}


////////////////////////////////////////////////////
///	Class  : ARM_Niederreiter
///	Routine: Copy constructor and operator=
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_Niederreiter::ARM_Niederreiter(const ARM_Niederreiter& rhs)
:	ARM_QuasiRandom( rhs), itsFirstUse(rhs.itsFirstUse)
{}

ARM_Niederreiter& ARM_Niederreiter::operator=(const ARM_Niederreiter& rhs )
{
	if(this!=&rhs)
	{
		ARM_QuasiRandom::operator=( rhs );
		itsFirstUse = rhs.itsFirstUse;
	}
	return *this;
}




////////////////////////////////////////////////////
///	Class  : ARM_Niederreiter
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_Niederreiter::~ARM_Niederreiter()
{}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_Niederreiter
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_Niederreiter::Clone() const
{
	return new ARM_Niederreiter(*this);
}


CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----