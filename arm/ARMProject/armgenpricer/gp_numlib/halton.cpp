/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file halton.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#include "gpnumlib/halton.h"
#include "gpnumlib/primenb.h"
#include "gpbase/ostringstream.h"
#include "expt.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////	ARM_Halton   ///////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_Halton
///	Routine: constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_Halton::ARM_Halton(int firstSimulations)
:	ARM_QuasiRandom(firstSimulations), 
	itsCountor(1), 
	itsBase(1)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_Halton
///	Routine: initBase
///	Returns: 
///	Action : init the base vector
////////////////////////////////////////////////////
void ARM_Halton::initBase( size_t dim )
{
	itsBase.resize(dim);

	for( size_t i=0; i<dim; ++i )
		itsBase[i] = ARM_PrimeNb::Prime(i+1);
}

////////////////////////////////////////////////////
///	Class  : ARM_Halton
///	Routine: DrawAll
///	Returns: 
///	Action : get the next random number
////////////////////////////////////////////////////
void ARM_Halton::DrawAll()
{
	int internalCountor,digit;
	double base_inv;
	
	
	for( size_t  i=0; i<itsDim; i++ )
	{
		internalCountor		= itsCountor;
		base_inv			= 1.0/((double) itsBase[i]);
		itsCurrentValues[itsCurrentBucketPos][itsCurrentPos][i] = 0.0;
		
		while( internalCountor!= 0 )
		{
			digit = internalCountor % itsBase[i];
			itsCurrentValues[itsCurrentBucketPos][itsCurrentPos][i] += ((double) digit) * base_inv;
			base_inv /= ((double) itsBase[i]);
			internalCountor /= itsBase[i];
		}
	}
	++itsCountor;
}


////////////////////////////////////////////////////
///	Class  : ARM_Halton
///	Routine: toString
///	Returns: 
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_Halton::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << "Halton Sequence with dim " << itsDim;
	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_Halton
///	Routine: SetDim
///	Returns: 
///	Action : reset the dimension!
////////////////////////////////////////////////////
void ARM_Halton::SetDim( size_t dim )
{
	itsDim		= dim;
	itsCountor	= 1;
	initBase(dim);
}


////////////////////////////////////////////////////
///	Class  : ARM_Halton
///	Routine: Copy constructor and operator=
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_Halton::ARM_Halton(const ARM_Halton& rhs)
:	ARM_QuasiRandom( rhs), itsCountor( rhs.itsCountor ), itsBase( rhs.itsBase )
{}


ARM_Halton& ARM_Halton::operator=(const ARM_Halton& rhs )
{
	if(this!=&rhs)
	{
		ARM_QuasiRandom::operator=( rhs );
		itsCountor	= rhs.itsCountor;
		itsBase		= rhs.itsBase;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_Halton
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_Halton::~ARM_Halton()
{}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_Halton
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_Halton::Clone() const
{
	return new ARM_Halton(*this);
}

CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----