/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file hammersley.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#include "firsttoinc.h"
#include "gpbase/removeidentifiedwarning.h"
#include "gpnumlib/hammersley.h"
#include "gpnumlib/primenb.h"
#include "gpbase/ostringstream.h"
#include "gpbase/gpvector.h"
#include <glob/expt.h>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////	ARM_Hammersley   ///////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_Hammersley
///	Routine: initBase
///	Returns: 
///	Action : init the base vector
////////////////////////////////////////////////////
void ARM_Hammersley::initBase( size_t dim )
{
	itsBase.resize(dim);
	for( size_t i=0; i<dim; ++i )
		itsBase[i] = ARM_PrimeNb::Prime(i);
}


////////////////////////////////////////////////////
///	Class  : ARM_QuasiRandom
///	Routine: constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_Hammersley::ARM_Hammersley(int firstSimulations)
:	ARM_QuasiRandom(firstSimulations), 
	itsCountor(1), 
	itsBase(1),
	itsTotalNbOfPoints(-1)
{
}

	


////////////////////////////////////////////////////
///	Class  : ARM_Hammersley
///	Routine: DrawAll
///	Returns: 
///	Action : get the next random number
////////////////////////////////////////////////////
void ARM_Hammersley::DrawAll()
{
	int internalCountor,digit;
	double base_inv;

	if( itsTotalNbOfPoints == -1 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "total nb of points not initialised!" );

	itsCurrentValues[itsCurrentBucketPos][itsCurrentPos][0] = ( (double) i_modp ( itsCountor, itsTotalNbOfPoints + 1 ) ) / ( (double) itsTotalNbOfPoints );

	for ( size_t i = 1; i<itsDim; ++i )
	{
		internalCountor	= itsCountor;
		base_inv		= 1.0 / ((double) itsBase[i]);
		itsCurrentValues[itsCurrentBucketPos][itsCurrentPos][i] = 0.0;
		
		while ( internalCountor != 0 )
		{
			digit = internalCountor % itsBase[i];
			itsCurrentValues[itsCurrentBucketPos][itsCurrentPos][i] += ( (double) digit ) * base_inv;
			base_inv = base_inv / ( (double) itsBase[i] );
			internalCountor = internalCountor / itsBase[i];
		}
	}
	++itsCountor;
}


////////////////////////////////////////////////////
///	Class  : ARM_Hammersley
///	Routine: toString
///	Returns: 
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_Hammersley::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << "hammersley Sequence with dim " << itsDim;
	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_Hammersley
///	Routine: SetDim
///	Returns: 
///	Action : reset the dimension!
////////////////////////////////////////////////////
void ARM_Hammersley::SetDim( size_t dim )
{
	initBase(dim);
	itsDim		= dim;
	itsCountor	= 1;
}


////////////////////////////////////////////////////
///	Class  : ARM_Hammersley
///	Routine: Copy constructor and operator=
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_Hammersley::ARM_Hammersley(const ARM_Hammersley& rhs)
:	ARM_QuasiRandom( rhs), 
	itsCountor( rhs.itsCountor ), 
	itsTotalNbOfPoints( rhs.itsTotalNbOfPoints ),
	itsBase( rhs.itsBase )

{}


ARM_Hammersley& ARM_Hammersley::operator=(const ARM_Hammersley& rhs )
{
	if(this!=&rhs)
	{
		ARM_QuasiRandom::operator=( rhs );
		itsCountor			= rhs.itsCountor;
		itsTotalNbOfPoints	= rhs.itsTotalNbOfPoints;
		itsBase				= rhs.itsBase;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_Hammersley
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_Hammersley::~ARM_Hammersley()
{}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_Hammersley
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_Hammersley::Clone() const
{
	return new ARM_Hammersley(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_Hammersley
///	Routine: reset
///	Returns: 
///	Action : reset the generator
////////////////////////////////////////////////////
void ARM_Hammersley::reset( size_t dim, const ARM_GP_T_Vector<size_t>& NbOfPointsList, size_t factorsNb )
{
	size_t NbOfPoints = 0;

	for(size_t i = 0; i < NbOfPointsList.size(); ++i)
		NbOfPoints += NbOfPointsList[i];

	itsTotalNbOfPoints = dim*NbOfPoints;
	ARM_QuasiRandom::reset(dim, NbOfPointsList, factorsNb);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_Hammersley
///	Routine: Clone
///	Returns: int
///	Action : returns the nonnegative remainder of integer division
//  Examples:
//        I         J     MOD  I_MODP   I_MODP Factorization
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
/////////////////////////////////////////////////////////////////
int ARM_Hammersley::i_modp ( int i, int j )
{
	if( j == 0 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "j == 0!" );
	int value = i % j;
	
	if ( value < 0 )
		value = value + abs ( j );

	return value;
}


CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----