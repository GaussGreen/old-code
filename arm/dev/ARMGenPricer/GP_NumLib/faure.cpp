/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file faure.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */


#include "gpbase/removeidentifiedwarning.h"
#include "gpnumlib/faure.h"
#include "gpnumlib/primenb.h"
#include "gpbase/ostringstream.h"
#include <glob/expt.h>
#include <cmath> /// for pow


CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////////////////
///               Gray Faure multidimentionnal sequence      ///
////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////
///
///	Class  : ARM_FaureSeq
///	Routine: constructor and DrawOne
///	Action : returns the smallest prime greater than or equal to N.
///  Purpose:
///    FAURE generates a new quasirandom Faure vector with each call.
///
///  Discussion:
///
///    This routine implements a method of H. Faure for computing
///    quasirandom numbers.  It is a merging and adaptation of
///    Bennett Fox's routines INFAUR and GOFAUR from ACM TOMS Algorithm 647.
///
///  Modified:
///    15 March 2003
///
///  Author:
///    John Burkardt
///
///  Reference:
///    H Faure,
///    Discrepance de suites associees a un systeme de numeration
///    (en dimension s),
///    Acta Arithmetica,
///    Volume XLI, 1982, pages 337-351, especially page 342.
///
///    Bennett Fox,
///    Algorithm 647:
///    Implementation and Relative Efficiency of Quasirandom 
///    Sequence Generators,
///    ACM Transactions on Mathematical Software,
///    Volume 12, Number 4, pages 362-376, 1986.
///
///  Parameters:
///
///    Input, int DIM_NUM, the spatial dimension, which should be
///    at least 2.
///
///    Input/output, int itsSeed, the seed, which can be used to index
///    the values.  On first call, set the input value of SEED to be 0
///    or negative.  The routine will automatically initialize data,
///    and set SEED to a new value.  Thereafter, to compute successive
///    entries of the sequence, simply call again without changing
///    SEED.  On the first call, if SEED is negative, it will be set
///    to a positive value that "skips over" an early part of the sequence
///    (This is recommended for better results).
///
///    Output, double QUASI(DIMEN), the next quasirandom vector.
//////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_QuasiRandom
///	Routine: constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_FaureSeq::ARM_FaureSeq(int firstSimulations)
:	ARM_QuasiRandom(firstSimulations), 
	itsSeed(-1), 
	itsCoeff(NULL), 
	itsYtemp(NULL), 
	itsHisum_save(-1), 
	itsHisum(0)
{
	Init();
}



	
////////////////////////////////////////////////////
///	Class  : ARM_FaureSeq
///	Routine: Init,DrawOne,DrawAll
///	Returns: 
///	Action : get the next random number
////////////////////////////////////////////////////
void ARM_FaureSeq::Init()
{
	itsQs		= ARM_PrimeNb::Prime_ge (itsDim);
	itsHisum_save	= -1;
	
	itsHisum = 3;
	itsSeed = ( int ) pow ( ( double ) itsQs, (( double ) ( itsHisum+1) ) ) - 1;
}



////////////////////////////////////////////////////
///	Class  : ARM_FaureSeq
///	Routine: DrawAll
///	Returns: 
///	Action : get the next random number
////////////////////////////////////////////////////
void ARM_FaureSeq::DrawAll()
{
	int i;
	int j;
	int k;
	int ktemp;
	int ltemp;
	int mtemp;

	double r;
	int ztemp;
	
	///  Negative or 0 seed indicates startup.
	if( itsSeed <= 0 )
	{
		Init();
	}
	///  Determine HISUM for the seed.
	else
	{
		itsHisum = 0;
		i = itsSeed;
		
		while ( itsQs <= i )
		{
			i = i / itsQs;
			itsHisum = itsHisum + 1;
		}
		
	}
	///
	///  Is it necessary to recompute the coefficient table?
	///
	if ( itsHisum_save != itsHisum )
	{
		itsHisum_save = itsHisum;
		///
		///  Apparently, trying to use NEW to set up a 2 dimensional array
		///  opens new vistas of obscurity I don't care to deal with right now.
		///  So let's do something STUPID and make COEF a 1 dimensional array.
		///  Thank goodness for progess!
		///
		itsCoeff = vector<int>((itsHisum+1)*(itsHisum+1),0);
		itsYtemp = vector<int>(itsHisum+1);
		

		for ( i = 1; i <= itsHisum; ++i )
		{
			itsCoeff[i*(itsHisum+1)+0] = 1;
			itsCoeff[i*(itsHisum+1)+i] = 1;
		}
		
		
		for ( j = 1; j <= itsHisum; j++ )
		{
			for ( i = j+1; i <= itsHisum; i++ )
			{
				itsCoeff[i*(itsHisum+1)+j] = 
					( itsCoeff[(i-1)*(itsHisum+1)+j] + itsCoeff[(i-1)*(itsHisum+1)+j-1] ) % itsQs;
			}
		}
	}
	
	///
	///  Find QUASI(1) using the method of Faure.
	///  SEED has a representation in base QS of the form: 
	///    Sum ( 0 <= J <= HISUM ) YTEMP(J) * QS**J
	///  We now compute the YTEMP(J)'s.
	///
	ktemp = ( int ) pow ( ( double ) itsQs, ( ( double ) ( itsHisum+1 )) );
	ltemp = itsSeed;
	
	for ( i = itsHisum; 0 <= i; i-- )
	{
		ktemp		= ktemp / itsQs;
		mtemp		= ltemp % ktemp;
		itsYtemp[i] = ( ltemp - mtemp ) / ktemp;
		ltemp		= mtemp;
	}
	///
	///  QUASI(K) has the form
	///    Sum ( 0 <= J <= HISUM ) YTEMP(J) / QS**(J+1)
	///  Compute QUASI(1) using nested multiplication.
	///
	r = ( ( double ) itsYtemp[itsHisum] );
	for ( i = itsHisum-1; 0 <= i; i-- )
	{
		r = ( ( double ) itsYtemp[i] ) + r / ( ( double ) itsQs );
	}
	
	itsCurrentValues[itsCurrentBucketPos][itsCurrentPos][0] = r / ( ( double ) itsQs );
	///
	///  Find components QUASI(2:DIM_NUM) using the Faure method.
	///
	for ( k = 1; k < itsDim; k++ )
	{
		r = 1.0 / ( ( double ) itsQs );
		
		for ( j = 0; j <= itsHisum; j++ )
		{
			ztemp = 0.0;
			for ( i = j; i <= itsHisum; i++ )
			{
				ztemp = ztemp + itsYtemp[i] * itsCoeff[i*(itsHisum+1)+j];
			}
			///
			///  New YTEMP(J) is:
			///    Sum ( J <= I <= HISUM ) ( old itsYtemp(i) * binom(i,j) ) mod QS.
			///
			itsYtemp[j] = ztemp % itsQs;
			itsCurrentValues[itsCurrentBucketPos][itsCurrentPos][k] += ( ( double ) itsYtemp[j] ) * r;
			r = r / ( ( double ) itsQs );
			
		}
		
	}
	///
	///  Update SEED.
	///
	itsSeed = itsSeed + 1;
}


////////////////////////////////////////////////////
///	Class  : ARM_FaureSeq
///	Routine: toString
///	Returns: 
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_FaureSeq::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << "Faure Sequence with dim " << itsDim;
	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_FaureSeq
///	Routine: SetDim
///	Returns: 
///	Action : reset the dimension!
////////////////////////////////////////////////////
void ARM_FaureSeq::SetDim( size_t dim )
{
	itsCoeff		= vector<int>(0);
	itsYtemp		= vector<int>(0);
	itsHisum_save	= -1;
	itsHisum		= 0;
	Init();
}


////////////////////////////////////////////////////
///	Class  : ARM_FaureSeq
///	Routine: Copy constructor and operator=
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_FaureSeq::ARM_FaureSeq(const ARM_FaureSeq& rhs)
:	ARM_QuasiRandom( rhs), 
	itsQs(			rhs.itsQs			),
	itsSeed(		rhs.itsSeed			), 
	itsHisum_save(	rhs.itsHisum_save	), 
	itsHisum(		rhs.itsHisum		),
	itsCoeff(		rhs.itsCoeff		), 
	itsYtemp(		rhs.itsYtemp		) 
{}

ARM_FaureSeq& ARM_FaureSeq::operator=(const ARM_FaureSeq& rhs )
{
	if(this!=&rhs)
	{
		ARM_QuasiRandom::operator=( rhs );
		itsQs			= rhs.itsQs;
		itsSeed			= rhs.itsSeed;
		itsHisum_save	= rhs.itsHisum_save;
		itsHisum		= rhs.itsHisum;
		itsCoeff		= rhs.itsCoeff;
		itsYtemp		= rhs.itsYtemp;
	}

	return *this;
}



////////////////////////////////////////////////////
///	Class  : ARM_FaureSeq
///	Routine: Destructor, Clone
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_FaureSeq::~ARM_FaureSeq()
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_FaureSeq
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_FaureSeq::Clone() const
{
	return new ARM_FaureSeq(*this);
}

CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----



















