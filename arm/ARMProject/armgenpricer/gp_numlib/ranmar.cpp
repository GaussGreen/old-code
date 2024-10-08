/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file ranmar.cpp
 *  \brief 
 *	\author  A. Chaix
 *	\version 1.0
 *	\date October 2005
 */


#include "gpnumlib/ranmar.h"
#include "gpbase/gpmatrix.h"


CC_BEGIN_NAMESPACE( ARM )

static float u[98], c1, cd, cm;
static int i97, j97;
static bool testvar = false;

/// initialisation
void rmarin(int ij, int kl)
{
/*
C This is the initialization routine for the random number generator RANMAR()
C NOTE: The seed variables can have values between:    0 <= IJ <= 31328
C                                                      0 <= KL <= 30081
C The random number sequences created by these two seeds are of sufficient
C length to complete an entire calculation with. For example, if sveral
C different groups are working on different parts of the same calculation,
C each group could be assigned its own IJ seed. This would leave each group
C with 30000 choices for the second seed. That is to say, this random
C number generator can create 900 million different subsequences -- with
C each subsequence having a length of approximately 10^30.
C
C Use IJ = 1802 & KL = 9373 to test the random number generator. The
C subroutine RANMAR should be used to generate 20000 random numbers.
C Then display the next six random numbers generated multiplied by 4096*4096
C If the random number generator is working properly, the random numbers
C should be:
C           6533892.0  14220222.0  7275067.0
C           6172232.0  8354498.0   10633180.0
*/
	int i, j, k, l, ii, jj, m;
	float s, t;

	if (ij<0 || ij>31328 || kl<0 || kl>30081) {
	
		throw;
		/// puts("The first random number seed must have a value between 0 and 31328.");
		/// puts("The second seed must have a value between 0 and 30081.");
		/// exit(1);
	}

	i = (ij/177)%177 + 2;
	j = ij%177 + 2;
	k = (kl/169)%178 + 1;
	l = kl%169;

	for (ii=1; ii<=97; ii++) {
		s = 0.0;
		t = 0.5;
		for (jj=1; jj<=24; jj++) {
			m = (((i*j)%179)*k) % 179;
			i = j;
			j = k;
			k = m;
			l = (53*l + 1) % 169;
			if ((l*m)%64 >= 32) s += t;
			t *= 0.5;
		}
		u[ii] = s;
	}

	c1 = 362436.0 / 16777216.0;
	cd = 7654321.0 / 16777216.0;
	cm = 16777213.0 / 16777216.0;

	i97 = 97;
	j97 = 33;

	testvar = true;
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
//////////	ARM_RandUniform_Ranmar   ////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Ranmar
///	Routine: copy constructor, assignment operator
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_RandUniform_Ranmar::ARM_RandUniform_Ranmar( const ARM_RandUniform_Ranmar& rhs )
:	ARM_UniformGenerator( rhs ), ij(rhs.ij), kl(rhs.kl)
{}

ARM_RandUniform_Ranmar::ARM_RandUniform_Ranmar(int _ij, int _kl)
: ij(_ij),kl(_kl)
{
	/// initialisation
	rmarin(ij,kl);
}

ARM_RandUniform_Ranmar& ARM_RandUniform_Ranmar::operator =(const ARM_RandUniform_Ranmar& rhs )
{
	if( this != & rhs )
	{
		ARM_UniformGenerator::operator =( rhs );
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Ranmar
///	Routine: desctructor (necessary even for pure virtual case)
///	Action : draw a vector of uniform random numbers
////////////////////////////////////////////////////
ARM_RandUniform_Ranmar::~ARM_RandUniform_Ranmar()
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Ranmar
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_RandUniform_Ranmar::Clone() const
{
	return new ARM_RandUniform_Ranmar(*this);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Ranmar
///	Routine: reset
///	Returns: void
///	Action : initializes generator
/////////////////////////////////////////////////////////////////
void ARM_RandUniform_Ranmar::reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb )
{
	/// initialisation
	rmarin(ij,kl);
}


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Ranmar
///	Routine: DrawOne
///	Returns: double
///	Action : do a next seed and get next number!
///			based on MRGK5 random number generator
////////////////////////////////////////////////////
double ARM_RandUniform_Ranmar::DrawOne()
{
	float uni;
	uni = u[i97] - u[j97];
	if (uni < 0.0) uni += 1.0;
	u[i97] = uni;
	i97--;
	if (i97==0) i97 = 97;
	j97--;
	if (j97==0) j97 = 97;
	c1 -= cd;
	if (c1<0.0) c1 += cm;
	uni -= c1;
	if (uni<0.0) uni += 1.0;
	
	return uni;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_Ranmar
///	Routine: draw
///	Returns: void
///	Action : generates random numbers
/////////////////////////////////////////////////////////////////
void ARM_RandUniform_Ranmar::draw( ARM_GP_Matrix& Matrix)
{	
	float uni;

	if (testvar==false) 
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+": ARM_RandUniform_Ranmar:: generator has not been initialized");
	}

	ARM_GP_Matrix::iterator begin = Matrix.begin();
	ARM_GP_Matrix::iterator end   = Matrix.end();

	for( ARM_GP_Matrix::iterator iter = begin; iter!=end; ++iter )
	{
		uni = u[i97] - u[j97];
		if (uni < 0.0) uni += 1.0;
		u[i97] = uni;
		i97--;
		if (i97==0) i97 = 97;
		j97--;
		if (j97==0) j97 = 97;
		c1 -= cd;
		if (c1<0.0) c1 += cm;
		uni -= c1;
		if (uni<0.0) uni += 1.0;
		
		*iter = uni;
	}
	
}





CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----