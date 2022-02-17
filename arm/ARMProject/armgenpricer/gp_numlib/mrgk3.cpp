/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file mrgk3.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#include "gpnumlib/mrgk3.h"
#include <cmath>		/// for floor


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
////////////////////////////////////////////////////
//////////	ARM_RandUniform_MRGK3   ////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_MRGK3
///	Routine: DrawOne
///	Returns: double
///	Action : do a next seed and get next number!
///			based on MRGK3 random number generator
///			see premia
////////////////////////////////////////////////////
double ARM_RandUniform_MRGK3::DrawOne()
{
	static double M1= 4294967087.0;
	static double M2= 4294944443.0;
	static double A12= 1403580.0;
	static double A13N= 810728.0;
	static double A21= 527612.0;
	static double A23N= 1370589.0;
	static double NORM= 2.328306549295728e-10;
	static double x10, x11, x12, x20, x21, x22;
	long k;
	double p1, p2;
	
	/// First call to the sequence 
	if (itsFirstUse) 
    {
		/// Initialization
		x10=231458761.;
		x11=34125679.;
		x12=45678213.;
		x20=57964412.;
		x21=12365487.;
		x22=77221456.;
		itsFirstUse = false;
    }
	
	/// For each call to the sequence, computation of a new point
	/// First generator
	p1= A12*x11 - A13N*x10;
	k= (long)floor(p1/M1); /// TOCHECK
	p1-= k*M1;
	
	if(p1 < 0.0) 
		p1+= M1;
	
	x10= x11;
	x11= x12;
	x12= p1;
	
	/// Second generator
	p2= A21*x22 - A23N*x20;
	k= (long)floor(p2/M2); /// TOCHECK
	p2-= k*M2;
	
	if(p2 < 0.0) 
		p2+= M2;
	
	x20= x21;
	x21= x22;
	x22= p2;
	
	
	/// Combination of the two generators 
	if (p1< p2) 
		return (p1- p2+ M1)*NORM;
	else 
		return (p1- p2)*NORM;
}




/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_MRGK3
///	Routine: copy constructor, assignment operator
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_RandUniform_MRGK3::ARM_RandUniform_MRGK3( const ARM_RandUniform_MRGK3& rhs )
:	ARM_UniformGenerator( rhs ), 
	itsFirstUse( rhs.itsFirstUse)
{}


ARM_RandUniform_MRGK3& ARM_RandUniform_MRGK3::operator =(const ARM_RandUniform_MRGK3& rhs )
{
	if( this != & rhs )
	{
		ARM_UniformGenerator::operator =( rhs );
		itsFirstUse		= rhs.itsFirstUse;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_MRGK3
///	Routine: desctructor (necessary even for pure virtual case)
///	Action : draw a vector of uniform random numbers
////////////////////////////////////////////////////
ARM_RandUniform_MRGK3::~ARM_RandUniform_MRGK3()
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_MRGK3
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_RandUniform_MRGK3::Clone() const
{
	return new ARM_RandUniform_MRGK3(*this);
}




CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----