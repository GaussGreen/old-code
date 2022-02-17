/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file mrgk5.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */


#include "gpnumlib/mrgk5.h"
#include <cmath>		/// for floor


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
////////////////////////////////////////////////////
//////////	ARM_RandUniform_MRGK5   ////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_MRGK5
///	Routine: DrawOne
///	Returns: double
///	Action : do a next seed and get next number!
///			based on MRGK5 random number generator
////////////////////////////////////////////////////
double ARM_RandUniform_MRGK5::DrawOne()
{
	static const double M1  = 4294949027.0;
	static const double M2  = 4294934327.0;
	static const double A12 = 1154721.0;
	static const double A14 = 1739991.0;
	static const double A15N= 1108499.0;
	static const double A21 = 1776413.0;
	static const double A23 = 865203.0;
	static const double A25N= 1641052.0;
	static const double NORM= 2.3283163396834613e-10;
	
	static double x10, x11, x12, x13, x14, x20, x21, x22, x23, x24;
	long k;
	double p1, p2,sample;
	
	/// First call to the sequence 
	if (itsFirstUse) 
	{
		/// Initialization
		x10= 231458761.;
		x11= 34125679.;
		x12= 45678213.;
		x13= 7438902.;
		x14= 957345.;
		
		x20= 57964412.;
		x21= 12365487.;
		x22= 77221456.;
		x23= 816403.;
		x24= 8488912.;
		itsFirstUse = false;
	}

	/// For each call to the sequence, computation of a new point
	/// First generator with Schrage method
	p1= A12*x13 - A15N*x10;
	
	if(p1> 0.0) 
		p1-= A14*M1;
	
	p1+= A14*x11;
	k= (long)floor(p1/M1);
	p1-= k*M1;
	
	if(p1< 0.0) 
		p1+= M1;
	
	x10= x11;
	x11= x12;
	x12= x13; 
	x13= x14;
	x14= p1;
	
	/// Second generator with Schrage method
	p2= A21*x24 - A25N*x20;
	
	if(p2> 0.0)
		p2-= A23*M2;
	
	p2+= A23*x22;
	k= (long)floor(p2/M2);
	p2-= k*M2;
	
	if(p2< 0.0)
		p2+= M2;
	
	x20= x21;
	x21= x22;
	x22= x23;
	x23= x24;
	x24= p2;
	
	/// Combination of the two generators
	if (p1<=p2)
		sample= (p1- p2+ M1)*NORM;
	else
		sample=(p1- p2)*NORM;

	return sample;
}




/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_MRGK5
///	Routine: copy constructor, assignment operator
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_RandUniform_MRGK5::ARM_RandUniform_MRGK5( const ARM_RandUniform_MRGK5& rhs )
:	ARM_UniformGenerator( rhs ), 
	itsFirstUse( rhs.itsFirstUse)
{}


ARM_RandUniform_MRGK5& ARM_RandUniform_MRGK5::operator =(const ARM_RandUniform_MRGK5& rhs )
{
	if( this != & rhs )
	{
		ARM_UniformGenerator::operator =( rhs );
		itsFirstUse		= rhs.itsFirstUse;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_MRGK5
///	Routine: desctructor (necessary even for pure virtual case)
///	Action : draw a vector of uniform random numbers
////////////////////////////////////////////////////
ARM_RandUniform_MRGK5::~ARM_RandUniform_MRGK5()
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandUniform_MRGK5
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_RandUniform_MRGK5::Clone() const
{
	return new ARM_RandUniform_MRGK5(*this);
}




CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----