/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file ran1.h
 *  \brief General file for the numerical recipees ran1
 *	random nb generator
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPNUMLIB_RAN1_H
#define _INGPNUMLIB_RAN1_H

/// the ran1 algorithm of numerical recipee is the same as the one in parkmillershufflgen!

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpnumlib/uniform.h"

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////
/// \class ARM_RandUniform_NRRan1 
//////////////////////////////
class ARM_RandUniform_NRRan1 : public ARM_UniformGenerator 
{
public:
	ARM_RandUniform_NRRan1(long seed=-1);
	virtual string toString(const string& indent="", const string& nextIndent="") const{ return "Numerical Recipes Ran1 Uniform Random Nb Generator"; }
	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb );

	/// copy constructor, assignment operator, destructor
	ARM_RandUniform_NRRan1( const ARM_RandUniform_NRRan1& rhs );
	ARM_RandUniform_NRRan1& operator =( const ARM_RandUniform_NRRan1& rhs );
	virtual ~ARM_RandUniform_NRRan1();

	/// standard ARM support
	virtual ARM_Object* Clone() const;

private:
	long itsCurrentSeed;
	long itsOriginalSeed;
	virtual double DrawOne();
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
