/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file ran2.h
 *  \brief General file for the numerical recipees ran2
 *	random nb generator
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */

#ifndef _INGPNUMLIB_RAN2_H
#define _INGPNUMLIB_RAN2_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpnumlib/uniform.h"

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////
/// \class ARM_RandUniform_NRRan2 
//////////////////////////////
class ARM_RandUniform_NRRan2 : public ARM_UniformGenerator 
{
public:
	ARM_RandUniform_NRRan2(long seed=-1);
	virtual string toString(const string& indent="", const string& nextIndent="") const{ return "Numerical Recipes Ran2 Uniform Random Nb Generator"; }
	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb );

	/// copy constructor, assignment operator, destructor
	ARM_RandUniform_NRRan2( const ARM_RandUniform_NRRan2& rhs );
	ARM_RandUniform_NRRan2& operator =( const ARM_RandUniform_NRRan2& rhs );
	virtual ~ARM_RandUniform_NRRan2();

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
