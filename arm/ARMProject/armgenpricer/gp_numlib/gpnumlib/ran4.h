/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file Ran4.h
 *  \brief General file for the numerical recipees Ran4
 *	random nb generator
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPNUMLIB_RAN4_H
#define _INGPNUMLIB_RAN4_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpnumlib/uniform.h"

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////
/// \class ARM_RandUniform_NRRan4 
//////////////////////////////
class ARM_RandUniform_NRRan4 : public ARM_UniformGenerator 
{
public:
	ARM_RandUniform_NRRan4(long seed=-1);
	virtual string toString(const string& indent="", const string& nextIndent="") const{ return "Numerical Recipes Ran4 Uniform Random Nb Generator"; }
	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb );

	/// copy constructor, assignment operator, destructor
	ARM_RandUniform_NRRan4( const ARM_RandUniform_NRRan4& rhs );
	ARM_RandUniform_NRRan4& operator =( const ARM_RandUniform_NRRan4& rhs );
	virtual ~ARM_RandUniform_NRRan4();

	/// standard ARM support
	virtual ARM_Object* Clone() const;

private:
	long itsCurrentSeed;
	long itsOriginalSeed;
	virtual double DrawOne();
	static void Psdes( unsigned long *lword, unsigned long *irword);
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
