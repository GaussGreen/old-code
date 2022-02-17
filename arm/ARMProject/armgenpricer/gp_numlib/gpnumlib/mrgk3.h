/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file mrgk3.h
 *  \brief General file for the MRGK3 algorithm
 *	random nb generator
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPNUMLIB_MRGK3_H
#define _INGPNUMLIB_MRGK3_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpnumlib/uniform.h"

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////
/// \class ARM_RandUniform_MRGK3
//////////////////////////////
class ARM_RandUniform_MRGK3 : public ARM_UniformGenerator 
{
public:
	ARM_RandUniform_MRGK3(): itsFirstUse(true) {}	/// no seed therefore... always the same sequence

	virtual string toString(const string& indent="", const string& nextIndent="") const{ return "MRGK3 Uniform Random Nb Generator"; }
	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb ){ itsFirstUse = true; }

	/// copy constructor, assignment operator, destructor
	ARM_RandUniform_MRGK3( const ARM_RandUniform_MRGK3& rhs );
	ARM_RandUniform_MRGK3& operator =( const ARM_RandUniform_MRGK3& rhs );
	virtual ~ARM_RandUniform_MRGK3();

	/// standard ARM support
	virtual ARM_Object* Clone() const;

private:
	bool itsFirstUse;
	virtual double DrawOne();
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
