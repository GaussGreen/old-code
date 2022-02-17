/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file mrgk5.h
 *  \brief General file for the MRGK5 algorithm
 *	random nb generator
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */

#ifndef _INGPNUMLIB_MRGK5_H
#define _INGPNUMLIB_MRGK5_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpnumlib/uniform.h"

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////
/// \class ARM_RandUniform_MRGK5 
/// same as Park Miller
//////////////////////////////
class ARM_RandUniform_MRGK5 : public ARM_UniformGenerator 
{
public:
	ARM_RandUniform_MRGK5(): itsFirstUse(true) {}	/// no seed therefore... always the same sequence

	virtual string toString(const string& indent="", const string& nextIndent="") const{ return "MRGK5 Uniform Random Nb Generator"; }
	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb ){ itsFirstUse = true; }

	/// copy constructor, assignment operator, destructor
	ARM_RandUniform_MRGK5( const ARM_RandUniform_MRGK5& rhs );
	ARM_RandUniform_MRGK5& operator =( const ARM_RandUniform_MRGK5& rhs );
	virtual ~ARM_RandUniform_MRGK5();

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



