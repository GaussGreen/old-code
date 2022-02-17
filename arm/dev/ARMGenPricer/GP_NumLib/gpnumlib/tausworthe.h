/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file tausworthe.h
 *  \brief General file for the Taus Worthe algorithm
 *	random nb generator
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPNUMLIB_TAUSWORTHE_H
#define _INGPNUMLIB_TAUSWORTHE_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpnumlib/uniform.h"

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////
/// \class ARM_RandUniform_Tausworthe 
//////////////////////////////
class ARM_RandUniform_Tausworthe : public ARM_UniformGenerator 
{
public:
	ARM_RandUniform_Tausworthe(): itsFirstUse(true) {}	/// no seed therefore... always the same sequence

	virtual string toString(const string& indent="", const string& nextIndent="") const{ return "Tausworthe Uniform Random Nb Generator"; }
	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPoints, size_t factorNb ){ itsFirstUse = true; }

	/// copy constructor, assignment operator, destructor
	ARM_RandUniform_Tausworthe( const ARM_RandUniform_Tausworthe& rhs );
	ARM_RandUniform_Tausworthe& operator =( const ARM_RandUniform_Tausworthe& rhs );
	virtual ~ARM_RandUniform_Tausworthe();

	/// standard ARM support
	virtual ARM_Object* Clone() const;

private:
	bool itsFirstUse;
	virtual double DrawOne();
	int Bit_random();
	unsigned long Random_word(int k);
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
