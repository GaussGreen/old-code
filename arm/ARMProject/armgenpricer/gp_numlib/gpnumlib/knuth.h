/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file Knuth.h
 *  \brief General file for the Knuth random gen algorithm
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPNUMLIB_KNUTH_H
#define _INGPNUMLIB_KNUTH_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpnumlib/uniform.h"

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////
/// \class ARM_RandUniform_Knuth
/// knuth is the same as NR ran3
//////////////////////////////
class ARM_RandUniform_Knuth : public ARM_UniformGenerator 
{
public:
	ARM_RandUniform_Knuth(): itsFirstUse(true) {}	/// no seed therefore... always the same sequence

	virtual string toString(const string& indent="", const string& nextIndent="") const{ return "Knuth Uniform Random Nb Generator"; }
	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb ){ itsFirstUse = true; }

	/// copy constructor, assignment operator, destructor
	ARM_RandUniform_Knuth( const ARM_RandUniform_Knuth& rhs );
	ARM_RandUniform_Knuth& operator =( const ARM_RandUniform_Knuth& rhs );
	virtual ~ARM_RandUniform_Knuth();

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
