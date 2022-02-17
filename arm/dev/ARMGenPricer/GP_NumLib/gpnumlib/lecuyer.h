/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file lecuyer.h
 *  \brief General file for the L'ECUYER
 *	random nb generator
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPNUMLIB_LECUYER_H
#define _INGPNUMLIB_LECUYER_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpnumlib/uniform.h"

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////
/// \class ARM_RandUniform_Lecuyer 
//////////////////////////////
class ARM_RandUniform_Lecuyer : public ARM_UniformGenerator 
{
public:
	ARM_RandUniform_Lecuyer(): itsFirstUse(true) {}	/// no seed therefore... always the same sequence

	virtual string toString(const string& indent="", const string& nextIndent="") const{ return "L'Ecuyer with Bayes & Durham shuffling Uniform Random Nb Generator"; }
	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb ){ itsFirstUse = true; }

	/// copy constructor, assignment operator, destructor
	ARM_RandUniform_Lecuyer( const ARM_RandUniform_Lecuyer& rhs );
	ARM_RandUniform_Lecuyer& operator =( const ARM_RandUniform_Lecuyer& rhs );
	virtual ~ARM_RandUniform_Lecuyer();

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
