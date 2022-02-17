/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file ranmar.h
 *  \brief General file for the RANMAR generator
 *	random nb generator
 *
 *	\author  A. Chaix
 *	\version 1.0
 *	\date October 2005
 */

#ifndef _INGPNUMLIB_RANMAR_H
#define _INGPNUMLIB_RANMAR_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpnumlib/uniform.h"

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////
/// \class ARM_RandUniform_Ranmar 
/// same as Park Miller
//////////////////////////////
class ARM_RandUniform_Ranmar : public ARM_UniformGenerator 
{	
private:
	int ij, kl;
public:
	ARM_RandUniform_Ranmar(int _ij = 1802, int _kl = 9373);

	virtual string toString(const string& indent="", const string& nextIndent="") const{ return "RANMAR Uniform Random Nb Generator"; }
	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb );
	

	/// copy constructor, assignment operator, destructor
	ARM_RandUniform_Ranmar( const ARM_RandUniform_Ranmar& rhs );
	ARM_RandUniform_Ranmar& operator =( const ARM_RandUniform_Ranmar& rhs );
	virtual ~ARM_RandUniform_Ranmar();

	/// core
	virtual	void draw( ARM_GP_Matrix& Matrix);
	
	/// standard ARM support
	virtual ARM_Object* Clone() const;

private:
	virtual double DrawOne();

};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/



