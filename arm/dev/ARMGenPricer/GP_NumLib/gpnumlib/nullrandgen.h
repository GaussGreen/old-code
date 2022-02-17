/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file mrgk5.h
 *  \brief General file for the Null Random Generator
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date December 2005
 */

#ifndef _INGPNUMLIB_NULLRANDGEN_H
#define _INGPNUMLIB_NULLRANDGEN_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpbase/assignop.h"
#include "gpnumlib/uniform.h"

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////
/// \class ARM_NullRandGen 
/// It will be used to tune the importance sampling
//////////////////////////////
class ARM_NullRandGen : public ARM_UniformGenerator 
{
public:
	ARM_NullRandGen(double seed) : itsSeed(seed) {}

	virtual string toString(const string& indent="", const string& nextIndent="") const{ return "Null Uniform Random Nb Generator"; }
	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb ) {}

	/// copy constructor, assignment operator, destructor
	ARM_NullRandGen( const ARM_NullRandGen& rhs );
	ASSIGN_OPERATOR(ARM_NullRandGen)
	virtual ~ARM_NullRandGen();

	/// standard ARM support
	virtual ARM_Object* Clone() const;

private:
	virtual double DrawOne();
	double itsSeed;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/



