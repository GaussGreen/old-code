/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file uniform.h
 *  \brief General file for the standard uniform 
 *	random nb generator
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */

#ifndef _INGPNUMLIB_UNIFORM_H
#define _INGPNUMLIB_UNIFORM_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpnumlib/random.h"

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////
/// \class ARM_UniformGenerator
///	just for 
//////////////////////////////
class ARM_UniformGenerator : public ARM_RandomGenerator
{
public:
	ARM_UniformGenerator();
	ARM_UniformGenerator( const ARM_UniformGenerator& );
	ARM_UniformGenerator& operator=( const ARM_UniformGenerator& );
	virtual ~ARM_UniformGenerator() = 0;

	virtual size_t dim() const { return 1; }
	virtual ARM_DistributionType GetDistributionType() const { return ARM_Uniform; }


private:
	virtual void ValidateResult( double result );
	virtual void ValidateDim( size_t dim );
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
