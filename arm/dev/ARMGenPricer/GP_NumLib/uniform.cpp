/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file uniform.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */


#include "gpnumlib/uniform.h"
#include "gpbase/ostringstream.h"
#include <glob/expt.h>
#include "gpbase/env.h"

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_UniformGenerator
///	Routine: ValidateResult
///	Returns: built object
///	Action : function to validate a result!
////////////////////////////////////////////////////
void ARM_UniformGenerator::ValidateResult( double result )
{
#if defined(__GP_STRICT_VALIDATION)
	if(result<-1.0 ||result>1.0)
	{
		CC_Ostringstream os;
		os  << "Uniform number generator: result " << result
			<< " should be between -1 and 1!";
	    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str());
	}
#endif
}


////////////////////////////////////////////////////
///	Class  : ARM_UniformGenerator
///	Routine: ValidateDim
///	Returns: 
///	Action : ability to be used with a certain dim!
////////////////////////////////////////////////////
void ARM_UniformGenerator::ValidateDim( size_t dim )
{}		/// nothing as it is does not depend on dim



/////////////////////////////////////////////////////////////////
///	Class  : ARM_UniformGenerator
///	Routine: constructor, copy constructor, assignment operator
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////

ARM_UniformGenerator::ARM_UniformGenerator()
:	ARM_RandomGenerator()
{}


ARM_UniformGenerator::ARM_UniformGenerator( const ARM_UniformGenerator& rhs )
:	ARM_RandomGenerator( rhs )
{}


ARM_UniformGenerator& ARM_UniformGenerator::operator=(const ARM_UniformGenerator& rhs )
{
	if( this != & rhs )
	{
		ARM_RandomGenerator::operator =( rhs );
	}
	return *this;
}



////////////////////////////////////////////////////
///	Class  : ARM_UniformGenerator
///	Routine: desctructor (necessary even for pure virtual case)
///	Action : draw a vector of uniform random numbers
////////////////////////////////////////////////////
ARM_UniformGenerator::~ARM_UniformGenerator()
{}



CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----
