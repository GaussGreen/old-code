/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file manipulator.h
 *
 *  \brief General file for the manipulator for a composite generator
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _INGPNUMLIB_MANIPULATOR_H
#define _INGPNUMLIB_MANIPULATOR_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpbase/functor.h"		/// for the def of the template for ARM_VoidToDbleFunctor
#include "gpbase/typedef.h"		/// for the def of the template for ARM_VoidToDbleFunctor
#include "typedef.h"
#include "gpbase/gpvector.h"

#include <vector>
CC_USING_NS(std,vector)

#include <string>
CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////
/// functors for the generation of
/// normal variables from uniform
/// pure abstract class!
/// no need to create the copy and assignment operator
/// as there is nothing to copy... and there should
///	be nothing to copy!
////////////////////////////////////////////////

struct ARM_ManipulatorMethod : public ARM_VoidToDbleFunctor
{
 	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb );
	virtual void reset(const ARM_GP_T_Vector<size_t>& nbOfPointsList, const ARM_GP_T_Vector<size_t>& factorNb);

	virtual string toString(const string& indent="", const string& nextIndent="") const = 0;
	ARM_ManipulatorMethod(const ARM_RandomGeneratorPtr&);
	ARM_ManipulatorMethod(const ARM_ManipulatorMethod&);

	virtual ARM_ManipulatorMethod* Clone() const						= 0;
	virtual ~ARM_ManipulatorMethod()									= 0;
	virtual ARM_RandomGeneratorPtr GetBaseGen() const {return itsRandomGen;};

protected:
	ARM_RandomGeneratorPtr itsRandomGen;

	virtual void resetwork(const ARM_GP_T_Vector<size_t>& nbOfPointsList, const ARM_GP_T_Vector<size_t>& factorNb){};
};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
