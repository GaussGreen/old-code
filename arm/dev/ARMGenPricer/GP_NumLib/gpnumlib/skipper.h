/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file skipper.h
 *
 *  \brief Generator Skipper
 *	\author  R. GUILLEMOT
 *	\version 1.0
 *	\date November 2005
 */

#ifndef _INGPNUMLIB_SKIPPER_H
#define _INGPNUMLIB_SKIPPER_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "manipulator.h"
#include "gpbase/assignop.h"
#include "typedef.h"

#include <string>
CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////
/// \struct ARM_SkipperFunction
///	\brief skip random number greater than a stddev value
////////////////////////////////////////////////

struct ARM_SkipperFunction : public ARM_ManipulatorMethod
{
public:
	ARM_SkipperFunction(const ARM_RandomGeneratorPtr&, double nbStdDev = 4.0);
	ARM_SkipperFunction(const ARM_SkipperFunction& rhs);
	ASSIGN_OPERATOR(ARM_SkipperFunction)
	virtual ARM_ManipulatorMethod* Clone() const;
	virtual ~ARM_SkipperFunction();

//	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb );
	virtual double operator()() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const;

//	ARM_RandomGeneratorPtr GetBaseGen() const {return itsRandomGen;}

private:
//	ARM_RandomGeneratorPtr itsRandomGen;
	double itsNbStdDev;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
