/*
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_calculator.h,v $
 * Revision 1.1  2004/03/02 15:08:43  ebenhamou
 * Initial version
 *
 */

#ifndef ARMLOCAL_GP_TARNANDCO_CALCULATORS_H
#define ARMLOCAL_GP_TARNANDCO_CALCULATORS_H

#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_local_gp_genericaddin.h"
#include "ARM_result.h"
#include <GP_Infra\gpinfra\gramfunctorarg.h>
#include <gpbase/gpvector.h>
#include <gpbase/gpmatrix.h>

using ARM::ARM_GramFctorArg;

class ARM_ReferenceValue;

/****************************************************************************
						TARN FX CALCULATOR
*****************************************************************************/

class ARM_TARNFXCalculator_CreateFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_TARNFXCalculator_CreateFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};
class ARM_TARNFXCalculator_InitFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_TARNFXCalculator_InitFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

/****************************************************************************
						INDIAN TARN FX CALCULATOR
*****************************************************************************/

class ARM_TARNCalculatorIndian_CreateFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_TARNCalculatorIndian_CreateFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};


#endif

//-----------------------------------------------------------------------------
/*---- End of file ----*/
