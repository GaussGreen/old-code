/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file transposer.h
 *
 *  \brief Conditional Normal Inverse Generator
 *	\author  R. GUILLEMOT
 *	\version 1.0
 *	\date November 2005
 */

#ifndef _INGPNUMLIB_CONDITIONAL_H
#define _INGPNUMLIB_CONDITIONAL_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpbase/assignop.h"
#include "gpbase/gpvector.h"
#include "invcum.h"
#include "typedef.h"

#include <string>
CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////
/// \struct ARM_CondInvCumFunction
///	\brief use of the conditional inverse of the cumulative function
////////////////////////////////////////////////

struct ARM_CondInvCumFunction : public ARM_InvCumFunction
{
public:
	ARM_CondInvCumFunction(const ARM_RandomGeneratorPtr&, invCumAlgoType algoType = INV_ERF_MORRO, double nbStdDev = 4.0);
	ARM_CondInvCumFunction(const ARM_CondInvCumFunction& rhs);
	ASSIGN_OPERATOR(ARM_CondInvCumFunction)
	virtual ARM_ManipulatorMethod* Clone() const;
	virtual ~ARM_CondInvCumFunction();

	virtual double operator()() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const;

	ARM_RandomGeneratorPtr GetBaseGen() const {return itsRandomGen;}

private:
	double itsNbStdDev;
	double itsMinProba;
	double itsMaxProba;
	double itsProbaCorrect;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
