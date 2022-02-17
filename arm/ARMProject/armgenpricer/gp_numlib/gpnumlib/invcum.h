/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file invcum.h
 *
 *  \brief General file for the inverse of a cumulative function of a normal distribution
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _INGPNUMLIB_INVCUM_H
#define _INGPNUMLIB_INVCUM_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpbase/assignop.h"
#include "manipulator.h"
#include "typedef.h"
#include "normalinvcum.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////
/// \struct ARM_InvCumFunction
///	\brief use of the inverse of the cumulative function
////////////////////////////////////////////////
struct ARM_InvCumFunction: public ARM_ManipulatorMethod 
{
	typedef double (*invCumAlgoFuncPtr)(double);
	enum invCumAlgoType
	{
		INV_ERF_MORRO,
		INV_ERF
	};

	ARM_InvCumFunction(const ARM_RandomGeneratorPtr& randomGen, invCumAlgoType algoType = INV_ERF_MORRO );
	virtual double operator()() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const;
//	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb );
	ASSIGN_OPERATOR(ARM_InvCumFunction)
	/// copy constructor, assignment operator, destructor
	ARM_InvCumFunction( const ARM_InvCumFunction& rhs );
	
	virtual ~ARM_InvCumFunction();

	/// standard ARM support
	virtual ARM_ManipulatorMethod* Clone() const;

	/// accessor
//	virtual ARM_RandomGeneratorPtr GetBaseGen() const { return itsRandomGen;}

protected:
//	ARM_RandomGeneratorPtr itsRandomGen;
	invCumAlgoFuncPtr itsFunc;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
