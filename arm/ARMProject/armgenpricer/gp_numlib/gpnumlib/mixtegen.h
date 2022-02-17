/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file mixtegen.h
 *
 *  \brief Make a mixte generator: it uses the first generator fo the first n dimensions
 *  \ then it uses the second. Typicaly we mix a quasi random and a traditional generators.
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _INGPNUMLIB_MIXTEGEN_H
#define _INGPNUMLIB_MIXTEGEN_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpbase/assignop.h"
#include "gpbase/gpvector.h"
#include "manipulator.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////
/// \struct ARM_InvCumFunction
///	\brief use of the inverse of the cumulative function
////////////////////////////////////////////////
struct ARM_MixteGen : public ARM_ManipulatorMethod 
{
	ARM_MixteGen(
		const ARM_RandomGeneratorPtr& randomGen1,
		const ARM_RandomGeneratorPtr& randomGen2,
		int firstNbTimes,
		int firstNbDims );
	virtual double operator()() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual void reset(const ARM_GP_T_Vector<size_t>& nbOfPointsList,  const ARM_GP_T_Vector<size_t>& factorNb );
	ASSIGN_OPERATOR(ARM_MixteGen)
	/// copy constructor, assignment operator, destructor
	ARM_MixteGen( const ARM_MixteGen& rhs );
	
	virtual ~ARM_MixteGen();

	/// standard ARM support
	virtual ARM_ManipulatorMethod* Clone() const;

	/// accessor
	virtual ARM_RandomGeneratorPtr GetBaseGen2() const { return itsRandomGen2;}

protected:
	virtual void resetwork(const ARM_GP_T_Vector<size_t>& nbOfPointsList, const ARM_GP_T_Vector<size_t>& factorNb);

protected:
	//	ARM_RandomGeneratorPtr itsRandomGen1;
	ARM_RandomGeneratorPtr itsRandomGen2;
	int itsFirstNbTimes;
	int itsFirstNbDims;
	int itsNbTimes;
	ARM_GP_T_Vector<size_t> itsNbOfPointsList;
	ARM_GP_T_Vector<size_t> itsStepNbFactor;
	mutable size_t itsCountor;
	mutable size_t itsDimCountor;
	mutable size_t itsTimeCountor;
	mutable size_t itsBucketCountor;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
