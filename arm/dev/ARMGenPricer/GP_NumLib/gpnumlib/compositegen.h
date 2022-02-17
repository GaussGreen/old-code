/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file compositegen.h
 *
 *  \brief General file for the composite generator
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _INGPNUMLIB_COMPOSITEGEN_H
#define _INGPNUMLIB_COMPOSITEGEN_H

#include "gpbase/port.h"
#include "random.h"
#include "typedef.h"
#include "gpbase/env.h"
#include "gpbase/typedef.h"
#include "gpbase/functor.h"

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////
/// \base class ARM_NormalCompositeGen
///	a composite generator takes as an input a generator itself
///	it is usually a uniform generator...
/// uses the composite design pattern!
////////////////////////////////////////////////
class ARM_NormalCompositeGen: public ARM_RandomGenerator
{
public:
	enum ARM_DrawMethod { BoxMuller_Method, InvCumDistribution, InvCumDistributionFast, CondInvCumDistribution, Transposer, Skipper, MixteGen };

	ARM_NormalCompositeGen(
		const ARM_RandomGeneratorPtr& randomGen1 = ARM_RandomGeneratorPtr(0),
		const ARM_RandomGeneratorPtr& randomGen2 = ARM_RandomGeneratorPtr(0),
		ARM_DrawMethod method = InvCumDistribution,
		double nbStdDevs = 4.0,
		int firstNbTimes = 0,
		int firstNbDims = 0,
		int order = 0,
		int firstSimulations = 0);
	virtual string toString(const string& indent, const string& nextIndent) const;
	size_t dim() const { return 1; }

	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb );
	virtual void reset(const ARM_GP_T_Vector<size_t>& nbOfPointsList, const ARM_GP_T_Vector<size_t>& factorNb);

	virtual ARM_DistributionType GetDistributionType() const 
	{ return ARM_RandomGenerator::ARM_Normal; }

	/// copy constructor, assignment operator, destructor
	ARM_NormalCompositeGen( const ARM_NormalCompositeGen& rhs );
	ARM_NormalCompositeGen& operator =( const ARM_NormalCompositeGen& rhs );
	virtual ~ARM_NormalCompositeGen();

	/// standard ARM support
	virtual ARM_Object* Clone() const;

	/// for stdDev computation
	virtual ARM_MomentFuncPtr StdDevComputationFunc() const;


private:
	ARM_ManipulatorMethodPtr itsFunction;
	virtual double DrawOne();
	virtual void ValidateResult( double result ){}; /// nothing
	virtual void ValidateDim( size_t dim )	{};
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
