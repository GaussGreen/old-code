/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file antitheticgen.h
 *
 *  \brief General file for the antithetic random nb generator
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPNUMLIB_ANTITHETICGEN_H
#define _INGPNUMLIB_ANTITHETICGEN_H

#include "gpbase/port.h"
#include "random.h"
#include "typedef.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////
/// \base class ARM_AntitheticOneGen
///	a class to take the antithetic variate from 
////////////////////////////////////////////////
class ARM_AntitheticOneGen : public ARM_RandomGenerator
{
public:
	ARM_AntitheticOneGen( const ARM_RandomGeneratorPtr& randomGen  );
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual size_t dim() const	{ return itsDim; }
	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb );

	/// copy constructor, assignment operator, destructor
	ARM_AntitheticOneGen( const ARM_AntitheticOneGen& rhs );
	ARM_AntitheticOneGen& operator =( const ARM_AntitheticOneGen& rhs );
	virtual ~ARM_AntitheticOneGen();

	/// standard ARM support
	virtual ARM_Object* Clone() const;

	/// std Dev computation (the default behavior is to compute the std dev normally!)
	virtual ARM_MomentFuncPtr StdDevComputationFunc() const;

	/// for distribution type checking
	virtual ARM_DistributionType GetDistributionType() const;

private:
	ARM_RandomGeneratorPtr itsRandGen;
	size_t itsCurrentIndex;
	ARM_GP_Vector* itsPreviousValues;
	size_t itsDim;
	size_t itsFactorDim;
	bool itsIsQuasiRandom;
	virtual double DrawOne();
	virtual void ValidateResult( double result );
	virtual void ValidateDim( size_t dim );
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
