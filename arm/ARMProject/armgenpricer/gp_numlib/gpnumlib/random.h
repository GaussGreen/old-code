/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file random.h
 *
 *  \brief General file for the random nb generator
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _INGPNUMLIB_RANDOM_H
#define _INGPNUMLIB_RANDOM_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "typedef.h"
#include <string>
CC_USING_NS(std,string)

#include "gpbase/gplinalgtypedef.h"
#include "gpbase/rootobject.h"


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_AntitheticOneGen;

//////////////////////////////
/// \class base class ARM_RandomGenerator
///		abstract class for random number generators!
//////////////////////////////
class ARM_RandomGenerator : public ARM_RootObject
{
public:
	/// because of the design that makes composition of object
	/// we cannot know the type of an object by a dynamic cast
	/// hence the enum here!
	enum ARM_DistributionType
	{
		ARM_Uniform,
		ARM_Normal
	};

	/// interface for pointor like
	inline double operator*() { return itsCurrentValue; }
	inline void operator++() { itsCurrentValue = DrawOne(); }

	/// standard interface
	inline double operator()() { return DrawOne(); }
	std::vector<double>& operator()( size_t Number );
	void draw( std::vector<double>& Vec);
	virtual	void draw( ARM_GP_Matrix& Matrix);

	/// reset enables to reset the random number generator!
	/// gives the dimension and indicates the nb of points to generate!
	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& NbOfPointsList, size_t factorNb )	= 0;
	virtual void reset( const ARM_GP_T_Vector<size_t>& NbOfPointsLits, const ARM_GP_T_Vector<size_t>& factorNb);

	virtual size_t dim() const												= 0;

	/// constructor, copy constructor, assignment operator, destructor
	ARM_RandomGenerator();
	ARM_RandomGenerator( const ARM_RandomGenerator& rhs );
	ARM_RandomGenerator& operator=(const ARM_RandomGenerator& rhs );
	virtual ~ARM_RandomGenerator()						= 0;

	/// for stdDev computation
	virtual ARM_MomentFuncPtr StdDevComputationFunc() const;

	/// for compatibility reason for the MC Method
	virtual ARM_DistributionType GetDistributionType() const = 0;

	/// the root name of all random nb generator is ARM_RANDOMGEN
	virtual ARM_CLASS_NAME GetRootName() { return ARM_RANDOMGEN; }

private:
	/// declared as friend to allow easy access to the various method!
	friend class ARM_AntitheticOneGen;

	/// function to generate a number
	virtual double DrawOne()						= 0;
	virtual double DrawOne(int nbfactor)			{return DrawOne();}; // redéfini uniquement pour les quasi
	
	/// function to do a range check and a dimension check!
	virtual void ValidateResult( double result )	= 0;
	virtual void ValidateDim( size_t dim )			= 0;
	double itsCurrentValue;

protected:
	bool	IsThisQuasiRandom;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
