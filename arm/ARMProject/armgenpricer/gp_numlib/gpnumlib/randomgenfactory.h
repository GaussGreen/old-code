/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file genfactory.h
 *
 *  \brief General file to create factory class for random numbers
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _INGPNUMLIB_GENFACTORY_H
#define _INGPNUMLIB_GENFACTORY_H

#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T> class ARM_SingletonHolder;

class ARM_RandomGenerator;

struct ARM_RandGenFactoryImp
{
	enum BaseGenType
	{
		/// quasi random generator
		Faure,
		Halton,
		Hammersley,
		Niederreiter,
		Sobol,
		NewSobol,

		/// random generator
		Knuth,
		Lecuyer,
		Mersene,
		MerseneStd,
		MRGK3,
		MRGK5,
		NR_Ran1,
		NR_Ran2,
		NR_Ran3,
		NR_Ran4,
		ParkMiller,
		Ranmar,
		Tausworthe,
		Null,
		UnknownBaseGenAlgorithm
	};

	enum TransformAlgo
	{
		/// normal generator
		BoxMuller,
		InvNormCum,
		InvNormCumFast,
		CondInvNormCum,
		Transposer,
		Skipper,
		AntitheticOne,
		MixteGen,

		UnknownTransformAlgo
	};

	enum RandGenOrder
	{
		/// draw order
		BucketOrder,
		PathOrder
	};

	ARM_RandomGenerator* CreateRandGen( 
		BaseGenType baseGenType	= UnknownBaseGenAlgorithm,
		TransformAlgo talgo		= UnknownTransformAlgo,
		const ARM_RandomGeneratorPtr& baseGen1 = ARM_RandomGeneratorPtr(NULL),
		const ARM_RandomGeneratorPtr& baseGen2 = ARM_RandomGeneratorPtr(NULL),
		int	seed				= -1,
		int dim					= 1,
		int factorDim			= 1,
		const ARM_GP_T_Vector<size_t>& nbOfPointsList			= ARM_GP_T_Vector<size_t>(1,10),
		double nbStdDevs=4.0,
		int firstNbTimes=0,
		int firstNbDims=0,
		RandGenOrder order=BucketOrder,
		int firstSimulations=0);

	ARM_RandomGenerator* CreateSimpleRandGen(
		const string& genType1,
		const string& genType2,
		const string& algo1,
		const string& algo2,
		const string& randGenOrder,
		int firstNbTimes,
		int firstNbDims,
		bool isAntithetic
		);

	ARM_RandomGenerator* ARM_RandGenFactoryImp::CreateSimpleRandGen( 
		int genType1,
		int genType2,
		int algo1,
		int algo2,
		int randGenOrder,
		int firstNbTimes,
		int firstNbDims,
		bool isAntithetic
		);
private:
	/// to forbid client from using it except for the singleton holder
	ARM_RandGenFactoryImp() {};
	friend class ARM_SingletonHolder<ARM_RandGenFactoryImp>;
};

extern ARM_SingletonHolder<ARM_RandGenFactoryImp> ARM_RandGenFactory;

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
