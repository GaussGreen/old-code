/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file genfactory.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */


#include "gpnumlib/randomgenfactory.h"

#include "gpbase/singleton.h"

/// quasi monte carlo
#include "gpnumlib/faure.h"
#include "gpnumlib/halton.h"
#include "gpnumlib/hammersley.h"
#include "gpnumlib/niederreiter.h"
#include "gpnumlib/sobol.h"
#include "gpnumlib/newsobol.h"
#include "gpnumlib/quasirandom.h"

/// random gen
#include "gpnumlib/antitheticgen.h"
#include "gpnumlib/compositegen.h"
#include "gpnumlib/knuth.h"
#include "gpnumlib/lecuyer.h"
#include "gpnumlib/mersenetwister.h"
#include "gpnumlib/mrgk3.h"
#include "gpnumlib/mrgk5.h"
#include "gpnumlib/parkmiller.h"
#include "gpnumlib/ran1.h"
#include "gpnumlib/ran2.h"
#include "gpnumlib/ran3.h"
#include "gpnumlib/ran4.h"
#include "gpnumlib/ranmar.h"
#include "gpnumlib/tausworthe.h"
#include "gpnumlib/nullrandgen.h"
#include "gpnumlib/argconvdefault.h"

/// other
#include "expt.h"
#include "gpbase/singleton.h"


CC_BEGIN_NAMESPACE( ARM )

/// macro to factorise code!
#define UNKNOWN_TYPE( msg ) \
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );


////////////////////////////////////////////////////
///	Class  : ARM_RandGenFactoryImp
///	Routine: CreateRandGen 
///	Returns: 
///	Action : create the most general random generator
////////////////////////////////////////////////////
ARM_RandomGenerator* ARM_RandGenFactoryImp::CreateRandGen( 
	BaseGenType baseGenType,
	TransformAlgo talgo,
	const ARM_RandomGeneratorPtr& baseGen1,
	const ARM_RandomGeneratorPtr& baseGen2,
	int	seed, 
	int dim,
	int factorDim,
	const ARM_GP_T_Vector<size_t>& nbOfPathsList,
	double nbStdDevs,
	int firstNbTimes,
	int firstNbDims,
	RandGenOrder order,
	int firstSimulations)
{
	ARM_RandomGenerator* randGen = NULL;

	switch(baseGenType)
	{
		/// quasi random generator
	case Faure:
		randGen = new ARM_FaureSeq(firstSimulations);
		break;
	case Halton:
		randGen = new ARM_Halton(firstSimulations);
		break;
	case Hammersley:
		randGen = new ARM_Hammersley(firstSimulations);
		break;
	case Niederreiter:
		randGen = new ARM_Niederreiter(firstSimulations);
		break;
	case Sobol:
		randGen = new ARM_Sobol(firstSimulations);
		break;
	case NewSobol:
		randGen = new ARM_NewSobol(seed, firstSimulations);
		break;
		
		/// random generator
	case Knuth:
		randGen = new ARM_RandUniform_Knuth;
		break;
	case Lecuyer:
		randGen = new ARM_RandUniform_Lecuyer;
		break;
	case Mersene:
		randGen = new ARM_FastMersenneTwister( seed );
		break;
	case MerseneStd:
		randGen = new ARM_OriginalMersenneTwister( seed );
		break;
	case MRGK3:
		randGen = new ARM_RandUniform_MRGK3;
		break;
	case MRGK5:
		randGen = new ARM_RandUniform_MRGK5;
		break;
	case ParkMiller:
		randGen = new ARM_RandUniform_Parkmiller;
		break;
	case NR_Ran1:
		randGen = new ARM_RandUniform_NRRan1(seed);
		break;
	case NR_Ran2:
		randGen = new ARM_RandUniform_NRRan2(seed);
		break;
	case NR_Ran3:
		randGen = new ARM_RandUniform_NRRan3(seed);
		break;
	case NR_Ran4:
		randGen = new ARM_RandUniform_NRRan4(seed);
		break;
	case Ranmar:
		randGen = new ARM_RandUniform_Ranmar;
		break;
	case Tausworthe:
		randGen = new ARM_RandUniform_Tausworthe;
		break;

	// Null Generator (used to tune the importance sampling)

	case Null:
		randGen = new ARM_NullRandGen(0.5);
		break;

	case UnknownBaseGenAlgorithm:
		{
			switch(talgo)
			{
			case BoxMuller:
				randGen = new ARM_NormalCompositeGen( baseGen1, baseGen2, ARM_NormalCompositeGen::BoxMuller_Method );
				break;
			case InvNormCum:
				randGen = new ARM_NormalCompositeGen( baseGen1, baseGen2, ARM_NormalCompositeGen::InvCumDistribution );
				break;
			case InvNormCumFast:
				randGen = new ARM_NormalCompositeGen( baseGen1, baseGen2, ARM_NormalCompositeGen::InvCumDistributionFast );
				break;
			case CondInvNormCum:
				randGen = new ARM_NormalCompositeGen( baseGen1, baseGen2, ARM_NormalCompositeGen::CondInvCumDistribution, nbStdDevs );
				break;
			case Transposer:
				randGen = new ARM_NormalCompositeGen( baseGen1, baseGen2, ARM_NormalCompositeGen::Transposer, nbStdDevs, firstNbTimes, firstNbDims, order, firstSimulations );
				break;
			case Skipper:
				randGen = new ARM_NormalCompositeGen( baseGen1, baseGen2, ARM_NormalCompositeGen::Skipper, nbStdDevs );
				break;
			case MixteGen:
				randGen = new ARM_NormalCompositeGen( baseGen1, baseGen2, ARM_NormalCompositeGen::MixteGen, nbStdDevs, firstNbTimes, firstNbDims );
				break;
			case AntitheticOne:
				randGen = new ARM_AntitheticOneGen(baseGen1);
				break;
			default:
				UNKNOWN_TYPE( "Unknown transform algorithm type" )
			}
		}
		break;

	
	default:
		UNKNOWN_TYPE( "Unknown base generator type" )
	}

	randGen->reset(factorDim*dim,nbOfPathsList,factorDim);

	return randGen;
}

////////////////////////////////////////////////////
///	Class  : ARM_RandGenFactoryImp
///	Routine: CreateRandGen 
///	Returns: 
///	Action : shortcut to create a simple rand gen
////////////////////////////////////////////////////
ARM_RandomGenerator* ARM_RandGenFactoryImp::CreateSimpleRandGen( 
	const string& genType1,
	const string& genType2,
	const string& algo1,
	const string& algo2,
	const string& randGenOrder,
	int firstNbTimes,
	int firstNbDims,
	bool isAntithetic
	)
{
	int genType1Int = ARM_ArgConv_BaseGenAlgoType.GetNumber(genType1);
	int genType2Int = ARM_ArgConv_BaseGenAlgoType.GetNumber(genType2);
	int algo1Int = ARM_ArgConv_TransformAlgoType.GetNumber(algo1);
	int algo2Int = ARM_ArgConv_TransformAlgoType.GetNumber(algo2);
	int randGenOrderInt = ARM_ArgConv_RandGenOrder.GetNumber(randGenOrder);

	ARM_RandomGenerator* mixteRandGen = ARM_RandGenFactory.Instance()->CreateSimpleRandGen(
		genType1Int,
		genType2Int,
		algo1Int,
		algo2Int,
		randGenOrderInt,
		firstNbTimes,
		firstNbDims,
		isAntithetic);

	return mixteRandGen;
}

////////////////////////////////////////////////////
///	Class  : ARM_RandGenFactoryImp
///	Routine: CreateRandGen 
///	Returns: 
///	Action : shortcut to create a simple rand gen
////////////////////////////////////////////////////
ARM_RandomGenerator* ARM_RandGenFactoryImp::CreateSimpleRandGen( 
	int genType1,
	int genType2,
	int algo1,
	int algo2,
	int randGenOrder,
	int firstNbTimes,
	int firstNbDims,
	bool isAntithetic
	)
{
	/// Base Generator 1
	ARM_RandomGeneratorPtr  pBaseRandomGen1( ARM_RandGenFactory.Instance()->CreateRandGen( 
		(ARM_RandGenFactoryImp::BaseGenType) genType2,
		ARM_RandGenFactoryImp::UnknownTransformAlgo ) );

	/// Base Generator 2
	ARM_RandomGeneratorPtr  pBaseRandomGen2(ARM_RandGenFactory.Instance()->CreateRandGen( 
		(ARM_RandGenFactoryImp::BaseGenType) genType1,
		ARM_RandGenFactoryImp::UnknownTransformAlgo ) );

	ARM_RandomGeneratorPtr normRandGen1(NULL); 
	ARM_RandomGeneratorPtr normRandGen2(NULL); 


	// Normal Rand Gen
	// Inv Norm Cum for the quasi random
	normRandGen1 = ARM_RandomGeneratorPtr( ARM_RandGenFactory.Instance()->CreateRandGen( 
		ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
		(ARM_RandGenFactoryImp::TransformAlgo) algo2,
		pBaseRandomGen1 ) );

	// Box Muller for the pseudo random
	normRandGen2 = ARM_RandomGeneratorPtr( ARM_RandGenFactory.Instance()->CreateRandGen( 
		ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
		(ARM_RandGenFactoryImp::TransformAlgo) algo1,
		pBaseRandomGen2 ) );

	/// tranposed random gen
	ARM_RandomGeneratorPtr transpRandGen2( ARM_RandGenFactory.Instance()->CreateRandGen( 
		ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
		ARM_RandGenFactoryImp::Transposer,
		ARM_RandomGeneratorPtr(normRandGen2),
		ARM_RandomGeneratorPtr(NULL),
		-1,
		1,
		1,
		ARM_GP_T_Vector<size_t>(1,10),
		4.0,
		0,
		0,
		(ARM_RandGenFactoryImp::RandGenOrder) randGenOrder));

	/// antithetic variates!
	ARM_RandomGeneratorPtr antiRandGen1, antiRandGen2;
	if (isAntithetic)
	{
		antiRandGen1 = ARM_RandomGeneratorPtr( ARM_RandGenFactory.Instance()->CreateRandGen( 
				ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
				ARM_RandGenFactoryImp::AntitheticOne,
				normRandGen1 ) );

		antiRandGen2 = ARM_RandomGeneratorPtr( ARM_RandGenFactory.Instance()->CreateRandGen(
				ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
				ARM_RandGenFactoryImp::AntitheticOne,
				transpRandGen2 ) );
	}
	else 
	{
		antiRandGen1 = normRandGen1;
		antiRandGen2 = transpRandGen2;
	}

	ARM_RandomGenerator* mixteRandGen = ARM_RandGenFactory.Instance()->CreateRandGen(
		ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
		ARM_RandGenFactoryImp::MixteGen,
		ARM_RandomGeneratorPtr(antiRandGen1),
		ARM_RandomGeneratorPtr(antiRandGen2),
		-1,
		1,
		1,
		ARM_GP_T_Vector<size_t>(1,10),
		4.0,
		firstNbTimes,
		firstNbDims);

	return mixteRandGen;
}

#undef UNKNOWN_TYPE

ARM_SingletonHolder<ARM_RandGenFactoryImp> ARM_RandGenFactory;

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

