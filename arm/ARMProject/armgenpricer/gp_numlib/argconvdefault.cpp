/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file GP_NumLib/argconvdefault.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */
/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

#include "gpnumlib/argconvdefault.h"
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/odefunctions.h"
#include "gpnumlib/regression.h"

CC_BEGIN_NAMESPACE( ARM )

ARGConvTable BaseGenAlgoTypeTable[] = 
{
	/// name				/// number
	/// Quasi Random Nb
	"Faure",				ARM_RandGenFactoryImp::Faure,
	"Halton",				ARM_RandGenFactoryImp::Halton,
	"Hammersley",			ARM_RandGenFactoryImp::Hammersley,
	"Niederreiter",			ARM_RandGenFactoryImp::Niederreiter,
	"Sobol",				ARM_RandGenFactoryImp::Sobol,
	"NewSobol",				ARM_RandGenFactoryImp::NewSobol,

	//// WARNING : The New Sobol will be refered as the limit of QuasiRandom Generator.
	//// Please put the new random Generetor at the right position
	
	/// Standard Random Nb
	"Knuth",				ARM_RandGenFactoryImp::Knuth,
	"Lecuyer",			ARM_RandGenFactoryImp::Lecuyer,
	"Mersene",			ARM_RandGenFactoryImp::Mersene,
	"MerseneStd",		ARM_RandGenFactoryImp::MerseneStd,
	"MRGK3",			ARM_RandGenFactoryImp::MRGK3,
	"MRGK5",			ARM_RandGenFactoryImp::MRGK5,
	"NR_Ran1",			ARM_RandGenFactoryImp::NR_Ran1,
	"NR_Ran2",			ARM_RandGenFactoryImp::NR_Ran2,		
	"NR_Ran3",			ARM_RandGenFactoryImp::NR_Ran3,
	"NR_Ran4",			ARM_RandGenFactoryImp::NR_Ran4,
	"ParkMiller",		ARM_RandGenFactoryImp::ParkMiller,
	"Ranmar",			ARM_RandGenFactoryImp::Ranmar,
	"Tausworthe",		ARM_RandGenFactoryImp::Tausworthe,
	"Null",				ARM_RandGenFactoryImp::Null,
	"UnknownBaseGenAlgorithm",
						ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_BaseGenAlgoType( BaseGenAlgoTypeTable, "Base Generator Algorithm Type" );
const ARM_ArgConvReverse ARM_ArgConvReverse_BaseGenAlgoType( BaseGenAlgoTypeTable, "Base Generator Algorithm Type");

ARGConvTable TransformAlgoTypeTable[] = 
{
	/// name				/// number
	"BoxMuller",			ARM_RandGenFactoryImp::BoxMuller,
	"InvNormCum",			ARM_RandGenFactoryImp::InvNormCum,
	"InvnormCumFast",		ARM_RandGenFactoryImp::InvNormCumFast,
	"CondInvnormCum",		ARM_RandGenFactoryImp::CondInvNormCum,
	"AntitheticOne",		ARM_RandGenFactoryImp::AntitheticOne,
	"Transposer",			ARM_RandGenFactoryImp::Transposer,
	"Skipper",				ARM_RandGenFactoryImp::Skipper,
	"MixteGen",				ARM_RandGenFactoryImp::MixteGen,
	"UnknownTransformAlgo",	ARM_RandGenFactoryImp::UnknownTransformAlgo,
 	
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_TransformAlgoType( TransformAlgoTypeTable, "Tranform Random Gen Algo Type" );
const ARM_ArgConvReverse ARM_ArgConvReverse_TransformAlgoType( TransformAlgoTypeTable, "Tranform Random Gen Algo Type" );


ARGConvTable RandGenOrderTable[] = 
{
	/// name				/// number
	"BucketOrder",			ARM_RandGenFactoryImp::BucketOrder,
	"PathOrder",			ARM_RandGenFactoryImp::PathOrder,
 	
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

extern const ARM_ArgConv ARM_ArgConv_RandGenOrder( RandGenOrderTable, "RandGen Order Str -> Order" );
extern const ARM_ArgConvReverse ARM_ArgConvReverse_RandGenOrder( RandGenOrderTable, "RandGen Order -> Order" );

ARGConvTable ODESolverTable[] = 
{
	/// name				/// number
	"RK4Constant",			ARM_ODEFunc::RK4Constant,
	"RK5Adaptative",		ARM_ODEFunc::RK5Adaptative,
 	
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

extern const ARM_ArgConv ARM_ArgConv_ODESolverType( ODESolverTable, "Solver Type Name -> Order" );
extern const ARM_ArgConvReverse ARM_ArgConvReverse_ODESolverType( ODESolverTable, "Order -> Solver Type Name" );

ARGConvTable RegModeTable[] = 
{
	/// name				/// number
	"LS",					ARM_Regression::LS,
	"LOESS",				ARM_Regression::LOESS,
 	
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

extern const ARM_ArgConv ARM_ArgConv_RegMode( RegModeTable, "Reg Mode Str -> Type" );
extern const ARM_ArgConvReverse ARM_ArgConvReverse_RegMode( RegModeTable, "Reg Mode -> Str" );


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

