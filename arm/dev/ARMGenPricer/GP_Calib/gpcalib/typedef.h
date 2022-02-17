/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *
 *	\file typedef.h
 *  \brief global typedef of namespace ARM for the gpcalib project
 *	\author  EBenhamou
 *	\version 1.0
 *	\date June 2004
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPCALIB_TYPEDEF_H
#define _INGPCALIB_TYPEDEF_H

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

/// gpbase
#include "gpbase/port.h"

/// gpinfra
#include "gpinfra/typedef.h"

///gpbase
#include "gpbase/valuetype.h"

/// gpcalib
#include "enumcalibration.h"

#include <utility>
#include <vector>
CC_USING_NS(std,vector)
#include <string>
CC_USING_NS(std,string)

/// forward declaration in global namespace (ARM kernel)
class ARM_Portfolio;

/// Start ARM namespace
CC_BEGIN_NAMESPACE( ARM )

/// forward declaration in ARM namespace
class ARM_CalibMethod;
class FunctionToSolve;
class FunctionToSolveWithDerivative;
class WeightedSquareFunc;
class ARM_ModelFitterDes;
class ARM_DensityFunctor;
//class ARM_MixtureDensityFunctor;
class ARM_VanillaSecurityDensity;


/// struct
struct ARM_VanillaArg;
struct ARM_CalibMethodType;
struct ARM_ModelFitterSolverType;
struct ARM_ModelFitterOptimizerType;

typedef ARM_CalibMethodType::MethodType ARM_MethodType;
typedef ARM_ModelFitterSolverType::SolverType ARM_SolverType;
typedef ARM_ModelFitterOptimizerType::OptimizerType ARM_OptimizerType;


/// vector declaration
typedef vector< ARM_CalibMethod* >                      ARM_CalibMethodVector;
typedef vector< ARM_VanillaArg*>                        ARM_VanillaArgVector;
typedef vector< ARM_StdPortfolio* >						ARM_PortfolioVector;

/// reference counted pointor
typedef ARM_CountedPtr< ARM_VanillaArg >                ARM_VanillaArgPtr;
typedef ARM_CountedPtr< FunctionToSolve >               FunctionToSolvePtr;
typedef ARM_CountedPtr< FunctionToSolveWithDerivative > FunctionToSolveWithDerivativePtr;
typedef ARM_CountedPtr< WeightedSquareFunc >            WeightedSquareFuncPtr;
typedef	ARM_CountedPtr< ARM_CalibMethod >				ARM_CalibMethodPtr;
typedef	ARM_CountedPtr< ARM_StdPortfolio >			    ARM_StdPortfolioPtr;
typedef	ARM_CountedPtr< ARM_ModelFitterDes >			ARM_ModelFitterDesPtr;
typedef ARM_CountedPtr< ARM_DensityFunctor >			ARM_DensityFunctorPtr;
typedef ARM_CountedPtr< ARM_VanillaSecurityDensity >	ARM_VanillaSecDensityPtr;

/// vector of ptr

typedef vector< ARM_DensityFunctorPtr >					ARM_DensityFunctorPtrVector;
typedef vector< ARM_VanillaSecDensityPtr >				ARM_VanillaSecDensityPtrVector;
typedef vector< ARM_VanillaArgPtr >						ARM_VanillaArgPtrVector;
typedef vector< ARM_CalibMethodPtr >                    ARM_CalibMethodPtrVector;
typedef vector< ARM_StdPortfolioPtr >                   ARM_StdPortfolioPtrVector;

typedef CC_NS(std,pair)< vector<string>,vector<ARM_GP_VALUE_TYPE> > ARM_RowInfo;

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
