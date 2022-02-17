/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file typedef.h
 *  \brief global typedef of namespace ARM
 *	\author  J-M Prié
 *	\version 1.0
 *	\date November 2003
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPINFRA_TYPEDEF_H
#define _INGPINFRA_TYPEDEF_H

#include "gpbase/removeidentifiedwarning.h"
#include "enumInfra.h"
#include "gpbase/typedef.h"
#include "gpbase/port.h"
#include "gpbase/gplinalgtypedef.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/calibdirection.h"

#include <functional>
#include <vector>
CC_USING_NS(std,vector)
#include <string>
CC_USING_NS(std,string)

/// forward declaration in global namespace (ARM kernel)
class ARM_Object;
class ARM_Date;
class ARM_StdPortfolio;
class ARM_ZeroCurve;

/// forward declaration in the ARM_GP namespace
CC_BEGIN_NAMESPACE( ARM_GP )
	template<typename T,typename U=T, typename V=U> struct BinaryFunc;
CC_END_NAMESPACE()

/// typdef definition in the std namespace
typedef CC_NS( std, pointer_to_unary_function  )<double,double> pDbleUnaryFunc;
typedef CC_NS( std, pointer_to_binary_function )<double,double,double> pDbleBinaryFunc;
typedef CC_NS( std, binary_function )<double,double,double> STLDbleBinaryFunc;

/// Start ARM namespace
CC_BEGIN_NAMESPACE( ARM )

/// forward declaration in ARM namespace
class ARM_GenSecurity;
class ARM_DealDescription;
class ARM_GramFunctionArg;
class ARM_PricingModel;
class ARM_ModelParam;
class ARM_PricingStates;
class ARM_PayoffStates;
class ARM_TimeInfo;
struct ARM_AdditionalTimeInfo;
class ARM_ModelFitter;
class ARM_CellRefCall;
class ARM_ModelCall;
class ARM_PricingAdviser;
class ARM_PricingModelIR;
class ARM_ModelNameMap;
class ARM_GramFctorArgDict;
class ARM_CstManager;
class ARM_ExerciseBoundary;
class ARM_PricingStatesContext;
class ARM_Numeraire;
class ARM_NumMethod;
class ARM_InfCurv;
struct ARM_PricerInfo;
struct ARM_SwapRate;

/// struct
struct ARM_GramFctorArg;
struct ARM_ExpNode;
struct ARM_GramFctor;
struct ARM_MarketData_ManagerRep;
struct ARM_CalibrationTarget;
struct ARM_SamplerBase;
struct ARM_SchedulerBase;
struct ARM_ModelParamBump;
struct ARM_FXDigitType;

/// template struct
template <typename T, typename U> struct ARM_Pair;
typedef ARM_ModelParamType::ParamNb ARM_ParamType;
typedef ARM_ModelParamType::DataNb  ARM_DataType;
typedef ARM_CalibrationTarget::TargetType ARM_MktTargetType;
typedef ARM_ModelParamBump::BumpType  ARM_BumpParamType;
typedef ARM_FXDigitType::DigitType  ARM_DigitType;

/// vector declaration
typedef vector< ARM_GramFunctionArg >					    ARM_GramFunctionArgVector;
typedef vector< ARM_GramFctorArg >						    ARM_GramFctorArgVector;
typedef vector< ARM_PricingModel* >						    ARM_PricingModelVector;
typedef vector< ARM_ModelParam* >						    ARM_ModelParamVector;
typedef vector< ARM_ModelParamType::ParamNb >			    ARM_ParamTypeVector;

/// reference counted pointor
typedef ARM_CountedPtr< ARM_DealDescription >			    ARM_DealDescriptionPtr;
typedef ARM_CountedPtr< ARM_PricingModel >				    ARM_PricingModelPtr;
typedef ARM_CountedPtr< ARM_PricingStates >				    ARM_PricingStatesPtr;
typedef ARM_CountedPtr< ARM_PayoffStates >				    ARM_PayoffStatesPtr;
typedef ARM_CountedPtr< ARM_TimeInfo >					    ARM_TimeInfoPtr;
typedef ARM_CountedPtr< ARM_AdditionalTimeInfo >			ARM_AdditionalTimeInfoPtr;
typedef ARM_CountedPtr< ARM_ModelFitter >				    ARM_ModelFitterPtr;
typedef ARM_CountedPtr< ARM_ExpNode >						ARM_ExpNodePtr;
typedef ARM_CountedPtr< ARM_GramFctor >						ARM_GramFctorPtr;
typedef	ARM_CountedPtr< ARM_CellRefCall >					ARM_CellRefCallPtr;
typedef	ARM_CountedPtr< ARM_ModelCall >						ARM_ModelCallPtr;
typedef	ARM_CountedPtr< ARM_PricingAdviser >				ARM_PricingAdviserPtr;
typedef	ARM_CountedPtr< ARM_GenSecurity >					ARM_GenSecurityPtr;
typedef ARM_CountedPtr< ARM_MarketData_ManagerRep >		    ARM_MarketData_ManagerRepPtr;
typedef ARM_CountedPtr< ARM_StdPortfolio >					ARM_StdPortfolioPtr;
typedef ARM_CountedPtr< ARM_PricingModelIR >				ARM_PricingModelIRPtr;
typedef ARM_CountedPtr< ARM_ModelNameMap >					ARM_ModelNameMapPtr;
typedef ARM_CountedPtr< ARM_CstManager>						ARM_CstManagerPtr;
typedef ARM_CountedPtr< ARM_ExerciseBoundary>				ARM_ExerciseBoundaryPtr;
typedef ARM_CountedPtr< ARM_PricingStatesContext>			ARM_PricingStatesContextPtr;
typedef ARM_CountedPtr< ARM_ZeroCurve >						ARM_ZeroCurvePtr;
typedef ARM_CountedPtr< ARM_InfCurv >						ARM_InfCurvPtr;
typedef ARM_CountedPtr< ARM_Numeraire >						ARM_NumerairePtr;
typedef ARM_CountedPtr< ARM_NumMethod >						ARM_NumMethodPtr;
typedef ARM_CountedPtr< ARM_PricerInfo >					ARM_PricerInfoPtr;
typedef ARM_CountedPtr< ARM_ModelParam >					ARM_ModelParamPtr;
typedef ARM_CountedPtr< ARM_SwapRate >						ARM_SwapRatePtr;

typedef ARM_GP_T_Matrix<ARM_ExpNodePtr>						ARM_GP_NodeMatrix;
typedef	ARM_CountedPtr< ARM_GP_NodeMatrix >					ARM_GP_NodeMatrixPtr;

/// vector of counted ptr
typedef vector< ARM_TimeInfoPtr >							ARM_TimeInfoPtrVector;
typedef vector< ARM_PricingStatesContextPtr >				ARM_PricingStatesContextPtrVector;
typedef vector< ARM_AdditionalTimeInfoPtr >					ARM_AdditionalTimeInfoPtrVector;
typedef ARM_CountedPtr< ARM_PricingStatesContextPtrVector >	ARM_PricingStatesContextPtrVectorPtr; 
typedef vector< ARM_SwapRatePtr >							ARM_SwapRatePtrVector;
typedef vector< ARM_SwapRatePtrVector >						ARM_SwapRatePtrVectorVector;

/// templated class instances
typedef ARM_Pair< size_t,size_t>							ARM_RowColCoords;
typedef CC_NS( std, pair )< ARM_VectorPtr, ARM_AdditionalTimeInfoPtr> ARM_NodeInfo;
typedef CC_NS( std, pair )< ARM_RowColCoords, ARM_Date >	ARM_SharedNodeInfo;

// multi types dictionary
typedef ARM_GramFctorArgDict								ARM_MultiTypeDict;

// Sampler Base & Scheduler Base
typedef ARM_CountedPtr< ARM_SchedulerBase >				ARM_SchedulerBasePtr; 
typedef ARM_CountedPtr< ARM_SamplerBase >				ARM_SamplerBasePtr;

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
