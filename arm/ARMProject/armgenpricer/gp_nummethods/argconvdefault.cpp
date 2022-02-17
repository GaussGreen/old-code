/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: argconvddefault.cpp,v $
 * Revision 1.1  2004/03/25 16:39:43  ebenhamou
 * Initial revision
 *
 */


/*! \file GP_NumMethods/argconvdefault.cpp
 *
 *  \brief
 *
 *	\author  J-M Prié
 *	\version 1.0
 *	\date November 2004
 */

#include "gpnummethods/argconvdefault.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"
#include "gpnummethods/pdenumericalschemes.h"
#include "gpnummethods/pde3Fnumericalschemes.h"
#include "gpnummethods/impsampler.h"
#include "gpnummethods/pathscheme.h"
#include "gpnummethods/cfmethod.h"

CC_BEGIN_NAMESPACE( ARM )

ARGConvTable SamplerTable[] = 
{
	/// name				/// number
	"NormalCentredSampler",         ARM_SamplerBase::NormalCentred,
	"MeanRevertingSampler",			ARM_SamplerBase::MeanReverting,
	"DriftedMeanRevertingSampler",	ARM_SamplerBase::DriftedMeanReverting,
	"MarkovianDriftSampler",		ARM_SamplerBase::MarkovianDrift,
 	
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_SamplerType( SamplerTable, "Sampler Str -> Type" );
const ARM_ArgConvReverse ARM_ArgConvReverse_SamplerType( SamplerTable, "Sampler Type -> Str" );

ARGConvTable TruncatorTable[] = 
{
	/// name				    /// number
	"StdDevTruncator",			ARM_TruncatorBase::StandardDeviation,
	"ArrowDebreuTruncator",     ARM_TruncatorBase::ArrowDebreu,
 	
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_TruncatorType( TruncatorTable, "Truncator Str -> Type" );
const ARM_ArgConvReverse ARM_ArgConvReverse_TruncatorType( TruncatorTable, "Truncator Type -> Str" );

ARGConvTable SchedulerTable[] = 
{
	/// name				/// number
	"ConstantVarianceScheduler",                ARM_SchedulerBase::ConstantVariance,
	"ConstantVarianceMeanRevertingScheduler",   ARM_SchedulerBase::ConstantVarianceMeanReverting,
	"MultiRegimeScheduler",                     ARM_SchedulerBase::MultiRegime,
	"TimeStepPerYearScheduler",					ARM_SchedulerBase::TimeStepPerYear,
 	
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_SchedulerType( SchedulerTable, "Scheduler Str -> Type" );
const ARM_ArgConvReverse ARM_ArgConvReverse_SchedulerType( SchedulerTable, "Scheduler Type -> Str" );

ARGConvTable ReconnectorTable[] = 
{
	/// name				/// number
	"DoNothingReconnector", ARM_ReconnectorBase::DoNothing,
	"MeanReconnector",      ARM_ReconnectorBase::Mean,
	"VarianceReconnector",  ARM_ReconnectorBase::Variance,
 	
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_ReconnectorType( ReconnectorTable, "Reconnector Str -> Type" );
const ARM_ArgConvReverse ARM_ArgConvReverse_ReconnectorType( ReconnectorTable, "Reconnector Type -> Str" );

ARGConvTable SmootherTable[] = 
{
	/// name				/// number
	"DoNothingSmoother",    ARM_SmootherBase::DoNothing,
	"LinearSmoother",       ARM_SmootherBase::Linear,
	"QuadraticSmoother",    ARM_SmootherBase::Quadratic,
	"CubicSmoother",        ARM_SmootherBase::Cubic,
 	
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_SmootherType( SmootherTable, "Smoother Str -> Type" );
const ARM_ArgConvReverse ARM_ArgConvReverse_SmootherType( SmootherTable, "Smoother Type -> Str" );

ARGConvTable PDENumSchemeTable[] = 
{
	/// name				/// number
	"None",					ARM_PDENumericalScheme::None,
	"Explicit1F",           ARM_PDENumericalScheme::Explicit1F,
	"CN1F",	                ARM_PDENumericalScheme::CN1F,
	"Explicit2F",           ARM_PDENumericalScheme::Explicit2F,
	"ADI2F",	            ARM_PDENumericalScheme::ADI2F,
	"CS3F",					ARM_PDENumericalScheme::CS3F,
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

extern const ARM_ArgConv ARM_ArgConv_PDENumSchemeType( PDENumSchemeTable, "NumScheme Str -> Type" );
extern const ARM_ArgConvReverse ARM_ArgConvReverse_PDENumSchemeType( PDENumSchemeTable, "NumScheme Type -> Str" );

ARGConvTable PDEBoundConditionTable[] = 
{
	/// name				/// number
	"None",					ARM_PDE3FCraigSneydNumericalScheme::None,
	"Dirichlet",			ARM_PDE3FCraigSneydNumericalScheme::Dirichlet,
	"VonNeumann",	        ARM_PDE3FCraigSneydNumericalScheme::VonNeumann,
	"Belkheir",				ARM_PDE3FCraigSneydNumericalScheme::Belkheir,
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

extern const ARM_ArgConv ARM_ArgConv_PDEBoundConditionType( PDEBoundConditionTable, "BoundCondition Str -> Type" );
extern const ARM_ArgConvReverse ARM_ArgConvReverse_PDEBoundConditionType( PDEBoundConditionTable, "BoundCondition Type -> Str" );

ARGConvTable PDEGridTypeTable[] = 
{
	/// name				/// number
	"StdDev",				ARM_PDE3FCraigSneydNumericalScheme::StdDev,
	"Fixed",				ARM_PDE3FCraigSneydNumericalScheme::Fixed,
	"BiReg",				ARM_PDE3FCraigSneydNumericalScheme::BiReg,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

extern const ARM_ArgConv ARM_ArgConv_PDEGridType( PDEGridTypeTable, "PDEGridType Str -> Type" );
extern const ARM_ArgConvReverse ARM_ArgConvReverse_PDEGridType( PDEGridTypeTable, "PDEGridType Type -> Str" );

ARGConvTable ImpSamplerTypeTable[] = 
{
	/// name				/// number
	"Dummy",				ARM_ImpSampler::DummyImpSampler,
	"Proportional",         ARM_ImpSampler::PropImpSampler,
 	
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

extern const ARM_ArgConv ARM_ArgConv_ImpSamplerType( ImpSamplerTypeTable, "Imp Sampler Str -> Type" );
extern const ARM_ArgConvReverse ARM_ArgConvReverse_ImpSamplerType( ImpSamplerTypeTable, "Imp Sampler Type -> Str" );


ARGConvTable PathSchemeTypeTable[] = 
{
	/// name				/// number
	"Incremental",			ARM_PathScheme::Incremental,
	"BrownianBridge",       ARM_PathScheme::BrownianBridge,
	"IncAdaptative",		ARM_PathScheme::IncAdaptative,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

extern const ARM_ArgConv ARM_ArgConv_PathSchemeType( PathSchemeTypeTable, "Path Scheme Str -> Type" );
extern const ARM_ArgConvReverse ARM_ArgConvReverse_PathSchemeType( PathSchemeTypeTable, "Path Scheme Type -> Str" );


ARGConvTable CFmethodTypeTable[] = 
{
	/// name				/// number
	"ANALYTIC",				ARM_CFMethod::Analytic,
	"INTEGRAL",				ARM_CFMethod::Integral,
	"UNKNOWN",				ARM_CFMethod::Unknown,
 	
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

extern const ARM_ArgConv ARM_ArgConv_CFmethodType( CFmethodTypeTable, "CF method Str -> Type" );
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CFmethodType( CFmethodTypeTable, "CF method Type -> Str" );


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

