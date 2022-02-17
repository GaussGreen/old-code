/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file argconvdefault.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

#include "gpbase/curve.h"

#include "gpcalib/argconvdefault.h"
#include "gpcalib/calibmethod.h"
#include "gpcalib/modelfitterdes.h"
#include "gpcalib/basket.h"
#include "gpcalib/densityfunctors.h"


CC_BEGIN_NAMESPACE( ARM )



ARGConvTable CalibMethodTable[] =
{	
    /// Type Name       /// number
    "Bootstrap1D",      ARM_CalibMethodType::Bootstrap1D,
    "BootstrapND",      ARM_CalibMethodType::BootstrapND,
    "Optimize",         ARM_CalibMethodType::Optimize,
    "Optimize1D",       ARM_CalibMethodType::Optimize1D,
    "Numerical",		ARM_CalibMethodType::Numerical,
	"HW2FOnly",			ARM_CalibMethodType::HW2FOnly,
	"Unknown",          ARM_CalibMethodType::Unknown,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CalibMethod( CalibMethodTable, "CalibMethod" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CalibMethod( CalibMethodTable, "CalibMethod" );

ARGConvTable TargetFuncTable[] =
{	
    /// Type Name           /// number
    "PRICE_TAR",            ARM_CalibrationTarget::PriceTarget,
    "IMPVOL_TAR",			ARM_CalibrationTarget::ImpliedVolatilityTarget,
    "UNKNOWN_TAR",          ARM_CalibrationTarget::UnknownTarget,


	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_TargetFuncMethod(TargetFuncTable, "CalibrationTargetFunc" );
const ARM_ArgConvReverse ARM_ArgConvReverse_TargetFuncMethod(TargetFuncTable, "CalibrationTargetFunc" );

ARGConvTable CalibDirectionMethodTable[] =
{	
    /// Type Name       /// number
    "CalibDirection_Forward",       CalibDirection_Forward,
    "CalibDirection_Backward",      CalibDirection_Backward,
    "CalibDirection_None",          CalibDirection_None,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CalibDirectionMethod( CalibDirectionMethodTable, "CalibDirection" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CalibDirectionMethod( CalibDirectionMethodTable, "CalibDirection" );

ARGConvTable CalibDirection2DMethodTable[] =
{	
    /// Type Name       /// number
    "Forward",       CalibDirection_Forward,
    "Backward",      CalibDirection_Backward,
    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CalibDirection2DMethod( CalibDirection2DMethodTable, "CalibDirection2D" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CalibDirection2DMethod( CalibDirection2DMethodTable, "CalibDirection2D" );


ARGConvTable SolverTypeTable[] =
{	
    /// Type Name				/// number
    "STD_NR",					ARM_ModelFitterSolverType::NewtonRaphson,
	"NR_WITHRETRIAL",			ARM_ModelFitterSolverType::NewtonRaphsonWithRetrial,
    "NR_WITHSMOOTHING",			ARM_ModelFitterSolverType::SmoothNewthonRhaphson,
	"NR_NOTHROW",				ARM_ModelFitterSolverType::NewthonRhaphsonNoThrow,
	"NR_WITHDICHO",				ARM_ModelFitterSolverType::NewtonRaphsonWithDichotomy,
	"DICHOTOMY",				ARM_ModelFitterSolverType::Dichotomy,
	"NAG_SOLVER",				ARM_ModelFitterSolverType::NagSolver,
	"BRENT",					ARM_ModelFitterSolverType::Brent,
	"NONE",						ARM_ModelFitterSolverType::NoSolverType,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};
const ARM_ArgConv ARM_ArgConv_SolverTypeMethod( SolverTypeTable, "SolverType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_SolverTypeMethod( SolverTypeTable, "SolverType" );

ARGConvTable OptimiseTypeTable[] =
{	
    /// Type Name					/// number
    "bounds_no_deriv",			ARM_ModelFitterOptimizerType::bounds_no_deriv,
	"nlin_lsq",					ARM_ModelFitterOptimizerType::nlin_lsq,
    "lsq_deriv",				ARM_ModelFitterOptimizerType::lsq_deriv,
    "lsq_check_deriv",			ARM_ModelFitterOptimizerType::lsq_check_deriv,
    "OptimiseWithBrent",		ARM_ModelFitterOptimizerType::OptimiseWithBrent,
	"None",						ARM_ModelFitterOptimizerType::NoOptimizerType,

		/// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_OptimizerTypeMethod( OptimiseTypeTable, "OptimiseType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_OptimizerTypeMethod( OptimiseTypeTable, "OptimiseType" );

ARGConvTable BasketCalibrationTypeTable[] =
{	
    /// Type Name       /// number
    "BASKET",			ARM_BasketCalib::FULL,
	"BASKET_SIMPLE",	ARM_BasketCalib::SIMPLE,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_BasketCalibrationType( BasketCalibrationTypeTable, "BasketCalibrationType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_BasketCalibrationType( BasketCalibrationTypeTable, "BasketCalibrationType" );

ARGConvTable BasketCalibrationStrikeTable[] =
{	
    /// Strike Name       /// number
    "EQUIVALENT",		ARM_BasketCalib::EQUIVALENT,
	"ATM",				ARM_BasketCalib::ATM,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_BasketCalibrationStrike( BasketCalibrationStrikeTable, "BasketCalibrationStrike" );
const ARM_ArgConvReverse ARM_ArgConvReverse_BasketCalibrationStrike( BasketCalibrationStrikeTable, "BasketCalibrationStrike" );

ARGConvTable MoneyTypeTable[] =
{	
    /// Type Name       /// number
    "GAUSS",			ARM_SmileViewer::GAUSS,
	"BLACK",			ARM_SmileViewer::BLACK,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_MoneyType( MoneyTypeTable, "MoneyType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_MoneyType( MoneyTypeTable, "MoneyType" );

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

