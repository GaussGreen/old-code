/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gramfunctordef.h 
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */



#ifndef _INGPINFRA_GRAMFUNCTORDEF_H
#define _INGPINFRA_GRAMFUNCTORDEF_H

#include "gpbase/port.h"
#include "gramfunctor.h"


CC_BEGIN_NAMESPACE( ARM )


/// maths function

/// exponential 
extern ARM_GP_UnaryOp ExpVecFctor;

/// logarithm
extern ARM_GP_UnaryOp LogVecFctor;

/// Square Root
extern ARM_GP_UnaryOp SqrtVecFctor;

/// Abs functor
extern ARM_GP_UnaryOp AbsVecFctor;

/// Cumulative Normal functor
extern ARM_GP_UnaryOp CumNormFctor;

/// closed form call functor
extern ARM_CFCallFctor CFCallFctor;

/// closed form greek functor
extern ARM_CFGreekFctor	CFGreekFctor;

/// option on mepi funds call functor
extern ARM_CF_MepiCallFctor CF_MepiCallFctor;

/// option on mepi funds greek functor
extern ARM_CF_MepiGreekFctor CF_MepiGreekFctor;

/// call on mepi funds greek functor
extern ARM_NumericalCallOnMepiDeltaFctor NumericalCallOnMepiDeltaFctor;

/// option on  funds call functor
extern ARM_CF_StochVolCallFctor CF_StochVolCallFctor;

/// option on  funds greek functor
extern ARM_CF_StochVolGreekFctor CF_StochVolGreekFctor;

/// Pow functor
extern ARM_GP_PowVector PowFctor;

/// Minimum operator
extern ARM_GP_BinaryOpWithDates MinFctor;

/// Maximum operator
extern ARM_GP_BinaryOpWithDates MaxFctor;

/// If operator
extern ARM_GP_IfVector IfFctor;

extern ARM_GP_SumSerieVector SumSerieFctor;

// Statistics Average operator
extern ARM_GP_StatAverageVector StatAverageFctor;

// Statistics Standard Deviation operator
extern ARM_GP_StatStdDevVector StatStdDevFctor;

/// DCF operator
extern ARM_GP_DCF DCFFctor;

/// General model functions
extern ARM_GP_MFactorFctor MFactorFctor;

/// interest rate functions

/// DF (discount factor) operator
extern ARM_GP_DF DFFctor;

/// Libor operator
extern ARM_GP_Libor LiborFctor;

/// Annuity operator
extern ARM_GP_Annuity AnnuityFctor;

/// Swap Rate operator
extern ARM_GP_SwapRate SwapRateFctor;

/// Spread operator
extern ARM_GP_Spread SpreadFctor;

/// Swap Functor
extern ARM_GP_Swap SwapFctor;

/// MaxRate Functor
extern ARM_GP_MaxRate MaxRateFctor;

/// Implied Vol Functor
extern ARM_GP_ImpliedVol ImpliedVolFctor;

/// BasisSwap Functor
extern ARM_GP_BasisSwap BasisSwapFctor;

/// Cpi Functor
extern ARM_GP_CPI CPIFctor;

/// Inflation Swap Rate Functor
extern ARM_GP_InfSwapRate InfSwapRateFctor;

/// Inflation Swap Functor
extern ARM_GP_InfSwap InfSwapFctor;

/// PV Functor
extern ARM_GP_PV PVFctor;

/// Unpay Functor
extern ARM_GP_UnPay UnPayFctor;

// Exercise Functor
extern ARM_GP_Exercise ExerciseFctor;

// Frontier Functor
extern ARM_GP_Frontier FrontierFctor;

// Trigger Functor
extern ARM_GP_Trigger TriggerFctor;

/// Caplet Functor
extern ARM_GP_CapDigital_Common CapletFctor;

/// Digital Functor
extern ARM_GP_CapDigital_Common DigitalFctor;

/// Swaption Functor
extern ARM_GP_Swaption SwaptionFctor;

/// SumOption Functor
extern ARM_GP_SumOption SumOptionFctor;

/// Cap Functor
extern ARM_GP_Cap CapFctor;

/// Spread Option Functor
extern ARM_GP_SpreadOption SpreadOptionFctor;

/// Corridor Functor
extern ARM_GP_Corridor CorridorFctor;

/// Double Digital Functor
extern ARM_GP_DoubleDigital DoubleDigitalFctor;

/// Price To Yield Functor
extern ARM_GP_PTYAndYTPFctor PTYFctor;

/// Yield To Price Functor
extern ARM_GP_PTYAndYTPFctor YTPFctor;


/// equity functions
extern ARM_GP_SpotFctor SpotFctor;

extern ARM_GP_FwdFctor FwdFctor;

extern ARM_GP_CallFctor CallFctor;

extern ARM_GP_GreekFctor GreekFctor;

extern ARM_GP_CallSpreadFctor CallSpreadFctor;

extern ARM_GP_EqDigitalFctor EqDigitalFctor;

extern ARM_GP_CallStripFctor CallStripFctor;

extern ARM_GP_RangeAccrualFctor RangeAccrualFctor;

CC_END_NAMESPACE()

#endif


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

