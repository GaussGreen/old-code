
#include "gpinfra/gramfunctordef.h"

/// gpbase
#include "gpbase/utilityport.h"
#include "gpbase/gpmatrix.h"

/// gpinfra
#include "gpinfra/pricingmodelir.h"

/// gpnumlib
#include "gpnumlib/gaussiananalytics.h"

/// gpclosedforms
#include "gpclosedforms/bondanalytic.h"

/// STL
#include <cmath>


CC_BEGIN_NAMESPACE( ARM )

/// exponential 
ARM_GP_UnaryOp ExpVecFctor( pDbleUnaryFunc( exp ), "Exp" );

/// logarithm
ARM_GP_UnaryOp LogVecFctor( pDbleUnaryFunc( log ), "Log" );

/// Square Root
ARM_GP_UnaryOp SqrtVecFctor( pDbleUnaryFunc( sqrt), "Sqrt" );

/// Abs function
ARM_GP_UnaryOp AbsVecFctor( pDbleUnaryFunc( fabs), "Abs" );

double cumulativeNormalFunction( double d )
{	return ARM_GaussianAnalytics::cdfNormal(d); }

/// Cumulative normal function
ARM_GP_UnaryOp CumNormFctor( pDbleUnaryFunc( cumulativeNormalFunction ), "CumNorm" );

/// Power operator
ARM_GP_PowVector PowFctor = ARM_GP_PowVector();

ARM_GP_SumSerieVector SumSerieFctor = ARM_GP_SumSerieVector();

// Statistics Average operator
ARM_GP_StatAverageVector StatAverageFctor = ARM_GP_StatAverageVector();

// Statistics Standard Deviation operator
ARM_GP_StatStdDevVector StatStdDevFctor = ARM_GP_StatStdDevVector();

/// If operator
ARM_GP_IfVector IfFctor = ARM_GP_IfVector();

/// Minimum operator
//FIXMEFRED: mig.vc8 : CC_Min -> CC_Min_c
ARM_GP_BinaryOpWithDates MinFctor( pDbleBinaryFunc( CC_Min_c<double> ), "Min" );

/// Maximum operator
//FIXMEFRED: mig.vc8 : CC_Max -> CC_Max_c
ARM_GP_BinaryOpWithDates MaxFctor( pDbleBinaryFunc( CC_Max_c<double> ), "Max" );

/// Day Count Fraction operator
ARM_GP_DCF DCFFctor = ARM_GP_DCF();

/// Discount Factor
ARM_GP_DF DFFctor = ARM_GP_DF();

/// Forward Rate
ARM_GP_Libor LiborFctor = ARM_GP_Libor();

/// Annuity key word
ARM_GP_Annuity AnnuityFctor = ARM_GP_Annuity();

/// Swap Rate
ARM_GP_SwapRate SwapRateFctor = ARM_GP_SwapRate();

/// Max Rate
ARM_GP_MaxRate MaxRateFctor = ARM_GP_MaxRate();

/// Implied Vol
ARM_GP_ImpliedVol ImpliedVolFctor = ARM_GP_ImpliedVol();

/// Spread
ARM_GP_Spread SpreadFctor = ARM_GP_Spread();

/// Swap
ARM_GP_Swap SwapFctor = ARM_GP_Swap();

/// Basis Swap 
ARM_GP_BasisSwap BasisSwapFctor = ARM_GP_BasisSwap();

/// CPI
ARM_GP_CPI CPIFctor = ARM_GP_CPI();

/// Inflation Swap Rate
ARM_GP_InfSwapRate InfSwapRateFctor = ARM_GP_InfSwapRate();

/// Inflation Swap Rate
ARM_GP_InfSwap InfSwapFctor = ARM_GP_InfSwap();

/// PV operator
ARM_GP_PV PVFctor = ARM_GP_PV() ;

/// UnPay operator
ARM_GP_UnPay UnPayFctor = ARM_GP_UnPay() ;

/// Exercise operator
ARM_GP_Exercise ExerciseFctor = ARM_GP_Exercise() ;

/// Frontier operator
ARM_GP_Frontier FrontierFctor = ARM_GP_Frontier() ;

/// Trigger operator
ARM_GP_Trigger TriggerFctor = ARM_GP_Trigger() ;

/// ModelFactor operator
ARM_GP_MFactorFctor MFactorFctor = ARM_GP_MFactorFctor();

/// Caplet Functor
ARM_GP_CapDigital_Common CapletFctor  = ARM_GP_CapDigital_Common( &ARM_PricingFunctionIR::VanillaCaplet,   "Caplet");

/// Digital Functor
ARM_GP_CapDigital_Common DigitalFctor = ARM_GP_CapDigital_Common( &ARM_PricingFunctionIR::VanillaDigital,  "Digital");

/// Swaption Functor
ARM_GP_Swaption SwaptionFctor = ARM_GP_Swaption();

// Sum Option
ARM_GP_SumOption SumOptionFctor = ARM_GP_SumOption();

/// Cap Functor
ARM_GP_Cap CapFctor = ARM_GP_Cap();

/// Spread Option Functor
ARM_GP_SpreadOption SpreadOptionFctor = ARM_GP_SpreadOption();

/// Corridor Functor
ARM_GP_Corridor CorridorFctor = ARM_GP_Corridor();

/// Double Digital
ARM_GP_DoubleDigital DoubleDigitalFctor = ARM_GP_DoubleDigital();

/// Price To Yield Functor
ARM_GP_PTYAndYTPFctor PTYFctor = ARM_GP_PTYAndYTPFctor( &ARM_BondAnalytics::PriceToYield, "PriceToYield" );

/// Yield To Price Functor
ARM_GP_PTYAndYTPFctor YTPFctor = ARM_GP_PTYAndYTPFctor( &ARM_BondAnalytics::YieldToPrice, "YieldToPrice" );



/// equity functions
ARM_GP_SpotFctor SpotFctor	= ARM_GP_SpotFctor();

ARM_GP_FwdFctor FwdFctor	= ARM_GP_FwdFctor();

ARM_GP_CallFctor CallFctor	= ARM_GP_CallFctor();

ARM_GP_GreekFctor GreekFctor= ARM_GP_GreekFctor();

ARM_GP_CallSpreadFctor CallSpreadFctor	= ARM_GP_CallSpreadFctor();

ARM_GP_EqDigitalFctor EqDigitalFctor = ARM_GP_EqDigitalFctor();

ARM_GP_CallStripFctor CallStripFctor = ARM_GP_CallStripFctor();

ARM_GP_RangeAccrualFctor RangeAccrualFctor = ARM_GP_RangeAccrualFctor();

/// closed form functors
ARM_CFCallFctor CFCallFctor = ARM_CFCallFctor ();

/// closed form greek functor
ARM_CFGreekFctor CFGreekFctor = ARM_CFGreekFctor ();

/// Option Mepi on Funds functors
ARM_CF_MepiCallFctor CF_MepiCallFctor = ARM_CF_MepiCallFctor ();

/// Option Mepi on Funds  greek functor
ARM_CF_MepiGreekFctor CF_MepiGreekFctor = ARM_CF_MepiGreekFctor ();

/// call on mepi funds greek functor
ARM_NumericalCallOnMepiDeltaFctor NumericalCallOnMepiDeltaFctor = ARM_NumericalCallOnMepiDeltaFctor();

/// Option  on Funds functors
ARM_CF_StochVolCallFctor CF_StochVolCallFctor = ARM_CF_StochVolCallFctor ();

/// Option  on Funds  greek functor
ARM_CF_StochVolGreekFctor CF_StochVolGreekFctor = ARM_CF_StochVolGreekFctor ();



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

