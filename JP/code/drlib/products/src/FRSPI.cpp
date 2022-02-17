//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRSPI.cpp
//
//   Description : Maths operations on rvalues
//
//   Date        : Aug 02
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/FRFunction.hpp"

DRLIB_BEGIN_NAMESPACE
#define FR_GET_VALUE(x) ((x)->func(x))
#define FR_GET_SIZE(x) ((x)->size(x))
/**  */
class FRSPI {

    static double SPILagAdj(FRFunction::Val* args) {
        enum {
            eAssetPrice = 0,
            eAssetAlloc,
            eRebalanceFlag,
            eLagSteps,
            eIndex
        };
        int index = args[eIndex].i;
        int lagSteps = args[eLagSteps].i;
        double lagAdj;
        
        if (index < lagSteps) {
            lagAdj = 0.;
        } else if (index == lagSteps) {
            double Ei_g = FR_GET_VALUE((*args[eAssetPrice].rValues)[0].dExp); // i-g == 0 
            double Ei = FR_GET_VALUE((*args[eAssetPrice].rValues)[index].dExp);
            double nEi_g = 
                FR_GET_VALUE((*args[eAssetAlloc].rValues)[0].dExp);// i-g == 0 
            lagAdj = nEi_g * (Ei - Ei_g);
        } else {
            bool Ri_g = FR_GET_VALUE((*args[eRebalanceFlag].
                                      rValues)[index-lagSteps].bExp);
            if (Ri_g) {
                double Ei_g = FR_GET_VALUE((*args[eAssetPrice].
                                            rValues)[index-lagSteps].dExp);
                double Ei = FR_GET_VALUE((*args[eAssetPrice].
                                          rValues)[index].dExp);
                double nEi_g = FR_GET_VALUE((*args[eAssetAlloc].
                                             rValues)[index-lagSteps].dExp);
                double nEi_g_1 =
                    FR_GET_VALUE((*args[eAssetAlloc].
                                  rValues)[index-lagSteps-1].dExp);
                lagAdj = (nEi_g - nEi_g_1) * (Ei - Ei_g);
            } else {
                lagAdj = 0.;
            }
        }
        return lagAdj;
    }

    static bool SPIReBalance(FRFunction::Val* args) {
        enum {
            eSustainableLoss = 0,
            eUnbalancedLoss,
            eCutoff,
            eOverExposureThreshold,
            eUnderExposureThreshold,
            eMaxLoss,
            eMinLoss,
            eIndex
        };
        int index = args[eIndex].i;
        bool reBal = false;
        bool cutoffi = FR_GET_VALUE((*args[eCutoff].rValues)[index].bExp);

        if (cutoffi) {
            if (index>0) {
                bool cutoffi_1 =
                    FR_GET_VALUE((*args[eCutoff].rValues)[index-1].bExp);
                reBal = !cutoffi_1;
            } else {
                reBal = false;
            }
        } else {
            double SL = FR_GET_VALUE(args[eSustainableLoss].dExp);
            double UL = FR_GET_VALUE(args[eUnbalancedLoss].dExp);
            double OT = args[eOverExposureThreshold].d;
            double UT = args[eUnderExposureThreshold].d;
            double MaxLoss = args[eMaxLoss].d;
            double MinLoss = args[eMinLoss].d;
            reBal = 
                ( (UL > Maths::min(SL, MaxLoss) * (1. + OT)) &&
                  (UL > MinLoss) ) ||
                ( (UL < SL * (1. - UT)) && (UL < MaxLoss) ) ||
                ( Maths::isZero(SL) && (UL > MinLoss) );
        }
        return reBal;
    }

    static double SPIUnbalancedExposure(FRFunction::Val* args) {
        enum {
            eAssetPriceArray = 0,
            eAssetAllocPrevArray, 
            eBondPrice, 
            eBondAllocPrev, 
            eLoanBalancePrev, 
            eBasket
        };
        FRIfaces::IRValueDoubleArray::RT* rtEi = 
            args[eAssetPriceArray].dExpArray; // for ease
        FRIfaces::IRValueDoubleArray::RT* rtnEi_1 = 
            args[eAssetAllocPrevArray].dExpArray; // for ease
        int numAssets = FR_GET_SIZE(rtEi);
        if (numAssets != FR_GET_SIZE(rtnEi_1)){
            throw ModelException("FRSPI::SPIUnbalancedExposure",
                                 "Length of arrays must be equal");
        }
        double Zi = FR_GET_VALUE(args[eBondPrice].dExp);
        double nZi_1 = FR_GET_VALUE(args[eBondAllocPrev].dExp);
        double LBi_1 = FR_GET_VALUE(args[eLoanBalancePrev].dExp);
        double B = FR_GET_VALUE(args[eBasket].dExp);
        double sum = 0.0;

        FRIfaces::IRValueDoubleArray::TGetValue* iGetValue = rtEi->func;
        FRIfaces::IRValueDoubleArray::TGetValue* ni_1GetValue = rtnEi_1->func;
        for(int iAsset = 0; iAsset < numAssets; iAsset++) {
            sum += ni_1GetValue(rtnEi_1, iAsset) * iGetValue(rtEi, iAsset);
        }
        double unbalExp = sum / ( nZi_1 * Zi + sum );
        unbalExp *= 1. + LBi_1 / B;
        
        return unbalExp;
    }

    /** faster version of SPIUnbalancedExposure for 1D - avoid array cost */
    static double SPIUnbalancedExposure1D(FRFunction::Val* args) {
        enum {
            eAssetPrice = 0,
            eAssetAllocPrev, 
            eBondPrice, 
            eBondAllocPrev, 
            eLoanBalancePrev, 
            eBasket
        };

        double E = FR_GET_VALUE(args[eAssetPrice].dExp);
        double nE = FR_GET_VALUE(args[eAssetAllocPrev].dExp);
        double Zi = FR_GET_VALUE(args[eBondPrice].dExp);
        double nZi_1 = FR_GET_VALUE(args[eBondAllocPrev].dExp);
        double LBi_1 = FR_GET_VALUE(args[eLoanBalancePrev].dExp);
        double B = FR_GET_VALUE(args[eBasket].dExp);
        double sum = nE * E;
        double unbalExp = sum / ( nZi_1 * Zi + sum );
        unbalExp *= 1. + LBi_1 / B;
        
        return unbalExp;
    }

    static double SPIUnbalancedCrash(FRFunction::Val* args) {
        enum {
            eAssetPriceArray = 0,
            eAssetAllocPrevArray, 
            eCrashSizeArray
        };
        // for ease ...
        FRIfaces::IRValueDoubleArray::RT* rtEi = 
            args[eAssetPriceArray].dExpArray;
        FRIfaces::IRValueDoubleArray::RT* rtnEi_1 = 
            args[eAssetAllocPrevArray].dExpArray;
        FRIfaces::IRValueDoubleArray::RT* rtC = 
            args[eCrashSizeArray].dExpArray;
        int numAssets = FR_GET_SIZE(rtEi);
        if (numAssets != FR_GET_SIZE(rtnEi_1) || numAssets !=FR_GET_SIZE(rtC)){
            throw ModelException("SPIUnbalancedExposure1D", "Arrays must be "
                                 "of equal length");
        }
        FRIfaces::IRValueDoubleArray::TGetValue* EiGetValue = rtEi->func;
        FRIfaces::IRValueDoubleArray::TGetValue* nEi_1GetValue = rtnEi_1->func;
        FRIfaces::IRValueDoubleArray::TGetValue* cGetValue = rtC->func;
        double             num = 0.0;
        double             den = 0.0;
        for(int iAsset = 0; iAsset < numAssets; iAsset++) {
            double nEi_1 = nEi_1GetValue(rtnEi_1, iAsset);
            double Ei = EiGetValue(rtEi, iAsset);
            double product = Ei * nEi_1;
            den += product;
            num += product * cGetValue(rtC, iAsset);
        }
        return (num/den);
    }

// not yet - needs some work on infrastructure
#define LATER 0
#if LATER
    static double SPIFeeAmount(FRFunction::Val* args) {
        return 0.;
    }

    static double SPILoanCharge(FRFunction::Val* args) {
        return 0.;
    }
#endif

    static double SPIBasketValue(FRFunction::Val* args) {
        enum {
            eAssetPriceArray = 0,
            eAssetAllocPrevArray, 
            eBondPrice,
            eBondAllocPrev,
            eDeductionsArray
        };
        // for ease ...
        FRIfaces::IRValueDoubleArray::RT* rtEi = 
            args[eAssetPriceArray].dExpArray;
        FRIfaces::IRValueDoubleArray::RT* rtnEi_1 = 
            args[eAssetAllocPrevArray].dExpArray;
        FRIfaces::IRValueDoubleArray::RT* rtX = 
            args[eDeductionsArray].dExpArray;
        int numAssets = FR_GET_SIZE(rtEi);
        if (numAssets != FR_GET_SIZE(rtnEi_1)){
            throw ModelException("FRSPI::SPIBasketValue",
                                 "Length of arrays must be equal");
        }
        FRIfaces::IRValueDoubleArray::TGetValue* iGetValue = rtEi->func;
        FRIfaces::IRValueDoubleArray::TGetValue* ni_1GetValue = rtnEi_1->func;
        FRIfaces::IRValueDoubleArray::TGetValue* XGetValue = rtX->func;
        double             Zi = FR_GET_VALUE(args[eBondPrice].dExp);
        double             nZi_1 = FR_GET_VALUE(args[eBondAllocPrev].dExp);
        int                i;
        double             B = 0.;

        for(i=0; i<numAssets; i++) {
            B += ni_1GetValue(rtnEi_1, i) * iGetValue(rtEi, i);
        }
        int XSize = FR_GET_SIZE(rtX);
        for(i = 0; i < XSize; i++) {
            B -= XGetValue(rtX, i);
        }
        B += nZi_1 * Zi;

        return B;
    }

public:
    static void registerFuncs(){
        { // SPILagAdj
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleType,   
                FRIfaces::doubleType,   
                FRIfaces::boolType,        
                FRIfaces::intType};
            static const char* paramNames[] = {
                "AssetPrice", 
                "AssetAlloc", 
                "RebalanceFlag", 
                "LagSteps"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::arrayExpression, 
                FRFunction::arrayExpression, 
                FRFunction::arrayExpression, 
                FRFunction::native}; // typically a const int

            FRFunction::registerGenericDoubleFunction("SPI_LAG_ADJ",
                                                      "Computes lag adjustment",
                                                      sizeof(types)/sizeof(FRIfaces::VarType),
                                                      types,
                                                      paramNames,
                                                      calcTypes,
                                                      SPILagAdj);
        }

        { // SPIReBalance
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleType,   
                FRIfaces::doubleType,   
                FRIfaces::boolType,        
                FRIfaces::doubleType,   
                FRIfaces::doubleType,   
                FRIfaces::doubleType,   
                FRIfaces::doubleType};
            static const char* paramNames[] = {
                "SustainableLoss", 
                "UnbalancedLoss", 
                "Cutoff", 
                "OverExposureThreshold", 
                "UnderExposureThreshold", 
                "MaxLoss", 
                "MinLoss"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::expression,
                FRFunction::expression,
                FRFunction::arrayExpression, 
                FRFunction::native,
                FRFunction::native,
                FRFunction::native,
                FRFunction::native};

            FRFunction::registerGenericBoolFunction("SPI_REBAL",
                                                    "Computes Rebalance flag",
                                                    sizeof(types)/sizeof(FRIfaces::VarType),
                                                    types,
                                                    paramNames,
                                                    calcTypes,
                                                    SPIReBalance);
        }
        { // SPIUnbalancedExposure
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleArrayType,   
                FRIfaces::doubleArrayType,   
                FRIfaces::doubleType,   
                FRIfaces::doubleType,   
                FRIfaces::doubleType,        
                FRIfaces::doubleType};
            static const char* paramNames[] = {
                "AssetPriceArray",  
                "AssetAllocPrevArray", 
                "BondPrice", 
                "BondAllocPrev", 
                "LoanBalancePrev", 
                "Basket"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::expression,
                FRFunction::expression,
                FRFunction::expression,
                FRFunction::expression,
                FRFunction::expression,
                FRFunction::expression};

            FRFunction::registerGenericDoubleFunction(
                "SPI_UNBAL_EXPO",
                "Computes unbalanced exposure",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes,
                SPIUnbalancedExposure);
        }
        { // SPIUnbalancedExposure1D
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleType,   
                FRIfaces::doubleType,   
                FRIfaces::doubleType,   
                FRIfaces::doubleType,   
                FRIfaces::doubleType,        
                FRIfaces::doubleType};
            static const char* paramNames[] = {
                "AssetPrice",  
                "AssetAllocPrev", 
                "BondPrice", 
                "BondAllocPrev", 
                "LoanBalancePrev", 
                "Basket"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::expression,
                FRFunction::expression,
                FRFunction::expression,
                FRFunction::expression,
                FRFunction::expression,
                FRFunction::expression};

            FRFunction::registerGenericDoubleFunction(
                "SPI_UNBAL_EXPO_1D",
                "Computes unbalanced exposure",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes,
                SPIUnbalancedExposure1D);
        }
        { // SPIUnbalancedCrash
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleArrayType,   
                FRIfaces::doubleArrayType,   
                FRIfaces::doubleArrayType};
            static const char* paramNames[] = {
                "AssetPriceArray", 
                "AssetAllocPrevArray", 
                "CrashSizeArray"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::expression,
                FRFunction::expression,
                FRFunction::expression};

            FRFunction::registerGenericDoubleFunction("SPI_UNBAL_CRASH",
                                                      "Computes unbalanced crash",
                                                      sizeof(types)/sizeof(FRIfaces::VarType),
                                                      types,
                                                      paramNames,
                                                      calcTypes,
                                                      SPIUnbalancedCrash);
        }
#if LATER
        { // SPIFeeAmount
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleArrayType,   
                FRIfaces::doubleArrayType,   
                FRIfaces::doubleArrayType,   
                FRIfaces::doubleType,   
                FRIfaces::doubleType,   
                FRIfaces::doubleType,   
                FRIfaces::doubleType,   
                FRIfaces::boolType,        
                FRIfaces::boolType};
            static const char* paramNames[] = {
                "AssetPriceArray", 
                "AssetAllocArray", 
                "AssetFeeRateArray",
                "BondPrice", 
                "BondAlloc", 
                "BondFeeRate",
                "DayCountDenominator",
                "FeePaymentFlag", 
                "calcAtEachStep"};

            FRFunction::registerGenericDoubleFunction("SPI_FEE_AMOUNT",
                                                      "",
                                                      sizeof(types)/sizeof(FRIfaces::VarType),
                                                      types,
                                                      params,
                                                      0,
                                                      SPIFeeAmount);
        }

        { // SPILoanCharge
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleType,   
                FRIfaces::doubleType,   
                FRIfaces::doubleType,   
                FRIfaces::boolType,   
                FRIfaces::boolType,
                XXX};
            static const char* paramNames[] = {
                "LoanBalance", 
                "DayCountDenominator", 
                "LoanSpread",
                "DeductionFlag", 
                "PotentialRebalanceFlag", 
                "Needs something for the rate calc?"};

            FRFunction::registerGenericDoubleFunction("SPI_LOAN_CHARGE",
                                                      "",
                                                      sizeof(types)/sizeof(FRIfaces::VarType),
                                                      types,
                                                      params,
                                                      0,
                                                      SPILoanCharge);
        }
#endif

        { // SPIBasketValue
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleArrayType,   
                FRIfaces::doubleArrayType,   
                FRIfaces::doubleType,   
                FRIfaces::doubleType,   
                FRIfaces::doubleArrayType};
            static const char* paramNames[] = {
                "AssetPriceArray", 
                "AssetAllocPrevArray", 
                "BondPrice", 
                "BondAllocPrev", 
                "DeductionsArray"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::expression,
                FRFunction::expression,
                FRFunction::expression,
                FRFunction::expression,
                FRFunction::expression};

            FRFunction::registerGenericDoubleFunction("SPI_BASKET_VALUE",
                                                      "Computes basket value",
                                                      sizeof(types)/sizeof(FRIfaces::VarType),
                                                      types,
                                                      paramNames,
                                                      calcTypes,
                                                      SPIBasketValue);
        }
    }
};

bool loadFRSPI() {
    FRSPI::registerFuncs();
    return true;
}
DRLIB_END_NAMESPACE


