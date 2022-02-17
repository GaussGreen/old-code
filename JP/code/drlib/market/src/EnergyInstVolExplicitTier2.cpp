//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyInstVolExplicitTier2.cpp
//
//   Description : Model volatility (not market vol) class for energy spreads
//
//   Author      : Lawrence Siu
//
//   Date        : May 01, 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include "edginc/EnergyInstVolExplicitTier2.hpp"


DRLIB_BEGIN_NAMESPACE

EnergyInstVolExplicitTier2::EnergyInstVolExplicitTier2():EnergyInstVolBase(TYPE)
{
}

EnergyInstVolExplicitTier2::~EnergyInstVolExplicitTier2()
{
}

void EnergyInstVolExplicitTier2::validatePop2Object()
{
    static const string method = "EnergyInstVolExplicitTier2::validatePop2Object()";

    EnergyInstVolBase::validatePop2Object();

    // do some checking...
    int numExpiries = fOptionExpiryDates.size();
    if ( numExpiries == 0 )
        ModelException(method,"No expiries for sigmas");

    fSigmas = DoubleMatrix(numExpiries, 1);
    for (int i=0; i<numExpiries; i++)
    {
        fSigmas[i][0] = sigmas[i];
    }
}


/*
int EnergyInstVolExplicitTier2::calibrate()
{
    static const string method = "EnergyInstVolExplicitTier2::calibrate()";

    int numCalibrated = deriveRatios();
    if ( numCalibrated > 1)
        fMaxPricingDate = fFutureMaturityDates[numCalibrated-1];

    return numCalibrated;
}


int EnergyInstVolExplicitTier2::deriveRatios()
{
    static const string method = "EnergyInstVolExplicitTier2::deriveRatios()";

    // Clear outputs
    fRatios.clear();

    // using local variables
    DateTimeArray fOptionExpiryDates = fFutureMaturityDates;
    DateTime fBaseDate = getBaseDate();
    double alpha = fAlphas[0];
    double beta = fAlphas[1];
    int numBenchmarks = fSigmas.numRows();

    // Check maturity dates so calibration is still possible even if the first future contract is expiring
    int isExpiryDate = (fOptionExpiryDates[0]==baseDate) ? 1 : 0;
    if (isExpiryDate)
    {
        fRatios.push_back(1.0);
    }
    
    // Loop through all points on the curve iteratively
    int i;
    double instRatio;
    double ratioFactor, maturityNowFactor, maturityNow;
    double ratioFactorS, maturityNowFactorS;
    double v1, v1Bar, v2Bar, v2;

    ratioFactor = exp(alpha*fOptionExpiryDates[isExpiryDate].daysDiff(fBaseDate)/365.0);
    ratioFactorS = exp(beta*fOptionExpiryDates[isExpiryDate].daysDiff(fBaseDate)/365.0);
    
    v1 = fSigmas[0][0];
    v1Bar = fSigmas[0][1];
    v2Bar = fSigmas[0][2];
    v2 = fSigmas[0][3];

    for (i=isExpiryDate; i<numBenchmarks; ++i)
    {
        maturityNow = fOptionExpiryDates[i].daysDiff(fBaseDate)/365.0;
        maturityNowFactor = exp(alpha*maturityNow);
        maturityNowFactorS = exp(beta*maturityNow);
    
        // inst vol ratio to first contract
        instRatio = getInstVolRatio(alpha,
                                    beta,
                                    v1Bar,
                                    v1,
                                    v2Bar,
                                    v2,
                                    maturityNowFactor,
                                    maturityNowFactorS,
                                    ratioFactor,
                                    ratioFactorS);
        
        fRatios.push_back(instRatio);
    }
    fMaxPricingDate = i==0 ? fBaseDate : fOptionExpiryDates[i-1] ;

    return i;
}*/

class EnergyInstVolExplicitTier2Helper
{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(EnergyInstVolExplicitTier2, clazz);
        SUPERCLASS(EnergyInstVolBase);
        EMPTY_SHELL_METHOD(defaultInstVolExplicitTier2);
        FIELD(sigmas, "Sigma array");
           
        Addin::registerConstructor("ENERGY_INST_VOL_EXPLICIT_TIER2",
                                   Addin::MARKET,
                                   "Creates a handle to an Energy Spred Instantaneous Vol Model",
                                   EnergyInstVolExplicitTier2::TYPE);
    }

    static IObject* defaultInstVolExplicitTier2()
    {
        return new EnergyInstVolExplicitTier2();
    }

};

CClassConstSP const EnergyInstVolExplicitTier2::TYPE =CClass::registerClassLoadMethod("EnergyInstVolExplicitTier2", 
                  typeid(EnergyInstVolExplicitTier2), EnergyInstVolExplicitTier2Helper::load);

// definition of TYPE for MarketWrapper template class
template <> CClassConstSP const EnergyInstVolExplicitTier2Wrapper::TYPE = 
       CClass::registerClassLoadMethod("EnergyInstVolExplicitTier2Wrapper", typeid(EnergyInstVolExplicitTier2Wrapper),load);


DRLIB_END_NAMESPACE
