//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyInstVolExplicitRegular.cpp
//
//   Description : Energy 2 factor vol surface. Based on drcommodityvolsurfacei.cpp
//                    and drcommodityInstaneousvolsurface.cpp in FXLIB.
//                 
//
//   Author      : Sean Chen
//
//   Date        : June 01, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include "edginc/EnergyInstVolExplicitRegular.hpp"


DRLIB_BEGIN_NAMESPACE

EnergyInstVolExplicitRegular::EnergyInstVolExplicitRegular() :
    EnergyInstVolBase(TYPE),
    x(-999.999),
    y(-999.999),
    z(-999.999)
{}

EnergyInstVolExplicitRegular::~EnergyInstVolExplicitRegular()
{}

void EnergyInstVolExplicitRegular::validatePop2Object()
{
    static const string method = "EnergyInstVolExplicitRegular::validatePop2Object()";

    EnergyInstVolBase::validatePop2Object();

    // do some checking...
    if (fOptionExpiryDates.size() == 0)
        ModelException(method, "Empty option expiry array!");

    if (sigma1s.size() == 0)
        ModelException(method, "Empty sigma1Bars array!");

    if (sigma1s.size() != fOptionExpiryDates.size())
        ModelException(method, "Number of sigma1Bars and option expiries mismatch!");

    if (fabs(rho) > 1.)
        ModelException(method, "Factor rho must lies within [-1,1]");

    if (w1 < 0.)
        ModelException(method, "Negative w1 is not allowed!");

    if (w2 < 0.)
        ModelException(method, "Negative w2 is not allowed!");

    if (!Maths::equals(w1 + w2, 1.))
        ModelException(method, "w1 and w2 do not sum up to 1!");

    // Convert inputs to QSRM form:
    convertInputs();

    fSigmas = DoubleMatrix(sigma1s.size(), 4);
    for (int i = 0; i < sigma1s.size(); ++i)
    {
        fSigmas[i][0] = sigma1s[i]; // sigma1s[i];
        fSigmas[i][1] = sigma1s[i]/x; // sigma1Bars[i];
        fSigmas[i][2] = sigma1s[i]/x*z; // sigma2Bars[i];
        fSigmas[i][3] = sigma1s[i]/x*y; // sigma2s[i];
    }

    // set future maturities = option expires
    fFutureMaturityDates = fOptionExpiryDates;    
}

void EnergyInstVolExplicitRegular::convertInputs()
{
    x = computeX();
    y = 1.0;
    z = computeZ();
}

double EnergyInstVolExplicitRegular::computeZ()
{
    return sqrt(1. - rho*rho)/rho;
}

double EnergyInstVolExplicitRegular::computeX()
{
    return w1/w2/rho;
}

/*
int EnergyInstVolExplicitRegular::calibrate()
{
    static const string method = "EnergyInstVolExplicitRegular::calibrate()";

    int numCalibrated = deriveRatios();
    if ( numCalibrated > 1)
        fMaxPricingDate = fFutureMaturityDates[numCalibrated-1];

    return numCalibrated;
}


int EnergyInstVolExplicitRegular::deriveRatios()
{
    static const string method = "EnergyInstVolExplicitRegular::deriveRatios()";

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
}
*/

class EnergyInstVolExplicitRegularHelper
{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(EnergyInstVolExplicitRegular, clazz);
        SUPERCLASS(EnergyInstVolBase);
        EMPTY_SHELL_METHOD(defaultInstVolExplicitRegular);
        FIELD(sigma1s, "Sigma1 array");
        //FIELD(sigma2s, "Sigma2 Array");
        //FIELD(sigma1Bars, "sigam1Bar array");
        //FIELD(sigma2Bars, "sigma2bar array");
        FIELD(w1, "Sampras Oil model const: weight for random factor1");
        FIELD(w2, "Sampras Oil model const: weight for random factor2");
        FIELD(rho, "Corr between factor 1 and 2 in the Sampras Oil Model");
        FIELD(x, "parameter x");
        FIELD_MAKE_TRANSIENT(x);
        FIELD(z, "parameter z");
        FIELD_MAKE_TRANSIENT(z);

        Addin::registerConstructor("ENERGY_INST_VOL_EXPLICIT_REGULAR",
                                   Addin::MARKET,
                                   "Creates a handle to an Energy Instantaneous Vol Model",
                                   EnergyInstVolExplicitRegular::TYPE);
    }

    static IObject* defaultInstVolExplicitRegular()
    {
        return new EnergyInstVolExplicitRegular();
    }

};

CClassConstSP const EnergyInstVolExplicitRegular::TYPE =CClass::registerClassLoadMethod("EnergyInstVolExplicitRegular", 
                  typeid(EnergyInstVolExplicitRegular), EnergyInstVolExplicitRegularHelper::load);

// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(EnergyInstVolExplicitRegularWrapper);

DRLIB_END_NAMESPACE
