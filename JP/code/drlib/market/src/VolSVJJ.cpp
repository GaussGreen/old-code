//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolSVJJ.cpp
//
//   Description : Heston + StockOnlyJump (Merton) + VolOnlyJump + CommonStockVolJump
//
//   Date        : 20 Nov 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolSVJJ.hpp"
#include "edginc/Maths.hpp"
#include "edginc/VolSV.hpp"
#include "edginc/VolSVJ.hpp"
#include "edginc/VolSVCJ.hpp"

DRLIB_BEGIN_NAMESPACE


VolSVJJ::VolSVJJ(const string& name,
                 double initialVol,
                 double meanVol,
                 double meanReversRate,
                 double volVol,
                 double correlation,
                 double stockCrashRate,
                 double stockCrashSizeMean,
                 double stockCrashSizeUncertainty,
                 double volCrashRate,
                 double volCrashSizeMean,
                 double commonCrashRate,
                 double commonStockCrashSizeMean,
                 double commonStockCrashSizeUncertainty,
                 double commonVolCrashSizeMean,
                 double stockVolCrashSizeCorrelation):
VolBaseParam(TYPE, name), initialVol(initialVol), meanVol(meanVol), meanReversRate(meanReversRate),
volVol(volVol), correlation(correlation), stockCrashRate(stockCrashRate), stockCrashSizeMean(stockCrashSizeMean),
stockCrashSizeUncertainty(stockCrashSizeUncertainty), volCrashRate(volCrashRate), volCrashSizeMean(volCrashSizeMean),
commonCrashRate(commonCrashRate), commonStockCrashSizeMean(commonStockCrashSizeMean), 
commonStockCrashSizeUncertainty(commonStockCrashSizeUncertainty), commonVolCrashSizeMean(commonVolCrashSizeMean), 
stockVolCrashSizeCorrelation(stockVolCrashSizeCorrelation) {
    validatePop2Object();
}


VolSVSP VolSVJJ::convert(VolSV* p) const {
    VolSVSP volSV;
    
    bool noStockJumps = 
        Maths::isZero(stockCrashRate) || 
        ( Maths::isZero(stockCrashSizeMean) && Maths::isZero(stockCrashSizeUncertainty) );

    bool noVolJumps = 
        Maths::isZero(volCrashRate) || Maths::isZero(volCrashSizeMean);

    bool noCommonJumps =
        Maths::isZero(commonCrashRate) || 
        ( Maths::isZero(commonStockCrashSizeMean) && Maths::isZero(commonStockCrashSizeUncertainty) && Maths::isZero(commonVolCrashSizeMean));

    
    if(noStockJumps && noVolJumps && noCommonJumps) {
        // If no jumps at all then convert to VolSV
        volSV = VolSVSP(new VolSV(getName(), initialVol, meanVol, meanReversRate, volVol, correlation));
    }
    
    return volSV;
}


VolSVJSP VolSVJJ::convert(VolSVJ* p) const {
    VolSVJSP volSVJ;

    bool noStockJumps = 
        Maths::isZero(stockCrashRate) || 
        ( Maths::isZero(stockCrashSizeMean) && Maths::isZero(stockCrashSizeUncertainty) );

    bool noVolJumps = 
        Maths::isZero(volCrashRate) || Maths::isZero(volCrashSizeMean);

    bool noCommonJumps =
        Maths::isZero(commonCrashRate) || 
        ( Maths::isZero(commonStockCrashSizeMean) && Maths::isZero(commonStockCrashSizeUncertainty) && Maths::isZero(commonVolCrashSizeMean));

    // Cannot support vol jumps in SVJ
    if(noVolJumps) {
        if(noCommonJumps) {
            // Either stock jumps only or no jumps at all
            volSVJ = VolSVJSP(new VolSVJ(getName(), initialVol, meanVol, meanReversRate, volVol, correlation, 
                stockCrashRate, stockCrashSizeMean, stockCrashSizeUncertainty, 0.0));
        } else {
            if(noStockJumps && Maths::isZero(commonVolCrashSizeMean)) {
                // Common jumps only and zero common vol jumps
                volSVJ = VolSVJSP(new VolSVJ(getName(), initialVol, meanVol, meanReversRate, volVol, correlation, 
                    commonCrashRate, commonStockCrashSizeMean, commonStockCrashSizeUncertainty, 0.0));
            }
        }
    }
    
    return volSVJ;
}


VolSVCJSP VolSVJJ::convert(VolSVCJ* p) const {
    VolSVCJSP volSVCJ;
    
    bool noStockJumps = 
        Maths::isZero(stockCrashRate) || 
        ( Maths::isZero(stockCrashSizeMean) && Maths::isZero(stockCrashSizeUncertainty) );
    
    bool noVolJumps = 
        Maths::isZero(volCrashRate) || Maths::isZero(volCrashSizeMean);

    bool noCommonJumps =
        Maths::isZero(commonCrashRate) || 
        ( Maths::isZero(commonStockCrashSizeMean) && Maths::isZero(commonStockCrashSizeUncertainty) && Maths::isZero(commonVolCrashSizeMean));

    if(noStockJumps && noVolJumps) {
        // If no idiosyncratic jumps then convert to VolSVCJ
        volSVCJ = VolSVCJSP(new VolSVCJ(
            getName(), initialVol, correlation, volVol, meanVol, meanReversRate,
            commonCrashRate, commonStockCrashSizeMean, commonStockCrashSizeUncertainty, commonVolCrashSizeMean, 
            stockVolCrashSizeCorrelation, false));
    } else if(noCommonJumps) {
        // Support only spot or only vol jumps
        if(noVolJumps) {
            // Only Spot jumps or not even that so can be modelled as SVCJ
            volSVCJ = VolSVCJSP(new VolSVCJ(
                getName(), initialVol, correlation, volVol, meanVol, meanReversRate,
                stockCrashRate, stockCrashSizeMean, stockCrashSizeUncertainty, 0.0, 
                0.0, false));
        } else if(noStockJumps) {
            // Only vol jumps or not even that so can be modelled as SVCJ
            volSVCJ = VolSVCJSP(new VolSVCJ(
                getName(), initialVol, correlation, volVol, meanVol, meanReversRate,
                volCrashRate, 0.0, 0.0, volCrashSizeMean, 0.0, false));
        }
    }
    
    return volSVCJ;
}



CVolParam* VolSVJJ::createVolParam() const{
    return new SVJJVolParam();
}


IObject* VolSVJJ::defaultCtor(){
    return new VolSVJJ();
}


void VolSVJJ::SVJJVolParam::ComputeImpVol(const CVolBase*          vol,
                                          const CLatticeDouble&    strikes,
                                          const DateTimeArray&     maturities,
                                          CLatticeDouble&          impV) const{
    // turn the vol into what we must have
    const VolSVJJ* myVol = static_cast<const VolSVJJ *>(vol);
    // then just pass through the parameterised vol
    myVol->ComputeImpVol(strikes, maturities, impV);
}

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolSVJJ::SVJJVolParam::spotVolSurfaceFromStrikes(
    const CVolBase*       vol,
    const CDoubleArray&   strikes) const{
    // turn the vol into what we must have
    const VolSVJJ* myVol = 
        static_cast<const VolSVJJ *>(vol);
    // then just pass through the parameterised vol
    return myVol->spotVolSurfaceFromStrikes(strikes);
}

void VolSVJJ::validatePop2Object(){
    try {
        Calibrator::IAdjustable::checkRange(this);
        // build diffusion + stockCrash + volCrash + commonCrash
        update();
    }
    catch(exception& e){
        throw ModelException(e, "VolSVJJ::validatePop2Object");
    }
}

void VolSVJJ::ComputeImpVol(const CLatticeDouble&      strikes,
                            const DateTimeArray&       maturities,
                            CLatticeDouble&            impV) const {
    static const string routine("VolSVJJ::ComputeImpVol");
    throw ModelException(routine, "Not supported");
    if ((maturities.size() != strikes.size()) ||
        (maturities.size() != impV.size())) {
        throw ModelException(routine, "Size mismatch between strikes ("+ 
                             Format::toString(strikes.size()) +
                             "), maturities ("+ 
                             Format::toString(maturities.size())+
                             ") and impV ("+ 
                             Format::toString(impV.size())+ ")");
    }
    
    for (int iMat = 0; iMat < maturities.size(); iMat++) {
        if (strikes[iMat].size() != impV[iMat].size()){
            throw ModelException(routine, "Size mismatch between strikes"
                                 " & maturities for Mat " +
                                 maturities[iMat].toString() +
                                 " (n "+ Format::toString(iMat) + ")");
        }
        for (int iStrike = 0; iStrike < strikes[iMat].size(); iStrike ++) {
            impV[iMat][iStrike] = 0.0;
        }
    }
}

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolSVJJ::spotVolSurfaceFromStrikes(
    const CDoubleArray&   strikes) const{
    static const string routine("VolSVJJ::spotVolSurfaceFromStrikes");
    try{
        throw ModelException(routine, "Not supported");

        const VolSurface* backbone = getBackboneSurface();
        const DateTimeArray& dates = backbone->getDates();
        CDoubleMatrix matrix(strikes.size(),
                             dates.size());
        for (int iStrike = 0; iStrike < strikes.size(); iStrike++) {
            for (int iMat = 0; iMat < dates.size(); iMat++) {
                matrix[iStrike][iMat] = 0.0;
            }
        }
        /** for performance need constructor that takes in
            cached values (to do) */
        VolSurface* volSurf = 
            new VolSurface(getName(),
                           timeMetric.get(),
                           strikes,
                           matrix,
                           backbone->getExpiries().get(),
                           baseDate);
        return volSurf;
    } catch (exception& e){
        throw ModelException(e, routine);
    }
} 

/** Invoked when Class is 'loaded' */
void VolSVJJ::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VolSVJJ, clazz);
    SUPERCLASS(VolBaseParam);
//    IMPLEMENTS(IVolatilityBS);
//    IMPLEMENTS(IVolatilityDVF);
    IMPLEMENTS(Calibrator::IAdjustable);
    IMPLEMENTS(IDynamicsParameter);
    IMPLEMENTS(MarketDataConvert::IConvert<VolSV>);
    IMPLEMENTS(MarketDataConvert::IConvert<VolSVJ>);
    IMPLEMENTS(MarketDataConvert::IConvert<VolSVCJ>);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(initialVol, "initialVol");
    FIELD(meanVol, "meanVol");
    FIELD(meanReversRate, "meanReversRate");
    FIELD(volVol, "volVol");
    FIELD(correlation, "correlation");
    FIELD(stockCrashRate, "stockCrashRate");
    FIELD(stockCrashSizeMean, "stockCrashSizeMean");
    FIELD(stockCrashSizeUncertainty, "stockCrashSizeUncertainty");
    FIELD(volCrashRate, "volCrashRate");
    FIELD(volCrashSizeMean, "volCrashSizeMean");
    FIELD(commonCrashRate, "commonCrashRate");
    FIELD(commonStockCrashSizeMean, "commonStockCrashSizeMean");
    FIELD(commonStockCrashSizeUncertainty, "commonStockCrashSizeUncertainty");
    FIELD(commonVolCrashSizeMean, "commonVolCrashSizeMean");
    FIELD(stockVolCrashSizeCorrelation, "stockVolCrashSizeCorrelation");

    // transient        
    FIELD(baseDate, "");
    FIELD_MAKE_TRANSIENT(baseDate);
    FIELD(diffusion, "");
    FIELD_MAKE_TRANSIENT(diffusion);
    FIELD(stockCrash, "");
    FIELD_MAKE_TRANSIENT(stockCrash);
    FIELD(volCrash, "");
    FIELD_MAKE_TRANSIENT(volCrash);
    FIELD(commonCrash, "");
    FIELD_MAKE_TRANSIENT(commonCrash);
    // add our fields and their ranges to central list
    Calibrator::IAdjustable::registerField(
        clazz, "initialVol",
        new Range(Heston::RangeDef::initialVol));
    Calibrator::IAdjustable::registerField(
        clazz, "meanVol", 
        new Range(Heston::RangeDef::meanVol));
    Calibrator::IAdjustable::registerField(
        clazz, "meanReversRate",
        new Range(Heston::RangeDef::meanReversRate));
    Calibrator::IAdjustable::registerField(
        clazz, "volVol", 
        new Range(Heston::RangeDef::volVol));
    Calibrator::IAdjustable::registerField(
        clazz, "correlation", 
        new Range(Heston::RangeDef::correlation));
    Calibrator::IAdjustable::registerField(
        clazz, "stockCrashRate",
        new Range(MertonCrash::RangeDef::crashRate));
    Calibrator::IAdjustable::registerField(
        clazz, "stockCrashSizeMean", 
        new Range(MertonCrash::RangeDef::crashSizeMean));
    Calibrator::IAdjustable::registerField(
        clazz, "stockCrashSizeUncertainty", 
        new Range(MertonCrash::RangeDef::crashSizeUncertainty));
    Calibrator::IAdjustable::registerField(
        clazz, "volCrashRate",
        new Range(Heston::VolCrash::RangeDef::crashRate));
    Calibrator::IAdjustable::registerField(
        clazz, "volCrashSizeMean", 
        new Range(Heston::VolCrash::RangeDef::crashSizeMean));
    Calibrator::IAdjustable::registerField(
        clazz, "commonCrashRate",
        new Range(Heston::CommonCrash::RangeDef::crashRate));
    Calibrator::IAdjustable::registerField(
        clazz, "commonStockCrashSizeMean", 
        new Range(Heston::CommonCrash::RangeDef::stockCrashSizeMean));
    Calibrator::IAdjustable::registerField(
        clazz, "commonStockCrashSizeUncertainty", 
        new Range(Heston::CommonCrash::RangeDef::stockCrashSizeUncertainty));
    Calibrator::IAdjustable::registerField(
        clazz, "commonVolCrashSizeMean", 
        new Range(Heston::CommonCrash::RangeDef::volCrashSizeMean));
    Calibrator::IAdjustable::registerField(
        clazz, "stockVolCrashSizeCorrelation", 
        new Range(Heston::CommonCrash::RangeDef::crashCorrelation));
}

VolSVJJ::VolSVJJ(): 
VolBaseParam(TYPE),
initialVol(Heston::DefaultVal::initialVol),
correlation(Heston::DefaultVal::correlation),
volVol(Heston::DefaultVal::volVol),
meanVol(Heston::DefaultVal::meanVol),
meanReversRate(Heston::DefaultVal::meanReversRate),
stockCrashRate(MertonCrash::DefaultVal::crashRate),
stockCrashSizeMean(MertonCrash::DefaultVal::crashSizeMean),
stockCrashSizeUncertainty(MertonCrash::DefaultVal::crashSizeUncertainty),
volCrashRate(Heston::VolCrash::DefaultVal::crashRate),
volCrashSizeMean(Heston::VolCrash::DefaultVal::crashSizeMean),
commonCrashRate(Heston::CommonCrash::DefaultVal::crashRate),
commonStockCrashSizeMean(Heston::CommonCrash::DefaultVal::stockCrashSizeMean),
commonStockCrashSizeUncertainty(Heston::CommonCrash::DefaultVal::stockCrashSizeUncertainty),
commonVolCrashSizeMean(Heston::CommonCrash::DefaultVal::volCrashSizeMean),
stockVolCrashSizeCorrelation(Heston::CommonCrash::DefaultVal::crashCorrelation){}

/*** Build the parameterised vol and cache any values **/
void VolSVJJ::buildCache() {
    const VolSurface* backbone = getBackboneSurface();
    baseDate = backbone->getBaseDate();
}

//////////// class FourierProcessSVJJ ///////////////////

/* Started Log Return */
Complex VolSVJJ::scalelessCumulant(const StFourierProcessLogRtn& process,
                                   const StFourierProductLogRtn& product, 
                                   const Complex& z, 
                                   const DateTime& matDate) const{
    static const string method = "VolSVJJ::scalelessCumulant";
    try{
        double tau = timeMetric->yearFrac(baseDate,
                                          matDate);

        Complex alpha, beta;
        calcJointLapAlphaBeta(tau,
                              z,
                              0.0,
                              0.0,
                              alpha,
                              beta);
        return (alpha + beta * Maths::square(initialVol));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Fwd starting Log Return  */
Complex VolSVJJ::scalelessCumulant(const FwdStFourierProcessLogRtn& process,
                                   const FwdStFourierProductLogRtn& product, 
                                   const Complex& z, 
                                   const DateTime& matDate) const{
    static const string method = "VolSVJJ::scalelessCumulant";
    try{
        /* Time fractions from today till start of option and then from
           start of option till maturity */
        double t = timeMetric->yearFrac(baseDate,
                                        product.getStartDate());
        double tau = timeMetric->yearFrac(product.getStartDate(),
                                          matDate);

        /* Calculate alpha, beta component from start date till maturity */
        Complex alpha_tT, beta_tT;
        calcJointLapAlphaBeta(tau,
                              z,
                              0.0,
                              0.0,
                              alpha_tT,
                              beta_tT);
        /* Then, from today till start date */
        Complex alpha_t, beta_t;
        calcJointLapAlphaBeta(t,
                              0.0,
                              beta_tT,
                              0.0,
                              alpha_t,
                              beta_t);
        return (alpha_t + beta_t * Maths::square(initialVol) + alpha_tT);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Started integrated variance */
Complex VolSVJJ::cumulant(const StFourierProcessIntVar& process,
                          const StFourierProductIntVar& product, 
                          const Complex z, 
                          const DateTime& matDate) const {    
    static const string method = "VolSVJJ::cumulant";
    try{
        double tau = timeMetric->yearFrac(baseDate,
                                          matDate);

        Complex alpha, beta;
        calcJointLapAlphaBeta(tau,
                              0.0,
                              0.0,
                              z,
                              alpha,
                              beta);
        return (alpha + beta * Maths::square(initialVol));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
    
/* Forward Starting integrated variance */
Complex VolSVJJ::cumulant(const FwdStFourierProcessIntVar& process,
                          const FwdStFourierProductIntVar& product, 
                          const Complex z, 
                          const DateTime& matDate) const {        
    static const string method = "VolSVJJ::cumulant";
    try{
        /* Time fractions from today till start of option and then from
           start of option till maturity */
        double t = timeMetric->yearFrac(baseDate,
                                        product.getStartDate());
        double tau = timeMetric->yearFrac(product.getStartDate(),
                                          matDate);

        /* Calculate alpha, beta component from start date till maturity */
        Complex alpha_tT, beta_tT;          
        calcJointLapAlphaBeta(tau,
                              0.0,
                              0.0,
                              z,
                              alpha_tT,
                              beta_tT);

        /* Then, from today till start date */
        Complex alpha_t, beta_t;
        calcJointLapAlphaBeta(t,
                              0.0,
                              beta_tT,
                              0.0,
                              alpha_t,
                              beta_t);
        
        return (alpha_t + beta_t * Maths::square(initialVol) + alpha_tT);

    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/** Calculates the components alpha and beta that appear in 
    the exponent of the time-t joint Bilateral Laplace 
    of X_T = (Y_T, V_T, I_T) where Y_T = ln(S_T / F(0, T)) is the dimension-less
    log spot at time T, V_T is the instantaneous variance at time T
    and I_T is the instantaneous variance from 0 to time T
    The Bilateral Laplace transform is of the form
        exp(alpha(tau, u1, u2, u3) 
            + u1 * Y_t 
            + beta(tau, u1, u2, u3) * V_t 
            + u3 * I_t)
    where tau = T - t and the complex numbers u1, u2, u3 are the frequencies 
    wrt Y_T, V_T and I_T, respectively. */
void VolSVJJ::calcJointLapAlphaBeta(double         tau,    // tau = T - t (in years)
                                    const Complex& u1,
                                    const Complex& u2,
                                    const Complex& u3,
                                    Complex&       alpha,
                                    Complex&       beta) const{
    static const string method = "VolSVJJ::calcJointLapAlphaBeta";

    try{    // in case Complex operations fail for some reason
        /* Calculate Heston contribution */
        diffusion->calcAlphaBeta(tau,
                                 u1,
                                 u2,
                                 u3,
                                 alpha,
                                 beta);

        /* Calculate stock-only jump contribution to alpha */
        alpha += stockCrash->calcAlpha(tau, u1);

        /* Calculate vol-only jump contribution to alpha */
        alpha += volCrash->calcAlpha(tau,
                                     u1,
                                     u2,
                                     u3);

        /* Calculate common jump contribution to alpha */
        alpha += commonCrash->calcAlpha(tau,
                                        u1,
                                        u2,
                                        u3);
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

CClassConstSP const VolSVJJ::TYPE =
CClass::registerClassLoadMethod("VolSVJJ", typeid(VolSVJJ), load);

CClassConstSP const VolSVJJ::SVJJVolParam::TYPE =
CClass::registerClassLoadMethod("VolSVJJ::SVJJVolParam", typeid(SVJJVolParam), load);

/** Needed for IAdjustable interface. Returns market data name for vol */
string VolSVJJ::getName() const{
    return CVolBaseParamSurface::getName();
}

/** Called after adjustments have been made to fields (eg calibrator) */
void VolSVJJ::fieldsUpdated(const CFieldArray& fields){
    update();
}

/** Called after adjustments have been made to fields */
void VolSVJJ::update(){
    // refresh diffusion
    diffusion = HestonSP(new Heston(initialVol,
                                    meanVol,
                                    meanReversRate,
                                    volVol,
                                    correlation));
    // refresh stockCrash
    stockCrash = MertonCrashSP(new MertonCrash(stockCrashRate,
                                               stockCrashSizeMean,
                                               stockCrashSizeUncertainty));
    // refresh volCrash
    volCrash = Heston::VolCrashSP(new Heston::VolCrash(volCrashRate,
                                                       volCrashSizeMean,
                                                       diffusion));

    // refresh commonCrash
    commonCrash = Heston::CommonCrashSP(new Heston::CommonCrash(commonCrashRate,
                                                                commonStockCrashSizeMean,
                                                                commonStockCrashSizeUncertainty,
                                                                commonVolCrashSizeMean,
                                                                stockVolCrashSizeCorrelation,
                                                                diffusion));
}

template<>
const CClassConstSP VolSVJJConvert::TYPE = 
CClass::registerInterfaceLoadMethod("MarketDataConvert::IConvert<VolSVJJ>", typeid(VolSVJJConvert), load);

/* external symbol to allow class to be forced to be linked in */
bool VolSVJJLinkIn(){
    return (VolSVJJ::TYPE != 0) && (VolSVJJConvert::TYPE != 0);
}

DRLIB_END_NAMESPACE
