//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolStochGarf.cpp
//
//   Description : CommonStockVolDiff + CommonStockVolJump
//
//   Date        : 20 Apr 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolSVJ.hpp"
#include "edginc/VolProcessedStochGarf.hpp"

DRLIB_BEGIN_NAMESPACE

// constructor 
VolStochGarf::VolStochGarf(): 
    VolBaseParam(TYPE),
    alpha(0.0),
    beta(0.0),
    meanReversRate(1.0),
    volVol(Heston::DefaultVal::volVol),
    correlation(0.0),
    atmVol(0.1){}

// constructor with scaled atm vol level
//VolStochGarf::VolStochGarf(const VolStochGarf& vol, double scale): VolBaseParam(TYPE),
VolStochGarf::VolStochGarf(const VolStochGarf& vol, double scale): VolBaseParam(TYPE),
    alpha(vol.alpha),
    beta(vol.beta),
    meanReversRate(vol.meanReversRate),
    volVol(vol.volVol),
    correlation(vol.correlation){
    atmVol = getATMVol(scale);
}

void VolStochGarf::validatePop2Object(){
    try {
        Calibrator::IAdjustable::checkRange(this);
    }
    catch(exception& e){
        throw ModelException(e, "VolStochGarf::validatePop2Object");
    }
}

DEFINE_TEMPLATE_TYPE(VolStochGarfArray);


IVolProcessed* VolStochGarf::getProcessedVol(const CVolRequest* volRequest,
                               const CAsset*      asset) const {
    return new VolProcessedStochGarf(VolStochGarfSP(copy(this)));    

    //VolRequestRaw request;
    //CVolProcessedSP processedVol(asset->getProcessedVol(&request));
    //return new VolProcessedStochGarf((dynamic_cast<VolStochGarf*>(processedVol.get())));    
}

/*
IVolProcessed* VolStochGarf::getProcessedVol(double scale) const {
    return  new VolProcessedStochGarf(VolStochGarfSP(copy(this)), scale);

    //VolRequestRaw request;
    //CVolProcessedSP processedVol(asset->getProcessedVol(&request));
    //return new VolProcessedStochGarf((dynamic_cast<VolStochGarf*>(processedVol.get())));    
}
*/

void VolStochGarf::StochGarfVolParam::ComputeImpVol(const CVolBase*          vol,
                                          const CLatticeDouble&    strikes,
                                          const DateTimeArray&     maturities,
                                          CLatticeDouble&          impV) const{
    // turn the vol into what we must have
    const VolStochGarf* myVol = static_cast<const VolStochGarf *>(vol);
    // then just pass through the parameterised vol
    myVol->ComputeImpVol(strikes, maturities, impV);
}

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolStochGarf::StochGarfVolParam::spotVolSurfaceFromStrikes(
    const CVolBase*       vol,
    const CDoubleArray&   strikes) const{
    // turn the vol into what we must have
    const VolStochGarf* myVol = 
        static_cast<const VolStochGarf *>(vol);
    // then just pass through the parameterised vol
    return myVol->spotVolSurfaceFromStrikes(strikes);
}

// just return 0 as not supported.
void VolStochGarf::ComputeImpVol(const CLatticeDouble&      strikes,
                            const DateTimeArray&       maturities,
                            CLatticeDouble&            impV) const {
    static const string routine("VolStochGarf::ComputeImpVol");
//    throw ModelException(routine, "Not supported");
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

double VolStochGarf::StochGarfVolParam::computeLocVol(const CVolBase* vol, 
                                     const double spot,
                                     const DateTime& maturity) const{
    const VolStochGarf* myVol = static_cast<const VolStochGarf *>(vol);
    // then just pass through the parameterised vol
    return myVol->computeLocVol(spot, maturity);
};

double VolStochGarf::computeLocVol(const double spot, const DateTime date) const{
    double x = log(spot/strikeRef);
    return atmVol + alpha * x + beta * x * x ;
}

// return ATM vol level, by taking scale.  Scale should be -1, 0, 1 for 3 status
double VolStochGarf::getATMVol(const double scale){
    return atmVol * exp(scale * volVol);
}

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolStochGarf::spotVolSurfaceFromStrikes(
    const CDoubleArray&   strikes) const{
    static const string routine("VolStochGarf::spotVolSurfaceFromStrikes");
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
void VolStochGarf::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VolStochGarf, clazz);
    SUPERCLASS(VolBaseParam);
    IMPLEMENTS(IVolatilityBS);
    IMPLEMENTS(IVolatilityDVF);
    IMPLEMENTS(Calibrator::IAdjustable);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(alpha, "alpha : skew parameters for Local vol");
    FIELD(beta, "beta : convexity parameters for Local vol");
    FIELD(meanReversRate, "mean Reversion Rate");
    FIELD(volVol, "volVol");
    FIELD(correlation, "correlation");
    FIELD(atmVol, "ATM vol");
    FIELD(strikeRef, "strike level of ATM");
    // transient        
    FIELD(baseDate, "");
    FIELD_MAKE_TRANSIENT(baseDate);
    // add our fields and their ranges to central list
//    Calibrator::IAdjustable::registerField(
//        clazz, "alpha",
//        new Range(Heston::RangeDef::alpha));
//    Calibrator::IAdjustable::registerField(
//        clazz, "beta", 
//        new Range(Heston::RangeDef::beta));
//    Calibrator::IAdjustable::registerField(
//        clazz, "meanReversRate",
//        new Range(Heston::RangeDef::meanReversRate));
//    Calibrator::IAdjustable::registerField(
//        clazz, "volVol", 
//        new Range(Heston::RangeDef::volVol));
//    Calibrator::IAdjustable::registerField(
//        clazz, "correlation", 
//        new Range(Heston::RangeDef::correlation));

}

/** populate from market cache */
void VolStochGarf::getMarket(const IModel* model, const MarketData* market) {

    // call parent method first
    VolBaseParam::getMarket(model, market);
    
    // const VolSurface* backbone = getBackboneSurface();
    /** returns a constant reference to surface to be used for the backbone */
    VolSurfaceSP backbone = VolSurfaceSP(VolSurfaceSP::dynamicCast(
                market->GetData(getName(),VolSurface::TYPE)));
    backbone->getMarket(model,market);  // put holiday and so on into volsurface.
    baseDate = backbone->getBaseDate();
}


CClassConstSP const VolStochGarf::TYPE =
CClass::registerClassLoadMethod("VolStochGarf", typeid(VolStochGarf), load);

CClassConstSP const VolStochGarf::StochGarfVolParam::TYPE =
CClass::registerClassLoadMethod("VolStochGarf::StochGarfVolParam", typeid(StochGarfVolParam), load);

/** Needed for IAdjustable interface. Returns market data name for vol */
string VolStochGarf::getName() const{
    return CVolBaseParamSurface::getName();
}

/* external symbol to allow class to be forced to be linked in */
bool VolStochGarfLinkIn(){
    return (VolStochGarf::TYPE != 0);
}

DRLIB_END_NAMESPACE
