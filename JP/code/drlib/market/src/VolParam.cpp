//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolParam.cpp
//
//   Description : Interface for Parametrized Vol Surface
//
//   Author      : Regis Guichard
//
//   Date        : 02 Mai 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolBaseParamSurface.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

// CVolParam stuff -----------------------------------------------------------
CVolParam::~CVolParam(){}

CClassConstSP const CVolParam::TYPE = CClass::registerClassLoadMethod(
    "VolParam", typeid(CVolParam), load);

void CVolParam::load(CClassSP& clazz){
    REGISTER(CVolParam, clazz);
    SUPERCLASS(CObject);
}

CVolParam::CVolParam(const CClassConstSP& clazz) : CObject(clazz){}

/** Given a maturity date and a strike, returns the implied vol.
    The default implementation is just a wrapper around the CLattice version.*/
double CVolParam::ComputeImpVol(const CVolBase* vol,
                                double          strike,
                                const DateTime& maturity) const {
    CSliceDouble singlePoint(1);
    DateTimeArray maturities(1, maturity);
    CSliceDouble strikes(1);
    strikes[0] = strike;
    ComputeImpVol(vol, strikes, maturities , singlePoint);
    return singlePoint[0];
};

/** Computes forward starting implied vols. Note default implementation
    provided. */
void CVolParam::computeFwdStartImpVol(
    const CVolBase*          vol,
    const FwdStart&          fwdStart,
    const CLatticeDouble&    strikes,
    bool                     strikesRelToSpotAtStart,
    const DateTimeArray&     maturities,
    CLatticeDouble&          impV) const{
    static const string method ("CVolParam::computeFwdStartImpVol");
    static const double YEAR_FRAC_FOR_INSTANENOUS_VOL = 
        10.0/DateTime::TIME_IN_YEAR;
    try{
        if (maturities.empty()){
            return; // get out now
        }
        // get atm vols from start date to each maturity
        const DateTime& startDate = fwdStart.getStartDate();
        if (startDate.isGreater(maturities[0])){
            // we could but what is the correct behaviour?
            throw ModelException(method, "This method does not support "
                                 "adjustment to vols before instrument start "
                                 "date");
        }
        double spot =  fwdStart.getSpot();
        // want to calculate vol from start date to each maturity using
        // strike = spot. Must use ComputeImpVol() for consistency
        vector<int> arrayOfOne(maturities.size(), 1);
        CLatticeDouble latticeOfSpot(arrayOfOne);
        // populate with spot
        for (int s = 0; s < maturities.size(); s++){
            latticeOfSpot[s][0] = spot;
        }
        // space for results
        CLatticeDouble  fwdVols(arrayOfOne);
        // do calc - note this is from today so must adjust
        ComputeImpVol(vol, latticeOfSpot, maturities, fwdVols);
        // .. so now do adjustment - get vol to start date
        double volToStart = ComputeImpVol(vol, spot, startDate);
        // then 'subtract' this vol
        const TimeMetricConstSP& metric = fwdStart.getTimeMetric();
        const DateTime& baseDate = fwdStart.getValueDate();
        double yearsToStart = metric->yearFrac(baseDate, startDate);
        double varToStart =  Maths::square(volToStart) * yearsToStart;
        for (int m = 0; m < maturities.size(); m++){
            double yearFrac = metric->yearFrac(startDate, maturities[m]);
            if (!Maths::isPositive(yearFrac)){
                // instaneous fwd vol is wanted - a bit of a pain
                DateTime adjDate = 
                    metric->impliedTime(startDate, 
                                        YEAR_FRAC_FOR_INSTANENOUS_VOL,
                                        yearFrac);
                fwdVols[m][0] = ComputeImpVol(vol, spot, adjDate);
            }
            double varToMat = Maths::square(fwdVols[m][0]) * 
                (yearFrac+yearsToStart);
            if (varToStart > varToMat) {
                throw ModelException(method, "Negative variance between (" + 
                                     startDate.toString() + ") and (" +
                                     maturities[m].toString() + 
                                     ") at strike " + 
                                     Format::toString(spot) + 
                                     " for vol " + vol->getName());
            }
            fwdVols[m][0] = sqrt((varToMat-varToStart)/yearFrac);
        }
        // then do the same except from today to effective maturity
        DateTimeArray effectiveMats(maturities.size());
        for (int i = 0; i < effectiveMats.size(); i++){
            effectiveMats[i] = baseDate.add(maturities[i].subtract(startDate));
        }
        // space for results
        CLatticeDouble  spotVols(arrayOfOne);
        // do calc
        ComputeImpVol(vol, latticeOfSpot, effectiveMats, spotVols);

       /* then calculate normal vol from baseDate to effectiveMats for
           adjusted strike (note: this will be 'tweaked' appropriately by
           VolParamTweak) */
        double strikeScalingFactor = strikesRelToSpotAtStart? 
            spot: spot/fwdStart.getFwdAtStart();
        // scale all strkes by strikeScalingFactor - avoid copy
        CLatticeDouble&    adjStrikes = const_cast<CLatticeDouble&>(strikes);
        try{
            adjStrikes.scale(strikeScalingFactor);
            ComputeImpVol(vol, strikes, effectiveMats, impV);
        } catch (exception&){
             adjStrikes.scale(1/strikeScalingFactor);
             throw;
        }
        adjStrikes.scale(1/strikeScalingFactor);
        
        // then add on spread
        for (int j = 0; j < impV.size(); j++){
            CSliceDouble& impVol = impV[j];
            double  spread = fwdVols[j][0]-spotVols[j][0];
            for (int i = 0; i < impVol.size(); i++){
                impVol[i] += spread;
                if (!Maths::isPositive(impVol[i])){
                    impVol[i] = 0.0; // floor it to 0 !!!
                }
            }
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }
}
CVolParam::FwdStart::~FwdStart(){}

/** Checks we're not fwd starting. String parameter is for exception */
void CVolParam::FwdStart::checkFwdStarting(const string& method) const{
    if (!fwdStarting){
        throw ModelException(method,
                             "Data not available if not forward starting");
    }
}

/** Constructor for not forward starting */
CVolParam::FwdStart::FwdStart(): fwdStarting(false), 
    spot(0.0), fwdAtStart(0.0){}

/** Constructor if forward starting. Uses asset to get spot and
    fwd at start date */
CVolParam::FwdStart::FwdStart(const DateTime&          valueDate,
                              const DateTime&          startDate,
                              const TimeMetricConstSP& timeMetric,
                              const CAsset*            asset):
    fwdStarting(startDate.isGreater(valueDate)), 
    valueDate(valueDate), startDate(startDate), metric(timeMetric),
    spot(asset->getSpot()), fwdAtStart(asset->fwdValue(startDate)){}

/** Constructor - determines everything from volRequest and asset */
CVolParam::FwdStart::FwdStart(const DateTime&          valueDate,
                              const CVolRequestDVF*    volRequest,
                              const TimeMetricConstSP& timeMetric,
                              const CAsset*            asset):
    valueDate(valueDate), startDate(volRequest->getStartDate()), 
    metric(timeMetric), spot(asset->getSpot()), 
    fwdAtStart(asset->fwdValue(startDate)){
    fwdStarting = startDate.isGreater(valueDate);
}


/** Constructor with everything specified explicitly */
CVolParam::FwdStart::FwdStart(const DateTime&          valueDate,
                              const DateTime&          startDate,
                              const TimeMetricConstSP& timeMetric,
                              double                   spot,
                              double                   fwdAtStart):
    fwdStarting(startDate.isGreater(valueDate)), 
    valueDate(valueDate), startDate(startDate), metric(timeMetric),
    spot(spot), fwdAtStart(fwdAtStart){}

/** Returns the instrument's start date */
const DateTime& CVolParam::FwdStart::getStartDate() const{
    checkFwdStarting("CVolParam::FwdStart::getStartDate");
    return startDate;
}
    
/** Should a forward starting adjustment be made */
bool CVolParam::FwdStart::isFwdStarting() const{
    return fwdStarting;
}

/** Returns value date */
const DateTime&  CVolParam::FwdStart::getValueDate() const{
    checkFwdStarting("CVolParam::FwdStart::getValueDate");
    return valueDate;
}

/** Returns the time metric */
const TimeMetricConstSP& CVolParam::FwdStart::getTimeMetric() const{
    checkFwdStarting("CVolParam::FwdStart::getTimeMetric");
    return metric;
}

/** Returns the asset's (whose vol this is) spot price */
double CVolParam::FwdStart::getSpot() const{
    checkFwdStarting("CVolParam::FwdStart::getSpot");
    return spot;
}

/** Returns the asset's (whose vol this is) spot price at the start
    date */
double CVolParam::FwdStart::getFwdAtStart() const{
    checkFwdStarting("CVolParam::FwdStart::getFwdAtStart");
    return fwdAtStart;
}



DRLIB_END_NAMESPACE
