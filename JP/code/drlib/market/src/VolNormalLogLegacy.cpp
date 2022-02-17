//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolNormalLogLegacy.cpp
//
//   Description : Parametrised Vol Surface: Backbone + NormDist(ln(k))
//                  Vol(K, T) = BackBone(StrikeRef) - Factor * CutOffSpeed *
//                             Skew1Y90100(T) * 
//                                 (N( LN(K/StrikeRef)/CutOffSpeed) -0.5)
//      Where Factor is such that the 90 - 100 smile is exactly the input on 1Y
//            Skew1Y90100(T) varies with T as a varying power function
//
//   Author      : Mark A Robson
//
//   Date        : 22 Aug 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolBaseParamSurface.hpp"
#include "edginc/VolatilityDVF.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"

DRLIB_BEGIN_NAMESPACE

/** A parameterised CVolBase that supports BS and DVF views of the world */
class VolNormalLogLegacy: public CVolBaseParamSurface,
                          virtual public IVolatilityBS,
                          virtual public IVolatilityDVF {
private:
    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class VNLVolParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        // constructor
        VNLVolParam(): CVolParam(TYPE){}

        ~VNLVolParam(){}

        virtual void ComputeImpVol(const CVolBase*          vol,
                                   const CLatticeDouble&    strikes,
                                   const DateTimeArray&     maturities,
                                   CLatticeDouble&          impV) const{
            // turn the vol into what we must have
            const VolNormalLogLegacy* myVol = 
                static_cast<const VolNormalLogLegacy *>(vol);
            // then just pass through the parameterised vol
            myVol->computeImpVol(strikes, maturities, impV);
        }

        /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
        virtual VolSurface* spotVolSurfaceFromStrikes(
            const CVolBase*       vol,
            const CDoubleArray&   strikes) const{
            // turn the vol into what we must have
            const VolNormalLogLegacy* myVol = 
                static_cast<const VolNormalLogLegacy *>(vol);
            // then just pass through the parameterised vol
            return myVol->spotVolSurfaceFromStrikes(strikes);
        }

    private:
        static void load(CClassSP& clazz){
            REGISTER(VNLVolParam, clazz);
            SUPERCLASS(CVolParam);
        }
    };
public:
    friend class VNLVolParam;
    static CClassConstSP const TYPE;

    ~VolNormalLogLegacy(){}

    void validatePop2Object(){
        static const char routine[] = "VolNormalLogLegacy::validatePop2Object";
        if (! Maths::isPositive(cutOffSpeedStart) ||
            ! Maths::isPositive(cutOffSpeedEnd) ||
            ! Maths::isPositive(cutOffSpeed1Y)) {
            throw ModelException(routine, "cutOffSpeed must be positive");
        }
    }

    /** method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const{
        return new VNLVolParam();
    }
    
private:

    void computeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const {
        static const string routine("VolNormalLogLegacy::computeImpVol");
        if ((maturities.size() != strikes.size()) ||
            (maturities.size() != impV.size())) {
            throw ModelException(routine, "Size mismatch between strikes ("+ 
                                 Format::toString(strikes.size()) +
                                 "), maturities ("+ 
                                 Format::toString(maturities.size())+
                                 ") and impV ("+ 
                                 Format::toString(impV.size())+ ")");
        }  
        // calculate backbone vols
        DoubleArray backBoneVols(maturities.size());
        const VolSurface* backbone = getBackboneSurface();
        CVolProcessedBSSP  atmVolCurve(backbone->getProcessedVol(strikeRef));
        atmVolCurve->CalcVol(baseDate, maturities, CVolProcessedBS::fromFirst,
                             backBoneVols);
        for (int iMat = 0; iMat < maturities.size(); iMat++) {
            if (strikes[iMat].size() != impV[iMat].size()){
                throw ModelException(routine, "Size mismatch between strikes"
                                     " & maturities for Mat " +
                                     maturities[iMat].toString() +
                                     " (n "+ Format::toString(iMat) + ")");
            }
            // floor yearToDate with 1D to avoid blowing up as t --> 0
            double yearToDate= Maths::max(1.0/DateTime::DAYS_PER_YEAR,
                                          baseDate.yearFrac(maturities[iMat]));
            double powerToDate = calcPowerToDate(yearToDate);
            double tToP = pow(yearToDate, powerToDate);
            double skew = skew1Y90100/tToP;
            double cutOffSpeed = calcCutOffSpeed(yearToDate);
            double factor = 1.0/(N1(log(90.0/100.0)/cutOffSpeed)-0.5);  
            for (int iStrike = 0; iStrike < strikes[iMat].size(); iStrike ++) {
                double logStrike = log(strikes[iMat][iStrike]/strikeRef);
                logStrike /= cutOffSpeed;
                logStrike  = N1(logStrike) - 0.5;
                // NB don't multiply by cutOffSpeed since we skipped
                // it from the factor calculation
                // to do: should we floor the vol here? 
                impV[iMat][iStrike]  = factor * skew * logStrike;
                impV[iMat][iStrike] += backBoneVols[iMat]; 
            }
        }
    }

    
    /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
    VolSurface* spotVolSurfaceFromStrikes(
        const CDoubleArray&   strikes) const{
        static const string routine(
            "VolNormalLogLegacy::spotVolSurfaceFromStrikes");
        try{
            const VolSurface* backbone = getBackboneSurface();
            CDoubleMatrix matrix(strikes.size(), atmVols.size());
            // loop over strikes
            for (int iStrike = 0; iStrike < strikes.size(); iStrike ++) {
                double logStrike  = log(strikes[iStrike]/strikeRef);
                // loop over maturities
                for (int iMat = 0; iMat < atmVols.size(); iMat++) {
                    double thisLogStrike = logStrike/cutOffSpeedCache[iMat];
                    thisLogStrike  = N1(thisLogStrike) - 0.5;
                    // cutOffSpeed is taken care of by 'factor'
                    double vol = factorCache[iMat] * skewCache[iMat] * 
                        thisLogStrike + atmVols[iMat];
                    if (vol < volFloor){
                        vol = volFloor;
                    }
                    matrix[iStrike][iMat] = vol;
                }
            }
            /** for performance use special constructors */
            return (new VolSurface(backbone, strikes, matrix));
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    // registered fields
    double        strikeRef;
    double        skew1Y90100;       /* one year skew */
    double        cutOffSpeedStart;
    double        cutOffSpeed1Y;
    double        cutOffSpeedEnd;
    double        powerStart;        /* to bend Skew1Y */
    double        power1Y;           /* to bend Skew1Y */
    double        powerEnd;          /* to bend Skew1Y */
    double        volFloor;          /* optional floor when calculating vols */
    // transient fields (won't appear in dd interface)
    DateTime              baseDate;
    TimeMetricConstSP     timeMetric;
    CDoubleArray          atmVols;     // Vol at each Bench Date
    DoubleArray           factorCache; // cache of 'factor' for bm dates
    DoubleArray           skewCache;   //cache of 'skew' for bm dates
    DoubleArray           cutOffSpeedCache;//cache of cutOffSpeed for bm dates

    //// computes 'alpha' aka cutOffSpeed
    double calcCutOffSpeed(double years) const {
        double alpha;
        if (years > 1){
            alpha = (cutOffSpeed1Y - cutOffSpeedEnd)/years + cutOffSpeedEnd;
        } else {
            alpha = (cutOffSpeed1Y -cutOffSpeedStart)*years + cutOffSpeedStart;
        }
        return alpha;
    }

    //// computes 'power' factor
    double calcPowerToDate(double yearToDate) const{
        double powerToDate;
        if (yearToDate <= 1.0) {
            powerToDate  = powerStart * (1.0 - yearToDate);
            powerToDate += power1Y * yearToDate;
        } else {
            powerToDate  = power1Y  / yearToDate;
            powerToDate += powerEnd * (1.0 - 1.0/yearToDate);
        }
        return powerToDate;
    }
        
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        // NB not visible to EAS/spreadsheet
        REGISTER(VolNormalLogLegacy, clazz);
        SUPERCLASS(CVolBaseParamSurface);
        IMPLEMENTS(IVolatilityBS);
        IMPLEMENTS(IVolatilityDVF);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(strikeRef, "strikeRef");
        FIELD(skew1Y90100, "One year skew");
        FIELD(cutOffSpeedStart, "CutOff Speed at start");
        FIELD(cutOffSpeed1Y, "CutOff Speed at one year");
        FIELD(cutOffSpeedEnd, "CutOff Speed at end");
        FIELD(powerStart, "to bend Skew1Y");
        FIELD(power1Y, "to bend Skew1Y");
        FIELD(powerEnd, "to bend Skew1Y");
        FIELD(volFloor, "floor when calculating vols");
        FIELD_MAKE_OPTIONAL(volFloor);
        FIELD(timeMetric, "used to throw at the VolSurface constructor");
        FIELD_MAKE_TRANSIENT(timeMetric);
        FIELD(baseDate, "Base Date");
        FIELD_MAKE_TRANSIENT(baseDate);
        FIELD(atmVols, "Vol at each Benchmark  Date");
        FIELD_MAKE_TRANSIENT(atmVols);
        FIELD(factorCache, "");
        FIELD_MAKE_TRANSIENT(factorCache);
        FIELD(skewCache, "");
        FIELD_MAKE_TRANSIENT(skewCache);
        FIELD(cutOffSpeedCache, "");
        FIELD_MAKE_TRANSIENT(cutOffSpeedCache);
    }

    VolNormalLogLegacy(): CVolBaseParamSurface(TYPE),
                          strikeRef(0.0), skew1Y90100(0.0), 
                          cutOffSpeedStart(0.0),
                          cutOffSpeed1Y(0.0), cutOffSpeedEnd(0.0),
                          powerStart(0.0), power1Y(0.0), powerEnd(0.0),
                          volFloor(0.0) {
        // empty
    }

    static IObject* defaultCtor(){
        return new VolNormalLogLegacy();
    }

protected:
    /*** Build the parameterised vol and cache any values **/
    void buildCache() {
        const VolSurface* backbone = getBackboneSurface();
        const DateTimeArray&    dates = backbone->getDates();
        CVolProcessedBSSP       atmVolCurve = CVolProcessedBSSP(
            backbone->getProcessedVol(strikeRef));

        baseDate = backbone->getBaseDate();
        timeMetric = atmVolCurve->GetTimeMetric();
        atmVols   = CDoubleArray(dates.size());
        atmVolCurve->CalcVol(baseDate, dates, 
                             CVolProcessedBS::fromFirst ,atmVols);
        factorCache = CDoubleArray(dates.size());
        skewCache = CDoubleArray(dates.size());
        cutOffSpeedCache = CDoubleArray(dates.size());
        for (int i = 0; i < dates.size(); i++){
            double yearToDate = baseDate.yearFrac(dates[i]);// calendar time
            double powerToDate = calcPowerToDate(yearToDate);
            double tToP = pow(yearToDate, powerToDate);
            skewCache[i] = skew1Y90100/tToP;
            cutOffSpeedCache[i] = calcCutOffSpeed(yearToDate);
            factorCache[i] = 1.0/(N1(LOG90OVER100/cutOffSpeedCache[i])-0.5);
        }
    }
private:
    static const double LOG90OVER100;
};

const double VolNormalLogLegacy::LOG90OVER100 = log(90.0/100.0);

CClassConstSP const VolNormalLogLegacy::TYPE =
CClass::registerClassLoadMethod("VolNormalLogLegacy", 
                                typeid(VolNormalLogLegacy), load);

CClassConstSP const VolNormalLogLegacy::VNLVolParam::TYPE =
CClass::registerClassLoadMethod("VolNormalLogLegacy::VNLVolParam", 
                                typeid(VNLVolParam), load);

// external symbol to allow class to be forced to be linked in
bool VolNormalLogLegacyLinkIn(){
    return (VolNormalLogLegacy::TYPE != 0);
}

DRLIB_END_NAMESPACE
