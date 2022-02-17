//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolNormalLog.cpp
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
#include "edginc/Calibrator.hpp"

DRLIB_BEGIN_NAMESPACE

/** A parameterised CVolBase that supports BS and DVF views of the world */
class VolNormalLog: public CVolBaseParamSurface,
                    virtual public IVolatilityBS,
                    virtual public IVolatilityDVF,
                    virtual public DeltaSurface::IShift,
                    public virtual Calibrator::IAdjustable{
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
            const VolNormalLog* myVol = 
                static_cast<const VolNormalLog *>(vol);
            // then just pass through the parameterised vol
            myVol->computeImpVol(strikes, maturities, impV);
        }

        /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
        virtual VolSurface* spotVolSurfaceFromStrikes(
            const CVolBase*       vol,
            const CDoubleArray&   strikes) const{
            // turn the vol into what we must have
            const VolNormalLog* myVol = 
                static_cast<const VolNormalLog *>(vol);
            // then just pass through the parameterised vol
            return myVol->spotVolSurfaceFromStrikes(strikes);
        }

    private:
        static void load(CClassSP& clazz){
            REGISTER(VNLVolParam, clazz);
            SUPERCLASS(CVolParam);
            EMPTY_SHELL_METHOD(defaultCtor);
        }
        static IObject* defaultCtor(){
            return new VNLVolParam();
        }
    };
public:
    friend class VNLVolParam;
    static CClassConstSP const TYPE;

    ~VolNormalLog(){}

    void validatePop2Object(){
        try {
            Calibrator::IAdjustable::checkRange(this);
        }
        catch(exception& e){
            throw ModelException(e, "VolNormalLog::validatePop2Object",
                                 "Failed for vol (" + getName() + ")");
        }
    }

    /** method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const{
        return new VNLVolParam();
    }
    

    string sensName(DeltaSurface* shift) const{
        return getName();
    }

    bool sensShift(DeltaSurface* shift){
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)){
            // only bother if non zero 
            // what to do - spot moves to spot + shift
            // simply move strike ref k -> k (S+dS)/S
            double spot    = shift->getSpot();
            double newSpot = spot * (1.0 + shiftSize);

            strikeRef *= newSpot/spot;

            // get new backbone
            buildCache();
        }
        // DON'T tweak the backbone !!!
        return false; // no more tweaking required here
    }


private:

    void computeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const {
        static const string routine("VolNormalLog::computeImpVol");
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
                                          metric->yearFrac(baseDate,
                                                           maturities[iMat]));
            double powerToDate = calcPowerToDate(yearToDate);
            double tToP = pow(yearToDate, powerToDate);
            double skew = skew1Y90100/tToP;
            double cutOffSpeedUp = calcCutOffSpeed(yearToDate, true);
            double cutOffSpeedDown = calcCutOffSpeed(yearToDate, false);
            double factorDown = 1.0/(N1(log(90.0/100.0)/cutOffSpeedDown)-0.5);
            double factorUp = factorDown*cutOffSpeedUp/cutOffSpeedDown;
            for (int iStrike = 0; iStrike < strikes[iMat].size(); iStrike ++) {
                double strike = strikes[iMat][iStrike];
                double logStrike = log(strike/strikeRef);
                bool   upside = strike >= strikeRef;
                logStrike /= upside? cutOffSpeedUp: cutOffSpeedDown;
                logStrike  = N1(logStrike) - 0.5;
                // NB don't multiply by cutOffSpeed since we skipped
                // it from the factor calculation
                // to do: should we floor the vol here? 
                impV[iMat][iStrike]  = (upside? factorUp: factorDown) * 
                    skew * logStrike;
                impV[iMat][iStrike] += backBoneVols[iMat];
                if (impV[iMat][iStrike] < volFloor){
                    impV[iMat][iStrike] = volFloor;
                }
            }
        }
    }

    
    /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
    VolSurface* spotVolSurfaceFromStrikes(
        const CDoubleArray&   strikes) const{
        static const string routine(
            "VolNormalLog::spotVolSurfaceFromStrikes");
        try{
            const VolSurface* backbone = getBackboneSurface();
            CDoubleMatrix matrix(strikes.size(), atmVols.size());
            // loop over strikes
            for (int iStrike = 0; iStrike < strikes.size(); iStrike ++) {
                double logStrike  = log(strikes[iStrike]/strikeRef);
                bool upside = strikes[iStrike] >= strikeRef;
                // loop over maturities
                for (int iMat = 0; iMat < atmVols.size(); iMat++) {
                    double thisLogStrike = logStrike/
                        (upside? cutOffSpeedUpCache[iMat]: 
                         cutOffSpeedDownCache[iMat]);
                    thisLogStrike  = N1(thisLogStrike) - 0.5;
                    // cutOffSpeed is taken care of by 'factor'
                    // we adjust factorDown to ensure smoothness of derivative
                    // wrt to strike
                    double factor = upside?
                        factorCacheDown[iMat]*cutOffSpeedUpCache[iMat]/
                         cutOffSpeedDownCache[iMat]:
                        factorCacheDown[iMat];
                    double vol = factor * skewCache[iMat] * 
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
    double        cutOffSpeedUpStart;
    double        cutOffSpeedDownStart;
    double        cutOffSpeedUpEnd;
    double        cutOffSpeedDownEnd;
    double        cutOffSpeedEndPoint;
    double        powerStart;        /* to bend Skew1Y */
    double        power1Y;           /* to bend Skew1Y */
    double        powerEnd;          /* to bend Skew1Y */
    double        volFloor;          /* optional floor when calculating vols */
    // transient fields (won't appear in dd interface)
    DateTime              baseDate;
    TimeMetricConstSP     metric;
    CDoubleArray          atmVols;     // Vol at each Bench Date
    DoubleArray           factorCacheDown; /* cache of down 'factor'
                                              for bm dates */
    DoubleArray           skewCache;   //cache of 'skew' for bm dates
    DoubleArray           cutOffSpeedUpCache; /* cache of cutOffSpeed for 
                                                 bm dates */
    DoubleArray           cutOffSpeedDownCache;/* cache of cutOffSpeed for
                                                  bm dates */

    /** computes 'alpha' aka cutOffSpeed. If upside is true uses 
        cutOffSpeedUp etc */
    double calcCutOffSpeed(double years, bool upside) const {
        double timeScaling = sqrt(years/cutOffSpeedEndPoint);
        double alpha;
        if (upside){
            alpha = cutOffSpeedUpStart + 
                (cutOffSpeedUpEnd - cutOffSpeedUpStart) * timeScaling;
                
        } else {
            alpha = cutOffSpeedDownStart + 
                (cutOffSpeedDownEnd - cutOffSpeedDownStart) * timeScaling;
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
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VolNormalLog, clazz);
        SUPERCLASS(CVolBaseParamSurface);
        IMPLEMENTS(IVolatilityBS);
        IMPLEMENTS(IVolatilityDVF);
        IMPLEMENTS(DeltaSurface::IShift);
        IMPLEMENTS(Calibrator::IAdjustable);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(strikeRef, "strikeRef");
        FIELD(skew1Y90100, "One year skew");
        FIELD(cutOffSpeedUpStart, "Upside cutoff speed at start");
        FIELD(cutOffSpeedDownStart, "Downside cutoff speed at start");
        FIELD(cutOffSpeedUpEnd, "Upside cutOff Speed at end");
        FIELD(cutOffSpeedDownEnd, "Downside cutoff Speed at end");
        FIELD(cutOffSpeedEndPoint, "When cutOffSpeedEnd applies")
        FIELD_MAKE_OPTIONAL(cutOffSpeedEndPoint);
        FIELD(powerStart, "to bend Skew1Y");
        FIELD(power1Y, "to bend Skew1Y");
        FIELD(powerEnd, "to bend Skew1Y");
        FIELD(volFloor, "floor when calculating vols");
        FIELD_MAKE_OPTIONAL(volFloor);
        FIELD(metric, "used to throw at the VolSurface constructor");
        FIELD_MAKE_TRANSIENT(metric);
        FIELD(baseDate, "Base Date");
        FIELD_MAKE_TRANSIENT(baseDate);
        FIELD(atmVols, "Vol at each Benchmark  Date");
        FIELD_MAKE_TRANSIENT(atmVols);
        FIELD(factorCacheDown, "");
        FIELD_MAKE_TRANSIENT(factorCacheDown);
        FIELD(skewCache, "");
        FIELD_MAKE_TRANSIENT(skewCache);
        FIELD(cutOffSpeedUpCache, "");
        FIELD_MAKE_TRANSIENT(cutOffSpeedUpCache);
        FIELD(cutOffSpeedDownCache, "");
        FIELD_MAKE_TRANSIENT(cutOffSpeedDownCache);

        // add our fields and their ranges to central list
        Calibrator::IAdjustable::registerField(
            clazz, "skew1Y90100",
            new InfiniteRange());
        Calibrator::IAdjustable::registerField(
            clazz, "cutOffSpeedUpStart",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "cutOffSpeedDownStart",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "cutOffSpeedUpEnd",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "cutOffSpeedDownEnd",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "cutOffSpeedEndPoint",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "powerStart",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "power1Y",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "powerEnd",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "volFloor",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "strikeRef",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));

    }

#define DEFAULT_ALPHA_END_POINT 30.0
    VolNormalLog(): 
    CVolBaseParamSurface(TYPE),
    strikeRef(0.0), 
    skew1Y90100(0.02),
    cutOffSpeedUpStart(0.1), 
    cutOffSpeedDownStart(2.0),
    cutOffSpeedUpEnd(2.0), 
    cutOffSpeedDownEnd(0.5),
    cutOffSpeedEndPoint(DEFAULT_ALPHA_END_POINT ), // 30 years
    powerStart(0.5), 
    power1Y(0.5), 
    powerEnd(0.5),
    volFloor(0.0){}

    static IObject* defaultCtor(){
        return new VolNormalLog();
    }

protected:
    /*** Build the parameterised vol and cache any values **/
    void buildCache() {
        try{
            const VolSurface* backbone = getBackboneSurface();
            const DateTimeArray&    dates = backbone->getDates();
            CVolProcessedBSSP       atmVolCurve = CVolProcessedBSSP(
                backbone->getProcessedVol(strikeRef));

            baseDate = backbone->getBaseDate();
            metric   = atmVolCurve->GetTimeMetric();
            atmVols  = CDoubleArray(dates.size());
            atmVolCurve->CalcVol(baseDate, dates, 
                                 CVolProcessedBS::fromFirst ,atmVols);
            factorCacheDown = CDoubleArray(dates.size());
            skewCache = CDoubleArray(dates.size());
            cutOffSpeedUpCache = CDoubleArray(dates.size());
            cutOffSpeedDownCache = CDoubleArray(dates.size());
            update();
        }
        catch(exception& e){
            throw ModelException(e, "VolNormalLog::buildCache");
        }
    }

    /** Needed for IAdjustable interface. Returns market data name for vol */
    virtual string getName() const{
        return CVolBaseParamSurface::getName();
    }

    /** Called after adjustments have been made to fields (eg calibrator) */
    virtual void fieldsUpdated(const CFieldArray& fields){
        update();
    }

    /** Called after adjustments have been made to fields */
    void update(){
        try{
            const VolSurface* backbone = getBackboneSurface();
            const DateTimeArray& dates = backbone->getDates();
            for (int i = 0; i < dates.size(); i++){
                // use trading time
                double yearToDate = metric->yearFrac(baseDate, dates[i]);
                //double yearToDate = baseDate.yearFrac(dates[i]);// calendar time
                double powerToDate = calcPowerToDate(yearToDate);
                double tToP = pow(yearToDate, powerToDate);
                skewCache[i] = skew1Y90100/tToP;
                cutOffSpeedUpCache[i] = calcCutOffSpeed(yearToDate, true);
                cutOffSpeedDownCache[i] = calcCutOffSpeed(yearToDate, false);
                factorCacheDown[i] = 
                    1.0/(N1(LOG90OVER100/cutOffSpeedDownCache[i])-0.5);
            }
        }
        catch(exception& e){
            throw ModelException(e, "VolNormalLog::update");
        }
    }

private:
    static const double LOG90OVER100;
};

const double VolNormalLog::LOG90OVER100 = log(90.0/100.0);

CClassConstSP const VolNormalLog::TYPE =
CClass::registerClassLoadMethod("VolNormalLog", typeid(VolNormalLog), load);

CClassConstSP const VolNormalLog::VNLVolParam::TYPE =
CClass::registerClassLoadMethod("VolNormalLog::VNLVolParam", 
                                typeid(VNLVolParam), load);

// external symbol to allow class to be forced to be linked in
bool VolNormalLogLinkIn(){
    return (VolNormalLog::TYPE != 0);
}

DRLIB_END_NAMESPACE
