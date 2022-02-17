//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolNormalLogMod.cpp
//
//   Description : Modified VolNormalLog following NY testing.
//                 Began with version 1.22 of VolNormalLog.
//                 VolNormalLogMod is "free standing" in that it completely specifies vol surface without the 
//                 need to interpolate backbone from a surface.
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : April 16, 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolBaseParamSurface.hpp"
#include "edginc/VolatilityDVF.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include "edginc/NonPricingModel.hpp"

DRLIB_BEGIN_NAMESPACE

/** A parameterised CVolBase that supports BS and DVF views of the world */
class VolNormalLogMod: public CVolBaseParamSurface,
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
            const VolNormalLogMod* myVol = 
                static_cast<const VolNormalLogMod *>(vol);
            // then just pass through the parameterised vol
            myVol->computeImpVol(strikes, maturities, impV);
        }

        /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
        virtual VolSurface* spotVolSurfaceFromStrikes(
            const CVolBase*       vol,
            const CDoubleArray&   strikes) const{
            // turn the vol into what we must have
            const VolNormalLogMod* myVol = 
                static_cast<const VolNormalLogMod *>(vol);
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

    ~VolNormalLogMod(){}

    void validatePop2Object(){
        const string method("VolNormalLogMod::validatePop2Object");

        try {
            Calibrator::IAdjustable::checkRange(this);

            int t     = DateTime::timeConvert(timeOfDay);
            int numBM = expiryStr.size();
            expiries = ExpiryArraySP(new ExpiryArray(numBM));
            for (int j=0; j<numBM; j++)
                (*expiries)[j] = ExpirySP(new MaturityTimePeriod(expiryStr[j], t));

            if (numBM != strikeRefArr.size())
                throw ModelException(method, "strike ref expiries and strike ref arrayy are different in size.");

            if (numBM != volArr.size())
                throw ModelException(method, "strike ref expiries and vol array are different in size.");

            for (int i = 0; i < strikeRefArr.size(); i++)
            {
                if (!Maths::isPositive(strikeRefArr[i]))
                {
                    throw ModelException(method, "strike refs must be in (0.00, +oo).");
                }
            }

        }
        catch(exception& e){
            throw ModelException(e, method, "Failed for vol ("+getName()+")");
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

            for (int i=0; i<strikeRefArr.size(); i++)
                strikeRefArr[i] *= newSpot/spot;

            // get new backbone
            buildCache();
        }
        // DON'T tweak the backbone !!!
        return false; // no more tweaking required here
    }

    // this is slow and lazy way
    // to do: review this !!!
   static double interpStrikeRef(const DateTime& date,
                                  const DateTimeArray& dates, 
                                  const  DoubleArray& strikes){
        Schedule x(dates, strikes, Schedule::INTERP_LINEAR);
        DateTime dateToInterp = date;

        // make sure dateToInterp is bracketed by dates 
        if (dates[0].isGreater(dateToInterp))
        {
            dateToInterp = dates[0];
        }
        if (dates[dates.size()-1].isLess(dateToInterp))
        {
            dateToInterp = dates[dates.size()-1];
        }

        return x.interpolate(dateToInterp);
    }

private:

    void computeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const {
        static const string routine("VolNormalLogMod::computeImpVol");
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

        const DateTimeArray& dates = backbone->getDates();

        CVolProcessedBSSP  refVolCurve(backbone->getProcessedVol(strikeRefArr[0])); // strike does not matter as we have only one column
        refVolCurve->CalcVol(baseDate, maturities, CVolProcessedBS::fromFirst,
                             backBoneVols);
        DoubleArray backBoneVars(maturities.size());
        refVolCurve->CalcVar(baseDate, maturities, CVolProcessedBS::fromFirst,
                             backBoneVars);
        
        for (int iMat = 0; iMat < maturities.size(); iMat++) {
            if (strikes[iMat].size() != impV[iMat].size()){
                throw ModelException(routine, "Size mismatch between strikes"
                                     " & maturities for Mat " +
                                     maturities[iMat].toString() +
                                     " (n "+ Format::toString(iMat) + ")");
            }
            // floor yearToDate with 1D to avoid blowing up as t --> 0
            double yearToDate= timeMetric->yearFrac(baseDate, maturities[iMat]);
            yearToDate= Maths::max(1.0/DateTime::DAYS_PER_YEAR, yearToDate);
            //double yearToDate= Maths::max(1.0/DateTime::DAYS_PER_YEAR,
              //                            baseDate.yearFrac(maturities[iMat]));
            double powerToDate = calcPowerToDate(yearToDate);
            double tToP = pow(yearToDate, powerToDate);
            double skew = skew1Y90100/tToP;
            double cutOffSpeedUp = calcCutOffSpeed(yearToDate, true, backBoneVars[iMat]);
            double cutOffSpeedDown = calcCutOffSpeed(yearToDate, false, backBoneVars[iMat]);
            double factorDown = 1.0/(N1(log(90.0/100.0)/cutOffSpeedDown)-0.5);
            double factorUp = factorDown*cutOffSpeedUp/cutOffSpeedDown;
            for (int iStrike = 0; iStrike < strikes[iMat].size(); iStrike ++) {
                double strike = strikes[iMat][iStrike];
                double strikeRef = interpStrikeRef(maturities[iMat], dates, strikeRefArr);
                double logStrike = log(strike/strikeRef);
                bool   upside = strike >= strikeRef;
                logStrike /= upside? cutOffSpeedUp: cutOffSpeedDown;
				double distance = logStrike;
				double quadratic = N1(distance*distance) - 0.5;
                logStrike  = N1(logStrike) - 0.5;
                // NB don't multiply by cutOffSpeed since we skipped
                // it from the factor calculation
                // to do: should we floor the vol here? 
                impV[iMat][iStrike]  = (upside? factorUp: factorDown) * 
                    skew * logStrike;
                impV[iMat][iStrike] += backBoneVols[iMat];
				impV[iMat][iStrike] += quadFactor*quadratic*exp(-distance*distance);
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
            "VolNormalLogMod::spotVolSurfaceFromStrikes");
        try{
            const VolSurface* backbone = getBackboneSurface();
            CDoubleMatrix matrix(strikes.size(), atmVols.size());

            if (atmVols.size() != strikeRefArr.size())
                throw ModelException(routine, "atm vol and strike ref array have different size!");

            // loop over strikes
            for (int iStrike = 0; iStrike < strikes.size(); iStrike ++) {
                // loop over maturities
                for (int iMat = 0; iMat < atmVols.size(); iMat++) {
                    bool upside = strikes[iStrike] >= strikeRefArr[iMat];
                    double thisLogStrike = log(strikes[iStrike]/strikeRefArr[iMat])/
                        (upside? cutOffSpeedUpCache[iMat]: 
                         cutOffSpeedDownCache[iMat]);
					double distance = thisLogStrike;
					double quadratic = N1(distance*distance) - 0.5;
                    thisLogStrike  = N1(thisLogStrike) - 0.5;
                    // cutOffSpeed is taken care of by 'factor'
                    // we adjust factorDown to ensure smoothness of derivative
                    // wrt to strike
                    double factor = upside?
                        factorCacheDown[iMat]*cutOffSpeedUpCache[iMat]/
                         cutOffSpeedDownCache[iMat]:
                        factorCacheDown[iMat];
                    double vol = factor * skewCache[iMat] * 
                        thisLogStrike + atmVols[iMat] + quadFactor*quadratic*exp(-distance*distance);
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
    StringArray   expiryStr; // param vol can only take string array for now
    string        timeOfDay; // EOD or SOD
    ExpiryArraySP expiries; // transient
    DoubleArray   strikeRefArr;
    DoubleArray   volArr;
    double        skew1Y90100;       /* one year skew */
    double        cutOffStDevUpStart;
    double        cutOffStDevUpMid;
    double        cutOffStDevUpEnd;
    double        cutOffUpMin;
    double        cutOffUpMax;
    double        cutOffStDevDownStart;
    double        cutOffStDevDownMid;
    double        cutOffStDevDownEnd;
    double        cutOffDownMin;
    double        cutOffDownMax;
    double        cutOffSpeedMidPoint;
    double        cutOffSpeedEndPoint;
	double        quadFactor;

    double        powerStart;        /* to bend Skew1Y */
    double        powerMid;           /* to bend Skew1Y */
    double        powerEnd;          /* to bend Skew1Y */
    double        volFloor;          /* optional floor when calculating vols */
    // transient fields (won't appear in dd interface)
    DateTime              baseDate;
    TimeMetricConstSP     timeMetric;
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
    double calcCutOffSpeed(double years, bool upside, double var) const {
        double cutStDev;
        double cutSpeed;
        if (upside){ 
            if (years < cutOffSpeedMidPoint) {
                cutStDev = cutOffStDevUpStart + 
                    years/cutOffSpeedMidPoint*(cutOffStDevUpMid - cutOffStDevUpStart);
            } else if (years >= cutOffSpeedMidPoint && years < cutOffSpeedEndPoint) {
                 cutStDev = cutOffStDevUpMid + (years - cutOffSpeedMidPoint)/(cutOffSpeedEndPoint - cutOffSpeedMidPoint)
                     *(cutOffStDevUpEnd - cutOffStDevUpMid);
            } else {
                cutStDev = cutOffStDevUpEnd;
            }

            cutSpeed = cutStDev*sqrt(var);

            cutSpeed = Maths::max(cutSpeed, cutOffUpMin);
            cutSpeed = Maths::min(cutSpeed, cutOffUpMax);
            
        } else { 
           if (years < cutOffSpeedMidPoint) {
                cutStDev = cutOffStDevDownStart + 
                    years/cutOffSpeedMidPoint*(cutOffStDevDownMid - cutOffStDevDownStart);
            } else if (years >= cutOffSpeedMidPoint && years < cutOffSpeedEndPoint) {
                 cutStDev = cutOffStDevDownMid + (years - cutOffSpeedMidPoint)/(cutOffSpeedEndPoint - cutOffSpeedMidPoint)
                     *(cutOffStDevDownEnd - cutOffStDevDownMid);
            } else {
                cutStDev = cutOffStDevDownEnd;
            }

            cutSpeed = cutStDev*sqrt(var);

            cutSpeed = Maths::max(cutSpeed, cutOffDownMin);
            cutSpeed = Maths::min(cutSpeed, cutOffDownMax);
        }    
        return cutSpeed;
    }

    //// computes 'power' factor
    double calcPowerToDate(double yearToDate) const{
        double powerToDate;
        if (yearToDate <= cutOffSpeedMidPoint) {
        
            powerToDate  = powerStart + (powerMid - powerStart)* yearToDate/cutOffSpeedMidPoint;
        } else {
            powerToDate  = powerEnd + (powerMid - powerEnd) * cutOffSpeedMidPoint / yearToDate;
        }
        return powerToDate;
    }
        
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VolNormalLogMod, clazz);
        SUPERCLASS(CVolBaseParamSurface);
        IMPLEMENTS(IVolatilityBS);
        IMPLEMENTS(IVolatilityDVF);
        IMPLEMENTS(DeltaSurface::IShift);
        IMPLEMENTS(Calibrator::IAdjustable);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(expiryStr, "expiry string array");
        FIELD(timeOfDay, " time of day for vol");
        FIELD(expiries, "expiries");
        FIELD_MAKE_TRANSIENT(expiries);
        FIELD(strikeRefArr, "strikeRefArr");
        FIELD(volArr, "volArr");
        FIELD(skew1Y90100, "One year skew");
        FIELD(cutOffStDevUpStart, "asd");
        FIELD(cutOffStDevUpMid, "asd");
        FIELD(cutOffStDevUpEnd, "asd");
        FIELD(cutOffUpMin, "asd");
        FIELD(cutOffUpMax, "asd");
        FIELD(cutOffStDevDownStart, "asd");
        FIELD(cutOffStDevDownMid, "asd");
        FIELD(cutOffStDevDownEnd, "asd");
        FIELD(cutOffDownMin, "asd");
        FIELD(cutOffDownMax, "asd");
        FIELD(cutOffSpeedMidPoint, "When cutOffSpeedMid applies")
        FIELD(cutOffSpeedEndPoint, "When cutOffSpeedEnd applies")
        FIELD_MAKE_OPTIONAL(cutOffSpeedEndPoint);
		FIELD(quadFactor, "Extra Quadratic term");
        FIELD(powerStart, "to bend Skew1Y");
        FIELD(powerMid, "to bend Skew1Y");
        FIELD(powerEnd, "to bend Skew1Y");
        FIELD(volFloor, "floor when calculating vols");
        FIELD_MAKE_OPTIONAL(volFloor);
        FIELD(timeMetric, "used to throw at the VolSurface constructor");
        FIELD_MAKE_TRANSIENT(timeMetric);
        FIELD(baseDate, "Base Date");
        FIELD_MAKE_TRANSIENT(baseDate);
        FIELD(factorCacheDown, "");
        FIELD_MAKE_TRANSIENT(factorCacheDown);
        FIELD(skewCache, "");
        FIELD_MAKE_TRANSIENT(skewCache);
        FIELD(cutOffSpeedUpCache, "");
        FIELD_MAKE_TRANSIENT(cutOffSpeedUpCache);
        FIELD(cutOffSpeedDownCache, "");
        FIELD_MAKE_TRANSIENT(cutOffSpeedDownCache);
        FIELD(atmVols, "");
        FIELD_MAKE_TRANSIENT(atmVols);

        // add our fields and their ranges to central list
        Calibrator::IAdjustable::registerField(
            clazz, "volArr",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "skew1Y90100",
            new InfiniteRange());
        Calibrator::IAdjustable::registerField(
            clazz, "cutOffStDevUpStart",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
         Calibrator::IAdjustable::registerField(
            clazz, "cutOffStDevUpMid",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "cutOffStDevUpEnd",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "cutOffUpMin",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "cutOffUpMax",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "cutOffStDevDownStart",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "cutOffStDevDownMid",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "cutOffStDevDownEnd",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "cutOffDownMin",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "cutOffDownMax",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "cutOffSpeedMidPoint",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "cutOffSpeedEndPoint",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
		Calibrator::IAdjustable::registerField(
            clazz, "quadFactor",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "powerStart",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "powerMid",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "powerEnd",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "volFloor",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
    }

#define DEFAULT_ALPHA_END_POINT 30.0
    VolNormalLogMod(): 
    CVolBaseParamSurface(TYPE),
    skew1Y90100(0.02),

    cutOffStDevUpStart(3.),
    cutOffStDevUpMid(2.),
    cutOffStDevUpEnd(1.),
    cutOffUpMin(0.001),
    cutOffUpMax(10.),
    cutOffStDevDownStart(3.),
    cutOffStDevDownMid(2.),
    cutOffStDevDownEnd(3.),
    cutOffDownMin(0.001),
    cutOffDownMax(10.),

    cutOffSpeedMidPoint(5. ),
    cutOffSpeedEndPoint(DEFAULT_ALPHA_END_POINT ), // 30 years
	quadFactor(0.1),
    powerStart(0.5), 
    powerMid(0.5), 
    powerEnd(0.5),
    volFloor(0.0){}

    static IObject* defaultCtor(){
        return new VolNormalLogMod();
    }

protected:
    // overriden because backbone is different
    void getMarket(const IModel*     model, 
                   const MarketData* market){
        try{

            if (!volSurfaceForBackbone){
                NonPricingModel nonPricingModel;
                VolSurfaceSP surface(VolSurfaceSP::dynamicCast(
    	            market->GetData(&nonPricingModel, getName(),VolSurface::TYPE))); // vol surface is just for time metric etc. vol values not used
            
                DoubleArray volStrike(1, strikeRefArr[0]);// only one column, strike does not matter
                volSurfaceForBackbone = VolSurfaceSP(new VolSurface(getName(),
                                                                    surface->getTimeMetric().get(),
                                                                    volStrike,
                                                                    CDoubleMatrix(volArr),
                                                                    expiries.get(),
                                                                    surface->getBaseDate()));
            }
        } catch (exception& e){
            throw ModelException(e, "VolNormalLogMod::getMarket", "failed creating refVol backbone.");
        }
        // call base method
        CVolBaseParamSurface::getMarket(model, market);
    }

    /*** Build the parameterised vol and cache any values **/
    void buildCache() {
        try{

        // uncomment the next line if you get an access error when the imsl integration routine fails. 
        //   imsl_error_options(IMSL_ERROR_MSG_PATH, "C:\\", 0);

            const VolSurface* backbone = getBackboneSurface();
            const DateTimeArray&    dates = backbone->getDates();
            CVolProcessedBSSP       refVolCurve = CVolProcessedBSSP(
                backbone->getProcessedVol(strikeRefArr[0])); // strike does not matter as we have one vector

            baseDate = backbone->getBaseDate();
            timeMetric = refVolCurve->GetTimeMetric();
            atmVols   = CDoubleArray(dates.size());
            refVolCurve->CalcVol(baseDate, dates, 
                                 CVolProcessedBS::fromFirst ,atmVols);
            factorCacheDown = CDoubleArray(dates.size());
            skewCache = CDoubleArray(dates.size());
            cutOffSpeedUpCache = CDoubleArray(dates.size());
            cutOffSpeedDownCache = CDoubleArray(dates.size());
            update();
        }
        catch(exception& e){
            throw ModelException(e, "VolNormalLogMod::buildCache");
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
			// const VolSurface* backbone = getBackboneSurface();
            const VolSurface* prevBackbone = getBackboneSurface();
			DoubleArray oneStrike(1, strikeRefArr[0]);
			VolSurface backbone(prevBackbone,oneStrike,CDoubleMatrix(volArr));
            const DateTimeArray& dates = backbone.getDates();

            CVolProcessedBSSP  refVolCurve(backbone.getProcessedVol(strikeRefArr[0])); // only one vector
            DoubleArray backBoneVars(dates.size());
            refVolCurve->CalcVar(baseDate, dates, CVolProcessedBS::fromFirst,
                                 backBoneVars);
        
            for (int i = 0; i < dates.size(); i++){
                double yearToDate = timeMetric->yearFrac(baseDate, dates[i]);// calendar time
                //double yearToDate = baseDate.yearFrac(dates[i]);// calendar time
                double powerToDate = calcPowerToDate(yearToDate);
                double tToP = pow(yearToDate, powerToDate);
                skewCache[i] = skew1Y90100/tToP;
                cutOffSpeedUpCache[i] = calcCutOffSpeed(yearToDate, true, backBoneVars[i]);
                cutOffSpeedDownCache[i] = calcCutOffSpeed(yearToDate, false, backBoneVars[i]);
                factorCacheDown[i] = 
                    1.0/(N1(LOG90OVER100/cutOffSpeedDownCache[i])-0.5);
            }
        }
        catch(exception& e){
            throw ModelException(e, "VolNormalLogMod::update");
        }
    }

private:
    static const double LOG90OVER100;
};

const double VolNormalLogMod::LOG90OVER100 = log(90.0/100.0);

CClassConstSP const VolNormalLogMod::TYPE =
CClass::registerClassLoadMethod("VolNormalLogMod", typeid(VolNormalLogMod), load);

CClassConstSP const VolNormalLogMod::VNLVolParam::TYPE =
CClass::registerClassLoadMethod("VolNormalLogMod::VNLVolParam", 
                                typeid(VNLVolParam), load);

// external symbol to allow class to be forced to be linked in
bool VolNormalLogModLinkIn(){
    return (VolNormalLogMod::TYPE != 0);
}

DRLIB_END_NAMESPACE
