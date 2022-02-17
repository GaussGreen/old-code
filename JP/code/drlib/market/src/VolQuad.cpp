//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolQuad.cpp
//
//   Description : Parametrised Quadratic Vol Surface:
//                 The equation within the quadratic region is given by:
//                 vol(K,T) = atmVols+(Skew1Y*z)+(smile1Y*(1/2)z^2)+(cubic1Y*(1/6)zt^3)
//
//                 where z = (K-S)/S  , K = strike, S = strike reference
//                       zt = (K-S)/S*cubic
//
//                 skew, smile and the strike cuttoff defining the quadratic region 
//                 vary with T as a varying power function
//
//   Author      : Stephen Hope
//
//   Date        : 19th April 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolBaseParamSurface.hpp"
#include "edginc/Maths.hpp"
#include "edginc/VolatilityDVF.hpp"
#include "edginc/Calibrator.hpp"

DRLIB_BEGIN_NAMESPACE

const double rootE = 1.64872127070013;

class VolQuad: public CVolBaseParamSurface,
               virtual public IVolatilityBS,
               virtual public IVolatilityDVF,
               public virtual Calibrator::IAdjustable,
               virtual public DeltaSurface::IShift {
private:
    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class VQVolParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        // constructor
        VQVolParam(): CVolParam(TYPE){}

        virtual void ComputeImpVol(const CVolBase*          vol,
                                   const CLatticeDouble&    strikes,
                                   const DateTimeArray&     maturities,
                                   CLatticeDouble&          impV)const{
            // turn the vol into what we must have
            const VolQuad* myVol = static_cast<const VolQuad*>(vol);
            // then just pass through the parameterised vol
            myVol->computeImpVol(strikes, maturities, impV);
        }

        /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
        virtual VolSurface*       spotVolSurfaceFromStrikes(
            const CVolBase*       vol,
            const CDoubleArray&   strikes) const{
            // turn the vol into what we must have
            const VolQuad* myVol = static_cast<const VolQuad*>(vol);
            // then just pass through the parameterised vol
            return myVol->spotVolSurfaceFromStrikes(strikes);
        }

    private:
        static void load(CClassSP& clazz){
            REGISTER(VQVolParam, clazz);
            SUPERCLASS(CVolParam);
            EMPTY_SHELL_METHOD(defaultCtor);
        }
        static IObject* defaultCtor(){
            return new VQVolParam();
        }
    };
public:
    friend class VQVolParam;
    static CClassConstSP const TYPE;
    
    void validatePop2Object(){
        static const string method("VolQuad::validatePop2Object");
        try {
            if (Maths::isNegative(strikeRef)) {
                throw ModelException(method,
                                     "strikeRef (" + Format::toString(strikeRef) + ") must be >= 0");
            }

            if (Maths::isNegative(powerStart)) {
                throw ModelException(method,
                                     "powerStart (" + Format::toString(powerStart) + ") must be >= 0");
            }
            if (Maths::isNegative(power1Y)) {
                throw ModelException(method,
                                     "power1Y (" + Format::toString(power1Y) + ") must be >= 0");
            }
            if (Maths::isNegative(powerEnd)) {
                throw ModelException(method,
                                     "powerEnd (" + Format::toString(powerEnd) + ") must be >= 0");
            }

            // Need strikePctCutDown <= 1.0 && strikePctCutUp >= 1.0 so that
            // K- < strikeRef < K+ and so that f(strikeRef) = atmVols
            
            // strikePctCutDown <= 1.0
            if (!(  Maths::isPositive(strikePctCutDown) && 
                   !Maths::isNegative(1.0 - strikePctCutDown))) {
                throw ModelException(method,
                                     "strikePctCutDown (" + 
                                     Format::toString(strikePctCutDown) + 
                                     ") must be in (0.0, 1.0]");
            }

            // strikePctCutUp >= 1.0
            if (Maths::isNegative(strikePctCutUp - 1.0)) {
                throw ModelException(method,
                                     "strikePctCutUp (" + 
                                     Format::toString(strikePctCutUp) + 
                                     ") must be >= 1.0");
            }

            if (Maths::isNegative(smoothWidthDown)) {
                throw ModelException(method,
                                     "smoothWidthDown (" + Format::toString(smoothWidthDown) + ") must be >= 0");
            }
            if (Maths::isNegative(smoothWidthUp)) {
                throw ModelException(method,
                                     "smoothWidthUp (" + Format::toString(smoothWidthUp) + ") must be >= 0");
            }

            if (Maths::isNegative(cubic1Y)) {
                throw ModelException(method,
                                     "cubic1Y (" +
                                     Format::toString(cubic1Y) + 
                                     ") must be >= 0");
            }
            useCubic = !Maths::isZero(cubic1Y);

            if (useCubic) {
                // Ensure that strikePctCubicStart > 1.0 so that 
                // K- <= strikeRef <= K+, K3 so that f(strikeRef) = atmVols
                if (Maths::isNegative(strikePctCubicStart - 1.0)) {
                    throw ModelException(method,
                                         "strikePctCubicStart (" +
                                         Format::toString(strikePctCubicStart) + 
                                         ") must be > 1.0 if cubic is switched on");
                }

                if (Maths::isNegative(strikePctCubic1Y - strikePctCubicStart)) {
                    throw ModelException(method,
                                         "strikePctCubic1Y (" +
                                         Format::toString(strikePctCubic1Y) + 
                                         ") must be > strikePctCubicStart (" +
                                         Format::toString(strikePctCubicStart) + 
                                         ") if cubic is switched on");
                }
            }
            if (Maths::isNegative(volFloor)) {
                throw ModelException(method,
                                     "volFloor (" +
                                     Format::toString(volFloor) + 
                                     ") must be >= 0");
            }
        }
        catch (exception& e) {
            throw ModelException(e,method,"Failed for vol ("+getName()+")");
        }            
    }
    
    /** method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const{
        return new VQVolParam();
    }
    
    /** Needed for IAdjustable interface. Returns market data name for vol */
    virtual string getName() const{
        return CVolBaseParamSurface::getName();
    }

    virtual string sensName(DeltaSurface* shift) const{
        return getName();
    }

    virtual bool sensShift(DeltaSurface* shift){
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
        static const string routine("VolQuad::ComputeImpVol");
        try
        {
            if ((maturities.size() != strikes.size()) ||
                (maturities.size() != impV.size())) 
            {
                throw ModelException(routine, 
                                     "Size mismatch between strikes ("+ 
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

            for (int iMat = 0; iMat < maturities.size(); iMat++) 
            {
                if (strikes[iMat].size() != impV[iMat].size())
                {
                    throw ModelException(routine, 
                                         "Size mismatch between strikes"
                                         " & maturities for Mat " +
                                         maturities[iMat].toString() +
                                         " (n "+ Format::toString(iMat) + ")");
                }
                double vol;  // the vol we're ultimately trying to calculate

                // calculate the power function for this maturity
                double yearToDate= Maths::max(1.0/DateTime::DAYS_PER_YEAR,
                                              metric->yearFrac(baseDate,
                                                               maturities[iMat]));

                double powerToDate = calcPowerToDate(yearToDate);
                double tToP = pow(yearToDate, powerToDate);

                // define some 'power adjusted' variables
                double paSkew      = skew1Y      / tToP;
                double paCurvature = curvature1Y /tToP;
                double paCubic     = cubic1Y     /tToP;
                double paStrikePctCutDown = exp(log(strikePctCutDown)* tToP);
                double paStrikePctCutUp   = exp(log(strikePctCutUp)  * tToP);
                double paSmoothWidthDown = smoothWidthDown * tToP;
                double paSmoothWidthUp   = smoothWidthUp   * tToP;
                // note this one is linearly (not power) adjusted
                double taStrikeCubic =  strikePctCubicStart * (1-yearToDate)
                                     + (strikePctCubic1Y    *   yearToDate);

                for (int iStrike = 0; iStrike < strikes[iMat].size(); 
                     iStrike++) 
                {
                    double thisStrike = strikes[iMat][iStrike];
                    
                    vol = computeVol(thisStrike, paSkew, paCurvature,
                                     paStrikePctCutDown, paStrikePctCutUp,
                                     paSmoothWidthDown, paSmoothWidthUp,
                                     paCubic, taStrikeCubic, iMat, 
                                     backBoneVols, true);
                   
                    if (vol < volFloor){
                        vol = volFloor;
                    }

                    impV[iMat][iStrike] = vol;
                }
            }
        }
        catch (exception& e)
        {
            throw ModelException(e, routine);
        }
    }
    
    
    /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
    VolSurface* spotVolSurfaceFromStrikes(
        const CDoubleArray&   strikes) const{
        static const string routine("VolQuad::spotVolSurfaceFromStrikes");
        try
        {
            const VolSurface* backbone = getBackboneSurface();
            const DateTimeArray& dates = backbone->getDates();
            CDoubleMatrix matrix(strikes.size(),
                                 atmVols.size());
            
            for (int iMat = 0; iMat < atmVols.size(); iMat++) {

                double vol;  // the vol we're ultimately trying to calculate
            
                // calculate the power function for this maturity
                double yearToDate = metric->yearFrac(baseDate, dates[iMat]);
                double powerToDate = calcPowerToDate(yearToDate);
                double tToP = pow(yearToDate, powerToDate);
                
                // define some 'power adjusted' variables
                double paSkew      = skew1Y      / tToP;
                double paCurvature = curvature1Y /tToP;
                double paCubic     = cubic1Y     /tToP;
                double paStrikePctCutDown = exp(log(strikePctCutDown)* tToP);
                double paStrikePctCutUp   = exp(log(strikePctCutUp)  * tToP);
                double paSmoothWidthDown = smoothWidthDown * tToP;
                double paSmoothWidthUp   = smoothWidthUp   * tToP;
                // note this one is linearly (not power) adjusted
                double taStrikeCubic =  strikePctCubicStart * (1-yearToDate)
                                     + (strikePctCubic1Y    *   yearToDate);
               
                for (int iStrike = 0; iStrike < strikes.size(); iStrike ++) {
                    // calculation and variable set up from here is dependent on which 'region' this strike falls in
                    
                    double thisStrike = strikes[iStrike];

                    vol = computeVol(thisStrike, paSkew, paCurvature,
                                     paStrikePctCutDown, paStrikePctCutUp,
                                     paSmoothWidthDown, paSmoothWidthUp,
                                     paCubic, taStrikeCubic, iMat, atmVols, 
                                     false);

                    if (vol < volFloor)
                    {
                        vol = volFloor;
                    }
                    matrix[iStrike][iMat] = vol;
                }
            }
            
            /** for performance need constructor that takes in
                cached values (to do) */
            VolSurface* volSurf = new VolSurface(backbone, strikes, matrix);
                
            return volSurf;
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    double computeVol(double thisStrike, double paSkew, double paCurvature,
                      double paStrikePctCutDown, double paStrikePctCutUp,
                      double paSmoothWidthDown, double paSmoothWidthUp,
                      double paCubic, double taStrikeCubic, int iMat, 
                      const DoubleArray& backBoneVols, 
                      bool implied)const{

        static const string routine("VolQuad::computeVol");
        try
        {
            double vol;
            DoubleArray backBone;
            if (implied)
            {
                backBone = backBoneVols;
            }   
            else
            {
                backBone = atmVols;
            }

            if ((thisStrike >= strikeRef*paStrikePctCutDown) &&
                (thisStrike <= strikeRef*paStrikePctCutUp)) {
                double z = thisStrike / strikeRef - 1;
                vol = backBone[iMat] + (paSkew * z) + (paCurvature * z * z / 2);
                if (useCubic) {
                    double zt = Maths::max(thisStrike/(strikeRef*taStrikeCubic)-1, 0.0);
                    vol += (paCubic * zt * zt * zt / 6);
                }
            } else {
                double sign;
                double zcut;
                double ztcut;
                double swidth;
                double scut;
                // Above upper cut off
                if (thisStrike > strikeRef * paStrikePctCutUp)
                {
                    sign = 1;
                    scut =  strikeRef * paStrikePctCutUp;
                    zcut = (strikeRef * paStrikePctCutUp) / strikeRef - 1;
                    if (useCubic) {
                        ztcut = Maths::max((paStrikePctCutUp / taStrikeCubic)-1, 0.0);
                    }
                    swidth = paSmoothWidthUp;
                }
                // Below lower cut off
                else
                {
                    sign = -1;
                    scut =  strikeRef * paStrikePctCutDown;
                    zcut = (strikeRef * paStrikePctCutDown) / strikeRef - 1;
                    if (useCubic) {
                        ztcut = Maths::max((paStrikePctCutDown / taStrikeCubic)-1, 0.0);
                    }
                    swidth = paSmoothWidthDown;
                }

                double refOverCub = 1 / taStrikeCubic;
                double c = -sign*swidth*rootE*(paSkew + (paCurvature * zcut)+
                                              (useCubic?(paCubic * refOverCub *ztcut*ztcut/2):0.0));
                double b = -swidth*swidth*(paCurvature + 
                                          (useCubic?(paCubic*refOverCub*refOverCub*ztcut):0.0));
                double a = (backBoneVols[iMat]+(paSkew*zcut)+ (paCurvature*zcut*zcut/2)+
                            (useCubic?(paCubic*ztcut*ztcut*ztcut/6):0.0))
                    -(b+c/rootE);

                if (!Maths::isZero(swidth)) {
                    double w = (thisStrike-scut)/(swidth*strikeRef);
                    
                    vol = a+b*exp(-0.5*w*w)+c*exp(-0.5*(w+sign)*(w+sign));
                }
                else {
                    // if smoothing region is zero sized, then avoid exp(-infinity)
                    vol = a;
                }
            }
            return vol;
        }
        catch (exception& e)
        {
            throw ModelException(e, routine);
        }
    }    
    
    // registered fields
    double        strikeRef;
    double        skew1Y;                 /* one year skew */
    double        curvature1Y;
    double        cubic1Y;
    double        powerStart;             /* to bend Skew1Y */
    double        power1Y;                /* to bend Skew1Y */
    double        powerEnd;               /* to bend Skew1Y */
    double        strikePctCubic1Y;
    double        strikePctCubicStart;
    double        strikePctCutDown;
    double        strikePctCutUp;
    double        smoothWidthDown;
    double        smoothWidthUp;
    double        volFloor;          /* optional floor when calculating vols */

    // transient fields (won't appear in dd interface)
    DateTime            baseDate;
    TimeMetricConstSP   metric;
    CDoubleArray        atmVols;     // Vol at each Bench Date
    bool                useCubic;


    //// computes 'power' factor  Move this to common area ?
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
        REGISTER(VolQuad, clazz);
        SUPERCLASS(CVolBaseParamSurface);
        IMPLEMENTS(IVolatilityBS);
        IMPLEMENTS(IVolatilityDVF);
        IMPLEMENTS(Calibrator::IAdjustable);
        IMPLEMENTS(DeltaSurface::IShift);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(strikeRef, "strikeRef");
        FIELD(skew1Y, "One year skew");
        FIELD(curvature1Y, "One year curvature");
        FIELD(cubic1Y, "One year cubic");
        FIELD(powerStart, "to bend Skew1Y");
        FIELD(power1Y, "to bend Skew1Y");
        FIELD(powerEnd, "to bend Skew1Y");
        FIELD(strikePctCubic1Y, "strike percent cubic 1Y");
        FIELD(strikePctCubicStart, "strike percent cubic start");
        FIELD(strikePctCutDown, "strike percent cut down");
        FIELD(strikePctCutUp, "strike percent cut up");
        FIELD(smoothWidthDown, "smooth width down");
        FIELD(smoothWidthUp, "smooth width up");
        FIELD(volFloor, "floor when calculating vols");
        FIELD_MAKE_OPTIONAL(volFloor);
        FIELD(baseDate, "Base Date");
        FIELD_MAKE_TRANSIENT(baseDate);
        FIELD(metric, "used to throw at the VolSurface constructor");
        FIELD_MAKE_TRANSIENT(metric);
        FIELD(atmVols, "Vol at each Benchmark  Date");
        FIELD_MAKE_TRANSIENT(atmVols);
        FIELD(useCubic, "Internal");
        FIELD_MAKE_TRANSIENT(useCubic);

        // add our fields and their ranges to central list
        Calibrator::IAdjustable::registerField(
            clazz, "skew1Y",
            new Range(Infinity(Infinity::Minus), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "curvature1Y",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "cubic1Y",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "powerStart",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "power1Y",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "powerEnd",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "strikePctCubicStart",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "strikePctCubic1Y",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "strikePctCutDown",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "strikePctCutUp",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "smoothWidthDown",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "smoothWidthUp",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "volFloor",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
    }

    VolQuad(): CVolBaseParamSurface(TYPE),
               strikeRef(0.0), skew1Y(0.0), curvature1Y(0.0),
               cubic1Y(0.0), powerStart(0.0), power1Y(0.0), powerEnd(0.0),
               strikePctCubic1Y(0.0), strikePctCubicStart(0.0),
               strikePctCutDown(0.0), strikePctCutUp(0.0), smoothWidthDown(0.0),
               smoothWidthUp(0.0), volFloor(0.0), useCubic(true)
        {
            // empty
        }
    
    static IObject* defaultCtor(){
        return new VolQuad();
    }

protected:
    /*** Build the parameterised vol and cache any values **/
    void buildCache() {
        const VolSurface* backbone = getBackboneSurface();
        const DateTimeArray&    dates = backbone->getDates();
        CVolProcessedBSSP       atmVolCurve = CVolProcessedBSSP(
            backbone->getProcessedVol(strikeRef));

        baseDate = backbone->getBaseDate();
        metric = atmVolCurve->GetTimeMetric();
        atmVols   = CDoubleArray(dates.size());
        atmVolCurve->CalcVol(baseDate, dates, 
                             CVolProcessedBS::fromFirst ,atmVols);
#if 0
        // NB This code was the old caching code - hashed out when the
        // term structure to the power was added. So its reintroduction will
        // need careful analysis.
        scaleCutOff = log(90.0 / 100.0) / cutOffSpeedStart;
        scaleCutOff = N1(scaleCutOff) - 0.5;
        scaleCutOff = 1.0 / scaleCutOff;

        atmSlopes = CDoubleArray(dates.size());

        for (int i = 0; i <atmSlopes.size(); i++) {
            double yearToDate = baseDate.yearFrac(dates[i]); // calendar time
            double powerToDate;
            if (yearToDate <= 1.0) {
                powerToDate  = powerStart * (1.0 - yearToDate);
                powerToDate += power1Y * yearToDate;
            } else {
                powerToDate  = power1Y  / yearToDate;
                powerToDate += powerEnd * (1.0 - 1.0/yearToDate);
            }
            double tToP = pow(yearToDate, powerToDate);
            double thisScaleCutOff = scaleCutOff
            atmSlopes[i]  = scaleCutOff * skew1Y90100;
            atmSlopes[i] /= tToP;
        }
#endif
    }
};

CClassConstSP const VolQuad::TYPE =
CClass::registerClassLoadMethod("VolQuad", typeid(VolQuad), load);

CClassConstSP const VolQuad::VQVolParam::TYPE =
CClass::registerClassLoadMethod("VolQuad::VQVolParam", 
                                typeid(VQVolParam), load);

// external symbol to allow class to be forced to be linked in
bool VolQuadLinkIn(){
    return (VolQuad::TYPE != 0);
}

DRLIB_END_NAMESPACE
