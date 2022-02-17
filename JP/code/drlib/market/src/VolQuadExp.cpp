//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolQuadExp.cpp
//
//   Description : QuadExp param function
//
//   Date        : 02 April 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolBaseParamSurface.hpp"
#include "edginc/VolatilityDVF.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Calibrator.hpp"

DRLIB_BEGIN_NAMESPACE

class VolQuadExp: public CVolBaseParamSurface,
                  virtual public IVolatilityBS,
                  virtual public IVolatilityDVF,
                  public virtual Calibrator::IAdjustable,
                  virtual public DeltaSurface::IShift {
private:
    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class QuadExpVolParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        // constructor
        QuadExpVolParam(): CVolParam(TYPE){}

        
        /** Calculates implied vols at a single slice ie one maturity,
            many strikes (implementation of pure virtual function in
            CVolParam) */
        virtual void ComputeImpVol(const CVolBase*          vol,
                                   const CLatticeDouble&    strikes,
                                   const DateTimeArray&     maturities,
                                   CLatticeDouble&          impV) const{
            // turn the vol into what we must have
            const VolQuadExp* myVol = static_cast<const VolQuadExp *>(vol);
            // then just pass through the parameterised vol
            myVol->computeImpVol(strikes, maturities, impV);
        }

        /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
        virtual VolSurface* spotVolSurfaceFromStrikes(
            const CVolBase*       vol,
            const CDoubleArray&   strikes) const{
            // turn the vol into what we must have
            const VolQuadExp* myVol = static_cast<const VolQuadExp *>(vol);
            // then just pass through the parameterised vol
            return myVol->spotVolSurfaceFromStrikes(strikes);
        }


    private:
        static void load(CClassSP& clazz){
            REGISTER(QuadExpVolParam, clazz);
            SUPERCLASS(CVolParam);
            EMPTY_SHELL_METHOD(defaultCtor);
        }
        static IObject* defaultCtor(){
            return new QuadExpVolParam();
        }
    };
public:
    friend class QuadExpVolParam;
    static CClassConstSP const TYPE;

    void validatePop2Object(){
        // to do ...
        static const char routine[] = "VolQuadExp::validatePop2Object";
        if (false) {
            throw ModelException(routine, "cutOffSpeed must be positive");
        }
    }

    /** method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const{
        return new QuadExpVolParam();
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
    /** calculates single vol using supplied parameters. Has no concept of
        forward starting */
    double computeFittedVolStrike(double strike,
                                  double stdev,
                                  double middleQuad,
                                  double middleLinear,
                                  double lowerCutScaled,
                                  double upperCutScaled,
                                  double smileVol) const {
        static const double EXP_MINUS_ONE_HALF = 0.60653065971263300;
        static const double EXP_PLUS_ONE_HALF  = 1.0 / EXP_MINUS_ONE_HALF;
        
        double percentStrike = strike /strikeRef;
        
        /* if percentStrike < lowerCutScaled */  
        if (percentStrike < lowerCutScaled)
        {
            double lowerSigmaScaled = lowerSigmaUnscaled * stdev;
            double lowerLinearExp = (middleQuad * lowerCutScaled + 
                                     middleLinear) 
                * lowerSigmaScaled * EXP_PLUS_ONE_HALF;
            double lowerQuadExp = -middleQuad*Maths::square(lowerSigmaScaled);
            double lowerConst = (0.5 * middleQuad * lowerCutScaled + 
                                 middleLinear) * lowerCutScaled
                - (lowerQuadExp + lowerLinearExp * EXP_MINUS_ONE_HALF);
            
            double percentStrikeScaled = (percentStrike - lowerCutScaled) / 
                lowerSigmaScaled;
            
            smileVol += lowerConst + lowerQuadExp * 
                exp(-0.5 * Maths::square(percentStrikeScaled))+lowerLinearExp*
                exp(-0.5 * Maths::square(percentStrikeScaled - 1.0));
        }
        /* if lowerCutScaled <= percentStrike <= upperCutScaled */
        else if(percentStrike <= upperCutScaled)
        {
            smileVol += (0.5 * middleQuad * percentStrike + middleLinear) * 
                percentStrike;
        }
        /* if percentStrike > upperCutScaled */
        else
        {
            double upperSigmaScaled = upperSigmaUnscaled * stdev;
            double upperLinearExp = - (middleQuad * upperCutScaled +
                                       middleLinear)
                * upperSigmaScaled * EXP_PLUS_ONE_HALF;
            double upperQuadExp = -middleQuad*Maths::square(upperSigmaScaled);
            double upperConst = (0.5 * middleQuad * upperCutScaled + 
                                 middleLinear) * upperCutScaled
                - (upperQuadExp + upperLinearExp * EXP_MINUS_ONE_HALF);
            
            double percentStrikeScaled = (percentStrike - upperCutScaled) /
                upperSigmaScaled;
            
            smileVol += upperConst + upperQuadExp * 
                exp(-0.5 * Maths::square(percentStrikeScaled)) 
                + upperLinearExp *
                exp(-0.5 * Maths::square(percentStrikeScaled + 1.0));
        }
        
        return smileVol;
    }
    
    /** Calculates implied vols at a single slice ie one maturity,
        many strikes */
    void computeImpVol(const CSliceDouble&        strike,
                       const DateTime&            maturity,
                       CSliceDouble&              impV) const {
        
        static const string routine("VolQuadExp::computeImpVol");
        const VolSurface* backbone = getBackboneSurface();
        const DateTime& baseDate = backbone->getBaseDate();
        DateTime myMaturity = maturity; // default

        double yrsToMat = baseDate.yearFrac(myMaturity);
        if (Maths::isZero(yrsToMat))
        {
            double volAtStrikeRef = 
                atmVolCurve->CalcVol(baseDate, myMaturity.rollDate(1));
            for (int iStrike = 0; iStrike < strike.size(); iStrike ++){
                impV[iStrike] = volAtStrikeRef;
            }
        } else {  
            double volAtStrikeRef = atmVolCurve->CalcVol(baseDate, myMaturity);
            double stdev = volAtStrikeRef * sqrt(yrsToMat);
            double var  = Maths::square(stdev);
            
            double middleQuad   = quadUnscaled / var * volAtStrikeRef;
            double middleLinear = linearUnscaled / stdev * volAtStrikeRef - 
                middleQuad;
            
            double lowerCutScaled = exp(lowerCutUnscaled * stdev);
            double upperCutScaled = exp(upperCutUnscaled * stdev);
            
            double smileVol = volAtStrikeRef - 0.5 * middleQuad - middleLinear;
            
            for (int iStrike = 0; iStrike < strike.size(); iStrike++) { 
                impV[iStrike] = computeFittedVolStrike(strike[iStrike],
                                                       stdev,
                                                       middleQuad,
                                                       middleLinear,
                                                       lowerCutScaled,
                                                       upperCutScaled,
                                                       smileVol);
            }
        }
    }

    /** calculates implied vols across strikes and maturities */
    void computeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const {
        static const string routine("VolQuadExp::computeImpVol");
        if ((maturities.size() != strikes.size()) ||
            (maturities.size() != impV.size())) {
            throw ModelException(routine, "Size mismatch between strikes ("+ 
                                 Format::toString(strikes.size()) +
                                 "), maturities ("+ 
                                 Format::toString(maturities.size())+
                                 ") and impV ("+ 
                                 Format::toString(impV.size())+ ")");
        }  
        try{
            for (int iMat = 0; iMat < maturities.size(); iMat++) {
                computeImpVol(strikes[iMat], maturities[iMat], impV[iMat]);
            }
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    
    /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
    VolSurface* spotVolSurfaceFromStrikes(
        const CDoubleArray&   strikes) const{
        static const string routine(
            "VolQuadExp::spotVolSurfaceFromStrikes");
        try{
            const VolSurface* backbone = getBackboneSurface();
            const DateTime& baseDate = backbone->getBaseDate();
            const DateTimeArray& atmVolDates = backbone->getDates();
            if (!atmVolsGood){
                if (atmVols.size() != atmVolDates.size()){
                    atmVols = CDoubleArray(atmVolDates.size());
                }
                atmVolCurve->CalcVol(baseDate, atmVolDates, 
                                     CVolProcessedBS::fromFirst ,atmVols);
                atmVolsGood = true;
            }
            CDoubleMatrix matrix(strikes.size(), atmVols.size());
            for (int iMat = 0; iMat < atmVols.size(); iMat++) {
                double yrsToMat = baseDate.yearFrac(atmVolDates[iMat]);
                if (Maths::isZero(yrsToMat))
                {
                    for (int iStrike = 0; iStrike < strikes.size(); iStrike++){
                        matrix[iStrike][iMat] = atmVols[iMat];
                    }
                }
                else
                {
                    double stdev = atmVols[iMat] * sqrt(yrsToMat);
                    double var  = Maths::square(stdev);

                    double middleQuad   = quadUnscaled / var * atmVols[iMat];
                    double middleLinear = linearUnscaled /
                        stdev * atmVols[iMat] - middleQuad;

                    double lowerCutScaled = exp(lowerCutUnscaled * stdev);
                    double upperCutScaled = exp(upperCutUnscaled * stdev);

                    double smileVol = atmVols[iMat] - 0.5 *
                        middleQuad - middleLinear;

                    for (int iStrike = 0; iStrike < strikes.size(); iStrike++){
                        matrix[iStrike][iMat] =
                            computeFittedVolStrike(strikes[iStrike],stdev,
                                                   middleQuad,
                                                   middleLinear,lowerCutScaled,
                                                   upperCutScaled,
                                                   smileVol);
                    }
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

    // registered fields
    double      linearUnscaled;
    double      quadUnscaled;
    double      lowerSigmaUnscaled;
    double      upperSigmaUnscaled;
    double      lowerCutUnscaled;
    double      upperCutUnscaled;
    double      strikeRef;

    //dummy params for reading existing quadExp only
    int         useFixedStrikeRef;
    double      tweakStrikeUnscaled;
    double      tweakTimeUnscaled;
    double      probDensRatioMin;
    int         useTweakingForStrikeDerivs;
    int         useTweakingForTimeDerivs; 

    // transient fields (won't appear in dd interface)
    mutable bool         atmVolsGood; // do the atmVols need calculating? $unregistered
    mutable CDoubleArray atmVols;     // Vol at each bench mark date
    CVolProcessedBSSP    atmVolCurve;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VolQuadExp, clazz);
        SUPERCLASS(CVolBaseParamSurface);
        IMPLEMENTS(IVolatilityBS);
        IMPLEMENTS(IVolatilityDVF);
        IMPLEMENTS(Calibrator::IAdjustable);
        IMPLEMENTS(DeltaSurface::IShift);
        EMPTY_SHELL_METHOD(defaultCtor);

        FIELD(linearUnscaled, 
                     "Slope of vol at strikeRef. Not scaled by strike or vol");
        FIELD(quadUnscaled,
                     "Curvature of vol at strikeRef. "
                     "Not scaled by strike or vol.");        
        FIELD(lowerSigmaUnscaled, "Transition speed from parabola"
                     " to flat.");
        FIELD(upperSigmaUnscaled, "Transition speed from parabola "
                     "to flat.");        
        FIELD(lowerCutUnscaled, "Unscaled transition strike from "
                     "parabola to flat.");
        FIELD(upperCutUnscaled, "Unscaled transition strike from "
                     "parabola to flat.");        
        FIELD(strikeRef, "Reference strike.");
        
        //dummy, not used for reading existing params only
        FIELD(useFixedStrikeRef, "useFixedStrikeRef");
        FIELD_MAKE_OPTIONAL(useFixedStrikeRef);
        FIELD(tweakStrikeUnscaled, "tweakStrikeUnscaled");
        FIELD_MAKE_OPTIONAL(tweakStrikeUnscaled);
        FIELD(tweakTimeUnscaled, "tweakTimeUnscaled");
        FIELD_MAKE_OPTIONAL(tweakTimeUnscaled);
        FIELD(probDensRatioMin, "probDensRatioMin");
        FIELD_MAKE_OPTIONAL(probDensRatioMin);
        FIELD(useTweakingForStrikeDerivs, "useTweakingForStrikeDerivs");
        FIELD_MAKE_OPTIONAL(useTweakingForStrikeDerivs);
        FIELD(useTweakingForTimeDerivs, "useTweakingForTimeDerivs");
        FIELD_MAKE_OPTIONAL(useTweakingForTimeDerivs);

        // transient fields
        FIELD(atmVols, "Vol at each Benchmark  Date");
        FIELD_MAKE_TRANSIENT(atmVols);
        FIELD(atmVolCurve, "Vol at each Benchmark  Date");
        FIELD_MAKE_TRANSIENT(atmVolCurve);

        // add our fields and their ranges to central list
        Calibrator::IAdjustable::registerField(
            clazz, "linearUnscaled",
            new InfiniteRange());
        Calibrator::IAdjustable::registerField(
            clazz, "quadUnscaled",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
    }

    VolQuadExp(): 
        CVolBaseParamSurface(TYPE),
        linearUnscaled(-0.25),
        quadUnscaled(0.01),
        lowerSigmaUnscaled(1.0),
        upperSigmaUnscaled(2.5), 
        lowerCutUnscaled(-4.0),
        upperCutUnscaled(0.0),
        strikeRef(0.0),
        useFixedStrikeRef(0),
        tweakStrikeUnscaled(0.0),
        tweakTimeUnscaled(0.0),
        probDensRatioMin(0.0),
        useTweakingForStrikeDerivs(0),
        useTweakingForTimeDerivs(0), 
        atmVolsGood(false) {}

    static IObject* defaultCtor(){
        return new VolQuadExp();
    }

protected:
    /*** Build the parameterised vol and cache any values **/
    void buildCache() {
        /* avoid lengthy calculations here and avoid calculations that
           may fail if the vols are suspect (this routine is called
           for every shift - scenarios may have lots of shifts)*/
        const VolSurface* backbone = getBackboneSurface();
        atmVolCurve = CVolProcessedBSSP(backbone->getProcessedVol(strikeRef));
        atmVolsGood = false;
    }
};

CClassConstSP const VolQuadExp::TYPE =
CClass::registerClassLoadMethod("VolQuadExp", typeid(VolQuadExp), load);

CClassConstSP const VolQuadExp::QuadExpVolParam::TYPE =
CClass::registerClassLoadMethod("VolQuadExp::QuadExpVolParam", 
                                typeid(QuadExpVolParam), load);

// external symbol to allow class to be forced to be linked in
bool VolQuadExpLinkIn(){
    return (VolQuadExp::TYPE != 0);
}

DRLIB_END_NAMESPACE
