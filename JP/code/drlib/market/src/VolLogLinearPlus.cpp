//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolLogLinearPlus.cpp
//
//   Description : Yet Another Prototype Parameterized Volatility
//
//   Author      : Regis Guichard
//
//   Date        : 30 Oct 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolBaseParamSurface.hpp"
#include "edginc/VolatilityDVF.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Spline.hpp"
#include "edginc/LinearInterpolator.hpp"

DRLIB_BEGIN_NAMESPACE

class VolLogLinearPlusDebugInfoAddin;
/** A parameterised CVolBase that supports BS and DVF views of the world */
class VolLogLinearPlus: public CVolBaseParamSurface,
                        virtual public IVolatilityBS,
                        virtual public IVolatilityDVF,
                        virtual public DeltaSurface::IShift,
                        public virtual Calibrator::IAdjustable{
private:
    friend class VolLogLinearPlusDebugInfoAddin;

    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class VolLogLinearPlusVolParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        // constructor
        VolLogLinearPlusVolParam(): CVolParam(TYPE){}

        ~VolLogLinearPlusVolParam(){}

        virtual void ComputeImpVol(const CVolBase*          vol,
                                   const CLatticeDouble&    strikes,
                                   const DateTimeArray&     maturities,
                                   CLatticeDouble&          impV) const{
            // turn the vol into what we must have
            const VolLogLinearPlus* myVol = 
                static_cast<const VolLogLinearPlus *>(vol);
            // then just pass through the parameterised vol
            myVol->computeImpVol(strikes, maturities, impV);
        }

        /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
        virtual VolSurface* spotVolSurfaceFromStrikes(
            const CVolBase*       vol,
            const CDoubleArray&   strikes) const{
            // turn the vol into what we must have
            const VolLogLinearPlus* myVol = 
                static_cast<const VolLogLinearPlus *>(vol);
            // then just pass through the parameterised vol
            return myVol->spotVolSurfaceFromStrikes(strikes);
        }

    private:
        static void load(CClassSP& clazz){
            REGISTER(VolLogLinearPlusVolParam, clazz);
            SUPERCLASS(CVolParam);
            EMPTY_SHELL_METHOD(defaultCtor);
        }
        static IObject* defaultCtor(){
            return new VolLogLinearPlusVolParam();
        }
    };
public:
    friend class VolLogLinearPlusVolParam;
    static CClassConstSP const TYPE;

    ~VolLogLinearPlus(){}

    void validatePop2Object(){
        static const string method("VolLogLinearPlus::validatePop2Object()");
        try{
            Maths::checkPositive(strikeRef, "strikeRef");
            Maths::checkPositive(spreadPowerShortTerm, "spreadPowerShortTerm");
            Maths::checkPositive(spreadPowerLongTerm, "spreadPowerLongTerm");
            if (!Maths::isPositive(tailDeltaDown)
                || !Maths::isPositive(midDeltaDown - tailDeltaDown)
                || !Maths::isPositive(0.5 - midDeltaDown)){
                throw ModelException(method,
                                     "tailDeltaDown, midDeltaDown should be such that "
                                     "0.0 < tailDeltaDown < midDeltaDown < 0.5;"
                                     "\ngot "
                                     + Format::toString(tailDeltaDown)
                                     + " and "
                                     + Format::toString(midDeltaDown)
                                     + " respectively");
            }
            if (!Maths::isPositive(tailDeltaUp)
                || !Maths::isPositive(midDeltaUp - tailDeltaUp)
                || !Maths::isPositive(0.5 - midDeltaUp)){
                throw ModelException(method,
                                     "tailDeltaUp, midDeltaUp should be such that "
                                     "0.0 < tailDeltaUp < midDeltaUp < 0.5;"
                                     "\ngot "
                                     + Format::toString(tailDeltaUp)
                                     + " and "
                                     + Format::toString(midDeltaUp)
                                     + " respectively");
            }
            int nbSkewFactors = skewFactors.size();
            if (nbSkewFactors != skewFactorBMs.size()){
                throw ModelException(method,
                                     "the number of skew factors ("
                                     + Format::toString(nbSkewFactors)
                                     + ") should be equal to the number of skew factor benchmarks ("
                                     + Format::toString(skewFactorBMs.size())
                                     + ")");
            }
            skewFactorDates.resize(nbSkewFactors);
            skewFactorTradYears.resize(nbSkewFactors);
            // spreads
            int nbSpreadBMs = spreadBMs.size();
            if (nbSpreadBMs == 0){
                throw ModelException(method,
                                     "the number of spread benchmarks ("
                                     + Format::toString(nbSpreadBMs)
                                     + ") should be greater than zero");
            }
            if (nbSpreadBMs != tailSpreadsDown.size()){
                throw ModelException(method,
                                     "the number of spread benchmarks ("
                                     + Format::toString(nbSpreadBMs)
                                     + ") should be equal to the number of lower tail volatility spreads ("
                                     + Format::toString(tailSpreadsDown.size())
                                     + ")");
            }
            if (nbSpreadBMs != midSpreadsUp.size()){
                throw ModelException(method,
                                     "the number of spread benchmarks ("
                                     + Format::toString(nbSpreadBMs)
                                     + ") should be equal to the number of upper middle volatility spreads ("
                                     + Format::toString(midSpreadsUp.size())
                                     + ")");
            }
            if (nbSpreadBMs != tailSpreadsUp.size()){
                throw ModelException(method,
                                     "the number of spread benchmarks ("
                                     + Format::toString(nbSpreadBMs)
                                     + ") should be equal to the number of upper tail volatility spreads ("
                                     + Format::toString(tailSpreadsUp.size())
                                     + ")");
            }
            midSpreadsDown.resize(nbSpreadBMs);
            for (int iSpreadBM = 0; iSpreadBM < nbSpreadBMs; ++iSpreadBM){
                midSpreadsDown[iSpreadBM] = 0.0;
            }
            spreadDates.resize(nbSpreadBMs);
            spreadTradYears.resize(nbSpreadBMs);
            // tail speeds
            int nbTailSpeedBMs = tailSpeedBMs.size();
            if (nbTailSpeedBMs != tailSpeedsDown.size()){
                throw ModelException(method,
                                     "the number of tail speed benchmarks ("
                                     + Format::toString(nbTailSpeedBMs)
                                     + ") should be equal to the number of lower tail tail speeds ("
                                     + Format::toString(tailSpeedsDown.size())
                                     + ")");
            }
            if (nbTailSpeedBMs != tailSpeedsUp.size()){
                throw ModelException(method,
                                     "the number of tail speed benchmarks ("
                                     + Format::toString(nbTailSpeedBMs)
                                     + ") should be equal to the number of lower tail tail speeds ("
                                     + Format::toString(tailSpeedsUp.size())
                                     + ")");
            }
            tailSpeedDates.resize(nbTailSpeedBMs);
            tailSpeedTradYears.resize(nbTailSpeedBMs);
            // check speeds are positive
            Maths::checkPositive(tailSpeedsDown, "tailSpeedsDown");
            Maths::checkPositive(tailSpeedsUp, "tailSpeedsUp");
            Calibrator::IAdjustable::checkRange(this);
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    /** method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const{
        return new VolLogLinearPlusVolParam();
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

    // copied and pasted from VolSpline
    void computeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const {
        static const string routine("VolLogLinearPlus::computeImpVol");

        if ((maturities.size() != strikes.size()) ||
            (maturities.size() != impV.size())) {
            throw ModelException(routine, "Size mismatch between strikes ("+ 
                                 Format::toString(strikes.size()) +
                                 "), maturities ("+ 
                                 Format::toString(maturities.size())+
                                 ") and impV ("+ 
                                 Format::toString(impV.size())+ ")");
        }  
    
        int nbBackboneDates = backboneDates.size();
        int iBackboneDateStart = 0;
        for (int iMat = 0; iMat < maturities.size(); iMat++) {
            int nbStrikes = strikes[iMat].size();
            if (nbStrikes != impV[iMat].size()){
                throw ModelException(routine, "Size mismatch between strikes"
                                     " & maturities for Mat " +
                                     maturities[iMat].toString() +
                                     " (n "+ Format::toString(iMat) + ")");
            }
            int iBackboneDate = iBackboneDateStart;
            while (iBackboneDate < nbBackboneDates && backboneDates[iBackboneDate] < maturities[iMat]) {
                ++iBackboneDate;
            }
            iBackboneDateStart = Maths::max(0, iBackboneDate - 1);
            int iBackboneDateEnd = Maths::min(iBackboneDate, nbBackboneDates - 1);
            CDoubleMatrixSP matrixsp(spotVolMatrixFromStrikes(&strikes[iMat][0],
                                                              nbStrikes,
                                                              iBackboneDateStart,
                                                              iBackboneDateEnd));
            CDoubleMatrix& matrix = *matrixsp;
            if (iBackboneDateEnd == iBackboneDateStart){
                int iStrike = 0;
                for (; iStrike < nbStrikes; ++iStrike){
                    impV[iMat][iStrike] = matrix[iStrike][0];
                }
            }
            else{
                double yearFrac = timeMetric->yearFrac(baseDate, maturities[iMat]);
                int iStrike = 0;
                for (; iStrike < nbStrikes; ++iStrike){
                    double varStart = Maths::square(matrix[iStrike][0] * backboneSqrTradYears[iBackboneDateStart]);
                    double varEnd = Maths::square(matrix[iStrike][1] * backboneSqrTradYears[iBackboneDateEnd]);
                    double yearFracDiff = backboneTradYears[iBackboneDateEnd] - backboneTradYears[iBackboneDateStart];
                    double slope = (varEnd - varStart) / yearFracDiff;
                    double var = varStart + slope * (yearFrac - backboneTradYears[iBackboneDateStart]);
                    if (Maths::isNegative(var)){
                        throw ModelException(routine,
                                             "Negative variance at maturity date " + 
                                             maturities[iMat].toString() + 
                                             " and at strike " +
                                             Format::toString(strikes[iMat][iStrike]));
                    }
                    else if(Maths::isZero(var)){
                        impV[iMat][iStrike] = 0.0;
                    } 
                    else{
                        impV[iMat][iStrike] = sqrt(var / yearFrac);
                    }
                }
            }
        }
    }

    // converts an absolute strike level into a level of moneyness
    double calcMoneyness(double strike,
                         double atmVol,
                         double tradTime) const{        
#if 1
        return ((strike / strikeRef - 1.0) / (atmVol * sqrt(tradTime)));
#else
        return (log(strike / strikeRef) / (atmVol * sqrt(tradTime)));
#endif
    }

    CDoubleMatrixSP spotVolMatrixFromStrikes(
            const double*         strikes,
            int                   strikes_size,
            int                   iBackboneDateStart,
            int                   iBackboneDateEnd) const{
        static const string routine("VolLogLinearPlus::spotVolMatrixFromStrikes");
        try{
            int nbBackboneDates = backboneDates.size();
            if (iBackboneDateStart < 0
                || iBackboneDateStart > iBackboneDateEnd
                || iBackboneDateEnd >= nbBackboneDates){
                throw ModelException(routine,
                                     Format::toString("Must have 0 <= iBackboneDateStart <= iBackboneDateEnd < nbBackboneDates;\n"
                                                      "got %ld, %ld and %ld respectively",
                                                      iBackboneDateStart,
                                                      iBackboneDateEnd,
                                                      nbBackboneDates));
            }

            CDoubleMatrixSP matrixsp(new CDoubleMatrix(strikes_size,
                                                       iBackboneDateEnd - iBackboneDateStart + 1));
            CDoubleMatrix& matrix = *matrixsp;
            int iMat = 0;
            int iBackboneDate = iBackboneDateStart;
            for (; iBackboneDate <= iBackboneDateEnd; iBackboneDate++, iMat++){
                double currTradYear = backboneTradYears[iBackboneDate];
                double atmVol = atmVols[iBackboneDate];
                int iStrike = 0;
                for (; iStrike < strikes_size; ++iStrike) {
                    // calculate the moneyness
                    double z = calcMoneyness(strikes[iStrike],
                                             atmVol,
                                             currTradYear);
                    double logVolSpread;
                    // lower tail
                    if (Maths::isNegative(z - gridMoneynesses[iBackboneDate][TAIL_DWN])){
                        double sigma = tailSpeeds[iBackboneDate][DWN];
                        double x = (z - gridMoneynesses[iBackboneDate][TAIL_DWN]) / sigma;
                        double a = extraCoeffs[iBackboneDate][DWN][0];
                        double b = extraCoeffs[iBackboneDate][DWN][1];
                        double c = extraCoeffs[iBackboneDate][DWN][2];
                        logVolSpread = a + b * exp(-0.5 * x * x) + c * exp(-0.5 * Maths::square(x - 1.0));
                    }
                    // upper tail
                    else if (Maths::isPositive(z - gridMoneynesses[iBackboneDate][TAIL_UP])){
                        double sigma = tailSpeeds[iBackboneDate][UP];
                        double x = (z - gridMoneynesses[iBackboneDate][TAIL_UP]) / sigma;
                        double a = extraCoeffs[iBackboneDate][UP][0];
                        double b = extraCoeffs[iBackboneDate][UP][1];
                        double c = extraCoeffs[iBackboneDate][UP][2];
                        logVolSpread = a + b * exp(-0.5 * x * x) + c * exp(-0.5 * Maths::square(x + 1.0));
                    }
                    // interpolation region
                    else{
                        logVolSpread = logVolSpreadInterPolys[iBackboneDate]->value(z);
                    }
                    matrix[iStrike][iMat] = atmVol * exp(logVolSpread);
                }
            }
            return matrixsp;
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
    VolSurface* spotVolSurfaceFromStrikes(
        const CDoubleArray&   strikes) const{
        static const string method("VolLogLinearPlus::spotVolSurfaceFromStrikes");
        try{
            CDoubleMatrixSP matrix(spotVolMatrixFromStrikes(&strikes[0],
                                                            strikes.size(),
                                                            0,
                                                            backboneDates.size() - 1));
            const VolSurface* backbone = getBackboneSurface();
            /** for performance use special constructors */
            return (new VolSurface(backbone, strikes, *matrix));
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    enum CutoffIndices{
        DWN = 0,
        UP,
        NB_CUTOFFS
    };

    enum StrikesIndices{
        TAIL_DWN = 0,
        MID_DWN,
        ATM,
        MID_UP,
        TAIL_UP,
        NB_STRIKES
    };

    // registered fields
    // 1Y specifics
    double        strikeRef;
    string        tenorRef;
    string        tenorTime;       // optional; applies to all BMs
    double        skew;
    double        midDeltaDown;
    double        midDeltaUp;
    double        tailDeltaDown;
    double        tailDeltaUp;
    double        spreadPowerShortTerm;
    double        spreadPowerLongTerm;
    // Term Structure
    StringArray   skewFactorBMs;        // polymorphic types not supported by IMS
    DoubleArray   skewFactors;
    StringArray   spreadBMs;
    DoubleArray   tailSpreadsDown;
    DoubleArray   midSpreadsUp;
    DoubleArray   tailSpreadsUp;
    StringArray   tailSpeedBMs;
    DoubleArray   tailSpeedsDown;
    DoubleArray   tailSpeedsUp;

    // transient fields (won't appear in dd interface)
    DateTime              baseDate;
    TimeMetricConstSP     timeMetric;
    CVolProcessedBSSP     atmVolCurve;
    DoubleArray           atmVols;              // vol at each bench date
    DateTimeArray         backboneDates;        // bench dates
    DoubleArray           backboneTradYears;    // trad years at each bench date
    DoubleArray           backboneSqrTradYears;    // trad year sqrt at each bench date
    // 1Y point specifics
    DateTime              date1Y;               // date at 1Y tenor $unregistered
    double                tradYear1Y;           // year frac to 1Y tenor
    // Term Structure
    DoubleArray           skewFactorsUsed;
    DoubleArray           midSpreadsDown;
    DoubleArrayArray      gridStrikes;          // 5 strikes
    DoubleArrayArray      gridMoneynesses;      // 5 moneynesses
    DoubleArrayArray      gridVols;             // 5 vol levels
    InterpolantArray      logVolSpreadInterPolys;   // spline interpolating polynomial
    DateTimeArray         skewFactorDates;
    DoubleArray           skewFactorTradYears;
    DateTimeArray         spreadDates;
    DoubleArray           spreadTradYears;
    DateTimeArray         tailSpeedDates;
    DoubleArray           tailSpeedTradYears;
    DoubleArrayArray      tailSpeeds;
    DoubleArrayArrayArray extraCoeffs;          // extrapolation coefficients

    static const double   exphalf;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VolLogLinearPlus, clazz);
        SUPERCLASS(CVolBaseParamSurface);
        IMPLEMENTS(IVolatilityBS);
        IMPLEMENTS(IVolatilityDVF);
        IMPLEMENTS(DeltaSurface::IShift);
        IMPLEMENTS(Calibrator::IAdjustable);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(strikeRef, "Reference strike");
        FIELD(tenorRef, "Reference tenor");
        FIELD(tenorTime, "Time of Day for Tenors");
        FIELD_MAKE_OPTIONAL(tenorTime);
        FIELD(skew, "Reference skew");
        FIELD(tailDeltaDown, "Reference lower tail delta");
        FIELD(midDeltaDown, "Reference lower middle delta");
        FIELD(midDeltaUp, "Reference upper middle delta");
        FIELD(tailDeltaUp, "Reference upper tail delta");
        FIELD(tailSpeedBMs, "Tail speed benchmarks (eg, '6M', '1Y')");
        FIELD(tailSpeedsDown, "Lower tail speed levels");
        FIELD(tailSpeedsUp, "Upper tail speed levels");
        FIELD(skewFactorBMs, "Skew power benchmarks (eg, '6M', '1Y')");
        FIELD(skewFactors, "Skew powers");
        FIELD(spreadBMs, "Volatility Spread benchmarks (eg, '6M', '1Y')");
        FIELD(tailSpreadsDown, "Reference lower tail volatility spread");
        FIELD(midSpreadsUp, "Reference upper middle volatility spread");
        FIELD(tailSpreadsUp, "Reference upper tail volatility spread");
        FIELD(spreadPowerShortTerm, "Short-term spread propagation power");
        FIELD(spreadPowerLongTerm, "Long-term spread propagation power");

        // transient
        FIELD(backboneDates, "");
        FIELD_MAKE_TRANSIENT(backboneDates);
        FIELD(timeMetric, "");
        FIELD_MAKE_TRANSIENT(timeMetric);
        FIELD(baseDate, "");
        FIELD_MAKE_TRANSIENT(baseDate);
        FIELD(atmVolCurve, "");
        FIELD_MAKE_TRANSIENT(atmVolCurve);
        FIELD(atmVols, "");
        FIELD_MAKE_TRANSIENT(atmVols);
        FIELD(backboneTradYears, "");
        FIELD_MAKE_TRANSIENT(backboneTradYears);
        FIELD(backboneSqrTradYears, "");
        FIELD_MAKE_TRANSIENT(backboneSqrTradYears);
        FIELD(tradYear1Y, "");
        FIELD_MAKE_TRANSIENT(tradYear1Y);        
        FIELD(midSpreadsDown, "");
        FIELD_MAKE_TRANSIENT(midSpreadsDown);
        FIELD(skewFactorsUsed, "");
        FIELD_MAKE_TRANSIENT(skewFactorsUsed);
        FIELD(gridStrikes, "");
        FIELD_MAKE_TRANSIENT(gridStrikes);
        FIELD(gridMoneynesses, "");
        FIELD_MAKE_TRANSIENT(gridMoneynesses);
        FIELD(gridVols, "");
        FIELD_MAKE_TRANSIENT(gridVols);
        FIELD(logVolSpreadInterPolys, "");
        FIELD_MAKE_TRANSIENT(logVolSpreadInterPolys);
        FIELD(extraCoeffs, "");
        FIELD_MAKE_TRANSIENT(extraCoeffs);
        FIELD(tailSpeeds, "");
        FIELD_MAKE_TRANSIENT(tailSpeeds);
        FIELD(skewFactorDates, "");
        FIELD_MAKE_TRANSIENT(skewFactorDates);
        FIELD(skewFactorTradYears, "");
        FIELD_MAKE_TRANSIENT(skewFactorTradYears);
        FIELD(spreadDates, "");
        FIELD_MAKE_TRANSIENT(spreadDates);
        FIELD(spreadTradYears, "");
        FIELD_MAKE_TRANSIENT(spreadTradYears);
        FIELD(tailSpeedDates, "");
        FIELD_MAKE_TRANSIENT(tailSpeedDates);
        FIELD(tailSpeedTradYears, "");
        FIELD_MAKE_TRANSIENT(tailSpeedTradYears);

        // add our fields and their ranges to central list
        Calibrator::IAdjustable::registerField(
            clazz, "skew",
            new InfiniteRange());
        Calibrator::IAdjustable::registerField(
            clazz, "skewFactors",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "midSpreadsUp",
            new InfiniteRange());
        Calibrator::IAdjustable::registerField(
            clazz, "tailSpreadsDown",
            new InfiniteRange());
        Calibrator::IAdjustable::registerField(
            clazz, "tailSpreadsUp",
            new InfiniteRange());
    }

    VolLogLinearPlus(): 
    CVolBaseParamSurface(TYPE),
    strikeRef(0.0), 
    tenorRef("1Y"),
    tenorTime("EOD"),
    skew(-0.02),
    midDeltaDown(0.25),
    midDeltaUp(0.25),
    tailDeltaDown(0.05),
    tailDeltaUp(0.05),
    spreadPowerShortTerm(0.5),
    spreadPowerLongTerm(0.5),
    tradYear1Y(0.0){}

    static IObject* defaultCtor(){
        return new VolLogLinearPlus();
    }

protected:
    /*** Build the parameterised vol and cache any values **/
    void buildCache() {
        static const string method("VolLogLinearPlus::buildCache");
        try{
            const VolSurface* backbone = getBackboneSurface();
            backboneDates = backbone->getDates();
            atmVolCurve = CVolProcessedBSSP(
                backbone->getProcessedVol(strikeRef));
            baseDate = backbone->getBaseDate();
            timeMetric = atmVolCurve->GetTimeMetric();
            int nbDates = backboneDates.size();
            atmVols.resize(nbDates);
            atmVolCurve->CalcVol(baseDate, 
                                 backboneDates, 
                                 CVolProcessedBS::fromFirst,
                                 atmVols);
            backboneTradYears.resize(nbDates);
            backboneSqrTradYears.resize(nbDates);
            gridStrikes.resize(nbDates);
            gridMoneynesses.resize(nbDates);
            gridVols.resize(nbDates);
            logVolSpreadInterPolys.resize(nbDates);
            extraCoeffs.resize(nbDates);
            tailSpeeds.resize(nbDates);
            for (int iDate = 0; iDate < atmVols.size(); ++iDate){
                backboneTradYears[iDate] = timeMetric->yearFrac(baseDate, backboneDates[iDate]);
                backboneSqrTradYears[iDate] = sqrt(backboneTradYears[iDate]);
                gridStrikes[iDate].resize(NB_STRIKES);
                gridMoneynesses[iDate].resize(NB_STRIKES);
                gridVols[iDate].resize(NB_STRIKES);
                extraCoeffs[iDate].resize(NB_CUTOFFS);
                for (int iCutoff = 0; iCutoff < NB_CUTOFFS; ++iCutoff){
                    extraCoeffs[iDate][iCutoff].resize(3);
                }
                tailSpeeds[iDate].resize(NB_CUTOFFS);
            }
            // 1Y point
            int tenorTime2 = DateTime::timeConvert(tenorTime);            
            MaturityPeriod period1Y(tenorRef);
            date1Y = DateTime(period1Y.toDate(baseDate).getDate(), tenorTime2);
            tradYear1Y = timeMetric->yearFrac(baseDate, date1Y);
            // skew factor BMs
            // first check the BMs are increasing
            // and that if the 1y ref is among the BMs then 
            // the corresponding skew factor == 1
            bool date1YFound = false;
            int iSkewFactorBM = 0;
            for (; iSkewFactorBM < skewFactorBMs.size(); ++iSkewFactorBM){
                MaturityPeriod period(skewFactorBMs[iSkewFactorBM]);
                skewFactorDates[iSkewFactorBM] = DateTime(period.toDate(baseDate).getDate(), tenorTime2);
                if (iSkewFactorBM > 0
                    && skewFactorDates[iSkewFactorBM - 1] >= skewFactorDates[iSkewFactorBM]){
                    throw ModelException(method,
                                         "The skew factor benchmarks must be increasing.\nCompare the "
                                         + Format::toString(iSkewFactorBM)
                                         + "-th benchmark date ("
                                         + skewFactorDates[iSkewFactorBM - 1].toString()
                                         + ") with the "
                                         + Format::toString(iSkewFactorBM + 1)
                                         + "-th benchmark date ("
                                         + skewFactorDates[iSkewFactorBM].toString()
                                         + ")");
                }
                if (skewFactorDates[iSkewFactorBM] == date1Y){
                    if (!Maths::isZero(skewFactors[iSkewFactorBM] - 1.0)){
                        throw ModelException(method,
                                             "the "
                                             + Format::toString(iSkewFactorBM + 1)
                                             + "-th skew factor ("
                                             + Format::toString(skewFactors[iSkewFactorBM])
                                             + ") should be equal to 1.0,\nas the "
                                             + Format::toString(iSkewFactorBM + 1)
                                             + "-th skew factor benchmark date is equal to the reference date ("
                                             + skewFactorDates[iSkewFactorBM].toString()
                                             + ")");
                    }
                    date1YFound = true;
                }
                skewFactorTradYears[iSkewFactorBM] 
                    = timeMetric->yearFrac(baseDate, skewFactorDates[iSkewFactorBM]);
            }
            // if the 1y ref date hasn't been found, need to insert it
            if (!date1YFound){
                DateTimeArray skewFactorDatesNew(skewFactorBMs.size() + 1);
                DoubleArray skewFactorTradYearsNew(skewFactorBMs.size() + 1);
                skewFactorsUsed.resize(skewFactorBMs.size() + 1);
                iSkewFactorBM = 0;
                int iNewSkewFactorBM = 0;
                while (iSkewFactorBM < skewFactorBMs.size()
                       && skewFactorDates[iSkewFactorBM] < date1Y){
                    skewFactorDatesNew[iNewSkewFactorBM] = skewFactorDates[iSkewFactorBM];
                    skewFactorTradYearsNew[iNewSkewFactorBM] = skewFactorTradYears[iSkewFactorBM];
                    skewFactorsUsed[iNewSkewFactorBM] = skewFactors[iSkewFactorBM];
                    ++iSkewFactorBM;
                    ++iNewSkewFactorBM;
                }
                skewFactorDatesNew[iNewSkewFactorBM] = date1Y;
                skewFactorTradYearsNew[iNewSkewFactorBM] = tradYear1Y;
                skewFactorsUsed[iNewSkewFactorBM] = 1.0;
                ++iNewSkewFactorBM;
                while (iSkewFactorBM < skewFactorBMs.size()){
                    skewFactorDatesNew[iNewSkewFactorBM] = skewFactorDates[iSkewFactorBM];
                    skewFactorTradYearsNew[iNewSkewFactorBM] = skewFactorTradYears[iSkewFactorBM];
                    skewFactorsUsed[iNewSkewFactorBM] = skewFactors[iSkewFactorBM];
                    ++iSkewFactorBM;
                    ++iNewSkewFactorBM;
                }
                skewFactorDates = skewFactorDatesNew;
                skewFactorTradYears = skewFactorTradYearsNew;
            }
            else{
                skewFactorsUsed = skewFactors;
            }
            // spread BMs
            // calculate mean reversion levels along the way
            // mean reversion levels are calculated as "half-lifes"
            int iSpreadBM = 0;
            for (; iSpreadBM < spreadBMs.size(); ++iSpreadBM){
                MaturityPeriod period(spreadBMs[iSpreadBM]);
                spreadDates[iSpreadBM] = DateTime(period.toDate(baseDate).getDate(), tenorTime2);
                if (iSpreadBM > 0
                    && spreadDates[iSpreadBM - 1] >= spreadDates[iSpreadBM]){
                    throw ModelException(method,
                                         "The volatility spread benchmarks must be increasing.\nCompare the "
                                         + Format::toString(iSpreadBM)
                                         + "-th benchmark date ("
                                         + spreadDates[iSpreadBM - 1].toString()
                                         + ") with the "
                                         + Format::toString(iSpreadBM + 1)
                                         + "-th benchmark date ("
                                         + spreadDates[iSpreadBM].toString()
                                         + ")");
                }
                spreadTradYears[iSpreadBM] 
                    = timeMetric->yearFrac(baseDate, spreadDates[iSpreadBM]);
            }
            // skew factor BMs
            int iTailSpeedBM = 0;
            for (; iTailSpeedBM < tailSpeedBMs.size(); ++iTailSpeedBM){
                MaturityPeriod period(tailSpeedBMs[iTailSpeedBM]);
                tailSpeedDates[iTailSpeedBM] = DateTime(period.toDate(baseDate).getDate(), tenorTime2);
                if (iTailSpeedBM > 0
                    && tailSpeedDates[iTailSpeedBM - 1] >= tailSpeedDates[iTailSpeedBM]){
                    throw ModelException(method,
                                         "The tail speed benchmarks must be increasing.\nCompare the "
                                         + Format::toString(iTailSpeedBM)
                                         + "-th benchmark date ("
                                         + tailSpeedDates[iTailSpeedBM - 1].toString()
                                         + ") with the "
                                         + Format::toString(iTailSpeedBM + 1)
                                         + "-th benchmark date ("
                                         + tailSpeedDates[iTailSpeedBM].toString()
                                         + ")");
                }
                tailSpeedTradYears[iTailSpeedBM] 
                    = timeMetric->yearFrac(baseDate, tailSpeedDates[iTailSpeedBM]);
            }
            update();
        }
        catch(exception& e){
            throw ModelException(e, method);
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
        static const string method("VolLogLinearPlus::update");
        try{
            // format the deltas into a vector
            DoubleArray deltas(NB_STRIKES);
            deltas[TAIL_DWN] = tailDeltaDown;
            deltas[MID_DWN] = midDeltaDown;
            deltas[MID_UP] = midDeltaUp;
            deltas[TAIL_UP] = tailDeltaUp;
            // will need to format the spreads and the log vol spreads into vectors, too
            DoubleArray spreads(NB_STRIKES);
            DoubleArray logSpreads(NB_STRIKES);
            // loop over all backbone dates
            int nbBackboneDates = backboneDates.size();
            int nbSkewFactorDates = skewFactorDates.size();
            int iSkewFactorDate = 0;
            int nbSpreadDates = spreadDates.size();
            int iSpreadDate = 0;
            int nbTailSpeedDates = tailSpeedBMs.size();
            int iTailSpeedDate = 0;
            for (int iBackboneDate = 0; iBackboneDate < nbBackboneDates; ++iBackboneDate){
                DateTime currDate = backboneDates[iBackboneDate];
                double currTradYear = backboneTradYears[iBackboneDate];
                // bracket current date with skew factor dates
                while (iSkewFactorDate < nbSkewFactorDates
                       && skewFactorDates[iSkewFactorDate] < currDate){
                    ++iSkewFactorDate;
                }
                // calculate current skew factor
                double currSkewFactor;
                // extrapolate flat if current date <= first skew factor date
                if (iSkewFactorDate == 0){
                    currSkewFactor = skewFactorsUsed[0];
                }
                // extrapolate flat if current date > last skew factor date
                else if (iSkewFactorDate == nbSkewFactorDates){
                    currSkewFactor = skewFactorsUsed[nbSkewFactorDates - 1];
                }
                // otherwise linear interpolation
                else{
                    currSkewFactor = skewFactorsUsed[iSkewFactorDate - 1] 
                                    + (skewFactorsUsed[iSkewFactorDate] - skewFactorsUsed[iSkewFactorDate - 1]) 
                                    / (skewFactorTradYears[iSkewFactorDate] - skewFactorTradYears[iSkewFactorDate - 1])
                                    * (currTradYear - skewFactorTradYears[iSkewFactorDate - 1]);
                }
                // calculate current skew
                double currSkew = skew / sqrt(currTradYear / tradYear1Y) * currSkewFactor;
                double currBeta = - currSkew / log(0.9);
                // bracket current date with spread dates
                while (iSpreadDate < nbSpreadDates
                       && spreadDates[iSpreadDate] < currDate){
                    ++iSpreadDate;
                }
                // calculate current spread
                // extrapolate using the short-term power if current date <= first spread date
                if (iSpreadDate == 0){
                    double effTradYear = pow(currTradYear / spreadTradYears[0], 
                                             - spreadPowerShortTerm);
                    spreads[TAIL_DWN] = tailSpreadsDown[0] * effTradYear;
                    spreads[MID_DWN] = midSpreadsDown[0] * effTradYear;
                    spreads[MID_UP] = midSpreadsUp[0] * effTradYear;
                    spreads[TAIL_UP] = tailSpreadsUp[0] * effTradYear;
                }
                // extrapolate using the long-term power if current date > last spread date
                else if (iSpreadDate == nbSpreadDates){
                    double effTradYear = pow(currTradYear / spreadTradYears[nbSpreadDates - 1], 
                                             - spreadPowerLongTerm);
                    spreads[TAIL_DWN] = tailSpreadsDown[nbSpreadDates - 1] * effTradYear;
                    spreads[MID_DWN] = midSpreadsDown[nbSpreadDates - 1] * effTradYear;
                    spreads[MID_UP] = midSpreadsUp[nbSpreadDates - 1] * effTradYear;
                    spreads[TAIL_UP] = tailSpreadsUp[nbSpreadDates - 1] * effTradYear;
                }
                // otherwise linear interpolation
                else{
                    spreads[TAIL_DWN] = tailSpreadsDown[iSpreadDate - 1] 
                                        + (tailSpreadsDown[iSpreadDate] - tailSpreadsDown[iSpreadDate - 1]) 
                                        / (spreadTradYears[iSpreadDate] - spreadTradYears[iSpreadDate - 1])
                                        * (currTradYear - spreadTradYears[iSpreadDate - 1]);
                    spreads[MID_DWN] = midSpreadsDown[iSpreadDate - 1] 
                                        + (midSpreadsDown[iSpreadDate] - midSpreadsDown[iSpreadDate - 1]) 
                                        / (spreadTradYears[iSpreadDate] - spreadTradYears[iSpreadDate - 1])
                                        * (currTradYear - spreadTradYears[iSpreadDate - 1]);
                    spreads[MID_UP] = midSpreadsUp[iSpreadDate - 1] 
                                        + (midSpreadsUp[iSpreadDate] - midSpreadsUp[iSpreadDate - 1]) 
                                        / (spreadTradYears[iSpreadDate] - spreadTradYears[iSpreadDate - 1])
                                        * (currTradYear - spreadTradYears[iSpreadDate - 1]);
                    spreads[TAIL_UP] = tailSpreadsUp[iSpreadDate - 1] 
                                        + (tailSpreadsUp[iSpreadDate] - tailSpreadsUp[iSpreadDate - 1]) 
                                        / (spreadTradYears[iSpreadDate] - spreadTradYears[iSpreadDate - 1])
                                        * (currTradYear - spreadTradYears[iSpreadDate - 1]);
                }
                // loop over all strikes
                // calculate strike corresponding to given delta
                // calculate corresponding volatility and log vol spread
                double currAtmVol = atmVols[iBackboneDate];
                double currLogAtmVol = log(currAtmVol);
                int iStrike = 0;
                for (; iStrike < NB_STRIKES; ++iStrike){
                    if (iStrike == ATM){
                        gridVols[iBackboneDate][iStrike] = currAtmVol;
                        gridStrikes[iBackboneDate][iStrike] = strikeRef;
                        logSpreads[iStrike] = 0.0;
                        continue;
                    }
                    double epsilon = (iStrike < ATM ? -1.0 : 1.0);
                    double d1 = epsilon * N1Inverse(deltas[iStrike]);
                    double a = currBeta * currBeta * currTradYear;
                    double currVol = currAtmVol + spreads[iStrike];
                    double c = currVol * currVol * currTradYear 
                               - 2.0 * d1 * currVol * sqrt(currTradYear);
                    double logstrike;
                    if (Maths::isZero(a)){
                        logstrike = 0.5 * c;
                    }
                    else{
                        double bdash = currVol * currBeta * currTradYear 
                                       - 1.0 - d1 * currBeta * sqrt(currTradYear);
                        double deltadash = bdash * bdash - a * c;
                        if (Maths::isNegative(deltadash)){
                            throw ModelException(method,
                                                 "It is not possible to back up a strike level "
                                                 "from the volatility spread ("
                                                 + Format::toString(spreads[iStrike])
                                                 + ") \ngiven at the "
                                                 + Format::toString(iStrike + 1)
                                                 + "-th strike for a "
                                                 + Format::toString(deltas[iStrike])
                                                 + "-delta");
                        }
                        logstrike = (-bdash - sqrt(deltadash)) / a;
                    }
                    gridVols[iBackboneDate][iStrike] = currVol + currBeta * logstrike;
                    gridStrikes[iBackboneDate][iStrike] = strikeRef * exp(logstrike);
                    logSpreads[iStrike] = log(gridVols[iBackboneDate][iStrike]) - currLogAtmVol;
                }
                // calculate moneyness
                for (iStrike = 0; iStrike < NB_STRIKES; ++iStrike){
                    gridMoneynesses[iBackboneDate][iStrike] 
                        = calcMoneyness(gridStrikes[iBackboneDate][iStrike],
                                        currAtmVol,
                                        currTradYear);
                }
#if 1
                DoubleArray tailLogVolSpreadSlopes(2);
                int iCutoff = DWN;
                int shift = 0;
                for (iStrike = TAIL_DWN; iCutoff < NB_CUTOFFS; ++iCutoff, iStrike = TAIL_UP, shift = 2){
                    DoubleArray coeffs(3);
                    polcof(&gridMoneynesses[iBackboneDate][iStrike - shift], 
                           &logSpreads[iStrike - shift], 
                           2, 
                           &coeffs[0]);
                    tailLogVolSpreadSlopes[iCutoff] = 2.0 * coeffs[2]  * gridMoneynesses[iBackboneDate][iStrike]
                                                      + coeffs[1];
                }
                // spline the log vol spreads by imposing zero 2nd derivatives
                // at the boundaries
                NRSpline spliner(1, tailLogVolSpreadSlopes[DWN],
                                 1, tailLogVolSpreadSlopes[UP]);
                logVolSpreadInterPolys[iBackboneDate] 
                    = NRSpline::InterpolantSP::constCast(
                                    NRSpline::computeInterp(spliner,
                                                            gridMoneynesses[iBackboneDate],
                                                            logSpreads));
#else
                CubicShapePresSplineInterpECnd spliner(false);
                logVolSpreadInterPolys[iBackboneDate] 
                    = Interpolator::InterpolantSP::constCast(
                                spliner.computeInterp(gridMoneynesses[iBackboneDate],
                                                      logSpreads));
#endif
                // bracket current date with tail speed dates
                while (iTailSpeedDate < nbTailSpeedDates
                       && tailSpeedDates[iTailSpeedDate] < currDate){
                    ++iTailSpeedDate;
                }
                // extrapolate flat if current date <= first tail speed date
                if (iTailSpeedDate == 0){
                    tailSpeeds[iBackboneDate][DWN] = tailSpeedsDown[0];
                    tailSpeeds[iBackboneDate][UP] = tailSpeedsUp[0];
                }
                // extrapolate flat if current date > last tail speed date
                else if (iTailSpeedDate == nbTailSpeedDates){
                    tailSpeeds[iBackboneDate][DWN] = tailSpeedsDown[nbTailSpeedDates - 1];
                    tailSpeeds[iBackboneDate][UP] = tailSpeedsUp[nbTailSpeedDates - 1];
                }
                // otherwise linear interpolation
                else{
                    tailSpeeds[iBackboneDate][DWN] = tailSpeedsDown[iTailSpeedDate - 1] 
                                          + (tailSpeedsDown[iTailSpeedDate] - tailSpeedsDown[iTailSpeedDate - 1]) 
                                          / (tailSpeedTradYears[iTailSpeedDate] - tailSpeedTradYears[iTailSpeedDate - 1])
                                          * (currTradYear - tailSpeedTradYears[iTailSpeedDate - 1]);
                    tailSpeeds[iBackboneDate][UP] = tailSpeedsUp[iTailSpeedDate - 1] 
                                          + (tailSpeedsUp[iTailSpeedDate] - tailSpeedsUp[iTailSpeedDate - 1]) 
                                          / (tailSpeedTradYears[iTailSpeedDate] - tailSpeedTradYears[iTailSpeedDate - 1])
                                          * (currTradYear - tailSpeedTradYears[iTailSpeedDate - 1]);
                }
                // calculate extrapolation coefs   
                for (iCutoff = DWN, iStrike = TAIL_DWN; iCutoff < NB_CUTOFFS; ++iCutoff, iStrike = TAIL_UP){
                    double dLogVolSpread_dZ 
                        = logVolSpreadInterPolys[iBackboneDate]->value(gridMoneynesses[iBackboneDate][iStrike], 1);
                    double d2LogVolSpread_dZ2 
                        = logVolSpreadInterPolys[iBackboneDate]->value(gridMoneynesses[iBackboneDate][iStrike], 2);
                    double sigma = tailSpeeds[iBackboneDate][iCutoff];
                    double dLogVolSpread_dX = sigma * dLogVolSpread_dZ;
                    double d2LogVolSpread_dX2 = sigma * sigma * d2LogVolSpread_dZ2;
                    double epsilon = (iCutoff == UP ? 1.0 : -1.0);
                    extraCoeffs[iBackboneDate][iCutoff][0] = logSpreads[iStrike]
                                                             + epsilon * dLogVolSpread_dX 
                                                             + d2LogVolSpread_dX2;
                    extraCoeffs[iBackboneDate][iCutoff][1] = - d2LogVolSpread_dX2;
                    extraCoeffs[iBackboneDate][iCutoff][2] = - epsilon * exphalf * dLogVolSpread_dX;
                    // check that the extrapolation is monotonous (ie that "b" and "c" coeffs
                    // are of the same sign)
                    // if not, sacrifice the continuity of the second derivative
                    // by killing the "b" exp term
                    if (Maths::isNegative(extraCoeffs[iBackboneDate][iCutoff][1] 
                                          * extraCoeffs[iBackboneDate][iCutoff][2])){
                        d2LogVolSpread_dX2 = 0.0;
                        extraCoeffs[iBackboneDate][iCutoff][0] = logSpreads[iStrike]
                                                                 + epsilon * dLogVolSpread_dX 
                                                                 + d2LogVolSpread_dX2;
                        extraCoeffs[iBackboneDate][iCutoff][1] = - d2LogVolSpread_dX2;
                        extraCoeffs[iBackboneDate][iCutoff][2] = - epsilon * exphalf * dLogVolSpread_dX;
                    }
                }
            }
        }
        catch(exception& e){
            throw ModelException(e, "VolLogLinearPlus::update");
        }
    }
};

const double VolLogLinearPlus::exphalf = exp(0.5);

CClassConstSP const VolLogLinearPlus::TYPE =
CClass::registerClassLoadMethod("VolLogLinearPlus", typeid(VolLogLinearPlus), load);

CClassConstSP const VolLogLinearPlus::VolLogLinearPlusVolParam::TYPE =
CClass::registerClassLoadMethod("VolLogLinearPlus::VolLogLinearPlusVolParam", 
                                typeid(VolLogLinearPlusVolParam), load);

// external symbol to allow class to be forced to be linked in
bool VolLogLinearPlusLinkIn(){
    return (VolLogLinearPlus::TYPE != 0);
}

/* Debug Info Addin */
class VolLogLinearPlusDebugInfoAddin: public CObject{
public:
    static CClassConstSP const TYPE;

private:
    CMarketDataSP market;
    string        name;
    
    static IObjectSP getInfo(VolLogLinearPlusDebugInfoAddin* params) {
        static const string routine = "VolLogLinearPlusDebugInfoAddin::prices";
        try {
            // create a non pricing model in order to get the market data
            MarketDataFetcherSP mdf(new MarketDataFetcherLN("VolLogLinearPlus"));
            NonPricingModel dummyModel(mdf);    
            MarketObjectSP mo(params->market->GetData(&dummyModel,
                                                      params->name,
                                                      CClass::forName("VolLogLinearPlus")));
            VolLogLinearPlus& vol = dynamic_cast<VolLogLinearPlus&>(*mo);
            ObjectArraySP info(new ObjectArray(2));
            (*info)[0] = DoubleArrayArraySP(copy(&vol.gridStrikes));
            (*info)[1] = DoubleArrayArraySP(copy(&vol.gridVols));
            return info;
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    VolLogLinearPlusDebugInfoAddin():
    CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VolLogLinearPlusDebugInfoAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultVolLogLinearPlusDebugInfoAddin);
        FIELD(market, "market");
        FIELD(name, "name");

        Addin::registerClassObjectMethod("VOLLOGLINEARPLUS_DEBUGINFO_GET",
                                         Addin::MARKET,
                                         "Get a VolLogLinearPlus's debug info",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)getInfo);
    }

    static IObject* defaultVolLogLinearPlusDebugInfoAddin(){
        return new VolLogLinearPlusDebugInfoAddin();
    }
};

CClassConstSP const VolLogLinearPlusDebugInfoAddin::TYPE =
CClass::registerClassLoadMethod("VolLogLinearPlusDebugInfoAddin", typeid(VolLogLinearPlusDebugInfoAddin), load);

DRLIB_END_NAMESPACE
