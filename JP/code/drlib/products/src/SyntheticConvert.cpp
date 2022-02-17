//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SyntheticConvert.cpp
//
//   Description : Synthetic Convertible 
//
//   Author      : Bruno O Melka
//
//   Date        : 01 Jun 2005
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/SVGenIRFloatRate.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"

DRLIB_BEGIN_NAMESPACE

class SyntheticConvert: public CInstrument,
                public virtual LastSensDate,
                virtual public Theta::Shift,
                public virtual IMCIntoProduct{
public:
    static CClassConstSP const TYPE; 

    virtual void validatePop2Object(){
        static const string method = "SyntheticConvert::validatePop2Object";
    }

    virtual void Validate() {
        static const string method = "SyntheticConvert::Validate";
        try {
            if (monitorDates[0].isLess(valueDate)) {
                throw ModelException(method,"first date of swap has to be spot or forward but not in the past");
            }
        }    
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    DateTime endDate(const Sensitivity* sensControl) const {
        return monitorDates[monitorDates.size() - 1];
    }

    DateTime getValueDate() const {
        return valueDate;
    }

    // roll through time 
    bool sensShift(Theta* theta){
        // roll today 
        valueDate = theta->rollDate(valueDate);
        return true; // continue to tweak components which implement Theta
    }
    
    /** Get the asset and discount market data */
    void GetMarket( const IModel*          model, 
                    const CMarketDataSP    market)
    {
        market->GetReferenceDate(valueDate);
        discount.getData(model, market);

        CAsset::getAssetMarketData(model, market.get(), ccyTreatment, discount, asset);
        ICDSParSpreads::getMarketData(model, market.get(), discount.getName(), credit);

        couponCurve = YieldCurveSP(discount.getSP().clone());
        couponCurve->setProjectionCurve(); // use 'growth' zc

    }

    /** Implementation of MonteCarlo::IntoProduct interface */
    IMCProduct* createProduct(const MonteCarlo* model) const; // see below


    /** Returns the name of the instrument's discount currency */
    string discountYieldCurveName() const {
        return discount.getName();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(SyntheticConvert, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultSyntheticConvert);
        FIELD(ccyTreatment,"Currency Treatment");
        FIELD(fixedDCC, "day count convention of swap");
        FIELD(fixedFrequency, "frequency of payment compared to monitoring");
        FIELD(floatFrequency, "frequency of payment compared to monitoring");
        FIELD(fixedCoupon, "fixed coupon");
        FIELD(strike, "strike of the option");
        FIELD(adjustment, "adjustment to forward to look like option");
        FIELD(recovery, "recovery");
        FIELD(swapNotional, "notional");
        FIELD(monitorDates, "fixing dates of the swap");
        FIELD(discount,"discount curve");
        FIELD_NO_DESC(couponCurve);
        FIELD_MAKE_TRANSIENT(couponCurve);
        FIELD(valueDate,"valuation Date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(cashPhys, "T=PV wrt ytm, F=PV wrt zcurve");
        FIELD(asset, "asset");
        FIELD(credit, "credit");

    }
     
    static IObject* defaultSyntheticConvert(){
        return new SyntheticConvert();
    }

private:

    SyntheticConvert():CInstrument(TYPE) {}; 
    SyntheticConvert(const SyntheticConvert& rhs);
    SyntheticConvert& operator=(const SyntheticConvert& rhs);

    class                   MC;
    friend class            MC;

    string                  ccyTreatment;           // None(vanilla), protected or struck
    string                  fixedDCC;               // day count convention of rate
    int                     fixedFrequency;         // compared to monitoring
    int                     floatFrequency;         // compared to monitoring
    double                  fixedCoupon;
    double                  strike;
    double                  adjustment;             // for pricing early optionality
    double                  recovery;
    double                  swapNotional;           // notional
    DateTimeArray           monitorDates;           // refix dates of the swap
    YieldCurveWrapper       discount;
    YieldCurveSP            couponCurve;            // transient - set to 'projection' curve    
    DateTime                valueDate;
    bool                    cashPhys;               // T=PV wrt ytm, F=PV wrt zcurve
    CAssetWrapper           asset;
    ICDSParSpreadsWrapper   credit;
};

CClassConstSP const SyntheticConvert::TYPE = CClass::registerClassLoadMethod(
    "SyntheticConvert", typeid(SyntheticConvert), SyntheticConvert::load);

/* MC product class for SyntheticConvert */
class SyntheticConvert::MC : public MCProductClient{

private:

    SVGenSpot               ::IStateVarSP  assetSV;                // asset state variable
        SVExpectedDiscFactorSP                   rateSV;                 // floating rate state variable
    SVSurvivalDiscFactorSP  survivalProbaSV;        // survival probability state variable
    SVDiscFactorSP                   discountSV;             // discounting dtate variable

    SVGenSpotSP                            assetGen;               // generator for asset
    SVGenExpectedDiscFactorSP              rateGen;                // generator for floating rate
    SVGenSurvivalDiscFactorSP              survivalProbaGen;       // generator for survival probablilty
    SVGenDiscFactorSP                      discountGen;            // generator for discounting
    
    double                              notional;               // from instrument
    double                              strike;
    double                              adjustment;             // for pricing early optionality
    double                              recovery;
    int                                 floatFrequency;
    int                                 fixedFrequency;
    DoubleArray                         fixedCoupons;
    
protected:
    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
        assetSV             =   assetGen->getSpotSV(newPathGen);
        rateSV              =   rateGen->getSVExpectedDiscFactor(newPathGen);
        survivalProbaSV     =   survivalProbaGen->getSVSurvivalDiscFactor(newPathGen);
        discountSV          =   discountGen->getSVDiscFactor(newPathGen);
    }

public:
    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        svCollector->append(assetGen.get()); 
        svCollector->append(rateGen.get());
        svCollector->append(survivalProbaGen.get());
        svCollector->append(discountGen.get());
    }

    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    MC(const SyntheticConvert*          inst,
       const SimSeriesSP&       simSeries,
       InstrumentSettlementSP   instSettle):
                MCProductClient(inst->asset.get(),
                                inst->valueDate,
                                inst->discount.get(),
                                IRefLevelSP(IRefLevel::Util::makeZero(inst->valueDate)), // fix!
                                simSeries, // fix!
                                IPastValuesSP(IPastValues::Util::makeTrivial(inst->valueDate, 0.0)), // fix
                                instSettle.get(),
                                inst->monitorDates[inst->monitorDates.size() - 1]), 
                notional(inst->swapNotional),
                strike(inst->strike),
                adjustment(inst->adjustment),
                recovery(inst->recovery),
                floatFrequency(inst->floatFrequency),
                fixedFrequency(inst->fixedFrequency),
                fixedCoupons((inst->monitorDates.size() - 1) / fixedFrequency),
                assetGen(new SVGenSpot(1, inst->monitorDates)),
                survivalProbaGen(new SVGenSurvivalDiscFactor(inst->valueDate,
                                                inst->credit.getSP(),
                                                inst->monitorDates)),
                discountGen(new SVGenDiscFactor(inst->valueDate,
                                                inst->discount.getSP(),
                                                inst->monitorDates))   {

        DayCountConventionSP dcc(DayCountConventionFactory::make(inst->fixedDCC));
        DateTimeArray floatPayDates((inst->monitorDates.size() - 1) / floatFrequency);
        DateTimeArray fixedPayDates((inst->monitorDates.size() - 1) / fixedFrequency);
        
        for (int i = 0; i < floatPayDates.size(); i++) {
            floatPayDates[i] = inst->monitorDates[(i+1)  * floatFrequency];
        }
        
        for (int j = 0; j < fixedPayDates.size(); j++) {
            fixedPayDates[j] = inst->monitorDates[(j+1) * fixedFrequency];
            fixedCoupons[j] = inst->fixedCoupon * dcc->years((j == 0)? inst->monitorDates[0]:fixedPayDates[j-1], fixedPayDates[j]);
        }   

        rateGen = SVGenExpectedDiscFactorSP(new SVGenExpectedDiscFactor(inst->monitorDates[0],
                                    inst->monitorDates[0],
                                    inst->couponCurve,
                                    floatPayDates,
                                    false)); // don't compute log
    }


    /** Called within the simulation loop */
    virtual void payoff(const IPathGenerator*  pathGen,
                        IMCPrices&                prices) {
        
        const SVPath& pathDf      = discountSV->path();
        const SVPath& pathSp      = survivalProbaSV->path();
        const SVPath& pathRate    = rateSV->path();
        const SVPath& pathAsset   = assetSV->path(0);

        int endIdx = pathDf.end();
        int iStep;

        double refSpot = pathAsset[0];
        double riskySum = 0.0;
        double prevZeroToCoupon = 1.0;
        double prevProba = 1.0;
        double prevSumCoupons = 0.0;
        double payOffSum = 0.0;


        DoubleArray riskyCoupon(endIdx - 1, 0.0);
        DoubleArray sumCoupons(endIdx - 1, 0.0);

        // pre-compute the value of future coupons and total sum in order to get renaining swap.
        for (int iSum = 1; iSum < endIdx; iSum++) {
            sumCoupons[iSum - 1] = prevSumCoupons;
            if ((iSum % floatFrequency) == 0) {
                double zeroToCoupon = pathRate[iSum / floatFrequency - 1];
                double floatCoupon = (prevZeroToCoupon / zeroToCoupon - 1) * pathDf[iSum];
                sumCoupons[iSum - 1] += floatCoupon;
                riskyCoupon[iSum - 1] += floatCoupon * pathSp[iSum];
                prevZeroToCoupon = zeroToCoupon;
            }
            if ((iSum % fixedFrequency) == 0) {
                double fixedCoupon = fixedCoupons[iSum / fixedFrequency - 1] * pathDf[iSum];
                sumCoupons[iSum - 1] -= fixedCoupon;
                riskyCoupon[iSum - 1] -= fixedCoupon * pathSp[iSum];
            }
            riskyCoupon[iSum - 1] += (1 - recovery) * pathDf[iSum] * (prevProba - pathSp[iSum]);
            riskySum += riskyCoupon[iSum - 1];
            prevSumCoupons = sumCoupons[iSum - 1];
            prevProba = pathSp[iSum];
        }

        // check for early exercise opportunities
        // adjustment allows to price early optionality.
        for (iStep = 1; iStep < endIdx; iStep++) {
            payOffSum = sumCoupons[iStep - 1];
            riskySum -= riskyCoupon[iStep - 1];
            double remainingSum = riskySum / (pathDf[iStep] * pathSp[iStep]);
            if (remainingSum + (pathAsset[iStep] / refSpot - strike) - adjustment > 0.0) {
                break;
            }
        }

        // compute value of option, whether at early exercise or at matutiry;
        int assetIdx = Maths::min(iStep, endIdx - 1);
        double option = Maths::max(0.0, (pathAsset[assetIdx] / refSpot - strike)) * pathDf[assetIdx];
        double myPayoff = payOffSum - option / strike; // divide by strike to get the equivalent number of shares.
        prices.add(notional * myPayoff); // finally scale by notional
        
    }

    virtual double pvFromPaymentDate() const{
        return 1.0; // already applied pv factor in payoff
    }

};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* SyntheticConvert::createProduct(const MonteCarlo* model) const {
    // the simSeries passed to the IMCProduct is redundant
    // well not quite. The MC wants to know the last sim date
    SimSeriesSP simSeries(new SimSeries(1)); /* create empty one */
    simSeries->addDates(DateTimeArray(1, monitorDates[monitorDates.size() - 1]));
    InstrumentSettlementSP instSettle(new CashSettlePeriod(0));
    return new MC(this, simSeries, instSettle);
}

// for class loading 
bool SyntheticConvertLoad() {
    return (SyntheticConvert::TYPE != 0);
}

DRLIB_END_NAMESPACE
