//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ForwardBasisIndexSeries.cpp
//
//   Description : IMCPrices a series of basis forwards.
//                 Primarily a test vehicle for SRM3
//
//   Date        : May 2006
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp"
#include "edginc/SVGenExpectedBasisForward.hpp"
#include "edginc/BootstrappedBasisIndexCurve.hpp"
#include "edginc/SVGenIRSwap.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/StubFactory.hpp"
#include <fstream>

DRLIB_BEGIN_NAMESPACE

/** ForwardOptionSeries product  - a strip of options on forwards.
    Derived from Generic1Factor just in case
    we ever want to book one - although several fields are a bit meaningless */
class ForwardBasisIndexSeries: public CInstrument, 
                               virtual public IMCIntoProduct{

private:
    /// fields ////////
    DateTime                valueDate;
    double                  notional;      // Size of deal 
    DateTimeArray           evalDates;     
    DateTimeArray           resetDates;    // when each forward starts
    BootstrappedBasisIndexCurveWrapper basisCurve;  // basis curve
    YieldCurveWrapper       discountCurve; // discounting curve (should be same as the underlying of basis)
    FXAssetWrapper          fxAsset;    // needed for foreign denominated basis
    bool                    doSquare;       // TODO: for testing only
    bool                    spreadOnly;
    bool                    dumpPrices; 
    string                  dumpFileName;
    bool                    doDeltaPricing;
    int                     nbDayLag;       // only use when doDeltaPricing = TRUE

public:
    static CClassConstSP const TYPE;
    friend class ForwardBasisIndexSeriesMC;

    // validation
    void validatePop2Object(){
        const string & method = "ForwardBasisIndexSeries::validatePop2Object";
        if (resetDates.empty()){
            throw ModelException(method, "No resetDates dates supplied");
        }
        if (evalDates.size() != resetDates.size()){
            throw ModelException(method, "The number of resetDates and evalDates dates"
                                 " must be the same");
        }
        if (valueDate > evalDates[0])
            throw ModelException(method, "resetDates cannot have dates in the past");
    }

    virtual void Validate(){
        const string & method = "ForwardBasisIndexSeries::validatePop2Object";
        // Check if fx is needed. If yes, check if it is provided
        if (discountYieldCurveName() != basisCurve->getRefCurve()->getName() &&
            fxAsset.isEmpty())
            throw ModelException(method, "Basis reference curve != product domestic "
            "discount curve.  Check if basis curve is foreign denominated.  If "
            "yes, make sure the corresponding fx asset was specified at the "
            "product level.");
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

private:
    /// helpers for MC -- creating generators ///
    //// Returns a SVGenDiscFactor for all the forwards
    SVGenDiscFactorSP createMCDF() const{
        return SVGenDiscFactorSP(new SVGenDiscFactor(valueDate, 
            discountCurve.getSP(), evalDates));
    }

    //// Returns a fx SVGenSpotSP for foreign denominated basis
    SVGenSpotSP createMCFX() const {
        return SVGenSpotSP(new SVGenSpot(1, evalDates));
    }

    //// Returns a MCExpectedDiscFactor for forward given by index idx
    SVGenExpectedDiscFactorSP createMCEDF(int idx) const{
        return SVGenExpectedDiscFactorSP(new SVGenExpectedDiscFactor(
            evalDates[idx], 
            evalDates[idx], 
            discountCurve.getSP(), 
            DateTimeArray(1, basisCurve->getRefPaymentDate(resetDates[idx])), 
            false));
    }

    //// Returns a MCIRSwap for the swap yield part of the basis coupon
    SVGenIRSwapSP createMCSwapYield(int idx) const {
        return createMCSwapYield(idx, 0); // no lag
    }
    SVGenIRSwapSP createMCSwapYield(int idx, int lag) const {
        return SVGenIRSwapSP(new SVGenIRSwap(
            basisCurve->getRefCurve().getSP(), // coupon curve
            basisCurve->getDiscountCurve().getSP(), /* Note: for basis, discount curve = reference curve */ // discount curve 
            evalDates[idx].rollDate(lag), // swap yield observation date
            resetDates[idx].rollDate(lag), // swap start date
            basisCurve->getRefPaymentDate(resetDates[idx].rollDate(lag)), // swap end date
            basisCurve->getBasisSwapLiborIvl()->toString(), // basisCurve->getRefSwapFixedIvl()->toString(), !?
            basisCurve->getBasisSwapLiborDcc()->toString(), // basisCurve->getRefSwapFixedDcc()->toString(), !?
            StubFactory::NONE, // stub type 
            false, // stub at end 
            basisCurve->getBadDayConvention()->toString(), // accural bad day conv
            basisCurve->getBadDayConvention()->toString(), // pay bad day conv
            HolidaySP((Holiday*)basisCurve->getHolidays().get()), 
            true )); // cash settle
    }

    //// Returns a MCForwardspread for forward given by index idx
    SVGenExpectedBasisFwdSpreadSP createMCEBF(int idx) const{
        return createMCEBF(idx, 0); // no lag
    }
    SVGenExpectedBasisFwdSpreadSP createMCEBF(int idx, int lag) const{
        return SVGenExpectedBasisFwdSpreadSP(new SVGenExpectedBasisFwdSpread(
            basisCurve.getSP(), 
            evalDates[idx].rollDate(lag), 
            DateTimeArray(1,resetDates[idx].rollDate(lag))));
    }


    ForwardBasisIndexSeries(): 
        CInstrument(TYPE), 
        notional(1.), 
        doSquare(false), 
        spreadOnly(false),
        dumpPrices(false), 
        dumpFileName("c:/temp/BasisPriceDump.txt"),
        doDeltaPricing(false),
        nbDayLag(1)
        {} // for reflection
    ForwardBasisIndexSeries(const ForwardBasisIndexSeries& rhs); // not implemented
    ForwardBasisIndexSeries& operator=(
        const ForwardBasisIndexSeries& rhs); // not implemented

    static IObject* defaultForwardBasisIndexSeries(){
        return new ForwardBasisIndexSeries();
    }

    /** Get the asset and discount market data */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market)
    {
        market->GetReferenceDate(valueDate);
        discountCurve.getData(model, market);
        basisCurve.getData(model, market);
        if (!fxAsset.isEmpty())
            fxAsset.getData(model, market);
    }

    virtual DateTime getValueDate()const
    {
        return valueDate;
    }

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const {
        return discountCurve.getName();
    }


    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ForwardBasisIndexSeries, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultForwardBasisIndexSeries);
        FIELD(valueDate,        "Valuation Date");
        FIELD(notional,         "Contract notional");
        FIELD_MAKE_OPTIONAL(notional);
        FIELD(basisCurve,        "Basis Curve Wrapper");
        FIELD(discountCurve,     "Discount Curve Wrapper");
        FIELD(fxAsset,           "fx asset of the foreign denominated basis curve");
        FIELD_MAKE_OPTIONAL(fxAsset);
        FIELD(evalDates,         "When each forward starts");
        FIELD(resetDates,        "When each forward ends");
        FIELD(doSquare, "Compute spread squared or not");
        FIELD_MAKE_OPTIONAL(doSquare);
        FIELD(spreadOnly, "Only report diffused spread");
        FIELD_MAKE_OPTIONAL(spreadOnly);
        FIELD(dumpPrices, "Dump simulated spreads & swap rates to text file");
        FIELD_MAKE_OPTIONAL(dumpPrices);
        FIELD(dumpFileName, "file name");
        FIELD_MAKE_OPTIONAL(dumpFileName);
        FIELD(doDeltaPricing, "If true, then output quantities as X(t+nbDayLag) - X(t)");
        FIELD_MAKE_OPTIONAL(doDeltaPricing);
        FIELD(nbDayLag, "The number of day lag in computing delta prices");
        FIELD_MAKE_OPTIONAL(nbDayLag);
    }
};

/* MC product class for super rainbow */
class ForwardBasisIndexSeriesMC : public MCProductClient {
    // a set of state variables for each option
    struct SV {
        vector<SVExpectedBasisFwdSpreadSP>  expFwd;   //!< expected fwd rate state variable
        vector<SVExpectedBasisFwdSpreadSP>  expFwdPlusLag; // for measuring delta spread
        vector<SVExpectedDiscFactorSP>      expDF;      //!< expected DF state var
        vector<SVGenIRSwap::IStateVarSP>    swapYield;   //!< swap yield state var
        vector<SVGenIRSwap::IStateVarSP>    swapYieldPlusLag; // for measuring delta yield
        SVDiscFactorSP                      df;         //!< df state variable
        MCPath::IStateVarSP                 fx;         // for foreign denominated basis
    };
    // a set of state variables generators for each set of state variables
    struct GenSV {
        vector<SVGenExpectedBasisFwdSpreadSP>  expFwd;   //!< Generator for expected fwd
        vector<SVGenExpectedBasisFwdSpreadSP>  expFwdPlusLag;
        vector<SVGenExpectedDiscFactorSP>      expDF;    //!< Gen for expected DF 
        vector<SVGenIRSwapSP>                  swapYield;    //!< Gen for swap yield 
        vector<SVGenIRSwapSP>                  swapYieldPlusLag;
        SVGenDiscFactorSP                      df;       //!< Generator for df
        SVGenSpotSP                            fx;      // for foreign denominated basis
    };
    const ForwardBasisIndexSeries*    inst;      //!< reference to original inst
    SV                                sv;        //!< state variables
    GenSV                             genSV;     //!< state variable generators
    vector<double>                    dcf;       //!< year fractions in accruals
    ofstream                          dumpFile;
    bool                              needFx;

private:
    /** Update our small collection of state variables */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen)
    {
        static const string routine("ForwardBasisIndexSeriesMC::pathGenUpdated");
        if (dynamic_cast<PastPathGen*>(newPathGen))
            return ;
        try {
            sv.df = genSV.df->getSVDiscFactor(newPathGen);
            if (needFx)
                sv.fx = genSV.fx->getSpotSV(newPathGen);

            sv.expFwd.clear();
            sv.expDF.clear();
            sv.swapYield.clear();
            for (size_t i = 0; i < genSV.expFwd.size(); ++i)
                sv.expFwd.push_back(genSV.expFwd[i]->getSVExpectedBasisFwdSpread(newPathGen));
            for (size_t i = 0; i < genSV.expDF.size(); ++i)
                sv.expDF.push_back(genSV.expDF[i]->getSVExpectedDiscFactor(newPathGen));
            for (size_t i = 0; i < genSV.swapYield.size(); ++i)
                sv.swapYield.push_back(genSV.swapYield[i]->getIRSwapSV(newPathGen));

            if (inst->doDeltaPricing) {
                sv.expFwdPlusLag.clear();
                sv.swapYieldPlusLag.clear();
                for (size_t i = 0; i < genSV.expFwd.size(); ++i)
                    sv.expFwdPlusLag.push_back(genSV.expFwdPlusLag[i]->getSVExpectedBasisFwdSpread(newPathGen));
                for (size_t i = 0; i < genSV.swapYield.size(); ++i)
                    sv.swapYieldPlusLag.push_back(genSV.swapYieldPlusLag[i]->getIRSwapSV(newPathGen));
            }
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    };

public:
    /** Ask for our small collection of state variables */
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        svCollector->append(genSV.df.get());
        if (needFx) 
            svCollector->append(genSV.fx.get());

        for (size_t i = 0; i < genSV.expFwd.size(); ++i)
            svCollector->append(genSV.expFwd[i].get());
        for (size_t i = 0; i < genSV.expDF.size(); ++i)
            svCollector->append(genSV.expDF[i].get());
        for (size_t i = 0; i < genSV.swapYield.size(); ++i)
            svCollector->append(genSV.swapYield[i].get());
        
        if (inst->doDeltaPricing) {
            for (size_t i = 0; i < genSV.expFwdPlusLag.size(); ++i)
                svCollector->append(genSV.expFwdPlusLag[i].get());
            for (size_t i = 0; i < genSV.swapYieldPlusLag.size(); ++i)
                svCollector->append(genSV.swapYieldPlusLag[i].get());            
        }
    }

    /** Need to call parent's constructor - a bit of a mess */
    ForwardBasisIndexSeriesMC(const ForwardBasisIndexSeries* inst,
                              SimSeriesSP                simSeries,
                              InstrumentSettlementSP     instSettle):
        MCProductClient(IMultiMarketFactors::asMulti(inst->discountCurve.getSP()).get(),
                        inst->valueDate,
                        inst->discountCurve.get(),
                        IRefLevelSP(IRefLevel::Util::makeZero(inst->valueDate)), // fix!
                        simSeries, // irrelevant - fix
                        IPastValuesSP(IPastValues::Util::makeTrivial(inst->valueDate, 0.0)), // fix
                        instSettle.get(),
                        inst->evalDates.back()),
        inst(inst), sv(), genSV(), needFx(!inst->fxAsset.isEmpty())
        {
            DayCountConventionConstSP basisDcc = inst->basisCurve->getBasisSwapBasisDcc();
            genSV.df = inst->createMCDF();
            if (needFx) 
                genSV.fx = inst->createMCFX();

            for (int i = 0; i < inst->resetDates.size(); ++i)
            {
                genSV.expFwd.push_back(inst->createMCEBF(i));
                genSV.expDF.push_back(inst->createMCEDF(i));
                genSV.swapYield.push_back(inst->createMCSwapYield(i));
                DateTime resetDate = inst->resetDates[i];
                DateTime payDate = inst->basisCurve->getRefPaymentDate(resetDate);
                dcf.push_back(basisDcc->years(resetDate, payDate)); //dcf.push_back(resetDate.yearFrac(payDate));

                if (inst->doDeltaPricing) {
                    genSV.expFwdPlusLag.push_back(inst->createMCEBF(i, inst->nbDayLag));
                    genSV.swapYieldPlusLag.push_back(inst->createMCSwapYield(i, inst->nbDayLag));
                }
            }
            
            if (inst->dumpPrices) {
                dumpFile.open(inst->dumpFileName.c_str()/*, ios_base::app*/);
                dumpFile.precision(10);
                dumpFile << "Idx:\tBasisSpreads:\tLn(SwapParYields):" << endl;
            }
        }

    /** Called within the simulation loop */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&             prices) 
    {
        double result, price = 0.0;
        double spread, swapParYield, swapAnnuity;
        double spreadPlusLag, swapParYieldPlusLag, swapAnnuityPlusLag;
        double spot2EvalDF, eval2PayoutDF;

        if (dynamic_cast<const PastPathGen*>(pathGen))
        {
            prices.add(0.0);
            return ;
        }

        for (size_t i = 0; i < sv.expFwd.size(); ++i) 
        {
            //  = EDF(T_eval, T_pay)*Expected_coupon(T_reset)*DCF(T_reset, T_pay)*DF(t_0, T_eval)
            
            //price += sv.expFwd[i]->firstSpread() * dcf[i] * sv.expDF[i]->firstDF() * sv.df->element(i);

            spot2EvalDF = sv.df->getDF(i);
            eval2PayoutDF = sv.expDF[i]->firstDF();
            spread = sv.expFwd[i]->firstSpread();
            sv.swapYield[i]->parYield(swapParYield, swapAnnuity);

            if (!inst->spreadOnly)
                result = (swapParYield - spread) * dcf[i] * eval2PayoutDF * spot2EvalDF; // basis coupon
            else
                result = spread;

            if (inst->doSquare)
                result *= result;

            price += result;

            if (inst->dumpPrices) {
                if (!inst->doDeltaPricing) {
                    dumpFile << "0" << "\t" << spread << "\t" << swapParYield;
                }
                else {
                    spreadPlusLag = sv.expFwdPlusLag[i]->firstSpread();
                    sv.swapYieldPlusLag[i]->parYield(swapParYieldPlusLag, swapAnnuityPlusLag);
                    dumpFile << "0" << "\t" << spreadPlusLag - spread << "\t" << swapParYieldPlusLag - swapParYield;
                }
            }
        }

        prices.add(price);

        if (inst->dumpPrices) {
            dumpFile << endl;
        }
    }

    virtual DateTimeArray getPastDates()
    {
        return DateTimeArray();
    }

};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* ForwardBasisIndexSeries::createProduct(const MonteCarlo* model) const {
     // the sim series here is irrelevant - to do: fix MCProductClient
    SimSeriesSP simSeries(new SimSeries(1));
    simSeries->addDates(DateTimeArray(1,
                                      resetDates.back())); // so getLastDate() works

    InstrumentSettlementSP instSettle(new CashSettlePeriod(0));
    return new ForwardBasisIndexSeriesMC(this, simSeries, instSettle);
}

CClassConstSP const ForwardBasisIndexSeries::TYPE = CClass::registerClassLoadMethod(
    "ForwardBasisIndexSeries", typeid(ForwardBasisIndexSeries), load);

// * for class loading (avoid having header file) */
bool ForwardBasisIndexSeriesLoad() {
    return (ForwardBasisIndexSeries::TYPE != 0);
}

DRLIB_END_NAMESPACE

