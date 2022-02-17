//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ASR.cpp
//
//   Description : ASR Products 
//
//   Author      : Bruno O Melka
//
//   Date        : 22 Nov 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/Maths.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/SVGenSpot.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/LatticeProdEDR.hpp"

#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/SampleList.hpp"

DRLIB_BEGIN_NAMESPACE

//**********************************************************************************************//
//*****   Override of Monte Carlo model that enables optimization for american exercise    *****//
//**********************************************************************************************//

class MonteCarloOpti : public MonteCarlo {
public:
    static CClassConstSP const TYPE;

    // registration, invoked when class is 'loaded'
    static void load(CClassSP& clazz){
        clazz->setPublic(); 
        REGISTER(MonteCarloOpti, clazz);
        SUPERCLASS(MonteCarlo);
        EMPTY_SHELL_METHOD(defaultMonteCarloOpti);
        FIELD(upperBound,"upper bound for optimization");
        FIELD_MAKE_OPTIONAL(upperBound);
        FIELD(lowerBound,"upper bound for optimization");
        FIELD_MAKE_OPTIONAL(lowerBound);
        FIELD(tolerance,"tolerance for optimization");
        FIELD_MAKE_OPTIONAL(tolerance);
        FIELD(numSubInterval,"parameter for optimization");
        FIELD_MAKE_OPTIONAL(numSubInterval);
        FIELD(optiNbIter,"number of iterations for optimization");
        FIELD_MAKE_TRANSIENT(optiNbIter);
        FIELD(isMax,"true if maximization");
        FIELD_MAKE_OPTIONAL(isMax);
        FIELD(doOpti,"true if optimization");
        FIELD_MAKE_OPTIONAL(doOpti);
    }

    static IObject* defaultMonteCarloOpti(){
        return new MonteCarloOpti();
    }

    // Override of the Price() function of MonteCarlo
    virtual void Price(CInstrument*  instrument,
                       Control*      control, 
                       CResults*     results);


    // constructor 
    MonteCarloOpti() : MonteCarlo (TYPE),
                        upperBound (1.0),
                        lowerBound(0.0),
                        tolerance (0.1),
                        numSubInterval (3),
                        optiNbIter (2),
                        isMax (true),
                        doOpti (false) {}

    // registered fields
    double    upperBound;
    double    lowerBound;
    double    tolerance;
    int        numSubInterval;
    int        optiNbIter;
    bool    isMax;
    bool    doOpti;    

protected:
    MonteCarloOpti(CClassConstSP clazz): MonteCarlo (clazz) {}

};

typedef smartPtr<MonteCarloOpti> MonteCarloOptiSP;

/** Override of the price function of MonteCarlo */
void MonteCarloOpti::Price(CInstrument*  instrument,
                       CControl*     control, 
                       CResults*     results)
{
    static const string method = "MonteCarloOpti::Price";
    try{
        auto_ptr<IMCProduct> prodAP(createProduct(instrument));
        prodAP->validate();
        // we call ourselves recursively
        bool isPricing = control->isPricing();
        if (isPricing){
            // instantiate MCPricing object
            createPricingObject(control, prodAP.get());
        }

        //********************      OPTIMIZATION PART       ********************//

        if (doOpti && isPricing) {

            // initialize parameters
            CControlSP ctrl(copy(control));
            ctrl->reset();
            CResultsSP rslts(new Results);
            int numDiv = (int)(log((upperBound - lowerBound) / tolerance) / log((double)numSubInterval / 2.0));
 
            // loop to get exercise boundary
            for (int i = 0; i < numDiv + 1; i++) {
                // optiNbIter = (int)((double)nbIter / (double)(numDiv + 1)) * (i + 1);
                // optiNbIter = optiNbIter - optiNbIter % (nbIter / nbSubSamples);
                optiNbIter = nbIter;
                MCPathGeneratorSP futurePathGenOpti(prodAP->price(this, ctrl.get(), rslts.get()));
            }
        }

        //**********************************************************************//

        doOpti = false; // turn off optimization
        MonteCarlo::Price(instrument, control, results);

    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Type registration */
CClassConstSP const MonteCarloOpti::TYPE =
    CClass::registerClassLoadMethod("MonteCarloOpti", typeid(MonteCarloOpti),
    MonteCarloOpti::load);


/**********************************************************************************************/
/*******************        ASR Instruments with 6 possible payoffs         *******************/
/**********************************************************************************************/

class ASR: public Generic1Factor,
                public virtual LastSensDate,
                public virtual FDModel::IIntoProduct,
                public virtual IMCIntoProduct {
public:

    static const string regular;
    static const string collared;
    static const string fixed;
    static const string holdback;
    static const string holdbackMod;
    static const string fixedCashAtMat;

    static CClassConstSP const TYPE; 
 
    virtual bool sensShift(Theta* shift){
        DateTime valueDateBeforeShift = valueDate;
        // value date is shifted by Generic1Factor::sensShift
        Generic1Factor::sensShift(shift);

        for (int i=0; i < monitorDates.size(); i++)
        {
            if( ( valueDateBeforeShift == assetHistory[i].date && Maths::isZero( assetHistory[i].amount ) ) ||
                ( valueDateBeforeShift < assetHistory[i].date && assetHistory[i].date <= valueDate ) )
            {
                assetHistory[i].amount = asset->getThetaSpotOnDate(shift, assetHistory[i].date);
            }
            else if( valueDate < assetHistory[i].date )
                break;
        }

        fillPastValues();
        
        return true;
    }

    /** checks inputs of the instrument */
    virtual void Validate() {
        static const string method = "ASR::Validate";
        try {
            if (type != regular && type != collared && type != fixed && type != holdback && type != holdbackMod && type != fixedCashAtMat) {
                throw ModelException(method,"type has to be either REGULAR, FIXEDDOLLAR, COLLARED, HOLDBACK, HOLDBACKMOD or FIXEDCASHATMAT");
            }
            if (monitorDates.size() != isDecisionDate.size()) {
                throw ModelException(method,"isDecisionDates must have same size as monitorDates");
            }
            if (monitorDates.size() != assetHistory.size()) {
                throw ModelException(method,"assetHistory must have same size as monitorDates");
            }
            if (!isDecisionDate[isDecisionDate.size() - 1]) {
                throw ModelException(method,"Last date must be decision date");
            }
            if ((type == regular || type == fixed || type == fixedCashAtMat) && (hedgeDays > 0)) {
                throw ModelException(method,"COLLARED, FIXEDDOLLAR and FIXEDCASHATMAT don't support hedge period");
            }
            if (type == collared || type == holdback || type == holdbackMod ) {
                int firstDecision = 0;
                while (!isDecisionDate[firstDecision]) {firstDecision++;}
                if (firstDecision <= hedgeDays) {
                    throw ModelException(method,"hedge period must end before exercizing starts");
                }
            }
            if ( (type != fixed  || type != fixedCashAtMat) && valueDate.getDate() > (startDate.getDate()) && Maths::isZero(initialSpot)) {
                throw ModelException(method,"When start date is before value date initialSpot must be provided");
            }

            //since we do division on 1 - spread in some cases
            if (isSpreadPerc && spread == 1.0){
                throw ModelException(method," spread should be smaller than 100% for percentage spread!");
            }
        }    

        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Returns last date */
    DateTime endDate(const Sensitivity* sensControl) const {
        return monitorDates[monitorDates.size() - 1];
    }

    /** Returns last date */
    DateTime getEndDate() const {
        return monitorDates[monitorDates.size() - 1];
    }

    /** Returns value date */
    DateTime getValueDate() const {
        return valueDate;
    }

    /** Get the asset and discount market data */
    void GetMarket(const IModel* model, const CMarketDataSP market) {

        market->GetReferenceDate(valueDate);
        CAsset::getAssetMarketData(model, market.get(), ccyTreatment, discount, asset);

        discount.getData(model, market);
        instSettle->getMarket(model, market.get());
        if(premiumSettle.get())
            premiumSettle->getMarket(model, market.get());

        /* keep it for use of assetHistoryWrapper later on
        if (valueDate > getStartDate()) {
           assetHistory.getData(model, market.get());
        }
        else {
            MarketObjectSP emptyAssetHist = MarketObjectSP(dynamic_cast<MarketObject *>(new AssetHistory("dummy",0)));
            assetHistory.setObject(emptyAssetHist);
        }
        CashFlowArraySP samples = CashFlowArraySP(new CashFlowArray(0));
        for (int i = 0; i < monitorDates.size(); ++i) {
            CashFlow dateOnly(monitorDates[i], 0);
            samples->push_back(dateOnly);
        }
        assetHistory.getSP()->getSamples(samples, valueDate);
        */

        fillPastValues();
    }

    /** Implementation of MonteCarlo::IntoProduct interface */
    IMCProduct* createProduct(const MonteCarlo* model) const;

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ASR, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultASR);
        FIELD(spread, "spread");
        FIELD(isSpreadPerc, "isSpreadPerc"); 
        FIELD_MAKE_OPTIONAL(isSpreadPerc);
        FIELD(ASRstrike, "ASRstrike"); 
        FIELD_MAKE_OPTIONAL(ASRstrike);
        FIELD(cash, "initial cash amount");
        FIELD_MAKE_OPTIONAL(cash);
        FIELD(initialShares, "initial number of shares");
        FIELD_MAKE_OPTIONAL(initialShares);
        FIELD(minCoeff, "determines minimum number of shares");
        FIELD_MAKE_OPTIONAL(minCoeff);
        FIELD(maxCoeff, "determines maximum number of shares");
        FIELD_MAKE_OPTIONAL(maxCoeff);
        FIELD(callStrike, "strike for the call");
        FIELD_MAKE_OPTIONAL(callStrike);
        FIELD(putStrike, "strike for the put");
        FIELD_MAKE_OPTIONAL(putStrike);
        FIELD(hedgeDays, "number of days in hedging period");
        FIELD_MAKE_OPTIONAL(hedgeDays);
        FIELD(monitorDates, "averaging schedule");
        FIELD_NO_DESC(pastValues);
        FIELD_MAKE_TRANSIENT(pastValues);
        FIELD(assetHistory, "Underlying sample history");
        FIELD(isDecisionDate, "averaging schedule");
        FIELD(type, "defines type of ASR products");
        FIELD_MAKE_OPTIONAL(type);
        FIELD(dimAvgOutFD, "number of slices in average dimension for FD");
        FIELD_MAKE_OPTIONAL(dimAvgOutFD);
        FIELD(dimAvgInFD, "dimension for hedging period in FD");
        FIELD_MAKE_OPTIONAL(dimAvgInFD);
        FIELD(stateStdev, "number of standard deviations for average state variables, default=3.5");
        FIELD_MAKE_OPTIONAL(stateStdev);
        FIELD(interpMethod, "interpolation method for state variable: QUADRATIC(default), LINEAR");
        FIELD_MAKE_OPTIONAL(interpMethod);
    }
     
    static IObject* defaultASR(){
        return new ASR();
    }

private:

    ASR() :
        Generic1Factor(TYPE),
        spread(0.0),
        isSpreadPerc(false),
        cash (1.0),
        initialShares(1.0), 
        minCoeff(1.0), 
        maxCoeff(1.0), 
        callStrike(0.0), 
        putStrike(0.0),            
        hedgeDays(0),
        type(regular),
        dimAvgOutFD(151),
        dimAvgInFD(1),
        ASRstrike(1.0), 
        stateStdev(3.5),
        interpMethod("QUADRATIC")  {}; 

    ASR(const ASR& rhs);
    ASR& operator=(const ASR& rhs);

    friend class ASRMC;
    friend class ASRMCSV;
    friend class ASRFD;
    friend class AvgState;

    void fillPastValues(){
        DoubleArray allValues(assetHistory.size());
        DateTimeArray allDates(assetHistory.size());
        for (int j = 0; j < assetHistory.size(); ++j) {
            allValues[j] = assetHistory[j].amount;
            allDates[j] = assetHistory[j].date;
        }
        DoubleMatrix allValuesAsMatrix(allValues);
        pastValues = IPastValuesSP(IPastValues::Util::makeSimple(allDates, allValuesAsMatrix));
    }

    double                  spread;            // spread
    bool                    isSpreadPerc;     //true: spread is % of Avg; False: spread is dollar amt
    double                  cash;            // cash amount received upfront for types 3 and 4
    double                  initialShares;    // initial number of shares delivered for type 4
    double                  minCoeff;        // determines minimum number of shares for type 4
    double                  maxCoeff;        // determines maximum number of shares for type 4
    double                  callStrike;        // strike for the call for type 2
    double                  putStrike;        // strike for the put for type 2
    int                     hedgeDays;        // number of days in hedging period for type 2 and 4
    DateTimeArray           monitorDates;    // averaging schedule
    IPastValuesSP           pastValues;        // all historic spots
    CashFlowArray           assetHistory;    // to be removed
    // AssetHistoryWrapper  assetHistory;    // Underlying sample history
    IntArray                isDecisionDate;    // record if dates are decision dates or not
    string                  type;            // determines type of ASR product
    int                     dimAvgOutFD;    // number of slices in average dimension for FD
    int                     dimAvgInFD;        // dimension for hedging period in FD
    double                  ASRstrike;      // Strike used only by REGULAR and COLLAR; 
                                                // given as PERCENTAGE IN BOTH FORWARD AND NONFORWARD STARTING  
    string                  interpMethod;   // interpolation method, 
    double                  stateStdev;    // number of standard deviations for state variables
};

CClassConstSP const ASR::TYPE = CClass::registerClassLoadMethod(
    "ASR", typeid(ASR), ASR::load);

const string ASR::regular    = "REGULAR";
const string ASR::collared    = "COLLARED";
const string ASR::fixed        = "FIXEDDOLLAR";
const string ASR::holdback    = "HOLDBACK";
const string ASR::holdbackMod    = "HOLDBACKMOD";
const string ASR::fixedCashAtMat = "FIXEDCASHATMAT";

/**********************************************************************************************/
/*********************           MC product class for ASR                **********************/
/**********************************************************************************************/

class ASRMC: public IMCProduct, virtual public IMCProductLN {

private:

    const ASR*            inst;
    double                initialSpot;
    double                ASRstrike; 
    double                hedgeAvgSoFar;
    double                hedgeSpot;
    double                averageSoFar;
    double                decisionCoeff;
    double                decisionRatio;
    int                   lastDecision;
    IntArray              decisionIndex;
    DoubleArray           pv;
    DoubleArrayArray      drifts;
    DoubleArrayArray      sumFwd;
    DoubleArray           priceArray;
    int                   numOpti;
    int                   numIter;

public:
    /** for the LogNormal path generator */
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGenerator,
                                    int                     iAsset) const {
        CVolRequestLNArray   reqarr(1); // one interp level/path per asset here
        reqarr[0] = CVolRequestLNSP(new ATMVolRequest());     
        return reqarr;
    }

    /** override the initial function. All discounting is done in payoff */
    double pvFromPaymentDate() const {
        return 1.0;
    }

    /** Use this opportunity to do any model driven initialisation
        of the instrument. e.g closed from barrier adjustments */
    virtual void initialiseLN(const MCPathGenerator* pathGenerator) const {}

    /** compute sum of drifts */
    void computeSumFwd(const DoubleArray& forwards, int num) {
        
        int idx = decisionIndex[num];
        // compute sum of drifts
        DoubleArray sumValues = DoubleArray(forwards.size() - idx - 1);
        DoubleArray fwd = DoubleArray(forwards.size() - idx - 1);
        double sum = 0.0;
        for (int i = 0; i < sumValues.size(); i++) {
            fwd[i] = forwards[idx+i+1] / forwards[idx];
            sum += fwd[i];
            sumValues[i] = sum;
        }
        // pick the ones for decicion dates
        drifts[idx].resize(decisionIndex.size() - num - 1);
        sumFwd[idx].resize(decisionIndex.size() - num - 1);
        for (int j = 0; j < sumFwd[idx].size(); j++) {
            drifts[idx][j] = fwd[decisionIndex[j + num + 1] - idx - 1];
            sumFwd[idx][j] = sumValues[decisionIndex[j + num + 1] - idx - 1];
        }
    }

    /** compute max condition for exercise checking */
    double getConditionMax(double average, double spot, double hedgeAvg, int num) {
        int idx = decisionIndex[num];
        double max;
        for (int i = 0; i < sumFwd[idx].size(); i++) {
            double expAverage = (average * (idx - inst->hedgeDays) + spot * sumFwd[idx][i]);
            expAverage /= (decisionIndex[num + i + 1] - inst->hedgeDays);
            double expPayoff = getPayoff(expAverage, hedgeAvg, decisionIndex[num + i + 1]);
            expPayoff /= (inst->type == ASR::fixed || inst->type == ASR::fixedCashAtMat || inst->type == ASR::holdback || inst->type == ASR::holdbackMod) ? drifts[idx][i] : 1.0;
            max = i == 0 ? expPayoff : Maths::max(max, expPayoff);
        }
        return max;
    }

    /** compute exercise payoff function */
    double getPayoff(double average, double hedgeAvg, int idx) {
        double payoff;
        payoff = inst->isSpreadPerc ? average * (1.0 - inst->spread): average - inst->spread;

        if (inst->type == ASR::collared) {
            double put = Maths::max(hedgeAvg * inst->putStrike - average, 0.0);
            double call = Maths::max(average - hedgeAvg * inst->callStrike, 0.0);             
            return (payoff - ASRstrike + put - call) * pv[idx]; 
        }
        else if (inst->type == ASR::fixed || inst->type == ASR::fixedCashAtMat || inst->type == ASR::holdback || inst->type == ASR::holdbackMod) {
            //return (average - inst->spread) / pv[idx];
            return payoff / pv[idx];
        }
        else {
            //return (average - ASRstrike - inst->spread) * pv[idx]; 
            return (payoff - ASRstrike) * pv[idx]; 
        }
    }

    /** compute final payoff function */
    double getFinalPayoff(double average, double hedgeAvg, double spot, double hedgeSpot, int idx) {
        
        double payoff;
        payoff = inst->isSpreadPerc ? average * (1.0 - inst->spread): average - inst->spread;
        const DateTimeArray& simDates = simSeries->getDates(0); //0 first asset

        if (inst->type == ASR::collared) {
            double put = Maths::max(hedgeAvg * inst->putStrike - average, 0.0);
            double call = Maths::max(average - hedgeAvg * inst->callStrike, 0.0);
            return (payoff - ASRstrike + put - call) * pv[idx]; 
        }
        else if (inst->type == ASR::fixed) {
            double myPayoff = 0.0;
            //if (simDates(idx).equals(inst->startDate, true)) {//initial cash flow happens at startDate
            if ( inst->monitorDates[0].isGreater(inst->valueDate)) {
                myPayoff = - inst->cash;
            }
            myPayoff  += inst->cash * spot * pv[idx] / payoff;
            return myPayoff;
        }
        else if (inst->type == ASR::fixedCashAtMat) {
            //return pv[idx] * ( inst->cash * spot / (average - inst->spread) - inst->cash);
            return pv[idx] * ( inst->cash * spot / (payoff) - inst->cash);
        }
        else if (inst->type == ASR::holdback || inst->type == ASR::holdbackMod) {
            double minShares = (inst->hedgeDays == 0) ? inst->cash / (inst->minCoeff * initialSpot) : inst->cash / (inst->minCoeff * hedgeAvg);
            double maxShares = (inst->hedgeDays == 0) ? inst->cash / (inst->maxCoeff * initialSpot) : inst->cash / (inst->maxCoeff * hedgeAvg);
            //double shares = inst->cash / (average - inst->spread);
            double shares = inst->cash / (payoff);

            double collar;
            double myPayoff = 0.0;
            if (inst->type == ASR::holdback){
                collar = Maths::max(shares - minShares, 0.0) - Maths::max(shares - maxShares, 0.0);
                collar *= spot * pv[idx];

//                if (simDates(idx).equals(inst->startDate, true)) {
//                    myPayoff = ((inst->hedgeDays == 0) ? minShares : inst->initialShares) * initialSpot * pv[0];
//                    myPayoff -= inst->cash *pv[0];
//                }

                if ( inst->monitorDates[0].isGreater(inst->valueDate)) {
                    myPayoff = ((inst->hedgeDays == 0) ? minShares : inst->initialShares) * initialSpot * pv[0];
                    myPayoff -= inst->cash * pv[0];
                }
            }
            else if (inst->type == ASR::holdbackMod){
               collar = Maths::max(shares - minShares, 0.0) - Maths::max(shares - maxShares, 0.0) - (inst->cash/initialSpot-minShares);
               collar *= spot * pv[idx];
            }

            if ( inst->monitorDates[inst->hedgeDays].isGreater(inst->valueDate)) {
                myPayoff += (inst->hedgeDays > 0) ? hedgeSpot * (minShares - inst->initialShares) * pv[inst->hedgeDays] : 0.0;
            }

            myPayoff +=  collar;
            return myPayoff;
        }
        else {
               //return (average - ASRstrike - inst->spread) * pv[idx]; 
               return (payoff - ASRstrike ) * pv[idx]; 
        }
    }
    
     /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    ASRMC(const ASR* inst, const SimSeriesSP& simSeries):
    IMCProduct(inst->asset.get(),
                inst->valueDate,
                inst->discount.get(),
                IRefLevelSP(IRefLevel::Util::makeZero(inst->valueDate)), // fix!
                simSeries, 
                inst->pastValues,
                inst->instSettle.get(),
                inst->monitorDates[inst->monitorDates.size() - 1]),
    inst(inst),
    //initialSpot((inst->valueDate.isGreater(inst->monitorDates[0])) ? inst->initialSpot : inst->asset->fwdValue(inst->monitorDates[0])),
    //ASRstrike((inst->monitorDates[0].isGreater(inst->valueDate)) ? (inst->asset->fwdValue(inst->monitorDates[0])*inst->ASRstrike) : (inst->initialSpot*inst->ASRstrike) ),

    //to review
    initialSpot((inst->valueDate.isGreater(inst->startDate)) ? inst->initialSpot : inst->asset->fwdValue(inst->startDate)),
    ASRstrike((inst->startDate.isGreater(inst->valueDate)) ? (inst->asset->fwdValue(inst->startDate)*inst->ASRstrike) : (inst->initialSpot*inst->ASRstrike) ),
    hedgeAvgSoFar(initialSpot),
    hedgeSpot(initialSpot),
    averageSoFar(0.0), 
    decisionCoeff(0.0),
    decisionRatio(1.0),
    pv(inst->isDecisionDate.size(), 1.0),
    drifts(inst->isDecisionDate.size()),
    sumFwd(inst->isDecisionDate.size()),
    priceArray(1, 0.0),
    numOpti(0),
    numIter(1) {

        // finds index of decision dates, compute pv for decision dates and forwards.
        DoubleArray forwards = inst->pastValues->getPastValues(inst->monitorDates, 0, inst->valueDate);
        int numPast = forwards.size();
        forwards.resize(inst->monitorDates.size());
        for (int i = 0; i < inst->isDecisionDate.size(); i++) {
            if (i > numPast - 1) {
                forwards[i] = inst->asset->fwdValue(inst->monitorDates[i]);
                pv[i] = inst->discount->pv(inst->valueDate, inst->monitorDates[i]);
            }
            if (inst->isDecisionDate[i]) {
                decisionIndex.push_back(i);
            }
        }
        
        // get last decision date
        if (decisionIndex.size() < 3) {            
            lastDecision = decisionIndex[0];
        }
        else {
            lastDecision = decisionIndex[decisionIndex.size() - 2];
            decisionRatio /= (double)(lastDecision - decisionIndex[0]);
        }

        /** Compute the expected sum of future drifts from a given point                        */ 
        /** where SumDrifts, Fwd and PV are computed in the MC constructor outside the loop        */    
        int num = 0;
        for (int k = 0; k < inst->isDecisionDate.size() - 1; k++) {
            if (inst->isDecisionDate[k]) {
                computeSumFwd(forwards, num);
                num++;
            }
        }
    }

    /** Called within the simulation loop */
    virtual void payoff(const IPathGenerator*    pathGen,
                        IMCPrices&                prices) {

        // compute averages
        double average = averageSoFar;
        double spot = initialSpot;
        double hedgeAvg = hedgeAvgSoFar;
        int num = 0;
        int iStep;

        for (iStep = pathGen->begin(0); iStep < pathGen->end(0); iStep++) {

            spot = pathGen->Path(0,0)[iStep];
            hedgeSpot = (iStep == inst->hedgeDays) ? spot : hedgeSpot;
            hedgeAvg = ((iStep <= inst->hedgeDays) && (iStep > 0)) ? (hedgeAvg * (iStep - 1) + spot) / iStep : hedgeAvg;
            average = (iStep > inst->hedgeDays) ? (average * (iStep - inst->hedgeDays - 1) + spot) / (iStep - inst->hedgeDays) : average;

            if (inst->isDecisionDate[iStep] && (iStep < inst->isDecisionDate.size() - 1) && !pathGen->doingPast()) {
                double scaling = decisionCoeff * (double)(lastDecision - iStep) * decisionRatio;
                double maxExpectedPayoff = getConditionMax(average, spot, hedgeAvg, num);
                double payoff = getPayoff(average, hedgeAvg, iStep);
                if (payoff > (maxExpectedPayoff * (1.0 + scaling))) {
                //if (average > (1.0 + decisionCoeff) * spot) {
                    break;
                }
                num++;
            }
        }

        if(pathGen->doingPast()) {
            averageSoFar = average;
            hedgeAvgSoFar = hedgeAvg;
        }

        // compute value of option.
        if (!pathGen->doingPast() || !hasFuture()) {
            double finalPayoff = getFinalPayoff(average, hedgeAvg, spot, hedgeSpot, decisionIndex[num]);
            prices.add(finalPayoff);
            priceArray[numOpti] += finalPayoff / numIter;
        }
    }

    /** IMCPrices the product. Returns the future path generator or null if there was no 'future' to price    */
    /** Apart from the OPTIMIZATION PART, this a copy of the price() function in MCProduct.                    */
    MCPathGeneratorSP price(MonteCarlo*       mcarlo,
                            Control*          control, 
                            Results*          results)
    {
        const static string routine("IMCProduct::price");
        MCPathGeneratorSP futurePathGenerator;
        try{

            bool isPricing = control->isPricing();
            OutputRequest* request;
            // get hold of the Pricing object
            MCPricing* pricing = mcarlo->getPricing().get();
            if (!pricing){ // internal error
                throw ModelException(routine, "Null pricing object");
            }
        
            // deal with historic dates
            IMCPricesSP prices; // note we pass this in empty and it might get filled
            MCPathGeneratorSP pastPathGenerator= handlePast(mcarlo, control, 
                                                             prices, true);       
            if (hasFuture()){ // sim needed
                // give payoff opportunity to establish future timeline before
                // the PathGen is configured for the simulation
                finaliseSimulationDates();
                // should we request the paths to be cached?
                bool cachePaths = isPricing && mcarlo->getCachedCachePaths() &&
                    !control->getSens()->empty();
                /* Some Pricing classes skip paths so the path generator needs
                   to provide 'random access' to the paths */
                int cacheMode = 
                    (pricing->needRandomAccessToPaths(mcarlo, this, control)?
                     IMCPathConfig::PATH_RANDOM_ACCESS: 0) |
                    (cachePaths? IMCPathConfig::PATH_CACHE: 0);
                // now get path generator for future dates
                futurePathGenerator = mcarlo->getPathConfig()->futurePathGenerator(
                    cacheMode, mcarlo->getNbIters(),
                    pastPathGenerator, this,
                    control, results, pricing->getTimeLine() );
                // switch doingThePast flag
                doingThePast = false;
                MCPathGenerator* pathGen = futurePathGenerator.get(); // for ease
                // tell the product that the generator has changed
                pathGenUpdated(pathGen);

                // where we store the results
                prices = pricing->createFuturePrices(mcarlo, this, 
                                                     control, futurePathGenerator);

                // If relevant allow configuration of product cache during tweaks
                if (!isPricing && mcarlo->getCachedCachePaths()) {
                    const IMultiFactors* multiFactor = getMultiFactors();
                    // then which greek we're doing
                    SensitivitySP sens(control->getCurrentSensitivity());
                    // getSensitiveAssets() does the real work
                    IntArray sensitivePhiAssets(multiFactor->getSensitiveAssets(
                        sens.get(), true)); // include phi etc
                    if (sensitivePhiAssets.empty()){
                        IntArray sensitiveAssets(multiFactor->getSensitiveAssets(
                            sens.get(), false)); // exclude phi etc
                        prices->configureCache(sensitiveAssets);
                    } 
                    else {
                        IntArray allAssets(multiFactor->NbAssets());
                        for (int i = 0; i < allAssets.size(); i++){
                            allAssets[i] = i;
                        }
                        prices->configureCache(allAssets);
                    }
                }    

                //********************      OPTIMIZATION PART       ********************//

                MonteCarloOpti* mcOpti = dynamic_cast<MonteCarloOpti*>(mcarlo);
                IMCPricesSP optPrices = prices;
                if (mcOpti) {
                    if (mcOpti->doOpti) {
                        numIter = mcOpti->optiNbIter;
                        priceArray = DoubleArray(mcOpti->numSubInterval + 1, 0.0);
                        double stepBound = (mcOpti->upperBound - mcOpti->lowerBound) / (double)(mcOpti->numSubInterval);
                        int optIdx;
                        double optPrice;
                        for (int idx = 0; idx < numIter; idx++) {
                            pathGen->generatePath(idx);
                            for (int i = 0; i < mcOpti->numSubInterval + 1; i++) {
                                numOpti = i;
                                decisionCoeff = mcOpti->lowerBound + stepBound * (double)(i);
                                payoff(pathGen, *optPrices);
                            }
                        }
                        for (int j = 0; j < mcOpti->numSubInterval + 1; j++) {
                            if ((j == 0) || (optPrice < priceArray[j] && mcOpti->isMax) || (optPrice > priceArray[j] && !mcOpti->isMax)) {
                                optPrice = priceArray[j];
                                optIdx = j;
                            }                            
                        }
                        mcOpti->upperBound = mcOpti->lowerBound + (optIdx + 1) * stepBound;
                        mcOpti->lowerBound = mcOpti->lowerBound + (optIdx - 1) * stepBound;
                    }
                    else {
                        decisionCoeff = (mcOpti->upperBound + mcOpti->lowerBound) / 2.0;
                        pricing->runSim(mcarlo, this, control, pathGen, *prices);
                    }
                }
                else {
                    pricing->runSim(mcarlo, this, control, pathGen, *prices);
                }

                //**********************************************************************//
            }
            // Pricing class is responsible for storing Price
            pricing->postResults(mcarlo, this, control, *prices,
                                 discount->getCcy(), results);
            // Allow product to store anything extra
            recordExtraOutput(control, results, *prices);
            if (isPricing){
                // save state of pathConfig if splitting into blocks
                if (mcarlo->getCachedCachePaths()){
                    mcarlo->initializePathConfigPostPricing();
                }
                mcarlo->setOrigPrices( prices ); // save for future use
                // do any events that we know how to do
                recordEvents(control, *prices, results);
            }
            if    ((isPricing) &&
                     (request = control->requestsOutput(OutputRequest::FWD_AT_MAT)) &&
                     (!request->getHasFinished())) {
                const DateTime& matDate = simSeries->getLastDate();
                if (matDate.isGreaterOrEqual(Today)){
                    /* do fwd @ mat by stock */
                    mfAsset->recordFwdAtMat(request, results, matDate);
                }
            }
        } catch (IPastValues::MissingSampleException& mse){
            mse.addMsg(routine);
            mse.setAssetName(mfAsset->getName(mse.getAssetIndex()));
            throw;
        } catch (ModelException& e){
            // check to see if MissingSampleException wrapped in a ModelException
            e.addMsg(routine);
            IPastValues::MissingSampleException* mse = 
                IPastValues::MissingSampleException::getInstance(e);
            if (mse){
                mse->setAssetName(mfAsset->getName(mse->getAssetIndex()));
            }
            throw;
        } catch (exception& e){
            throw ModelException(e, routine);
        }
        return futurePathGenerator;
    }

};


/**********************************************************************************************/
/***************       MC product class for ASR with State Variable            ****************/
/**********************************************************************************************/

class ASRMCSV: public MCProductClient, virtual public IMCProductLN {

private:

    SVGenSpot::IStateVarSP    assetSV;        // asset state variable
    SVGenSpotSP               assetGen;       // generator for asset

    const ASR*             inst;
    double                 initialSpot;
    double                 ASRstrike; 
    double                 hedgeAvgSoFar;
    double                 hedgeSpot;
    double                 averageSoFar;
    double                 decisionCoeff;
    double                 decisionRatio;
    int                    lastDecision;
    IntArray               decisionIndex;
    DoubleArray            pv;
    DoubleArrayArray       drifts;
    DoubleArrayArray       sumFwd;
    DoubleArray            priceArray;
    int                    numOpti;
    int                    numIter;

protected:
    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
        assetSV    = assetGen->getSpotSV(newPathGen);
    }


public:
    /** for the LogNormal path generator */
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGenerator,
                                    int                     iAsset) const {
        CVolRequestLNArray   reqarr(1); // one interp level/path per asset here
        reqarr[0] = CVolRequestLNSP(new ATMVolRequest());     
        return reqarr;
    }

    /** override the initial function. All discounting is done in payoff */
    double pvFromPaymentDate() const {
        return 1.0;
    }

    /** Use this opportunity to do any model driven initialisation
        of the instrument. e.g closed from barrier adjustments */
    virtual void initialiseLN(const MCPathGenerator* pathGenerator) const {}

    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        svCollector->append(assetGen.get()); 
    }

    /** compute sum of drifts */
    void computeSumFwd(const DoubleArray& forwards, int num) {
        
        int idx = decisionIndex[num];
        // compute sum of drifts
        DoubleArray sumValues = DoubleArray(forwards.size() - idx - 1);
        DoubleArray fwd = DoubleArray(forwards.size() - idx - 1);
        double sum = 0.0;
        for (int i = 0; i < sumValues.size(); i++) {
            fwd[i] = forwards[idx+i+1] / forwards[idx];
            sum += fwd[i];
            sumValues[i] = sum;
        }
        // pick the ones for decicion dates
        drifts[idx].resize(decisionIndex.size() - num - 1);
        sumFwd[idx].resize(decisionIndex.size() - num - 1);
        for (int j = 0; j < sumFwd[idx].size(); j++) {
            drifts[idx][j] = fwd[decisionIndex[j + num + 1] - idx - 1];
            sumFwd[idx][j] = sumValues[decisionIndex[j + num + 1] - idx - 1];
        }
    }

    /** compute max condition for exercise checking */
    double getConditionMax(double average, double spot, double hedgeAvg, int num) {
        int idx = decisionIndex[num];
        double max;
        for (int i = 0; i < sumFwd[idx].size(); i++) {
            double expAverage = (average * (idx - inst->hedgeDays) + spot * sumFwd[idx][i]);
            expAverage /= (decisionIndex[num + i + 1] - inst->hedgeDays);
            double expPayoff = getPayoff(expAverage, hedgeAvg, decisionIndex[num + i + 1]);
            expPayoff /= (inst->type == ASR::fixed || inst->type == ASR::fixedCashAtMat || inst->type == ASR::holdback || inst->type == ASR::holdbackMod) ? drifts[idx][i] : 1.0;
            max = i == 0 ? expPayoff : Maths::max(max, expPayoff);
        }
        return max;
    }

    /** compute exercise payoff function */
    double getPayoff(double average, double hedgeAvg, int idx) {

        double payoff;
        payoff = inst->isSpreadPerc ? average * (1.0 - inst->spread): average - inst->spread;

        if (inst->type == ASR::collared) {
            double put = Maths::max(hedgeAvg * inst->putStrike - average, 0.0);
            double call = Maths::max(average - hedgeAvg * inst->callStrike, 0.0);
            return (payoff  - ASRstrike + put - call) * pv[idx]; 
        }
        else if (inst->type == ASR::fixed || inst->type == ASR::fixedCashAtMat || inst->type == ASR::holdback || inst->type == ASR::holdbackMod) {
            return (payoff) / pv[idx];
        }
        else {
             return (payoff - ASRstrike ) * pv[idx];
        }
    }

    /** compute final payoff function */
    double getFinalPayoff(double average, double hedgeAvg, double spot, double hedgeSpot, int idx) {

        double payoff;
        payoff = inst->isSpreadPerc ? average * (1.0 - inst->spread): average - inst->spread;
        const DateTimeArray& simDates = simSeries->getDates(0); //0 first asset

        if (inst->type == ASR::collared) {
            double put = Maths::max(hedgeAvg * inst->putStrike - average, 0.0);
            double call = Maths::max(average - hedgeAvg * inst->callStrike, 0.0);
            return (payoff  - ASRstrike + put - call) * pv[idx]; 
        }
        else if (inst->type == ASR::fixed) {
            double myPayoff = 0.0;
//            if (simDates(idx).equals(inst->startDate, true)) {//initial cash flow happens at startDate
            if ( inst->monitorDates[0].isGreater(inst->valueDate)) {
                myPayoff = - inst->cash;
            }
            myPayoff  += inst->cash * spot * pv[idx] / payoff;
            return myPayoff;
        }
        else if (inst->type == ASR::fixedCashAtMat) {
            return pv[idx] * ( inst->cash * spot  / (payoff) - inst->cash);
        }
        else if (inst->type == ASR::holdback || inst->type == ASR::holdbackMod) {
            double minShares = (inst->hedgeDays == 0) ? inst->cash / (inst->minCoeff * initialSpot) : inst->cash / (inst->minCoeff * hedgeAvg);
            double maxShares = (inst->hedgeDays == 0) ? inst->cash / (inst->maxCoeff * initialSpot) : inst->cash / (inst->maxCoeff * hedgeAvg);
            //double shares = inst->cash / (average - inst->spread);
            double shares = inst->cash / (payoff);

            double collar;
            double myPayoff = 0.0;
            if (inst->type == ASR::holdback){
                collar = Maths::max(shares - minShares, 0.0) - Maths::max(shares - maxShares, 0.0);
                collar *= spot * pv[idx];

//                if (simDates(idx).equals(inst->startDate, true)) {
//                if ( inst->monitorDates[0].isGreater(inst->valueDate)) {
//                    myPayoff = ((inst->hedgeDays == 0) ? minShares : inst->initialShares) * initialSpot * pv[0];
//                    myPayoff -= inst->cash *pv[0];
//                }

                if ( inst->monitorDates[0].isGreater(inst->valueDate)) {
                    myPayoff = ((inst->hedgeDays == 0) ? minShares : inst->initialShares) * initialSpot * pv[0];
                    myPayoff -= inst->cash *pv[0];
                }
            }
            else if (inst->type == ASR::holdbackMod){
               collar = Maths::max(shares - minShares, 0.0) - Maths::max(shares - maxShares, 0.0) - (inst->cash/initialSpot-minShares);
               collar *= spot * pv[idx];
            }

            if (inst->monitorDates[inst->hedgeDays].isGreater(inst->valueDate)) {
                myPayoff += (inst->hedgeDays > 0) ? hedgeSpot * (minShares - inst->initialShares) * pv[inst->hedgeDays] : 0.0;
            }
            myPayoff +=  collar;
            return myPayoff;
        }
        else {
            //return (average - ASRstrike - inst->spread) * pv[idx]; 
            return (payoff - ASRstrike ) * pv[idx]; 
        }
    }
    
     /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    ASRMCSV(const ASR* inst, const SimSeriesSP& simSeries):
    MCProductClient(inst->asset.get(),
                inst->valueDate,
                inst->discount.get(),
                IRefLevelSP(IRefLevel::Util::makeZero(inst->valueDate)), // fix!
                simSeries, 
                inst->pastValues,
                inst->instSettle.get(),
                inst->monitorDates[inst->monitorDates.size() - 1]),
    inst(inst),
//    initialSpot((inst->monitorDates[0].isGreater(inst->valueDate)) ? inst->asset->fwdValue(inst->monitorDates[0]) : inst->initialSpot),
//    ASRstrike((inst->monitorDates[0].isGreater(inst->valueDate)) ? (inst->asset->fwdValue(inst->monitorDates[0])*inst->ASRstrike) : (inst->initialSpot*inst->ASRstrike) ), 
    //to review
    initialSpot((inst->valueDate.isGreater(inst->startDate)) ? inst->initialSpot : inst->asset->fwdValue(inst->startDate)),
    ASRstrike((inst->startDate.isGreater(inst->valueDate)) ? (inst->asset->fwdValue(inst->startDate)*inst->ASRstrike) : (inst->initialSpot*inst->ASRstrike) ),

    hedgeAvgSoFar(initialSpot),
    hedgeSpot(initialSpot),
    averageSoFar(0.0), 
    decisionCoeff(0.0),
    decisionRatio(1.0),
    pv(inst->isDecisionDate.size(), 1.0),
    drifts(inst->isDecisionDate.size()),
    sumFwd(inst->isDecisionDate.size()),
    priceArray(1, 0.0),
    numOpti(0),
    numIter(1),
    assetGen(new SVGenSpot(1, inst->monitorDates)) {

        // finds index of decision dates, compute pv for decision dates and forwards.
        DoubleArray forwards = inst->pastValues->getPastValues(inst->monitorDates, 0, inst->valueDate);
        int numPast = forwards.size();
        forwards.resize(inst->monitorDates.size());
        for (int i = 0; i < inst->isDecisionDate.size(); i++) {
            if (i > numPast - 1) {
                forwards[i] = inst->asset->fwdValue(inst->monitorDates[i]);
                pv[i] = inst->discount->pv(inst->valueDate, inst->monitorDates[i]);
            }
            if (inst->isDecisionDate[i]) {
                decisionIndex.push_back(i);
            }
        }
        
        // get last decision date
        if (decisionIndex.size() < 3) {            
            lastDecision = decisionIndex[0];
        }
        else {
            lastDecision = decisionIndex[decisionIndex.size() - 2];
            decisionRatio /= (double)(lastDecision - decisionIndex[0]);
        }

        /** Compute the expected sum of future drifts from a given point                        */ 
        /** where SumDrifts, Fwd and PV are computed in the MC constructor outside the loop        */    
        int num = 0;
        for (int k = 0; k < inst->isDecisionDate.size() - 1; k++) {
            if (inst->isDecisionDate[k]) {
                computeSumFwd(forwards, num);
                num++;
            }
        }
    }

    /** Called within the simulation loop */
    virtual void payoff(const IPathGenerator*    pathGen,
                        IMCPrices&                prices) {

        const SVPath& pathAsset = assetSV->path(0);

        // compute averages
        double average = averageSoFar;
        double spot = initialSpot;
        double hedgeAvg = hedgeAvgSoFar;
        int num = 0;
        int iStep;

        for (iStep = pathAsset.begin(); iStep < pathAsset.end(); iStep++) {

            spot = pathAsset[iStep];
            hedgeSpot = (iStep == inst->hedgeDays) ? spot : hedgeSpot;
            hedgeAvg = ((iStep <= inst->hedgeDays) && (iStep > 0)) ? (hedgeAvg * (iStep - 1) + spot) / iStep : hedgeAvg;
            average = (iStep > inst->hedgeDays) ? (average * (iStep - inst->hedgeDays - 1) + spot) / (iStep - inst->hedgeDays) : average;

            if (inst->isDecisionDate[iStep] && (iStep < inst->isDecisionDate.size() - 1) && !pathGen->doingPast()) {
                double scaling = decisionCoeff * (double)(lastDecision - iStep) * decisionRatio;
                double maxExpectedPayoff = getConditionMax(average, spot, hedgeAvg, num);
                double payoff = getPayoff(average, hedgeAvg, iStep);
                if (payoff > (maxExpectedPayoff * (1.0 + scaling))) {
                    break;
                }
                num++;
            }
        }

        if(pathGen->doingPast()) {
            averageSoFar = average;
            hedgeAvgSoFar = hedgeAvg;
        }

        // compute value of option.
        if (!pathGen->doingPast() || !hasFuture()) {
            double finalPayoff = getFinalPayoff(average, hedgeAvg, spot, hedgeSpot, decisionIndex[num]);
            prices.add(finalPayoff);
            priceArray[numOpti] += finalPayoff / numIter;
        }
    }

    /** IMCPrices the product. Returns the future path generator or null if there was no 'future' to price    */
    /** Apart from the OPTIMIZATION PART, this a copy of the price() function in MCProduct.                    */
    MCPathGeneratorSP price(MonteCarlo*       mcarlo,
                            Control*          control, 
                            Results*          results)
    {
        const static string routine("IMCProduct::price");
        MCPathGeneratorSP futurePathGenerator;
        try{

            bool isPricing = control->isPricing();
            OutputRequest* request;
            // get hold of the Pricing object
            MCPricing* pricing = mcarlo->getPricing().get();
            if (!pricing){ // internal error
                throw ModelException(routine, "Null pricing object");
            }
        
            // deal with historic dates
            IMCPricesSP prices; // note we pass this in empty and it might get filled
            MCPathGeneratorSP pastPathGenerator= handlePast(mcarlo, control, 
                                                             prices, true);
            if (hasFuture()){ // sim needed
                // give payoff opportunity to establish future timeline before
                // the PathGen is configured for the simulation
                finaliseSimulationDates();
                // should we request the paths to be cached?
                bool cachePaths = isPricing && mcarlo->getCachedCachePaths() &&
                    !control->getSens()->empty();
                /* Some Pricing classes skip paths so the path generator needs
                   to provide 'random access' to the paths */
                int cacheMode = 
                    (pricing->needRandomAccessToPaths(mcarlo, this, control)?
                     IMCPathConfig::PATH_RANDOM_ACCESS: 0) |
                    (cachePaths? IMCPathConfig::PATH_CACHE: 0);
                // now get path generator for future dates
                futurePathGenerator = mcarlo->getPathConfig()->futurePathGenerator(
                    cacheMode, mcarlo->getNbIters(),
                    pastPathGenerator, this,
                    control, results, pricing->getTimeLine() );
                // switch doingThePast flag
                doingThePast = false;
                MCPathGenerator* pathGen = futurePathGenerator.get(); // for ease
                // tell the product that the generator has changed
                pathGenUpdated(pathGen);

                // where we store the results
                prices = pricing->createFuturePrices(mcarlo, this, 
                                                     control, futurePathGenerator);

                // If relevant allow configuration of product cache during tweaks
                if (!isPricing && mcarlo->getCachedCachePaths()) {
                    const IMultiFactors* multiFactor = getMultiFactors();
                    // then which greek we're doing
                    SensitivitySP sens(control->getCurrentSensitivity());
                    // getSensitiveAssets() does the real work
                    IntArray sensitivePhiAssets(multiFactor->getSensitiveAssets(
                        sens.get(), true)); // include phi etc
                    if (sensitivePhiAssets.empty()){
                        IntArray sensitiveAssets(multiFactor->getSensitiveAssets(
                            sens.get(), false)); // exclude phi etc
                        prices->configureCache(sensitiveAssets);
                    } 
                    else {
                        IntArray allAssets(multiFactor->NbAssets());
                        for (int i = 0; i < allAssets.size(); i++) {
                            allAssets[i] = i;
                        }
                        prices->configureCache(allAssets);
                    }
                }    

                //********************      OPTIMIZATION PART       ********************//

                MonteCarloOpti* mcOpti = dynamic_cast<MonteCarloOpti*>(mcarlo);
                IMCPricesSP optPrices = prices;
                if (mcOpti) {
                    if (mcOpti->doOpti) {
                        numIter = mcOpti->optiNbIter;
                        priceArray = DoubleArray(mcOpti->numSubInterval + 1, 0.0);
                        double stepBound = (mcOpti->upperBound - mcOpti->lowerBound) / (double)(mcOpti->numSubInterval);
                        int optIdx;
                        double optPrice;
                        for (int idx = 0; idx < numIter; idx++) {
                            pathGen->generatePath(idx);
                            for (int i = 0; i < mcOpti->numSubInterval + 1; i++) {
                                numOpti = i;
                                decisionCoeff = mcOpti->lowerBound + stepBound * (double)(i);
                                payoff(pathGen, *optPrices);
                            }
                        }
                        for (int j = 0; j < mcOpti->numSubInterval + 1; j++) {
                            if ((j == 0) || (optPrice < priceArray[j] && mcOpti->isMax) || (optPrice > priceArray[j] && !mcOpti->isMax)) {
                                optPrice = priceArray[j];
                                optIdx = j;
                            }                            
                        }
                        mcOpti->upperBound = mcOpti->lowerBound + (optIdx + 1) * stepBound;
                        mcOpti->lowerBound = mcOpti->lowerBound + (optIdx - 1) * stepBound;
                    }
                    else {
                        decisionCoeff = (mcOpti->upperBound + mcOpti->lowerBound) / 2.0;
                        pricing->runSim(mcarlo, this, control, pathGen, *prices);
                    }
                }
                else {
                    pricing->runSim(mcarlo, this, control, pathGen, *prices);
                }

                //**********************************************************************//
            }
            // Pricing class is responsible for storing Price
            pricing->postResults(mcarlo, this, control, *prices,
                                 discount->getCcy(), results);
            // Allow product to store anything extra
            recordExtraOutput(control, results, *prices);
            if (isPricing){
                // save state of pathConfig if splitting into blocks
                if (mcarlo->getCachedCachePaths()){
                    mcarlo->initializePathConfigPostPricing();
                }
                mcarlo->setOrigPrices( prices ); // save for future use
                // do any events that we know how to do
                recordEvents(control, *prices, results);
            }
            if    ((isPricing) &&
                     (request = control->requestsOutput(OutputRequest::FWD_AT_MAT)) &&
                     (!request->getHasFinished())) {
                const DateTime& matDate = simSeries->getLastDate();
                if (matDate.isGreaterOrEqual(Today)){
                    /* do fwd @ mat by stock */
                    mfAsset->recordFwdAtMat(request, results, matDate);
                }
            }
        } catch (IPastValues::MissingSampleException& mse){
            mse.addMsg(routine);
            mse.setAssetName(mfAsset->getName(mse.getAssetIndex()));
            throw;
        } catch (ModelException& e){
            // check to see if MissingSampleException wrapped in a ModelException
            e.addMsg(routine);
            IPastValues::MissingSampleException* mse = 
                IPastValues::MissingSampleException::getInstance(e);
            if (mse){
                mse->setAssetName(mfAsset->getName(mse->getAssetIndex()));
            }
            throw;
        } catch (exception& e){
            throw ModelException(e, routine);
        }
        return futurePathGenerator;
    }

};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* ASR::createProduct(const MonteCarlo* model) const {
    // the simSeries passed to the IMCProduct is redundant
    // well not quite. The MC wants to know the last sim date
    SimSeriesSP simSeries(new SimSeries(1)); /* create empty one */
    simSeries->addDates(monitorDates);
    if (model->stateVarUsed()){
        return new ASRMCSV(this, simSeries);
    }
    return new ASRMC(this, simSeries);
}


/**********************************************************************************************/
/*********************           FD product class for ASR                **********************/
/**********************************************************************************************/
/************************ state variables operation *********************/
// shorter name for convenience
typedef TreeSliceLayer::StateSupport::InterpMethod AVG_INTERP;
// average state
class ASRFD;
class AvgState : public TreeSliceLayer::StateSupport{
public:
    // helper function
    static double calcAvg(const ASR* inst){
        DoubleArray pastClosings = inst->pastValues->getPastValues(inst->monitorDates, 0, inst->valueDate);
        double avg = 0.0;
        for (int i = 0; i < pastClosings.size(); i++) {
            avg += pastClosings[i];
        }
        if (pastClosings.size() > 0)
            avg /= (double)pastClosings.size();

        return avg;
    }

    AvgState(const ASRFD* prod, const ASR* inst, FDProductSP payoffIndex, int dim, AVG_INTERP method) :
        StateSupport( "AvgState", calcAvg(inst), method ),
        prod(prod),
        inst(inst),
        payoffIndex(payoffIndex),
        dim(dim),
        sampleCount(inst->monitorDates.size()) {}

    /** state variable support using slice layers, TreeSliceLayer::StateSupport:: */
    // update state contents
    void update(int step, const TreeSlice& s);
    // set grid levels
    void setGrid( int step, bool init = false );
    // compute grid updates at a reset step. This is for each slice node level.
    virtual void transition(
        const vector< const TreeSlice * > & srcSlices,
        const vector< double > & gridIn,
        vector< double > & gridOut ) const;

    const ASRFD*    prod;
    const ASR*      inst;
    FDProductSP     payoffIndex;
    int             dim;
    double          avgSoFar;
    int             firstResetStep;
    BoolArray       isResetStep;
    int             sampleCount;
};
typedef refCountPtr<AvgState > AvgStateSP;

/************************ end state variables operation *********************/
/////////// FDProduct class ////////
class ASRFD : public LatticeProdEDR, 
              virtual public IFDProductLN {
private:
    friend class AvgState;

    const ASR*        inst;
    BoolArray         stepExercise;
    bool              isFwd;
    double            initialSpot;
    double            ASRstrike; 
    
    // price slice
    TreeSliceSP value;
    
    // shorter for convenience
    AVG_INTERP interpMethod;

public: 
    AvgStateSP avgOutState;
//!!!    AvgStateSP avgInState;

    /** constructor */
    ASRFD(const ASR* inst, FDModel* m) : 
        LatticeProdEDR(m), 
        inst(inst)
//        isFwd(inst->monitorDates[0].isGreater(inst->valueDate)),
//        initialSpot(isFwd ? inst->asset->fwdValue(inst->monitorDates[0]) : inst->initialSpot),
//        ASRstrike(initialSpot*inst->ASRstrike)
    {
        //to review
        isFwd = inst->startDate.isGreater(inst->valueDate);
        initialSpot= isFwd ? inst->asset->fwdValue(inst->startDate) : inst->initialSpot;
        ASRstrike = initialSpot*inst->ASRstrike;

        //validate for ASRFD with hedge period
        if (inst->hedgeDays > 0 ) {            
            throw ModelException("ASRFD: ", "ASRFD can't handle heddge period for now.");
        }

        if( tree1f ){
            throw ModelException("ASRFD: ", "tree1f can't handle ASRFD yet, please use FD1D!");
        }

        if( tree1f )
            tree1f->setDiscountCurve(inst->discount.getSP());

        payoffIndex = model->createProduct(IProdCreatorSP(new
            IndexSpecEQ(inst->asset.getName(), inst->asset, inst->ccyTreatment)));

        if (inst->interpMethod == "QUADRATIC")
            interpMethod = TreeSliceLayer::StateSupport::INTERP_QUADRATIC;
        else if (inst->interpMethod == "LINEAR")
            interpMethod = TreeSliceLayer::StateSupport::INTERP_LINEAR;
        else
            throw ModelException("ASRFD: ", "unknown interpMethod: " + inst->interpMethod);

        // state variable support, in and out average should be different !!!
        avgOutState = AvgStateSP(new AvgState(this, inst, payoffIndex, inst->dimAvgOutFD, interpMethod));
//!!!        avgInState = AvgStateSP(new AvgState(this, inst, payoffIndex, inst->dimAvgInFD, interpMethod));
    }

    /** vanilla, quanto or struck */
    virtual string getCcyTreatment() const {return inst->ccyTreatment;}
    
    /** initialisation, called ONCE only before initModel() for each new model instance */
    virtual void init(Control*  control) const;

    /** set flags for averaging time steps,
        return the first step of the monitorDates after valueDate */
    int setStepFlags(BoolArray& stepFlags, const DateTimeArray& stepDates, const DateTimeArray& monitorDates);

    /** initialising and setting product variables this is called per pricing call before each pricing */
    virtual void initProd();

    /** isInitValue == true, payoff at T for backward or value at t=0 for fwd induction
        isInitValue == false, payoff boundary condition, for KO, early exercise etc. */
    virtual void update(int& step, FDProduct::UpdateType type);

    /** scale for notional */ 
    double scalePremium(const double& fairValue, YieldCurveConstSP disc);

    /** output prices and any additional results */
    virtual void recordOutput(Control* control, YieldCurveConstSP disc, Results* results);

    /** local method, product payoff method at maturity */
    void prod_BWD_T(const TreeSlice & s, int step);

    /** local method, this is payoff boundary condition, for KO, early exercise etc */
    void prod_BWD(const TreeSlice & s, int step);
    
    /** Returns start date */
    virtual DateTime getStartDate() const {
//        return isFwd ? inst->monitorDates[0] : inst->valueDate;
//to review
        return isFwd ? inst->startDate : inst->valueDate;
    }

    /** for the LogNormal model */
    CVolRequestLNSP getVolInterp(int iAsset) const {
        // get strike and maturity date from instrument
        DateTime matDate = inst->getEndDate();
        //double volStrike = (isFwd) ? (ASRstrike / inst->asset->fwdValue(inst->monitorDates[0])) : ASRstrike;
        //to review
        double volStrike = (isFwd) ? (ASRstrike / inst->asset->fwdValue(inst->startDate)) : ASRstrike;

        DateTime imntStartDate = getStartDate();
        CVolRequestLNSP   reqarr;         
        reqarr = CVolRequestLNSP(new LinearStrikeVolRequest(volStrike, imntStartDate, matDate, isFwd));        
        return reqarr;
    }
};

void ASRFD::update(int& step, FDProduct::UpdateType type)
{
    const TreeSlice & s = payoffIndex->getValue(step);

    if( type == FDProduct::BWD_T )
    {
        // set starting grid
        avgOutState->setGrid( step, true );
//!!!     avgInState->setGrid( step, true );

        prod_BWD_T(s, step);
    }
    else if( type == FDProduct::BWD )
        prod_BWD(s, step);

    avgOutState->update(step, s);

    //avgInState->update(step, s);
    }

/** initialise tree1f - allow product customisation 
    must not init product variables here, use initProd() instead */
void ASRFD::init(CControl* control) const{

    static const string method = "ASRFD::init()";
    try {
        if (inst->dimAvgOutFD < 2) {
            throw ModelException(method,"dimAvgOutFD must ba at least 2 when pricing with tree");
        }
        if ((inst->hedgeDays > 0) && (inst->dimAvgInFD < 1)) {
            throw ModelException(method,"dimAvgInFD must be at least 1 when pricing with tree and hedge period");
        }
        // customize tree parameters here and set up the tree
        DateTimeArray segDates;
        segDates.resize(2);
        segDates[0] = getStartDate();
        segDates[1] = inst->getEndDate();
        IntArray density(1,1);

        // all averaging dates are copied to critical dates
        DateTimeArraySP critDates(new DateTimeArray(inst->monitorDates));
        // add critical dates
        model->addCritDates(*critDates);

        // prepare model set up
        model->initSegments(segDates, density);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** set flags for averaging time steps,
    return the first step of the monitorDates after valueDate */
int ASRFD::setStepFlags(BoolArray& stepFlags, const DateTimeArray& stepDates, const DateTimeArray& monitorDates) {

    static const string method = "ASRFD::setStepAverage";
    try {
        int iMonitor = 0;
        while (monitorDates[iMonitor] < stepDates[0]) {
            iMonitor++;
        }

        int firstStep = -1;

        for (int iStep = 0; iStep < stepDates.size(); iStep++) {
            if (stepDates[iStep].equals(monitorDates[iMonitor], true)) {
                stepFlags[iStep] = true;
                if( firstStep < 0 ) firstStep = iStep;
                iMonitor++;
            }
            else {
                stepFlags[iStep] = false;
            }
        }

        return firstStep;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** initialising and setting product variables
    this is called per pricing call before tree sweep call (after InitTree)  */
void ASRFD::initProd()
{
    static const string method = "ASRFD::initProd()";
    try {
        //validate against forward except for regular ASR
        //forward start here is defined as monitorDates[0] > value date, not including time
//        if ( (inst->type != ASR::regular) && (inst->monitorDates[0].getDate() > inst->valueDate.getDate())){
//            throw ModelException(method,"Only Regular ASR supports forward start.");
//        }

//to review
        //set default value for startDate = avg out date[0]
        //if  avg out date[0] > valueDate (not including time), user need to input startDate
        if ( (inst->type != ASR::regular) && (inst->startDate.getDate() > inst->valueDate.getDate())){
            throw ModelException(method,"Only Regular ASR supports forward start, please check the startDate.");
        }

        int lastStep = model->getLastStep();
        stepExercise.resize(lastStep + 1);
        DateTimeArray exerciseDates;
        for (int i = 0; i < inst->isDecisionDate.size(); i++) {
            if (inst->isDecisionDate[i]) {
                exerciseDates.push_back(inst->monitorDates[i]);
            }
        }
        setStepFlags(stepExercise, model->getDates(), exerciseDates);

        model->registerStateSupport( avgOutState.get() );

        startDEV( value = model->createSlice( inst->discount->getName() ) ); 
        // for avgOut
        avgOutState->isResetStep.resize(lastStep + 1);
        avgOutState->firstResetStep =
            setStepFlags(avgOutState->isResetStep, model->getDates(), inst->monitorDates);
        if(model->getDate(0) == inst->valueDate) // value date is already counted as history
            avgOutState->isResetStep[0] = false;
        //!!! to do: for avg in state
        //!!! avgInState->isResetStep.resize(lastStep + 1);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** premium scaling */
double ASRFD::scalePremium(const double& fairValue, YieldCurveConstSP disc) {

    double fwdAtStart = 0.0;
    double fwdStartDF = 1.0;
    if (isFwd) {
        fwdAtStart = initialSpot;
        fwdStartDF = disc->pv(inst->valueDate, model->getDate(0));
    }
    double scalingFactor = InstrumentUtil::scalePremium(
                                    inst->oneContract,
                                    isFwd,
                                    inst->notional,
                                    fwdAtStart,
                                    initialSpot);
    {
        if (inst->type == ASR::regular || inst->type == ASR::collared)
            return (fairValue*scalingFactor*fwdStartDF);
        else
            return (-fairValue*scalingFactor*fwdStartDF);
    }
}

/** output results */
void ASRFD::recordOutput(Control* control, YieldCurveConstSP disc, Results* results) {

    double price = model->getPrice0( *value );
    double scaledPrice = scalePremium( price, disc );
    results->storePrice( scaledPrice, disc->getCcy() );

    // take care of additional outputs
    if (control && control->isPricing()) {

        DateTime matDate = inst->getEndDate();

        OutputRequest* request = NULL;
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTime paymentDate = inst->instSettle->settles(matDate, inst->asset.get());
            DateTimeArray date(1, paymentDate);
            OutputRequestUtil::recordPaymentDates(control,results,&date); 
        }
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            // if we're past the last sample then we know the payoff
            // given fv (fv = pv*c) => back out final payment as we have
            // all the data to hand
            //DateTime end = avgOutState->getLastDate();
            if (inst->valueDate.isGreater(matDate)) {
                DateTime paymentDate = inst->instSettle->settles(matDate, inst->asset.get());
                double pv = disc->pv(paymentDate);
                CashFlow cf(paymentDate, scaledPrice/pv);
                CashFlowArray cfl(1, cf);
                OutputRequestUtil::recordKnownCashflows(control,
                                                        results,
                                                        disc->getCcy(),
                                                        &cfl);   
            }
        }


        double        indVol;
        // calculate indicative vol
        if (matDate.isGreater(inst->valueDate)) {
            DateTime imntStartDate = getStartDate();

            // get vol request
            double volStrike  = ASRstrike; 
            LinearStrikeVolRequest volRequest(volStrike, imntStartDate, matDate, isFwd);

            try {
                // interpolate the vol
                CVolProcessedSP  vol(inst->asset->getProcessedVol(&volRequest));
                // cast to the BS vol we're expecting
                CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);
                // this should never happen if our get market data has worked properly
                if (!volBS) {
                    throw ModelException("ASRFD::recordOutput", "No Black Scholes Vol");
                }
                // calculate the indicative vol
                indVol = volBS->CalcVol(imntStartDate, matDate);
            }
            catch (exception& ) {
                indVol = 0.0;
            }
        }
        else {
            indVol = 0.0;
        }
    }
}

/** local method, product payoff method at maturity */
void ASRFD::prod_BWD_T(const TreeSlice & s, int step)
{
    static const string method = "ASRFD::prod_BWD_T";
    try {
        double settlementPV = inst->instSettle->pvAdjust(model->getDate(step), inst->discount.get(), inst->asset.get());

		const TreeSlice & outGrid = avgOutState->getGridSlice();
//!!!		const TreeSlice & inGrid = avgInState->getGridSlice();

        // compute payoff at maturity for each type
        if (inst->type == ASR::regular) {
            if (inst->isSpreadPerc){
                *value = settlementPV * (outGrid * ( 1  - inst->spread )- ASRstrike); 
            }else{
                *value = settlementPV * (outGrid - inst->spread - ASRstrike ); 
            }
        }
        else if (inst->type == ASR::collared) {
            if (inst->hedgeDays == 0){
                if (inst->isSpreadPerc){
                    *value = settlementPV * (outGrid *( 1.0 - inst->spread) - ASRstrike 
                            + smax(initialSpot * inst->putStrike - outGrid, 0.0) // put
                            - smax(outGrid - initialSpot * inst->callStrike, 0.0)); // call 
                }else{
                    *value = settlementPV * (outGrid - inst->spread - ASRstrike 
                            + smax(initialSpot * inst->putStrike - outGrid, 0.0) // put
                            - smax(outGrid - initialSpot * inst->callStrike, 0.0)); // call 
                }
            }
                    
/*//!!!
            else
                *value = settlementPV * (outGrid - initialSpot - inst->spread
                        + smax(inGrid * inst->putStrike - outGrid, 0.0) // put
                        - smax(outGrid - inGrid * inst->callStrike, 0.0)); // call
*/
        }
        else if (inst->type == ASR::fixed) {
            if (inst->isSpreadPerc){
                *value =
                    settlementPV * s * inst->cash / (outGrid * (1.0 - inst->spread));
            }else{
                *value =
                    settlementPV * s * inst->cash / (outGrid - inst->spread);
            }
        }
        else if (inst->type == ASR::fixedCashAtMat) {
            if (inst->isSpreadPerc){
                *value =
                    settlementPV * (s * inst->cash / (outGrid * (1.0 - inst->spread)) - inst->cash);
            }else{
                *value =
                    settlementPV * (s * inst->cash / (outGrid - inst->spread) - inst->cash);
            }
        }
        else if (inst->type == ASR::holdback || inst->type == ASR::holdbackMod) {
            double minShares = inst->cash / (inst->minCoeff * initialSpot);
            double maxShares = inst->cash / (inst->maxCoeff * initialSpot);

            // for performance

            if (inst->type == ASR::holdback){ 
                if (inst->isSpreadPerc){
                    #define shares ( inst->cash / ( outGrid * (1.0 - inst->spread ) ))
                    *value =
                        settlementPV * s * (smax(shares - minShares, 0.0) - smax(shares - maxShares, 0.0));
                    #undef shares
                }else{
                    #define shares ( inst->cash / ( outGrid - inst->spread ) )
                    *value =
                        settlementPV * s * (smax(shares - minShares, 0.0) - smax(shares - maxShares, 0.0));
                    #undef shares
                }
            }
            else if (inst->type == ASR::holdbackMod){
                if (inst->isSpreadPerc){
                    #define shares ( inst->cash / ( outGrid * (1.0 - inst->spread ) ))
                    *value =
                        settlementPV * s * (smax(shares - minShares, 0.0) - smax(shares - maxShares, 0.0) - (inst->cash/initialSpot - minShares));
                    #undef shares
                }else{
                    #define shares ( inst->cash / ( outGrid - inst->spread ) )
                    *value =
                        settlementPV * s * (smax(shares - minShares, 0.0) - smax(shares - maxShares, 0.0) - (inst->cash/initialSpot - minShares));
                    #undef shares
                }
            }
        }
    } 
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** local method, product payoff method at steps earlier than maturity */
void ASRFD::prod_BWD(const TreeSlice & s, int step)
{
    static const string method = "ASRFD::prod_BWD";
    try {
        double settlementPV = inst->instSettle->pvAdjust(model->getDate(step), inst->discount.get(), inst->asset.get());

		const TreeSlice & outGrid = avgOutState->getGridSlice();
//!!!		const TreeSlice & inGrid = avgInState->getGridSlice();

        // check for avg out

//        if (avgOutState->isResetStep[step])
//        {
//            int dimHedge = inst->dimAvgInFD;
//
//            // for holdback type add shares at the end of hedge period
//            // for holdbackMod type: shares and cash net out to 0 at the end of hedge period 
//            if (inst->type == ASR::holdback) {
//                /*
//                if (countAvg == inst->hedgeDays && inst->hedgeDays != 0) {
//                    for (int j = pStart; j <= pEnd; ++j) {
//                        double minShares = inst->cash / (inst->minCoeff * hedgeAvg[j / dimHedge]);
//                        price[j] = price[j] + s * (minShares - inst->initialShares);
//                    }
//                }
//                */
//                
//                if (avgOutState->sampleCount == 1 && (inst->monitorDates[0].isGreater(inst->valueDate))) {
//                    double initial = (inst->hedgeDays == 0) ? inst->cash / inst->minCoeff : initialSpot * inst->initialShares;
//                    initial -= inst->cash ;
//                    *value += initial;
//                }
//            }    
//            else if (inst->type == ASR::fixed) { 
//                if (avgOutState->sampleCount == 1 && (inst->monitorDates[0].isGreater(inst->valueDate))) {             
//                    *value -= inst->cash;
//                }
//            }
//
//        }

        // for holdback type add shares at the start Date
        // for holdbackMod type: shares and cash net out to 0 at the start Date 

        if (model->getDate(step).equals(inst->startDate, true)) {
//          &&  (inst->valueDate.equals(inst->startDate, true)) ){ // initial amt is already counted as history
            if (inst->type == ASR::holdback) {
                double initial = (inst->hedgeDays == 0) ? inst->cash / inst->minCoeff : initialSpot * inst->initialShares;
                initial -= inst->cash ;
                *value += initial;
            }
            else if (inst->type == ASR::fixed) { 
                *value -= inst->cash;
            }
        }

        // computes exercise price for american options
        if( ! stepExercise[step] )
            return;

        if (inst->type == ASR::regular) {
            if (inst->isSpreadPerc){
                *value = smax(*value, settlementPV * (outGrid * (1.0 - inst->spread) - ASRstrike)); 
            }else{
                *value = smax(*value, settlementPV * (outGrid  - inst->spread - ASRstrike)); 
            }
        }
        else if (inst->type == ASR::collared) {
            if (inst->hedgeDays == 0){
                if (inst->isSpreadPerc){
                    *value =smax(*value, 
                            settlementPV * (outGrid * (1.0 - inst->spread ) - ASRstrike 
                            + smax(initialSpot * inst->putStrike - outGrid, 0.0) // put
                            - smax(outGrid - initialSpot * inst->callStrike, 0.0))); // call 
                }else{
                    *value =smax(*value, 
                            settlementPV * (outGrid - inst->spread - ASRstrike 
                            + smax(initialSpot * inst->putStrike - outGrid, 0.0) // put
                            - smax(outGrid - initialSpot * inst->callStrike, 0.0))); // call 
                }
            }

/*//!!!
            else
                *value =smax(*value, 
                        settlementPV * (outGrid - initialSpot - inst->spread
                        + smax(inGrid * inst->putStrike - outGrid, 0.0) // put
                        - smax(outGrid - inGrid * inst->callStrike, 0.0))); // call
*/
        }
        // below two cases we take smin because the holder wants to minimise the value
        else if (inst->type == ASR::fixed) {
            if (inst->isSpreadPerc){
                *value = smin(*value,
                    settlementPV * s * inst->cash / (outGrid * (1.0 - inst->spread)));
            }else{
                *value = smin(*value,
                    settlementPV * s * inst->cash / (outGrid - inst->spread));
            }
        }
        else if (inst->type == ASR::fixedCashAtMat) {
            if (inst->isSpreadPerc){
                *value = smin(*value,
                    settlementPV * (s * inst->cash / (outGrid* (1.0 - inst->spread)) - inst->cash));

            }else{
                *value = smin(*value,
                    settlementPV * (s * inst->cash / (outGrid - inst->spread) - inst->cash));
            }
        }
        else if (inst->type == ASR::holdback || inst->type == ASR::holdbackMod) {
            double minShares = inst->cash / (inst->minCoeff * initialSpot);
            double maxShares = inst->cash / (inst->maxCoeff * initialSpot);

            // for performance
            if (inst->type == ASR::holdback){
                if (inst->isSpreadPerc){
                    #define shares ( inst->cash / ( outGrid * (1.0 - inst->spread ) ) )
                    *value = smin(*value,
                        settlementPV * s * (smax(shares - minShares, 0.0) - smax(shares - maxShares, 0.0)));
                    #undef shares
                }else{
                    #define shares ( inst->cash / ( outGrid - inst->spread ) )
                    *value = smin(*value,
                        settlementPV * s * (smax(shares - minShares, 0.0) - smax(shares - maxShares, 0.0)));
                    #undef shares
                }
            }
            else if (inst->type == ASR::holdbackMod){
                if (inst->isSpreadPerc){
                    #define shares ( inst->cash / ( outGrid *(1.0- inst->spread) ) )
                    *value = smin(*value,
                        settlementPV * s * (smax(shares - minShares, 0.0) - smax(shares - maxShares, 0.0) - (inst->cash/initialSpot - minShares)));
                    #undef shares
                }else{
                    #define shares ( inst->cash / ( outGrid - inst->spread ) )
                    *value = smin(*value,
                        settlementPV * s * (smax(shares - minShares, 0.0) - smax(shares - maxShares, 0.0) - (inst->cash/initialSpot - minShares)));
                    #undef shares
                }
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** create a fd payoff product */
FDProductSP ASR::createProduct(FDModel* model) const {
    return FDProductSP(new ASRFD(this, model));
}

/************************ state operation methods *********************/
/////////////// avg out operations ////////////
void AvgState::setGrid( int step, bool init )
{
    if( init || step != firstResetStep )
    {
        // computing average range using stateStdev
        DateTimeArray sampleDates(inst->monitorDates);
        sampleDates.resize(sampleCount);
        DoubleArray values = inst->pastValues->getPastValues(sampleDates, 0, inst->valueDate);
        DoubleArray weights(sampleCount, 1.0);
        values.resize(sampleCount);
        SampleList samples(sampleDates, values, weights);
        // compute avg varince
        ATMVolRequest volReq;
        CVolProcessedBSSP  volBS(inst->asset->getProcessedVol(&volReq));
        double avgVariance = samples.averageVariance(volBS.get(),
                                                    inst->valueDate,
                                                    true);
        double fwd = samples.expectedAverage(inst->asset.get(), inst->valueDate);
        double maxAvg = fwd*exp(-0.5*avgVariance + inst->stateStdev*sqrt(avgVariance));
        double minAvg = fwd*exp(-0.5*avgVariance - inst->stateStdev*sqrt(avgVariance));
      
        int n = dim; // *sqrt((double)sampleCount/inst->monitorDates.size());
        n = Maths::max(n, 2);
        // build a log array for average
        double increment = exp((log(maxAvg) - log(minAvg))/(n - 1.0));

        currGrid.resize(n);
        currGrid[0] = minAvg;
        for (int j = 1; j < n; ++j) {
            currGrid[j] = currGrid[j-1] * increment;
        }
    }
    else
    {
        // collapse grid to dim=1
        currGrid.resize(1); // IMPORTANT: currGrid.resize(1, todayValue) does not work because resize down
        currGrid[0] = todayValue;
    }

    if( init )
        prevGrid = currGrid;

    populateGrid( init ? step : step - 1 );
}

// compute grid transition through a reset/monitoring
void AvgState::transition(
    const vector< const TreeSlice * > & srcSlices,
    const vector< double > & gridIn,
    vector< double > & gridOut ) const
{
    int nbGrid = gridOut.size();
    ASSERT( nbGrid == gridIn.size() );

    double s = srcSlices[ 0 ]->calc();

    for( int i = 0; i < nbGrid; ++i )
        gridOut[ i ] = ( gridIn[ i ] * ( sampleCount - 1 ) + s ) / sampleCount;
}

// interpolation etc.
void AvgState::update(int step, const TreeSlice& s)
{
    if( isResetStep[step] )
    {
        setGrid( step );
        vector< const TreeSlice * > srcSlices( 1, &s );
        interpolate( step, srcSlices );
        --sampleCount;
    }
}

/** for class loading */ 
bool ASRLoad() {
    return (ASR::TYPE != 0);
}

DRLIB_END_NAMESPACE
