//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CorpAct.cpp
//
//   Description : CorpAct instrument
//
//   Author      : Qing Hou
//
//   Date        : 05 Feb 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Average.hpp"
#include "edginc/AtMaturity.hpp"
#include "edginc/Black.hpp"

DRLIB_BEGIN_NAMESPACE

/** CorpAct instrument, also called Risk Arb instrument
 */
#ifdef EDR_RISKARB_USE_GENERICNFBASE
class CorpAct: public GenericNFBase, 
#else
    class CorpAct: public GenericNFactor,
#endif
                   public CClosedFormLN::IIntoProduct 
{
public:
    static CClassConstSP const TYPE;
 
    /** input data validation */
    void validatePop2Object();

    /** Implementation of ClosedForm::IntoProduct interface */
    CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;
  
    DateTime endDate(const Sensitivity* sensControl) const;

    //// roll through time (setting historic values)
    bool sensShift(Theta* theta);

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
        returns true if it is dead (and priced), false if it is not dead */
    bool priceDeadInstrument(CControl* control, CResults* results) const;


    /** Validate instrument having aquired market data */
    void Validate(){
        // first call parent
#ifdef EDR_RISKARB_USE_GENERICNFBASE
        GenericNFBase::Validate();
#else
        GenericNFactor::Validate();
#endif
        static const string method("CorpAct::Validate");
        // check they're both CAssets
        const Asset *assetACQ = 
            dynamic_cast<const Asset*>(assets->getFactor(ACQ).get());
        const Asset *assetTAR =
            dynamic_cast<const Asset*>(assets->getFactor(TAR).get());
        if (!assetACQ || !assetTAR){
            throw ModelException(method, "Underlyings must both be assets");
        }
        if( assetACQ->getCcyTreatment() == CAsset::CCY_TREATMENT_PROTECTED ||
            assetTAR->getCcyTreatment() == CAsset::CCY_TREATMENT_PROTECTED )
            throw ModelException(method, "Do not allow ccy protection");
    }

private:
    friend class CorpActClosedForm;

    CorpAct();
    //CorpAct(const CorpAct& rhs);
    //CorpAct& operator=(const CorpAct& rhs);

    double cashFactor(int type) const;
        
    // for each contingenet value element
    double getContElemValue(int type, double strike) const;

    void getContingentValue(double &contingentValue, double &cashValue) const;

    double getTargetValue() const;

    void addOutputRequests(Control* control,
                           Results* results, double prorateCashValue) const;

    static void load(CClassSP& clazz);

    static IObject* defaultCorpAct(){
        return new CorpAct();
    }

    enum { ACQ=0, TAR=1 };

protected:

    DateTime                closingDate;

    // averaging dates and historical prices to be compatible to heritage EDG booking 
    // this supercedes refLevel and pastValues in GenericNFBase.
    CashFlowArraySP                     histClosingPrice;

    double                  spotAtClosingACQ;
    double                  spotAtClosingTAR;

    double                  proration;
    double                  COBUntendered;

    bool                    hasCashElection;
    double                  percentCash;
    double                  fixedCash;

    DoubleArraySP           multipliers;
    IntArraySP              payoffFlags;
    DoubleArraySP           strikes;
        
    // field not used for pricing
    bool                    hasProrationRisk;
    double          acqFailPrice;
    double          tarFailPrice;
    double          preDealSpread;
    double          initialSpread;
    DateTime        announceDate;

    bool                        accretingInterest;
    double          accretionRate;
    DateTime            accretionBaseDate;
    int             accretionFrequency;
    bool                        proratedDividend;
    int                             debug;
    string          riskArbType;

    // internal field
    int             nbElem;     // # of payoff elements
    SampleListSP    samples;        // get data from histClosingPrice
};

void CorpAct::validatePop2Object()
{
    static const string method("CorpAct::validatePop2Object");
#ifdef EDR_RISKARB_USE_GENERICNFBASE
    GenericNFBase::validatePop2Object();
#else
    GenericNFactor::validatePop2Object();
#endif
    if( assets->NbAssets() != 2 )
        throw ModelException(method, "Must have 2 asset for acquirer and target");

    // need to create samples from histClosingPrice. if does not exist, from refLevel/pastValues
    try {

        if( !!histClosingPrice && histClosingPrice->size() != 0 )
        {
            int i, nbSample = histClosingPrice->size();
            DateTimeArray dates(nbSample);
            DoubleArray values(nbSample);
            DoubleArray weights(nbSample);
            for(i=0; i<nbSample; i++)
            {
                dates[i] = (*histClosingPrice)[i].date;
                values[i] = (*histClosingPrice)[i].amount;
                weights[i] = 1.0/nbSample;
            }
            samples = SampleListSP(new SampleList(dates, values, weights));
        }
#ifdef EDR_RISKARB_USE_GENERICNFBASE
        else if (!!refLevel && refLevel->numDates(0) != 0 ) 
        {
            if( !pastValues )
                throw ModelException(method, "Past values list is empty");

            if( pastValues->getNumAssets() == 0 )
                throw ModelException(method, "Past values needs to be at least 1 asset");

            int i, nbSample = refLevel->numDates(0);
            const DateTimeArray &dates = refLevel->getAllDates();
            const DoubleArray &pastAmts = pastValues->getPastValues(dates, 0, valueDate); // past value for ACQ only
            int nbPastAmts = pastAmts.size();
            DoubleArray values(nbSample);
            DoubleArray weights(nbSample);
            for(i=0; i<nbSample; i++)
            {
                values[i] = (i<nbPastAmts?pastAmts[i]:0.0);
                weights[i] = 1.0/nbSample;
            }               
            samples = SampleListSP(new SampleList(dates, values, weights));
        }
        else
            throw ModelException(method, "Both histClosingPrice and refLevel are empty");
#else
        else
            throw ModelException(method, "HistClosingPrice is empty");
#endif

    }  catch (exception& e) {
        throw ModelException(e, method, "Failed to create sample list");
    }


    if( !!samples && closingDate < samples->getLastDate() )
        throw ModelException(method, "Closing must be >= last sample date");

    if( !multipliers )
    {
        if( !!payoffFlags || !!strikes )
            throw ModelException(method, "Multiplier, payoffFlags and strike arrays must be of same length");
    }
    else
    {
        // allow no payoff elements. ie. only cash vs target
        nbElem = multipliers->size();
        
        if( !payoffFlags || nbElem != payoffFlags->size() || 
            !strikes || nbElem != strikes->size() )
            throw ModelException(method, "Multiplier, payoffFlags and strike arrays must be of same length");
    
        int i;
        for(i=0; i<nbElem; i++)
        {
            if( (*payoffFlags)[i] != 0 && (*payoffFlags)[i] != 1 )
                throw ModelException(method, "Invalid payoff type. Must be 0 (stock) or 1(cash)");

            if( (*strikes)[i] < 0 )
                throw ModelException(method, "Strike can not be negative");

            if( !Maths::isZero( (*strikes)[i] ) && ( !samples || samples->numDates(0) == 0 ) )
                throw ModelException(method, "Sample must be provided if strike is nonzero");
        }
    }

    if( proration < 0 || proration > 1 )
        throw ModelException(method, "Proration must be between 0 and 1");

    if( hasCashElection && (percentCash < 0 || percentCash > 1 ) )
        throw ModelException(method, "Percent cash must be between 0 and 1 if hasCashElection");
}

/** when to stop tweaking */
DateTime CorpAct::endDate(const Sensitivity* sensControl) const {
    return closingDate;
}

bool CorpAct::sensShift(Theta* theta)
{
    const Asset *assetACQ = 
        &dynamic_cast<const Asset&>(*assets->getFactor(ACQ));
    const Asset *assetTAR =
        &dynamic_cast<const Asset&>(*assets->getFactor(TAR));

    // record valueDate before it changes, also closing prices if cross boundary.
    DateTime today = valueDate;
    if( today < closingDate && theta->rollDate(today) >= closingDate )
    {
        spotAtClosingACQ = assetACQ->fwdValue(closingDate);
        spotAtClosingTAR = assetTAR->fwdValue(closingDate);
    }

#ifdef EDR_RISKARB_USE_GENERICNFBASE
    GenericNFBase::sensShift(theta); // call parent's method. it shift value date
#else
    GenericNFactor::sensShift(theta); // call parent's method. it shift value date
#endif
    samples->roll(assetACQ, today, valueDate, false); // use fwd for theta roll
    return true; // continue to tweak components which implement Theta
}

double CorpAct::cashFactor(int type) const
{
    return (hasCashElection?(type==1?percentCash:(1-percentCash)):1.0);
}

/** price a dead instrument
    returns true if it is dead (and priced), false if it is not dead */
bool CorpAct::priceDeadInstrument(CControl* control, CResults* results) const
{
    static string method = "CorpAct::priceDeadInstrument";

    if( valueDate < closingDate ) return false;

    if (!Maths::isPositive(spotAtClosingACQ) || !Maths::isPositive(spotAtClosingTAR))
        throw ModelException(method, "Spot at closing has not been set");

    double value = 0;
    double cashValue = 0;

    // if we're past settlement, it's just worth nothing
    // if before settlement, calc intrinsic value
    const Asset *assetACQ = 
        &dynamic_cast<const Asset&>(*assets->getFactor(ACQ));
    if (valueDate < instSettle->settles(closingDate, assetACQ))
    {
        double contingentValue = fixedCash *(hasCashElection?percentCash:1.0);;
        cashValue = contingentValue;
        for(int i=0; i<nbElem; i++)
        {
            double valTmp;
            if( (*payoffFlags)[i] == 0 )
            {
                valTmp = spotAtClosingACQ;
                if( !Maths::isZero((*strikes)[i]) )
                    valTmp *= Maths::max(0.0, (1 - (*strikes)[i]/samples->sumToDate(valueDate)));
            }
            else
            {
                valTmp = Maths::max(0.0, samples->sumToDate(valueDate) -  (*strikes)[i]);
            }

            valTmp *= (*multipliers)[i] * cashFactor((*payoffFlags)[i]);
            contingentValue += valTmp;
            if( (*payoffFlags)[i] == 1 ) cashValue += valTmp;
        }
        
        double targetValue = spotAtClosingTAR;

        // (contingentValue - tarValue) * proration - COB * (1 - proration)
        value = (contingentValue - targetValue) * proration -
            COBUntendered * (1. - proration);
    }

    results->storePrice(value, discount->getCcy());
    addOutputRequests(control, results, cashValue * proration);

    return true;
}

double CorpAct::getTargetValue() const
{
    const Asset *assetTAR =
        &dynamic_cast<const Asset&>(*assets->getFactor(TAR));
    double fwd = assetTAR->fwdValue(closingDate);
    double disc = instSettle->pv(valueDate,
                                 closingDate, 
                                 discount.get(), 
                                 assetTAR);
    return fwd * disc;
}

double CorpAct::getContElemValue(int type, double strike) const
{
    double value = 0;
    const Asset *assetACQ = 
        &dynamic_cast<const Asset&>(*assets->getFactor(ACQ));

    switch(type)
    {
    case 0:
        double fwdX, disc;
        
        fwdX = assetACQ->fwdValue(closingDate);
        disc = instSettle->pv(valueDate,
                              closingDate, 
                              discount.get(), 
                              assetACQ);
        
        // simple case if strike is 0 or sample in the past
        if( Maths::isZero( strike ) )
        {
            value = fwdX * disc;
        }
        else if ( samples->countFutureSamples(valueDate) == 0 )
        {
            double fwdZ = samples->sumToDate(valueDate);
            value = fwdX * disc * Maths::max(0.0, 1.0 - strike/fwdZ);
        }
        else
        {
            // obtain volInterp from average and get interpolated vol
            CVolRequestLNSP volRequest(Average::volInterpSpot(valueDate,
                                                              valueDate,
                                                              closingDate,
                                                              false, // fwdStarting,
                                                              strike,
                                                              samples.get()));
            CVolProcessedBSSP  volBS(assetACQ->getProcessedVol(volRequest.get()));

            // fwd, variance and covariance
            double fwdZ = samples->expectedAverage(assetACQ, valueDate);
            double varZ = samples->averageVariance(volBS.get(), valueDate, true); // usePastWeight

            DateTimeArray dtArray(1);
            DoubleArray dbArray(1);
            dtArray[0] = closingDate;
            dbArray[0] = 1.0;
            SampleList closingSample(dtArray, dbArray, dbArray);
            double covarXZ = samples->averageCovariance(&closingSample, volBS.get(), valueDate);

            strike *= fwdX * exp(varZ - covarXZ) / fwdZ;

            value = disc * Black::price(true, // isCall
                                        fwdX, 
                                        strike,
                                        1.0,    // need to take disc outside since strike is "fwd value"
                                        varZ);
        }

        break;
    case 1:
        // if type 1, directly use avg pricing
        CClosedFormLN   modelLocal;
        CResults        resultLocal;
        SensitivityArraySP sensLocal(new SensitivityArray(0));
        OutputRequestArraySP reqLocal(new OutputRequestArray(0));
        CControl        controlLocal(sensLocal, reqLocal, false, "");
        AtMaturity      dummyPremiumSettle;
        AverageSP   avg(Average::makeAvgSpot(true, //isCall
                                             closingDate,
                                             strike,
                                             samples.get(),
                                             instSettle.get(),
                                             &dummyPremiumSettle,
                                             assetACQ,
                                             assetACQ->getCcyTreatment(),
                                             discount.get(),
                                             valueDate,
                                             false,       // fwdStarting
                                             valueDate,   // dummy valueDate
                                             true,        // oneContract
                                             1.0,         // dummy notional
                                             1.0));       // dummy initialSpot
        modelLocal.Price(avg.get(), &controlLocal, &resultLocal);
        value = resultLocal.retrievePrice();
        break;
    }

    return value;
}

void CorpAct::getContingentValue(double &contingentValue, double &cashValue) const
{        
    const Asset *assetACQ = 
        &dynamic_cast<const Asset&>(*assets->getFactor(ACQ));
    double discFactor = instSettle->pv(valueDate,
                                       closingDate, 
                                       discount.get(), 
                                       assetACQ);

    contingentValue = fixedCash * discFactor *(hasCashElection?percentCash:1.0);
    cashValue = contingentValue;
    for(int i=0; i<nbElem; i++)
    {
        try {
            double valTmp = getContElemValue((*payoffFlags)[i], (*strikes)[i]);
            valTmp *= (*multipliers)[i] * cashFactor((*payoffFlags)[i]);
            contingentValue += valTmp;
            if( (*payoffFlags)[i] == 1 ) cashValue += valTmp;
        } catch(exception& e) {
            throw ModelException(e, "CorpAct::getContingentValue", "Failed for " + Format::toString(i+1) + "th payoff");
        }
    }
    return;
}

void CorpAct::addOutputRequests(Control* control,
                                Results* results, double prorateCashValue) const
{
    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        const Asset *assetACQ = 
            &dynamic_cast<const Asset&>(*assets->getFactor(ACQ));

        // FWD_AT_MAT
#ifdef EDR_RISKARB_USE_GENERICNFBASE
        GenericNFBase::addRequests(control, 
                                   results,
                                   closingDate);
#else
        GenericNFactor::addRequests(control, 
                                    results,
                                    closingDate);
#endif

        // IND_VOL. using the first payoff
        OutputRequest* request = NULL;
        if (control->requestsOutput(OutputRequest::IND_VOL, request)) 
        {
            double indVol = 0.0;
            if( closingDate>valueDate )
            {
                // obtain volInterp from average and get interpolated vol
                // if strike=0, the AtmVolRequest is returned
                double strike = ((!strikes)||strikes->size()==0)?assetACQ->getSpot():(*strikes)[0];
                CVolRequestLNSP volRequest(Average::volInterpSpot(valueDate,
                                                                  valueDate,
                                                                  closingDate,
                                                                  false, // fwdStarting,
                                                                  strike,
                                                                  samples.get()));
                CVolProcessedBSSP  volBS(assetACQ->getProcessedVol(volRequest.get()));
                indVol = volBS->CalcVol(valueDate, closingDate);
            }
            results->storeRequestResult(request, indVol); 
        }

        // discount factor
        request = NULL;
        if (control->requestsOutput(OutputRequest::DISCOUNT_FACTOR, request)) 
        {
            double disc = 1.0;
            if(closingDate>valueDate)
            {
                disc = instSettle->pv(valueDate,
                                      closingDate, 
                                      discount.get(), 
                                      assetACQ);
            }
            results->storeRequestResult(request, disc); 
        }

        // cash value
        request = NULL;
        if (control->requestsOutput(OutputRequest::PRORATED_CASH_VALUE, request)) 
        {
            results->storeRequestResult(request, prorateCashValue); 
        }

    }
}

// for reflection
CorpAct::CorpAct(): 
#ifdef EDR_RISKARB_USE_GENERICNFBASE
    GenericNFBase(TYPE),
#else
    GenericNFactor(TYPE),
#endif
    spotAtClosingACQ(0), spotAtClosingTAR(0), proration(1), COBUntendered(0), 
    hasCashElection(false), percentCash(0), fixedCash(0), 
// fields not used for pricing
    hasProrationRisk(false), acqFailPrice(0), tarFailPrice(0), preDealSpread(0),
    initialSpread(0), accretingInterest(false), accretionRate(0), 
    accretionFrequency(0), proratedDividend(false), debug(0), nbElem(0)
{}

// Invoked when Class is 'loaded'
void CorpAct::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Risk arb instrument to price the spread between "
                          "target share value and acquirer payoff value");
    REGISTER(CorpAct, clazz);
#ifdef EDR_RISKARB_USE_GENERICNFBASE
    SUPERCLASS(GenericNFBase);
#else
    SUPERCLASS(GenericNFactor);
#endif
    IMPLEMENTS(CClosedFormLN::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultCorpAct);
    FIELD(closingDate,       "Deal closing date");
    FIELD(histClosingPrice,                 "Historical samples for averaging");
    FIELD_MAKE_OPTIONAL(histClosingPrice);
    FIELD(proration,         "Proration. Default 100%");
    FIELD_MAKE_OPTIONAL(proration);
    FIELD(hasCashElection,   "Has cash election. Default false");
    FIELD_MAKE_OPTIONAL(hasCashElection);
    FIELD(fixedCash,         "Fixed cash. Default 0");
    FIELD_MAKE_OPTIONAL(fixedCash);
    FIELD(percentCash,       "Percent of cash. Default 0");
    FIELD_MAKE_OPTIONAL(percentCash);
    FIELD(COBUntendered,     "Cost of borrow untendered");
    FIELD_MAKE_OPTIONAL(COBUntendered);
    FIELD(multipliers,              "Multipliers");
    FIELD(payoffFlags,              "Payoff types, 0 (stock) or 1 (cash)");
    FIELD(strikes,                  "Strikes");
    FIELD(spotAtClosingACQ,  "Spot at closing of acquirer");
    FIELD_MAKE_OPTIONAL(spotAtClosingACQ);
    FIELD(spotAtClosingTAR,  "Spot at closing of target");
    FIELD_MAKE_OPTIONAL(spotAtClosingTAR);
    //non pricing member
    FIELD(hasProrationRisk,  "Strike of each instrument in payoff");
    FIELD_MAKE_OPTIONAL(hasProrationRisk);
    FIELD(acqFailPrice,              "Acquirer price after deal fails");
    FIELD_MAKE_OPTIONAL(acqFailPrice);
    FIELD(tarFailPrice,              "Target price after deal fails");
    FIELD_MAKE_OPTIONAL(tarFailPrice);
    FIELD(preDealSpread,             "Another number for estimating loss if failure");
    FIELD_MAKE_OPTIONAL(preDealSpread);
    FIELD(initialSpread,             "Spread immediately following announcement");
    FIELD_MAKE_OPTIONAL(initialSpread);
    FIELD(announceDate,              "Strike of each instrument in payoff");
    FIELD_MAKE_OPTIONAL(announceDate);
    FIELD(accretingInterest, "does the cash amount accrete");
    FIELD_MAKE_OPTIONAL(accretingInterest);
    FIELD(accretionRate,             "annual accretion rate as percent for cash amount");
    FIELD_MAKE_OPTIONAL(accretionRate);
    FIELD(accretionBaseDate, "date from which accreting begins");
    FIELD_MAKE_OPTIONAL(accretionBaseDate);
    FIELD(accretionFrequency,"number of times a year accretion occurs");
    FIELD_MAKE_OPTIONAL(accretionFrequency);
    FIELD(proratedDividend,  "are you paid a prorated dividend amount");
    FIELD_MAKE_OPTIONAL(proratedDividend);
    FIELD(debug,                     "not used");
    FIELD_MAKE_OPTIONAL(debug);
    FIELD(riskArbType,               "not used for pricing");
    FIELD_MAKE_OPTIONAL(riskArbType);
    // internal member
    FIELD(samples,                                  "");
    FIELD_MAKE_TRANSIENT(samples);  // hide from interface
    FIELD(nbElem,            "");
    FIELD_MAKE_TRANSIENT(nbElem);   // hide from interface
};

CClassConstSP const CorpAct::TYPE = CClass::registerClassLoadMethod(
    "CorpAct", typeid(CorpAct), CorpAct::load);

// * for class loading (avoid having header file) */
bool CorpActLoad() {
    return (CorpAct::TYPE != 0);
}


/** private class */
class CorpActClosedForm: public CClosedFormLN::IProduct{
private:
    const CorpAct*  corpact; // a reference

public:
    CorpActClosedForm(const CorpAct* corpact): corpact(corpact){}

    void price(CClosedFormLN*   model,
               Control*        control, 
               CResults*       results) const;
};

void CorpActClosedForm::price(CClosedFormLN*   model,
                              Control*        control, 
                              CResults*       results) const
{
    static const string method = "CorpActClosedForm::price";
    try {
        double         premium;         // the fair value 

        // valueDate >= matDate is taken care of here
        if(corpact->priceDeadInstrument(control, results)){
            return; // dead instrument priced
        }
        
        double contingentValue, cashValue;
        corpact->getContingentValue(contingentValue, cashValue);

        double targetValue = corpact->getTargetValue();

        // (contingentValue - tarValue) * proration - COB * (1 - proration)
        premium = (contingentValue - targetValue) * corpact->proration -
            corpact->COBUntendered * (1. - corpact->proration);

        results->storePrice(premium, corpact->discount->getCcy());
        corpact->addOutputRequests(control, results, cashValue * corpact->proration);

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* CorpAct::createProduct(
    CClosedFormLN* model) const{

    return new CorpActClosedForm(this);
}

DRLIB_END_NAMESPACE

