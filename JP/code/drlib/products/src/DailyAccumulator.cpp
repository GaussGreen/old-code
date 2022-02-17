//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DailyAccumulator.cpp
//
//   Author      : Keith Law
//
//   Description : Accumulator using new Observationbuilder to build the dates on the fly, used for
//                 daily monitoring trades only.
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/FD1F.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/Tree1fLV.hpp"
#include "edginc/Tree1fLN.hpp"
#include "edginc/Tree1f.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/PhysicalDelivery.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"
#include "edginc/TreeSliceOper.hpp"
// for Datebuilder and Asset History 
#include "edginc/DateBuilder.hpp"
#include "edginc/VarSwapUtilities.hpp"
#include "edginc/ObservationBuilder.hpp"
#include "edginc/MarketObservable.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/PhysicalAndCashSettlement.hpp"


DRLIB_BEGIN_NAMESPACE

class DailyAccumulator : public Generic1Factor, 
                         virtual public FDModel::IIntoProduct,
                         virtual public LegalTerms::Shift,
                         virtual public LastSensDate
{
public:
    static CClassConstSP const TYPE;

    // override base implementation if required
    virtual void GetMarket(const IModel*          model, 
                           const CMarketDataSP    market);
    
    virtual void Validate();

    virtual void validatePop2Object();

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;

    // below are copied from Vanilla
    virtual DateTime getValueDate() const;

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
        returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

    virtual DateTime endDate(const Sensitivity* sensControl) const;

    void addOutputRequests(Control* control,
                           Results* results,
                           const double& fairValue,
                           const double& indVol) const;
    
    bool sensShift(Theta* shift);

    bool sensShift(LegalTerms* shift);

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel *     model);

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    double computeHistContracts(const DateTime& evalDate) const;

    CashFlowArraySP getKnownCashFlows() const;

    DateTimeArraySP getPaymentDates() const;

    PhysicalDeliveryArray getPhysicalDelivery() const;

    double decRound(double value, int deci) const //function for getting rounded strike, barrier etc. as user request
    {
        double mult = ::pow(10.0, deci);
        double product = value * mult;
        int med = Maths::round(product);
        product = (double)(med / mult);
        return product;
    }

    double productValue(double spot, double strike) const
    {
        double val = spot - strike;
        if( productAccumulated == 2 )
            val = Maths::max(0.0, val);
        else if( productAccumulated == 3 )
            val = Maths::max(0.0, -val);
        return val;
    }

    void populateSamples(const CAsset*               asset,
                         const ObservationSource*    assetHistorySource,
                         const DateTime&             valueDate,
                         const DateTimeArray&        obsDates,
                         const ObservationTypeArray& obsTypes,
                         DoubleArray&                obsSamples,
                         SamplingConventionSP        sampleRule=SamplingConventionSP()) {
    
    if (!sampleRule) {
        sampleRule.reset(new UnadjustedConvention());
    }

    static const string routine = "DailyAccumulator::populateSamples";
    try {
        // Fill in historic samples
        FixingTypeSP fixType(new AssetFixType(asset->getTrueName()));
        
        for (int i = 0; i < obsDates.size() && obsDates[i] <= valueDate; i++) {
            obsSamples[i] = asset->pastValue(obsDates[i],
                                              obsTypes[i].get(),
                                              assetHistorySource,
                                              fixType.get(),
                                              0,
                                              sampleRule.get());
            // Temp validation
            if (Maths::isZero(obsSamples[i])) {
                throw ModelException(routine, "Invalid sample (0.0) found for date " +
                                              obsDates[i].toString() + " for asset " +
                                              asset->getName());
            }
        }

        } catch(exception& e) {
        
            throw ModelException(e, routine);
        }
    }
                         
 
    /** return first index in sorted dates that is on or after aDate */
    int getIndexOnOrAfterDate(const DateTimeArray dates, const DateTime aDate) const //new
    {
        int idx = 0;
        while ( idx < dates.size() && dates[idx] < aDate ){
            idx ++;
        }
        return idx;
    }

    /** find pay date in PaymentDates */
   DateTime getPayDate(const DateTime& today) const
    { 
        int i = 0, j = 0;
        
        bool found = false;
        
        while (!found &&  i < monitorDates->size() )
        {
            if ( (*monitorDates)[i].equals(today) ) 
            {
                found = true;
            }

            i++;
        }

        if (found) // payment date is the first date equal to or after the supplied date
        {
            while(j < monitoringEndDates->size() && today > (*monitoringEndDates)[j]) j++; 

            return (*paymentDates)[j]; 
        }
        else
        {
            throw ModelException("DailyAccumulator::getPayDate", "date supplied is not in monitoring dates!");
        }
    }

    // Gives the payment date corresponding to the monitoring date that is equal to 
    // or the first after today            
    DateTime getNextPaymentDate(const DateTime& today) const
    {
        DateTime nextPayDate;

        bool found = false;
        int i = 0;
        int j = 0;
        while (!found && i < monitorDates->size())
        {
            if ( today <= (*monitorDates)[i] )
            {
                found = true;

                while((*monitorDates)[i] > (*monitoringEndDates)[j]) j++; 

                nextPayDate = (*paymentDates)[j];
            }
            i++;
        }
        if (!found) // in this case, we are after all monitoring dates so result here is last payment date
        {
            nextPayDate = paymentDates->back();
        }
    
        return nextPayDate;
    }

    // check whether historical hit has happened
    bool checkHit(const DateTime& evalDate) const
    {
        double pastSpotValue = 0.0;
        DateTime currentDate;
        double barLevel = legalBarrier * initialSpot;
        if(decimalRound>=0) barLevel = decRound(barLevel,decimalRound);

        int i = 0, j = 0;
        if (!isHit && isGlobalKO && payAtHit)
        {
            while ( i<monitorDates->size() && ((*monitorDates)[i] <= evalDate) &&(!isHit) ) 
            {
            
                currentDate = (*monitorDates)[i];

                pastSpotValue = (*obsSamples)[i];
            
                if( (isUpAndOut&& (pastSpotValue > barLevel * (1-FP_MIN))) || (!isUpAndOut && (pastSpotValue < barLevel * (1+FP_MIN))) )
                {            
                    isHit = true;
                    hitDate = currentDate;
                    spotAtHit = pastSpotValue;
                }

                i++;
            }
        }
        return isHit;
    }

    // computes historical contracts observed strictly prior to current moment evalDate
    // the boolean includeCurrentObservation is to include historical contract observed exactly at evalDate
    double getKoHistContracts(const DateTime& evalDate, const bool includeCurrentObservation =false) const
    {
        int i = 0, j = 0;
        double nbContracts = 0.0;
        double pastSpotValue = 0.0;
        double barLevel = legalBarrier * initialSpot;
        bool checkCondition;
        double iRebate = rebate * initialSpot;
        if(decimalRound>=0) barLevel = decRound(barLevel,decimalRound);

        KoRebate = 0.0;
        isKO = false;

        DateTime currentDate;
        DateTime settleDate;
        DateTime PayDate = getNextPaymentDate(evalDate);

        while ( i<monitorDates->size() && ((*monitorDates)[i] <= evalDate) && (!isKO) ) 
        {
            currentDate = (*monitorDates)[i];
            settleDate = getPayDate(currentDate);
 
            pastSpotValue = (*obsSamples)[i];

            // in case of global ko past monitoring period knock-outs are taken into account
            if(isGlobalKO)
            {   
                if( (isUpAndOut&& (pastSpotValue > barLevel * (1-FP_MIN))) || (!isUpAndOut && (pastSpotValue < barLevel * (1+FP_MIN))) )
                {            
                    isKO = true;
                    knockedOutDate = currentDate;
                }
            }
            
            // in any case check this monitoring period's past values
            if(!includeCurrentObservation)
            {
                checkCondition = (settleDate > evalDate) && (settleDate.equals(PayDate));
            }
            else // if includeCurrentObservation, we take a settle date happening NOW into account
            {
                checkCondition = (settleDate >= evalDate) && (settleDate.equals(PayDate));
            }
            
            if (checkCondition)
            {            
                // checking KO
                if( (isUpAndOut&& (pastSpotValue > barLevel * (1-FP_MIN))) || (!isUpAndOut && (pastSpotValue < barLevel * (1+FP_MIN))) )
                {            
                    isKO = true;
                    knockedOutDate = DateTime(currentDate);
                    iRebate = initialSpot* (rebate - (i+1)*rebate / (monitorDates->size()));
                    KoRebate = uniformRebate? (rebate* initialSpot) : iRebate;
                }
                if (!isKO)
                {
                    nbContracts +=  participation;

                    if( hasExtraParticipation() )
                    {
                        double extraStrike = strike * initialSpot;
                        if(decimalRound>=0) extraStrike = decRound(extraStrike,decimalRound);
                        if( (isUpAndOut && pastSpotValue<extraStrike ) || (!isUpAndOut && pastSpotValue>extraStrike) )
                            nbContracts +=  extraParticipation;
                    }
                }

                if(dayKoOnly) isKO = false;
            }
            i++;    
        }
        
        return nbContracts;
    }



    /* true if we are at the end of a monitoring period */
    bool isEndMonitoringPeriod(const DateTime& testDate) const 
    {
        int i = 0;
        bool isEndPeriod = false;    
        bool monitoringDateFound = false;

        // checking whether testDate is a monitoring date
        while ( i<monitorDates->size() && !monitoringDateFound )
        {
            if( testDate.equals((*monitorDates)[i]) )
            {
                monitoringDateFound = true;
            }
            i++;
        }

        // now looking whether we are at a payment date  
        if (monitoringDateFound)
        {

            for(int j = 0; j < monitoringEndDates->size() ; j++){ 
                if(testDate.equals((*monitoringEndDates)[j])) isEndPeriod = true;
            }

        }

        return isEndPeriod;
    }

private:
    friend class DailyAccumulatorHelper;
    friend class DailyAccumulatorFDProd;

protected:
    bool hasExtraParticipation() const
    {
        return (!Maths::isZero(extraParticipation));
    }

    static void load(CClassSP& clazz);
    DailyAccumulator();

    IDateBuilderSP          monitorEndEachPeriod;     // end date for each monitoring period
    IDateBuilderSP          paymentDateEachPeriod;            // one for each monitoring period
    CashFlowArray           sampleOverride;            // for overriding Asset History purpose
    double                  strike;                    
    double                  participation;              
    double                  barrierLevel;               
    double                  legalBarrier;        // legal barrier used for reporting only
    double                  rebate;                     // uniform if uniformRebate, linearly decrease to zero if Notes form
    bool                    uniformRebate;                    // rebate will be fixed if true, else (notes) rebate will linearly decrease to zero
    bool                    payAtHit;            // True: shares delivered at KO, else at end of monitor period
    IObservationBuilderSP   observationBuilder;     //!< Factory for observations
    string                  assetHistorySource;     //!< Designates asset history source
    PhysicalAndCashSettlementSP phyAndCashSettle;   // physical and cash settlement
    
    // extra participation available if spot reach the strike
    double                  extraParticipation;  

    bool                    isUpAndOut;    
    bool                    isGlobalKO;
    bool                    includeDivs;
    string                  productToAccumulate;    //"forward", "call", "put"
    
    mutable bool            isKO; // $unregistered
    mutable DateTime        knockedOutDate; // $unregistered
    mutable double          KoRebate; // $unregistered    
 
    bool                    dayKoOnly; // Only the KO date not accumulating
    int                     decimalRound; //number of decimal places, defaulted to '-1' for no rounding
   
    // Transient fields
    DateTimeArraySP         monitorDates;
    ObservationTypeArraySP  obsTypes;
    DoubleArraySP           obsSamples;
    ObservationSourceSP     assetHistorySourceObject; //!< object for asset history source
    DateTimeArraySP         monitoringEndDates;
    DateTimeArraySP         paymentDates;
    int                     productAccumulated;
    mutable bool            isHit; // for payAtHit purpose
    mutable DateTime        hitDate; 
    mutable double          spotAtHit;
    bool                    isContinuouslyMonitored;
};

typedef smartPtr<DailyAccumulator> DailyAccumulatorSP;

// helpers
class DailyAccumulatorHelper {
public:
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(DailyAccumulator, clazz);
        SUPERCLASS(Generic1Factor);
        EMPTY_SHELL_METHOD(defaultDailyAccumulator);
        // same as vanilla
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(ISensitiveStrikes);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(LegalTerms::Shift);
        FIELD(monitorEndEachPeriod, "End date of each monitoring period");
        FIELD(paymentDateEachPeriod, "Payment dates: Corresponding payment date for each Monitoring End Date");  
        FIELD(sampleOverride, "Past Sample Dates: For overriding Asset History Purpose");  
        FIELD_MAKE_OPTIONAL(sampleOverride);
        FIELD(strike, "Strike");
        FIELD(participation, "Participation");
        FIELD(barrierLevel,        "Barrier Level for Knock-Out Condition");
        FIELD(legalBarrier,        "Reported Barrier");
        FIELD(rebate,        "Rebate in case of Knock-Out");
        FIELD(uniformRebate,        "True: uniformRebate; False: Note form");
        FIELD(payAtHit, "True: Pay or deliver shares at KO; False: deliver at monitoring end");
        FIELD(isUpAndOut,        "Is product Up and Out ? If not, product is Down and Out");
        FIELD(isGlobalKO,   "Does a Knock-Out affect all the product ?");
        FIELD(includeDivs,   "Are intermediary dividends included ?");
        FIELD(productToAccumulate, "What we accumulated can be a forward, a call or a put");
        FIELD(extraParticipation, "Extra participation on a day if spot is above/below strike for downKO/upKO respectively");
        FIELD(assetHistorySource, "Asset history source");
        FIELD_MAKE_OPTIONAL(assetHistorySource);
        FIELD(observationBuilder, "Observation builder");
        FIELD(phyAndCashSettle, "Physical and Cash Settlement");

        FIELD(dayKoOnly, "if true, only the KO date is not accumulating, while the rest of period continues");  
        FIELD_MAKE_OPTIONAL(dayKoOnly);

        FIELD(decimalRound, "Round strike and barrier to the number of decimal places");  
        FIELD_MAKE_OPTIONAL(decimalRound);

        // Transient fields
        FIELD(monitorDates, "monitoring dates");
        FIELD_MAKE_TRANSIENT(monitorDates);    
        FIELD(obsTypes, "Observation types");
        FIELD_MAKE_TRANSIENT(obsTypes);    
        FIELD(obsSamples, "Observation samples");    
        FIELD_MAKE_TRANSIENT(obsSamples);    
        FIELD(assetHistorySourceObject, "Asset history source object");    
        FIELD_MAKE_TRANSIENT(assetHistorySourceObject);
        FIELD(monitoringEndDates, "Monitor End Dates each period");    
        FIELD_MAKE_TRANSIENT(monitoringEndDates);
        FIELD(paymentDates, "Payment Date each peiod");    
        FIELD_MAKE_TRANSIENT(paymentDates);
        FIELD(productAccumulated, "1-forward, 2-call, 3-put");    
        FIELD_MAKE_TRANSIENT(productAccumulated);
        FIELD(isHit, "knocked-out already and payAtHit");    
        FIELD_MAKE_TRANSIENT(isHit);
        FIELD(hitDate, "hit date");    
        FIELD_MAKE_TRANSIENT(hitDate);
        FIELD(spotAtHit, "stock price at hit date");    
        FIELD_MAKE_TRANSIENT(spotAtHit);
        FIELD(isContinuouslyMonitored,   "Is Barrier Continuously monitored ?");
        FIELD_MAKE_TRANSIENT(isContinuouslyMonitored);
    }

    static IObject* defaultDailyAccumulator(){
        return new DailyAccumulator();
    }
};


CClassConstSP const DailyAccumulator::TYPE = CClass::registerClassLoadMethod(
    "DailyAccumulator", typeid(DailyAccumulator), DailyAccumulatorHelper::load);


void DailyAccumulator::GetMarket(const IModel*         model, 
                            const CMarketDataSP    market)
{
    market->GetReferenceDate(valueDate);

    if (asset.usingCache() || !Tree1fLN::TYPE->isInstance(model))
    {// should always come in - just to cope with old regression convertion
        CAsset::getAssetMarketData(model, market.get(), ccyTreatment, 
                                   discount, asset);
    }

    discount.getData(model, market);
    instSettle->getMarket(model, market.get());

          // Grab samples from asset history
    populateSamples(asset.get(), 
            assetHistorySourceObject.get(), valueDate, *monitorDates, *obsTypes, *obsSamples);

    if (premiumSettle.get())
    {
        premiumSettle->getMarket(model, market.get());
    }
    
    // override asset history if there is any...
    if(!Maths::isZero(sampleOverride.size()))
        {
            int i = 0;
            int j = 0;
            for(i=0; i<sampleOverride.size(); i++)
            {
                bool overrideAH = false;
                while(!overrideAH && j<monitorDates->size() )
                {
                    if((*monitorDates)[j].equals(sampleOverride[i].date))
                    {
                        (*obsSamples)[j] = sampleOverride[i].amount;
                        overrideAH = true;
                    }
                    j++;
                }
            }
        }

    isHit = checkHit(valueDate);
}

bool DailyAccumulator::avoidVegaMatrix(const IModel* model)
{
    return CTree1fLV::TYPE->isInstance(model);
}

// constructor: initialize optional ks
DailyAccumulator::DailyAccumulator(): Generic1Factor(TYPE), 
        uniformRebate(true),
        payAtHit(false),
        isUpAndOut(true),
        isContinuouslyMonitored(false),
        isGlobalKO(false),
        includeDivs(false),
        productAccumulated(1),
        isKO(false),
        knockedOutDate(),
        isHit(false),
        hitDate(),
        spotAtHit(0.0),
        dayKoOnly(false),
        decimalRound(-1),
        assetHistorySource(IMarketObservable::DEFAULT_SOURCE){}

void DailyAccumulator::validatePop2Object(){
    static const string method = "DailyAccumulator::validatePop2Object";

        productAccumulated = (productToAccumulate=="forward")?1:(productToAccumulate=="call")?2:3;
        monitorDates = observationBuilder->dateList();
        obsTypes = observationBuilder->obsTypes();
        obsSamples = DoubleArraySP(new DoubleArray(monitorDates->size()));
        monitoringEndDates = monitorEndEachPeriod->dates();
        paymentDates = paymentDateEachPeriod->dates();
 
        if((isUpAndOut && (*obsTypes)[0]->toString()=="High") ||(!isUpAndOut && (*obsTypes)[0]->toString()=="Low")) isContinuouslyMonitored = true;

         // payAtHit is only used for isGlobalKO
        if(payAtHit && !isGlobalKO)
        throw ModelException(method, "payAtHit flag is only used for isGlobalKO");  


       // validate monitor end dates are in ascending order
	    for(int i = 0; i< (monitoringEndDates->size()-1); i++)
	    {
            if( (*monitoringEndDates)[i] >= (*monitoringEndDates)[i+1] )
            throw ModelException(method, "MonitoringEndDate ("+(*monitoringEndDates)[i+1].toString()+") need to be in ascending order");
	    }

        // validate each monitor end date is in monitorDates
	    for(int i = 0; i< monitoringEndDates->size(); i++)
	    {
            bool found = false;
            int j = 0;
            while(!found && j < monitorDates->size())
            {
                if((*monitorDates)[j].equals((*monitoringEndDates)[i])) found = true;
                j++;
            }
            if (found == false)
                throw ModelException(method, "MonitoringEndDate ("+(*monitoringEndDates)[i].toString()+") is not included in monitor dates");
            
	    }

        // the last monitoringEndDates need to be the same as last monitorDates
        if(!(monitorDates->back().equals(monitoringEndDates->back())))
            throw ModelException(method, "The last monitoringEndDate needs to be the same as the last monitorDates");

        // dayKoOnly cannot be used together with payAtHit, isGlobalKO or Notes form
        if(( dayKoOnly && payAtHit)||(dayKoOnly && isGlobalKO)||(dayKoOnly && !uniformRebate))
        throw ModelException(method, "dayKoOnly cannot be used together with payAtHit, isGlobalKO or Notes (non-uniformRebate) form");  

        // build the observation source object
        assetHistorySourceObject 
            = ObservationSourceSP(new ObservationSource(assetHistorySource));
           
        if (!(productToAccumulate=="forward"||productToAccumulate=="call"||productToAccumulate=="put"))
            throw ModelException(method, "productToAccumulate can only be foward, call or put");
 
        if( includeDivs && productToAccumulate != "forward" )
        throw ModelException(method, "IncludeDivs flag only allowed if accumulate fwd contract.");      

}

void DailyAccumulator::Validate()
{
    static const string method = "DailyAccumulator::Validate";
    // just check the things that aren't/cannot be checked in 
    // validatePop2Object
 
    // validate against struck until thoroughly tested and agreed methodology
    if (ccyTreatment == Asset::CCY_TREATMENT_STRUCK) {
        throw ModelException(method, "ccy struck type not supported.");
    }

    // validation against forward start
    if ( fwdStarting ) {
        throw ModelException(method, "Forward Starting feature not supported");
    }

    AssetUtil::assetCrossValidate(asset.get(),
                                  fwdStarting,
                                  startDate,
                                  valueDate,
                                  discount,
                                  this);  

    //validate that SpotAtStart is not zero or negative
    if ( !Maths::isPositive(initialSpot) )
    {
        throw ModelException(method,
                                     "Initial Spot value must be positive ");
    }

    DateTime start = valueDate;
    int numPastSampleDates = monitorDates->size() - start.numFutureDates(*monitorDates);

}

/** returns the current value date */
DateTime DailyAccumulator::getValueDate() const
{
   return valueDate;
}

//////////////////////////////////////////////////
//
// We call getKoHistContracts() many times which is unnecessary.
// ideally should preprocess once with all historical accumulation 
// pre-computed to be used for pricing and for KnowCF/PhysDelivery
//
//////////////////////////////////////////////////
CashFlowArraySP DailyAccumulator::getKnownCashFlows() const
{
    DateTime settleDate,pastPayDate,monitorDate, pastDate = valueDate;
    double strikes = strike * initialSpot;
    if(decimalRound>=0) strikes = decRound(strike * initialSpot,decimalRound);
    CashFlowArraySP cfl(new CashFlowArray(0));
       
    // return empty cfl if all monitor dates are in future
    if ( valueDate < monitorDates->front() ) return cfl;
    
    // computing the scaling factor
    double fwdStrt = fwdStarting?asset->fwdValue(startDate):0.0;
    double scalingFactor = InstrumentUtil::scalePremium(oneContract,fwdStarting,notional,fwdStrt,initialSpot);

    // prepare to later add div
    DateTimeArrayConstSP exDivDates;
    DoubleArrayConstSP exDivAmounts;
    if(includeDivs)
    {
        DividendListSP exDivs = AssetUtil::getAllDivsBetweenDates(asset, startDate, pastDate);       
        exDivDates = exDivs->getExDivDates();
        exDivAmounts = exDivs->getDivAmounts();
    }
    
    // find the latest monitor date, can be pastDate
    int iEnd = monitorDates->size()-1, i; 
    double usedSample = 0.0;
    while( iEnd >= 0 && (*monitorDates)[iEnd] > pastDate ) 
        iEnd--;
    if( iEnd < 0 )    // sanity check, should never happen
        throw ModelException("DailyAccumulator::getKnownCashFlows", "Internal error.");
    
    int j = 0;
    int k = 0;

    for(i=0; i<=iEnd; i++)
    {
        if (k < monitoringEndDates->size() && (*monitoringEndDates)[k] < (*monitorDates)[i]) k++;

        bool isEndPeriod = (i == (monitorDates->size()-1) || (*monitorDates)[i] == (*monitoringEndDates)[k] );
        
        // for cash settle, add periodic payment
        // if KO or end of monitor period, include current period accumulated contract.
        // otherwise, only include previous fully historical periods        
        if( (!instSettle->isPhysical()|| phyAndCashSettle->isPhysicalAndCash) &&
            ( isEndPeriod || ( iEnd==i && isKO )||( iEnd==i && isHit) ))
        {
            monitorDate = (isHit && ((*monitorDates)[i] > hitDate))?hitDate:(*monitorDates)[i];                       
            usedSample = (isHit && ((*monitorDates)[i] > hitDate))?spotAtHit:(*obsSamples)[i];
            
            settleDate =(payAtHit && isHit && iEnd==i)? instSettle->settles(hitDate, asset.get()):
                                    instSettle->settles((*paymentDates)[k], asset.get());

            double histContracts = getKoHistContracts(monitorDate,true);
            double prodVal = histContracts*productValue(usedSample, strikes);
            prodVal = phyAndCashSettle->isPhysicalAndCash? prodVal * phyAndCashSettle->cashWeight : prodVal;
            if( !Maths::isZero(prodVal) )
                cfl->push_back(CashFlow(settleDate, scalingFactor*prodVal));
        }
        
        // adding dividend payments as they become known
        if( includeDivs && exDivDates->size() > 0 )
        {
            DateTime monitorStart = (*monitorDates)[i];
            DateTime monitorEnd = isEndPeriod?monitorStart:(*monitorDates)[i+1];
            while( j < exDivDates->size() && 
                (*exDivDates)[j] >= monitorStart &&
                (*exDivDates)[j] < monitorEnd &&
                !Maths::isZero((*exDivAmounts)[j])) // if zero div do nothing
            {
                // contract does not include participation on ex div date
                DateTime divDate = DateTime( (*exDivDates)[j].getDate(), DateTime::START_OF_DAY_TIME);
                settleDate = instSettle->settles(divDate, asset.get());
                double divHistContracts = getKoHistContracts(divDate);
                double divAmount = divHistContracts*(*exDivAmounts)[j];
                if (!Maths::isZero(divAmount))
                    cfl->push_back(CashFlow(settleDate, scalingFactor*divAmount));
                j++;
            }
        }
    }
    
    return cfl;
}

PhysicalDeliveryArray DailyAccumulator::getPhysicalDelivery() const
{
    DateTime pastDate = valueDate, pastPayDate, monitorDate, settleDate;
    double usedSample = 0.0;
    double strikes = strike*initialSpot;
    if(decimalRound>=0) strikes = decRound(strike*initialSpot,decimalRound);
    PhysicalDeliveryArray pda;
  
    // return empty cfl if all monitor dates are in future
    if ( valueDate < monitorDates->front() ) return pda;
    
    // computing the scaling factor
    double fwdStrt = fwdStarting?asset->fwdValue(startDate):0.0;
    double scalingFactor = InstrumentUtil::scalePremium(oneContract,fwdStarting,notional,fwdStrt,initialSpot);
    
    // find the latest monitor date, can be pastDate
    int iEnd = monitorDates->size()-1; 
    while( iEnd >= 0 && (*monitorDates)[iEnd] > pastDate ) iEnd--;
    
    // sanity check, should never happen
    if( iEnd < 0 )
        throw ModelException("DailyAccumulator::getKnownCashFlows", "Internal error.");
    
    int k = 0;

    for(int i=0; i<=iEnd; i++)
    {
        if (k < monitoringEndDates->size() && (*monitoringEndDates)[k] < (*monitorDates)[i]) k++;
        bool isEndPeriod = (i == (monitorDates->size()-1) || (*monitorDates)[i] == (*monitoringEndDates)[k] );
        
        // for cash settle, add periodic payment
        // if KO or end of monitor period, include current period accumulated contract.
        // otherwise, only include previous fully historical periods        
        if( isEndPeriod || ( iEnd==i && isKO )||( iEnd==i && isHit) )
        {
            monitorDate = (isHit && ((*monitorDates)[i] > hitDate))?hitDate:(*monitorDates)[i];
            usedSample = (isHit && ((*monitorDates)[i] > hitDate))?spotAtHit:(*obsSamples)[i];

            pastPayDate = (payAtHit && isHit && iEnd==i)? hitDate:(*paymentDates)[k];

            settleDate = instSettle->settles(pastPayDate, asset.get());
            
            double histContracts = getKoHistContracts(monitorDate,true);
            histContracts = phyAndCashSettle->isPhysicalAndCash? histContracts * phyAndCashSettle->physicalWeight : histContracts;
            double spotAtMat = usedSample;
            
            if( !Maths::isZero(histContracts) &&
                (productAccumulated == 1 ||
                (productAccumulated == 2 && spotAtMat > strikes) ||
                (productAccumulated == 3 && spotAtMat < strikes) ))
            {
                if( productAccumulated == 3 ) histContracts = - histContracts;
                
                if(phyAndCashSettle->phyDeliveredAtSpot)
                    {
                    PhysicalDelivery delivery(scalingFactor*histContracts,usedSample, monitorDate, settleDate);
                    pda.push_back(delivery);
                    }
                else 
                    {
                    PhysicalDelivery delivery(scalingFactor*histContracts,strikes, monitorDate, settleDate);
                    pda.push_back(delivery);
                    }                
            }    
        }
        
    }
     
    return pda;
}

DateTimeArraySP DailyAccumulator::getPaymentDates() const
{
    DateTimeArraySP payDates(new DateTimeArray(paymentDates->size()));
    
    for(int i=0; i<paymentDates->size(); i++)
        (*payDates)[i] = (isHit && (hitDate < (*paymentDates)[i]))?
        instSettle->settles(hitDate, asset.get()):instSettle->settles((*paymentDates)[i], asset.get());

    return payDates;
}

/** price a dead instrument until settlement - exercised, expired, knocked out etc.
returns true if it is dead (and priced), false if it is not dead */
bool DailyAccumulator::priceDeadInstrument(CControl* control, CResults* results) const
{
    static const string method = "DailyAccumulator::priceDeadInstrument";

     // computing the scaling factor
    double fwdStrt = fwdStarting?asset->fwdValue(startDate):0.0;
    double scalingFactor = InstrumentUtil::scalePremium(oneContract,fwdStarting,notional,fwdStrt,initialSpot);
    double histContracts = 0.0;
    double prodVal = 0.0;
    double strikes = strike*initialSpot;
    if(decimalRound>=0) strikes = decRound(strike * initialSpot,decimalRound);

    try  
    {
        bool deadInstrument = false;
  
        DateTime settlementDate = paymentDates->back();
        DateTime settleDate = instSettle->settles(settlementDate, asset.get());
         
        if(instSettle->isPhysical()&& !phyAndCashSettle->isPhysicalAndCash){
            if(isKO||isHit||(valueDate >= settlementDate))
                {   
                    DateTime endDate = payAtHit? hitDate:getPayDate(knockedOutDate);
                    if (valueDate >= endDate)
                        {
                            results->storePrice(0.0, discount->getCcy());
                            addOutputRequests(control, results, 0.0, 0.0);
                            deadInstrument = true;
                        }
                }
        }
        else{ // cash settle or isPhysicalAndCash settle
     
                DateTime endDate = settlementDate;
                DateTime endSettleDate = settleDate;
                double spotAtMat = obsSamples->back();

                if (isKO||isHit){
                    endDate = payAtHit? hitDate:getPayDate(knockedOutDate);
                    endSettleDate = instSettle->settles(endDate, asset.get());
                    spotAtMat = spotAtHit;
                }
                
                if ((valueDate >= endDate) && (valueDate < endSettleDate)){
                    histContracts = getKoHistContracts(endDate,true);
                    prodVal = productValue(spotAtMat, strikes);
                    double value = histContracts*prodVal*scalingFactor;
                    value *= discount->pv(valueDate, endSettleDate);
                    // physical weight and cash weight does NOT neccessarily add up to 1, in phyDeliveredAtSpot case, the cash weight is always 100%...
                    // while physical weight determines the physical delivery portion
                    // if it is physical and cash settled, AND physical delivered at Spot, the remaining portion is (1 - physical weight)
                    value = phyAndCashSettle->isPhysicalAndCash?(phyAndCashSettle->phyDeliveredAtSpot? value*(1-phyAndCashSettle->physicalWeight): value* phyAndCashSettle->cashWeight):value;
                    results->storePrice(value, discount->getCcy());
                    addOutputRequests(control, results, 0.0, 0.0);
                    deadInstrument = true;
                }

                if (valueDate >= endSettleDate){ // for cash settle
                    results->storePrice(0.0, discount->getCcy());
                    addOutputRequests(control, results, 0.0, 0.0);
                    deadInstrument = true;
                }
        }
 
        // record KNOWN_CASHFLOWS
        if ( control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS) ) {
            CashFlowArraySP allCFs = getKnownCashFlows();

            OutputRequestUtil::recordKnownCashflows(control,
                results,
                discount->getCcy(),
                allCFs.get()); 
        }
        
        // PAYMENT_DATES
        if ( control->requestsOutput(OutputRequest::PAYMENT_DATES) ) {
            DateTimeArraySP dates = getPaymentDates();
            
            OutputRequestUtil::recordPaymentDates(control,results,dates.get()); 
        }
        
        // PHYSICAL_DELIVERY
        if (instSettle->isPhysical()) 
        {
            PhysicalDeliveryArray pda = getPhysicalDelivery();
            if ( control->requestsOutput(OutputRequest::PHYSICAL_DELIVERY) ) 
            {
                PhysicalDelivery::recordPhysicalDelivery(control,
                                                             results,
                                                             asset->getTrueName(),
                                                             &pda);
            }                
        }

        return deadInstrument;
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }          
}

/** Rolls the value date and sets initial spot if rolling over start date */
bool DailyAccumulator::sensShift(Theta* shift)
{    
    DateTime newDate = shift->rollDate(valueDate);

    DateTime matDate = paymentDates->back(); 

    /* also need to fill in historical levels for any samples that are in the
           future with respect to the value date but are in the past with respect
           to the new date */

    if (!obsSamples)
    {
        obsSamples = DoubleArraySP( new DoubleArray(monitorDates->size(), 0.0) );
    }
    else
    {
        obsSamples->resize(monitorDates->size());
    }


    DoubleArray wgts(monitorDates->size(), 1.0);
    SampleList monDatesValues(*monitorDates, *obsSamples, wgts);
    monDatesValues.roll(asset.get(),valueDate,newDate,true); 

    monitorDates = DateTimeArraySP( new DateTimeArray(monDatesValues.getDates()) );
    obsSamples = DoubleArraySP( new DoubleArray(monDatesValues.getValues()) );

    // roll today 
    valueDate = newDate;
    
    return true;
};

/** Satisfy LegalTerms::Shift interface */
bool DailyAccumulator::sensShift(LegalTerms* shift) {
    // Set the barriers for pricing equal to the economic barriers
    // legalBarrier is the same size as barrierLevel (checked in Validate())
    {     
    barrierLevel = legalBarrier;
    }

    return true;
}

/** when to stop tweaking */
DateTime DailyAccumulator::endDate(const Sensitivity* sensControl) const {
    DateTime matDate = monitoringEndDates->back();
    DateTime instEnd  = instSettle->settles(matDate, asset.get());
    DateTime assetEnd = asset->settleDate(matDate);
    DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
    return end;
}

void DailyAccumulator::addOutputRequests(Control* control,
                                 Results* results,
                                 const double& fairValue,
                                 const double& indVol) const
{
    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        DateTime matDate = paymentDates->back();
        // DELAY_PRICE
        InstrumentUtil::delayPriceHelper(control,
                                         results,
                                         fairValue,
                                         valueDate,
                                         discount.get(),
                                         asset.get(),
                                         premiumSettle.get());
        // IND_VOL
        InstrumentUtil::recordIndicativeVol(control,results,indVol);
        // FWD_AT_MAT
        try{
            InstrumentUtil::recordFwdAtMat(control,
                                           results,
                                           matDate,
                                           valueDate,
                                           asset.get());
        }
        catch(exception&)
        {// continue if fwd failed - this hapens now for quanto asset with CEVj vol
        }
    }
}

// product class
class DailyAccumulatorFDProd: public LatticeProdEDRIns
{
public:
    DailyAccumulatorFDProd( const DailyAccumulator * acc, FDModel * mdl ) :
        LatticeProdEDRIns( mdl, 1, acc->isGlobalKO ? 3 : 4 ),
        inst( acc )
    {
        // second: create spot payoff
        payoffIndex = model->createProduct( IProdCreatorSP( new
            IndexSpecEQ( inst->asset.getName(), inst->asset, inst->ccyTreatment ) ) );

        // retrieve div dates between value date and last payment date if necessary
        if(inst->includeDivs)
        {
            DividendListSP exDivs = AssetUtil::getAllDivsBetweenDates(inst->asset, 
                                                 inst->valueDate, 
                                                 inst->paymentDates->back());

            DateTimeArrayConstSP exDivDates = exDivs->getExDivDates();
            DoubleArrayConstSP   exDivAmounts = exDivs->getDivAmounts();
            // change div time to be able to insert into tree
            if (exDivDates->size() >0)
            {
                divDates.resize(exDivDates->size());
                divAmounts.resize(exDivDates->size());
                for (int j=0; j<divDates.size();j++)
                {
                    divDates[j] = DateTime( (*exDivDates)[j].getDate(), DateTime::START_OF_DAY_TIME);
                    divAmounts[j] = (*exDivAmounts)[j];
                }
            }
        }
    }

    /** vanilla, quanto or struck */
    virtual string getCcyTreatment() const
    {
        return inst->ccyTreatment;
    }

    /** ignore start date if not forward starting */
    virtual DateTime getStartDate() const
    {
        return inst->fwdStarting ? inst->startDate : inst->valueDate;
    }

    /** this sets up the timeline */
    virtual void init( CControl * control ) const
    {
        static const string method = "DailyAccumulatorFDProd::Init";
        try {
            int i;
            int j = 0;

            // all exe dates are critical dates
            DateTimeArray critDates;
            int num = inst->monitorDates->size();
            for (i=0; i<num; i++) {          
                critDates.push_back((*(inst->monitorDates))[i]);
                if ((*(inst->monitorDates))[i] <= (*(inst->monitoringEndDates))[j]){
                    critDates.push_back((*(inst->paymentDates))[j]);
                }
                if ((*(inst->monitorDates))[i] == (*(inst->monitoringEndDates))[j]) {j++;}
            }

            if( tree1f ) //here
            {
                // default to NODE_INSERTION smoothing
                if (tree1f->GetSmoothMethod() == CTree1f::DEFAULT) {
                    tree1f->SetSmoothMethod(CTree1f::NODE_INSERTION);
                }

                tree1f->NumOfPrice = numPrices;
                tree1f->NumOfInsertNode =
                    ( tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION ? numIns : 0 );

                // default to DollarDivTreament = false !!! to be removed after EAS/EIS can handle it correctly !!!
                tree1f->SetDivAmountTreatment(false);

                if( inst->fwdStarting )
                    tree1f->controlSameGridFwdStart(inst->ccyTreatment);        
            }

            // add critical dates
            model->addCritDates( critDates );
            model->addCritDates( divDates );

            // define start and end points of tree
            DateTimeArray segDates(2);
            segDates[0] = getStartDate();
            segDates[1] = inst->paymentDates->back();

            // timeline configuration
            // 'density factor' for timeline
            IntArray density( 1, 1 );

            // prepare timeline set up
            model->initSegments( segDates, density );
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** initialising and setting product variables */
    // this is called per pricing call before tree sweep call (after InitTree)
    virtual void initProd()
    {
        static const string method = "DailyAccumulatorFDProd::initProd";
        try {
            int i = 0;
            int j = 0;
            // timeline exists at this point
            int lastStep = model->getLastStep();

            initSlices( numPrices );
            initInsertNode();

            DateTime endDate(inst->valueDate.getDate()+364, inst->valueDate.getTime()); 

            stepIsMonitoringDate.resize(lastStep+1, false);
            stepIsKODate.resize(lastStep+1, false);
            stepIsPayDate.resize(lastStep+1, false);
            stepStrikeLevel.resize(lastStep+1, 0.0);
            stepParticipation.resize(lastStep+1, 0.0);
            stepBarrierLevel.resize(lastStep+1, 0.0);
            stepRebateLevel.resize(lastStep+1, 0.0);
            stepIsDividendDate.resize(lastStep+1, false);
            stepDividend.resize(lastStep+1, 0.0);
            stepExtraParticipation.resize(lastStep+1, 0.0);

            // past call dates are not taken into account

            DateTimeArraySP monitorDatesSP = inst->monitorDates;
            DateTimeArray dts = *monitorDatesSP;
             int monitorIdx = inst->getIndexOnOrAfterDate(dts, model->getDate(0));

             paymentDates.resize(0);
             barrierLevel = inst->barrierLevel * inst->initialSpot;
             rebate.resize(0);
             participation = inst->participation;
             strikes = inst->strike * inst->initialSpot;
             extraParticipation = inst->extraParticipation;
             legalBarrier = inst->legalBarrier * inst->initialSpot;

             if(inst->decimalRound>=0)
             {
                 barrierLevel = inst->decRound(barrierLevel,inst->decimalRound);
                 strikes = inst->decRound(strikes,inst->decimalRound);
                 legalBarrier = inst->decRound(legalBarrier,inst->decimalRound);
             }
                      
            
             for (i=0; i< inst->monitorDates->size(); i++) {          
                
                    paymentDates.push_back((*(inst->paymentDates))[j]);
                    iRebate = inst->uniformRebate? (inst->rebate * inst->initialSpot): (inst->rebate - (inst->rebate)*(i+1)/(inst->monitorDates->size()))*inst->initialSpot;
                    rebate.push_back(iRebate);
                  
                if ((*(inst->monitorDates))[i] > (*(inst->monitoringEndDates))[j]) {j++;}
            }

            for (i=0; i <= lastStep; i++) 
            {
                const DateTime treeDate = model->getDate(i);

                if (inst->isUpAndOut)
                {    
                    stepBarrierLevel[i] = 9999999999999.0 * inst->initialSpot;
                }
                else
                {
                    stepBarrierLevel[i] = 0.0;
                }
                
                if (monitorIdx < dts.size() && 
                    treeDate.equals( dts[monitorIdx],false ) )
                {
                    stepBarrierLevel[i] = barrierLevel;
                    stepRebateLevel[i] = rebate[monitorIdx];
                    stepParticipation[i] = participation;
                    stepIsKODate[i] = true; /* discrete monitoring uses modified continous monitoring barrier level
                                                 but the algorithmic treatment is the same */
                    stepStrikeLevel[i] = strikes;
                    if( inst->hasExtraParticipation() )
                        stepExtraParticipation[i] = extraParticipation;

                    if (treeDate.equals( dts[monitorIdx] ) ) 
                    // here date AND time must match for declaring the step to be a monitoring date
                    {
                        stepIsMonitoringDate[i] = true; 
                        monitorIdx ++;
                    }
                }            
            }   
            
            // past payment dates are not taken into account 
            int payIdx = 0;
            while ( (payIdx < paymentDates.size())
                        && (paymentDates[payIdx] < model->getDate(0)) )
            {
                payIdx ++;
            }

            for (i=0; i<=lastStep; i++) 
            {
                const DateTime treeDate = model->getDate(i);

                stepIsPayDate[i] = false;

                if ( payIdx < paymentDates.size() && 
                    treeDate.equals(paymentDates[payIdx]) )
                {
                    stepIsPayDate[i] = true; 
                    stepStrikeLevel[i] = strikes;
                    
                    //payIdx ++;
                    int k = payIdx;
                    
                    while ( (payIdx < paymentDates.size())
                        && (paymentDates[payIdx]).equals( paymentDates[k]))
                    {
                        payIdx ++;
                    }        
                }
            }

            // dealing with dividends if necessary
            if(inst->includeDivs)
            {
                    DividendListSP exDivs = AssetUtil::getAllDivsBetweenDates(inst->asset, 
                                                         inst->valueDate, 
                                                         paymentDates[paymentDates.size() - 1]);

                    DateTimeArrayConstSP exDivDates = exDivs->getExDivDates();
                    DoubleArrayConstSP   exDivAmounts = exDivs->getDivAmounts();
                
                    if (exDivDates->size() >0)
                    {
                            divDates.resize(exDivDates->size());
                            divAmounts.resize(exDivDates->size());
                            for (j=0; j<divDates.size();j++)
                            {
                                divDates[j] = DateTime( (*exDivDates)[j].getDate(), DateTime::START_OF_DAY_TIME);
                                divAmounts[j] = (*exDivAmounts)[j];
                            }
                    }
            }

            if(inst->includeDivs && divDates.size() > 0) 
            {
                // past exDividend dates are not taken into account
                int divIdx = 0;
                while ( divIdx < divDates.size() &&
                    ( divDates[divIdx] < model->getDate(0)) ){
                    divIdx ++;
                }
                
                for (i=0; i<=lastStep; i++) 
                {
                    const DateTime treeDate = model->getDate(i);

                    stepIsDividendDate[i] = false;
                    stepDividend[i] = 0.0;

                    if ( divIdx < divDates.size() && 
                        treeDate.equals( divDates[divIdx]) )
                    {
                        stepIsDividendDate[i] = true; 
                        stepDividend[i] = divAmounts[divIdx];
                        divIdx ++;
                    }
                }
            }

            //     computing the number of contracts in the current monitoring period
            //   this needs a definition of what a monitoring period is
            KoHistContracts = inst->getKoHistContracts(inst->valueDate);

        
            // compute knownCF value which is not already paid out and which is determined before today
            // physical delivery value is dropped on trade date, so no need to compute known value
            CashFlowArraySP knownCF = inst->getKnownCashFlows();
            knownCFValue = 0;
            if( !!knownCF && knownCF->size() != 0 )
            {
                DateTime valDate = inst->valueDate;
                DateTime cutoffDate = inst->instSettle->settles(valDate, inst->asset.get());
                for (i=0; i<knownCF->size(); i++)
                {
                    if( (*knownCF)[i].date < cutoffDate && (*knownCF)[i].date > inst->valueDate )
                        knownCFValue += (*knownCF)[i].amount * inst->discount->pv(inst->valueDate, (*knownCF)[i].date);
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // here just take care of scaling by notional and any additional 
    // discounting for fwd starting case
    double scalePremium(const double& fairValue, YieldCurveConstSP disc)
    {
        double fwdAtStart = 0.0;
        double fwdStartDF = 1.0;
        if (inst->fwdStarting)
        {
            fwdAtStart = inst->asset->fwdValue(inst->startDate);
            fwdStartDF = disc->pv(inst->valueDate, model->getDate(0));
        }

        double scalingFactor = InstrumentUtil::scalePremium(
            inst->oneContract,
            inst->fwdStarting,
            inst->notional,
            fwdAtStart,
            inst->initialSpot);

        return fairValue * scalingFactor * fwdStartDF + knownCFValue;
    }

    /** extra output requests */
    virtual void recordOutput(Control* control, YieldCurveConstSP disc, Results* results)
    {
        // get prices at t=0
        // save price
        double price = scalePremium(model->getPrice0( *slices[0] ), disc);
        results->storePrice(price, disc->getCcy());

        // take care of additional outputs
        if ( control && control->isPricing() )
        {
            DateTime       matDate = (*(inst->paymentDates))[(inst->paymentDates)->size() - 1]; 
            double         indVol;
            // calculate indicative vol
            try{
                if ( matDate.isGreater(inst->valueDate) )
                {

                    DateTime imntStartDate = inst->fwdStarting? 
                        inst->startDate: inst->valueDate;

                    // get vol request
                    CVolRequestConstSP lnVolRequest = GetLNRequest();

                    // interpolate the vol
                    CVolProcessedSP  vol(inst->asset->getProcessedVol(lnVolRequest.get()));
                    // cast to the type of vol we're expecting
                    CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);
                    // this should never happen if our get market data has worked properly
                    if (!vol){
                        throw ModelException("VanillaTreeFDProd::recordOutputRequests", 
                                             "No Black Scholes Vol");
                    }

                    // dealing with the barrier output
                    OutputRequest* request = NULL;
                    request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
                    if (request) 
                    {
                        // report ladder level today 
                    //    double spotScale = inst->initialSpot;
                        BarrierLevelArray levels;
                        double barlevel;
                        DateTime barDate;
                
                        // store all ladder level as upper barrier
                       
                        DateTimeArraySP monitorDates = inst->monitorDates;
                        double size = monitorDates->size();
                        double legalBarrier = inst->legalBarrier * inst->initialSpot;
                        if(inst->decimalRound>=0) legalBarrier = inst->decRound(legalBarrier,inst->decimalRound);

                        // use the barrierWindow routine 
                        DateTime upperDate = BarrierLevel::barrierWindow(inst->valueDate);


                        for (int i=0; i<size; i++)
                        {                             
                            barDate = (*monitorDates)[i];
                            barlevel = legalBarrier;
                          
                            // record only future barrier levels
                            if ((inst->valueDate < barDate) && (barDate < upperDate))
                            {
                                levels.push_back(BarrierLevel(inst->isUpAndOut,barDate, barlevel,inst->isContinuouslyMonitored));
                            }
                        }

                        OutputRequestUtil::recordBarrierLevels(control,
                                                               results,
                                                               inst->asset->getTrueName(),
                                                               &levels);
                    } 

                    // calculate the indicative vol
                    indVol = volBS->CalcVol(imntStartDate, matDate);
                }
                else
                {
                    indVol = 0.0;
                }
            }
            catch(exception&)
            {// continue if indicative vol fail - non BS vol
                indVol = 0.0;
            }

            inst->addOutputRequests(control,
                                    results,
                                    price,
                                    indVol);

        }
    }
    
    /** returns a vol request for log-normal vol */
    //for set up timeline only,  to be reviewed/removed */
    virtual CVolRequestConstSP GetLNRequest() const
    {
        // get strike and maturity date from instrument
        DateTime matDate   = inst->monitoringEndDates->back();
        double   volStrike = inst->strike * inst->initialSpot;

        CVolRequestConstSP volRequest(
            new LinearStrikeVolRequest(volStrike, getStartDate(), 
                                       matDate, inst->fwdStarting));
        return volRequest;
    }

    /** called before each update() */
    virtual void preCalc(int step)
    {
        static const string method("DailyAccumulatorFDProd::preCalc");
        try {

            if( tree1f )
            {
                int idx = tree1f->getSliceIndex(step);
                const DateTime treeDate = model->getDate(step);

                if ( ( treeDate <= paymentDates[paymentDates.size() - 1])) 
                                                    // && (stepIsKODate[step]) )
                                                  
                {
                    if (treeDate.equals(paymentDates[paymentDates.size() - 1]))
                    {
                        //KOlevel =  inst->barrierLevel * inst->initialSpot;
                        KOlevel =  barrierLevel;
                    }
                    else
                    {
                        KOlevel =  stepBarrierLevel[step];
                    }

                    // making a correction in case barrier monitoring is discrete
                    if (!inst->isContinuouslyMonitored) {
                        vector<double> vol;
                        // adjust barrier if needed
                        tree1f->GetStepVol(step, vol, &KOlevel, 0, 0); // get vol at barrier
                        Barrier::BarrierAdjustment(vol[0], inst->isUpAndOut, KOlevel);         
                    }

                    if (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION) 
                    { //assuming Lower Barrier is index 0 and Upper uses 1 !!!
                        tree1f->SetInsertNode(idx, 0, KOlevel, 0); // insert barrier level
                    }
                }   
                else
                {
                    if (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION)
                    { 
                        tree1f->SetInsertNode(idx, 0, -1.0, 0); // if not a barrier date, insert dummy;
                    }
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** isInitValue == true, payoff at T for backward or value at t=0 for fwd induction
        isInitValue == false, payoff boundary condition, for KO, early exercise etc. */
    virtual void update(int& step, FDProduct::UpdateType type)
    {
        // we assume just need one und level for spot here
        const TreeSlice & spot = payoffIndex->getValue( step );
        int bot, top;
        spot.getCalcRange( bot, top );

        const vector< TreeSliceSP > & price = slices;
        int pStart = 0, pEnd = price.size() - 1;

        if (type == FDProduct::BWD_T){
            prod_BWD_T( spot,
                        step,
                        bot,
                        top,
                        pStart, 
                        pEnd,
                        price);

            //insert nodes
            if (tree1f && tree1f->NumOfInsertNode>0)
            {
                prod_BWD_T( *insNodes,
                            step,
                            0,
                            tree1f->NumOfInsertNode-1,
                            pStart, 
                            pEnd,
                            *insPrices);
            }

        }
        else if(type == FDProduct::BWD){
            prod_BWD( spot,
                      step,
                      bot,
                      top,
                      pStart, 
                      pEnd,
                      price);

            //insert nodes
            if (tree1f && tree1f->NumOfInsertNode>0)
            {
                prod_BWD( *insNodes,
                          step,
                          0,
                          tree1f->NumOfInsertNode-1,
                          pStart, 
                          pEnd,
                          *insPrices);

            }
        }
    }

    /** product payoff method at maturity */
    void prod_BWD_T(
        const TreeSlice & spot,
        int step,
        int bot,
        int top,
        int pStart,
        int pEnd,
        const vector< TreeSliceSP > & price)
    {
        static const string method("DailyAccumulatorFDProd::prod_BWD_T");
        try {
            double * s = spot.getValues();
            const vector< double * > & p = getValues( price );

            int j;   

            // get strike and maturity date from instrument
            DateTime matDate = inst->monitoringEndDates->back();
            //double matStrike  = inst->strike * inst->initialSpot;
            double matStrike  = strikes;
            //DateTimeArraySP paymentDates = inst->paymentDates;
            const DateTime payDate = paymentDates.back();    
            //double matParticipation = inst->participation;
            double matParticipation = participation;

            nextPayDate = payDate;
            double barLevel = KOlevel; //inst->barrierLevel.back();
            double rebates = inst->uniformRebate?(inst->rebate * inst->initialSpot):0.0;
            int prod = inst->productAccumulated;

            // initializing: (*price[1]) is the forward contract payoff (will be useful for progressive style case)
            // (*price[0]) is the product value and (*price[3]) is an auxiliary vector
            // used to store product future value at strike dates

            (*price[1]) = (prod-2)*(prod-3)*0.5*(spot - matStrike) 
                     - (prod-1)*(prod-3)*smax(spot - matStrike, 0.)
                     + (prod-1)*(prod-2)*0.5*smax(matStrike - spot, 0.);
            // (*price[1]) is the  "contract for settlement purposes", and can be based on a forward, call or put
            
            
            (*price[2]) = (*price[1]); // (*price[2]) is the  "contract for monitoring purposes"
            (*price[0]) = 0.; // value will start to show up at the last monitoring dates
            if( !inst->isGlobalKO )
                (*price[3]) = 0.; // future remaining value

            if (matDate.equals(payDate))
            {
                for (j=bot; j<=top; j++)
                {
                    // checking if barrier is breached just now
                    if( (inst->isUpAndOut&& (s[j] > barLevel*(1-FP_MIN))) || (!inst->isUpAndOut && (s[j] < barLevel*(1+FP_MIN))) )
                    {        
                        p[0][j] = rebates;
                    }
                    else
                    {
                        p[0][j] = p[1][j]*matParticipation;
                        if( inst->hasExtraParticipation() )
                        {
                            if( (inst->isUpAndOut && s[j]<matStrike) || (!inst->isUpAndOut && s[j]>matStrike) )
                                p[0][j] += p[1][j] * inst->extraParticipation;
                        }
                    }
                }
            }
        }
        catch(exception& e) {
            throw ModelException(e, method);
        }
    }

    // this is payoff boundary condition, for KO, early exercise etc.
    void prod_BWD(
        const TreeSlice & spot,
        int step,
        int bot,
        int top,
        int pStart,
        int pEnd,
        const vector< TreeSliceSP > & price)
    {
        static const string method("DailyAccumulatorFDProd::prod_BWD");
        try {
            double * s = spot.getValues();
            const vector< double * > & p = getValues( price);

            int j;
            const DateTime treeDate = model->getDate(step);    

            double barLevel = KOlevel;
            double rebates = stepRebateLevel[step];
            double settlementPV;
            int prod = inst->productAccumulated;
    
    
            // if we are at a settlement date, (*price[1])[j] is computed
            if(stepIsPayDate[step])
            {
                nextPayDate = treeDate;
                (*price[1]) = (prod-2)*(prod-3)*0.5*(spot - stepStrikeLevel[step]) 
                         - (prod-1)*(prod-3)*smax(spot - stepStrikeLevel[step], 0.)
                         + (prod-1)*(prod-2)*0.5*smax(stepStrikeLevel[step] - spot, 0.);
            }

            // if we are at a dividend date, (*price[2])[j] incorporates the dividend payment
            if(inst->includeDivs && stepIsDividendDate[step])
            {
                (*price[1]) += stepDividend[step];
                (*price[2]) += stepDividend[step];
            }

            // if we are at the end of a monitoring period, several operations are performed:
            if( (step > 0) && stepIsMonitoringDate[step] )
            {
                if (inst->isEndMonitoringPeriod(treeDate))
                {
                    if(!inst->isGlobalKO)
                        (*price[3]) = (*price[0]); // (*price[3]) is updated if necessary
            
                    (*price[2]) = (*price[1]); /* we're getting into the previous monitoring date 
                    hence "monitoring" fwd reverts to "settlement" fwd as this one has already been updated */
                }
            }

            settlementPV = inst->discount->pv(treeDate,nextPayDate);
            
            if(stepIsKODate[step]) // if continuously monitored case
            {
                for (j=bot; j<=top; j++)
                {
                    if( (inst->isUpAndOut&& (s[j] > barLevel*(1-FP_MIN))) || (!inst->isUpAndOut && (s[j] < barLevel*(1+FP_MIN))) )
                    {            
                        if(inst->payAtHit) settlementPV = 1.0;
                        p[0][j] = (inst->dayKoOnly?p[0][j]:(inst->isGlobalKO?0.0:p[3][j])) + rebates*settlementPV;
                    }
                }
            }

            if(stepIsMonitoringDate[step] && step > 0)
            {
                for (j=bot; j<=top; j++)
                {
                    // adding payment of the forward contract maturing today and checking Knock-Out
                    if( (inst->isUpAndOut&& (s[j] >barLevel*(1-FP_MIN))) || (!inst->isUpAndOut && (s[j] <barLevel*(1+FP_MIN))) )
                    {   
                        if(inst->payAtHit) 
                        {
                            p[2][j] = inst->productValue(s[j], stepStrikeLevel[step]);
                            settlementPV = 1.0;
                        }
                        p[0][j] =  (inst->dayKoOnly?p[0][j]:(inst->isGlobalKO?0.0:p[3][j])) + rebates*settlementPV;
                    }
                    else
                    {
                        p[0][j] += p[2][j]*stepParticipation[step];
                    }
                }

                if( inst->hasExtraParticipation() )
                {
                    if( inst->isUpAndOut )
                    {
                        j=bot; 
                        while(j<=top && s[j] < stepStrikeLevel[step] )
                        {
                            p[0][j] += p[2][j]*stepExtraParticipation[step];
                            j++;
                        }
                    }
                    else
                    {
                        j=top; 
                        while(j>=bot && s[j] > stepStrikeLevel[step] )
                        {
                            p[0][j] += p[2][j]*stepExtraParticipation[step];
                            j--;
                        }
                    }
                } 
            } 
            
            // adding historical contracts
            if (step ==0)  
            {
                // case current forward is already Knocked-Out: KoRebate and KOFractionElapsed have already
                // been computed so total price can be determined as:
                if (inst->isKO)
                {
                    if( inst->isGlobalKO )
                        (*price[0]) = 0.;
                    else
                        (*price[0]) = (*price[3]);
                    
                    (*price[0]) += KoHistContracts*(*price[2]) + inst->KoRebate*settlementPV;
                }
                else
                {
                    (*price[0]) += KoHistContracts*(*price[2]);
                }

            }    
        }
        catch(exception& e) {
            throw ModelException(e, method);
        }
    }
    
private:
    const DailyAccumulator* inst;
    vector<bool>      stepIsMonitoringDate; // whether a step is a strike date
    vector<bool>      stepIsKODate; // whether a step is a KO observation date
    vector<bool>      stepIsPayDate; // whether a step is a payment date
    vector<bool>      stepIsDividendDate; // whether a step is a dividend payment date
    vector<double>    stepStrikeLevel; // level of strike on step
    vector<double>    stepParticipation; // level of participation on step
    vector<double>    stepBarrierLevel; // level of barrier on step
    vector<double>    stepRebateLevel; // level of rebate on step
    vector<double>    stepDividend; // amount of dividend paid
    vector<double>    stepExtraParticipation; // level of participation on step for extraParticipation

    double              KoHistContracts;
    double              KOlevel;
    DateTime            nextPayDate;
    DateTimeArray       divDates;
    DoubleArray         divAmounts;

    double              knownCFValue;

    DateTimeArray       paymentDates;
    double              barrierLevel;
    DoubleArray         rebate;
    double              participation;
    double              strikes;
    double              extraParticipation;
    double              legalBarrier;
    double              iRebate;
};

/** create a fd payoff tree - new for all fd/tree state variable interface */
FDProductSP DailyAccumulator::createProduct(FDModel* model) const
{
    return FDProductSP( new DailyAccumulatorFDProd(this, model) );
}

/** returns all strikes on the vol surface to which 
    this instrument is sensitive */
DoubleArraySP DailyAccumulator::getSensitiveStrikes(OutputNameConstSP outputName,
                                               const IModel*      model)
{
    DoubleArraySP sensStrikes = DoubleArraySP(new DoubleArray(0));

    if (avoidVegaMatrix(model)) {
        throw ModelException("DailyAccumulator::getSensitiveStrikes", 
                             "VEGA_MATRIX is not valid for this instrument");
    }

    // get start date for vol interpolation
    DateTime imntStartDate = fwdStarting ? startDate:valueDate;
    
    // get last exercise date in exercise schedule
    DateTime  maturityDate =  monitoringEndDates->back();

    double strikes = strike * initialSpot;

    // create a vol request object to be passed on
    LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(strikes, 
                                                                   imntStartDate, 
                                                                   maturityDate,
                                                                   fwdStarting));

    SensitiveStrikeDescriptor sensStrikeDesc;
    sensStrikeDesc.forwardOnly = false;

    asset->getSensitiveStrikes(volRequest.get(), outputName, 
                               sensStrikeDesc, sensStrikes);

    return sensStrikes;
}

extern bool DailyAccumulatorLoad()
{
    return true && DailyAccumulator::TYPE;
}

DRLIB_END_NAMESPACE
