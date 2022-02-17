//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Accumulator.cpp
//
//   Author      : Jay Blumenstein
//
//   Description : Accumulator option.
//                 A stream of forward contracts with a set of Knock-out conditions.
//
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
#include "edginc/SampleList.hpp"
#include "edginc/PhysicalDelivery.hpp"
#include "edginc/VegaParallel.hpp"

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"
#include "edginc/TreeSliceOper.hpp"

DRLIB_BEGIN_NAMESPACE

class Accumulator : public Generic1Factor, 
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

    double productValue(double spot, double strike) const
    {
        double val = spot - strike;
        if( productAccumulated == 2 )
            val = Maths::max(0.0, val);
        else if( productAccumulated == 3 )
            val = Maths::max(0.0, -val);
        return val;
    }

    DateTime getLastDate() const
    {
         return monitoringDatesAndStrikes.back().date;
    }

    /** return first index in sorted dates that is on or after aDate */
    int getIndexOnOrAfterDate(const DateTimeArray& dates, const DateTime& aDate) const
    {
        int idx = 0;
        while ( idx < dates.size() && dates[idx] < aDate ){
            idx ++;
        }
        return idx;
    }

    /** find pay date in PaymentDates */
    const DateTime& getPayDate(const DateTime& today) const
    { 
        int i = 0;
        bool found = false;

        DateTimeArraySP monitorDates = getDatesFromCFL(monitoringDatesAndStrikes);
        
        while (!found &&  i < monitorDates->size() )
        {
            if ( (*monitorDates)[i].equals(today) ) 
            {
                found = true;
            }

            i++;
        }

        if (found)
        {
            return (*paymentDates)[i - 1];
        }
        else
        {
            throw ModelException("Accumulator::getPayDate", "date supplied is not in monitoring dates!");
        }
    }

    // Gives the payment date corresponding to the monitoring date that is equal to 
    // or the first after today            
    DateTime getNextPaymentDate(const DateTime& today) const
    {
        DateTimeArraySP monitorDates = getDatesFromCFL(monitoringDatesAndStrikes);
        DateTime nextPayDate;

        bool found = false;
        int i = 0;
        while (!found && i < monitorDates->size())
        {
            if ( today <= (*monitorDates)[i] )
            {
                found = true;
                nextPayDate = (*paymentDates)[i];
            }
            i++;
        }
        if (!found) // in this case, we are after all monitoring dates so result here is last payment date
        {
            nextPayDate = (*paymentDates)[paymentDates->size()-1];
        }
    
        return nextPayDate;
    }

    // Gives the monitoring date corresponding to the monitoring date that is equal to 
    // or the first after today            
    DateTime getNextMonitoringDate(const DateTime& today) const
    {
        DateTimeArraySP monitorDates = getDatesFromCFL(monitoringDatesAndStrikes);
        DateTime nextMonitorDate;

        bool found = false;
        int i = 0;
        while (!found && i < monitorDates->size())
        {
            if ( today <= (*monitorDates)[i] )
            {
                found = true;
                nextMonitorDate = (*monitorDates)[i];
            }
            i++;
        }
        if (!found) // in this case, we are after all monitoring dates so result here is last monitoring date
        {
            nextMonitorDate = (*monitorDates)[monitorDates->size()-1];

        }
    
        return nextMonitorDate;
    }

    // computes historical contracts observed strictly prior to current moment evalDate
    // the boolean includeCurrentObservation is to include historical contract observed exactly at evalDate
    double getKoHistContracts(const DateTime& evalDate, const bool includeCurrentObservation =false) const
    {
        int i = 0;
        double nbContracts = 0.0;
        double pastSpotValue = 0.0;
        double barLevel;
        bool checkCondition;

        KoRebate = 0.0;
        isKO = false;

        DateTime currentDate;
        DateTime settleDate;
        DateTime PayDate = getNextPaymentDate(evalDate);
        DateTimeArraySP monitorDates = getDatesFromCFL(monitoringDatesAndStrikes);

        while ( i<monitorDates->size() && ((*monitorDates)[i] <= evalDate) && (!isKO) ) 
        {
            currentDate = (*monitorDates)[i];
            settleDate = getPayDate(currentDate);
            
            pastSpotValue = (*histMonSamples)[i];
            barLevel = barrierLevel[i];

            // in case of global ko past monitoring period knock-outs are taken into account
            if(isGlobalKO)
            {   
                if( (isUpAndOut&& (pastSpotValue >= barLevel)) || (!isUpAndOut && (pastSpotValue <= barLevel)) )
                {            
                    isKO = true;
                    knockedOutDate = DateTime(currentDate);
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
                if( (isUpAndOut&& (pastSpotValue >= barLevel)) || (!isUpAndOut && (pastSpotValue <= barLevel)) )
                {            
                    isKO = true;
                    knockedOutDate = DateTime(currentDate);
                    //KoDate = currentDate;
                    KoRebate = rebate[i];
                }
                if (!isKO)
                {
                    nbContracts +=  participation[i];

                    if( hasExtraParticipation() )
                    {
                        double extraStrike = monitoringDatesAndStrikes[i].amount;
                        if( (isUpAndOut && pastSpotValue<extraStrike ) || (!isUpAndOut && pastSpotValue>extraStrike) )
                            nbContracts +=  (*extraParticipation)[i];
                    }
                }            
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
        bool nextMonitorDateFound = false;
        DateTime nextMonitoringDate;

        DateTimeArraySP monitorDates = getDatesFromCFL(monitoringDatesAndStrikes);

        // checking whether testDate is a monitoring date
        while ( i<monitorDates->size() && !monitoringDateFound )
        {
            if( testDate.equals((*monitorDates)[i]) )
            {
                monitoringDateFound = true;
                if (i<monitorDates->size()-1)
                { 
                    nextMonitorDateFound = true;
                    nextMonitoringDate = (*monitorDates)[i+1];
                }            
            }
            i++;
        }

        /* part commented out to return false if date is not a monitoring period
        if (!monitoringDateFound)
        {
            throw ModelException("Accumulator::isEndMonitoringPeriod", "date entered is not a monitoring date");
        }*/

        // now looking whether we are at a payment date
        if (monitoringDateFound)
        {
            if ( nextMonitorDateFound && ((*paymentDates)[i-1] < (*paymentDates)[i]) )
            {
                isEndPeriod = true;
            }
            if ( !nextMonitorDateFound)
            {
                isEndPeriod = true;
            }        
        }

        return isEndPeriod;
    }

    DateTime getEndCurrentMonitoringPeriod(const DateTime& today) const
    {
        DateTimeArraySP monitorDates = getDatesFromCFL(monitoringDatesAndStrikes);
        DateTime endMonitoringDate;
        DateTime payDate;

        bool foundMonitor = false;
        bool foundEndPeriod = false;
        int i = 0;
        int j = 0;

        while (!foundMonitor && i < monitorDates->size())
        {
            if ( today.equals((*monitorDates)[i]))
            {
                foundMonitor = true;
                payDate = (*paymentDates)[i];
            }
            i++;
        }

        if(foundMonitor)
        {
            j = i;
            while (!foundEndPeriod && j < monitorDates->size())
            {
                if( j== monitorDates->size()-1)
                {
                    foundEndPeriod = true;
                    endMonitoringDate = (*monitorDates)[monitorDates->size()-1];
                }
                if( (*paymentDates)[j]> payDate)
                {
                    foundEndPeriod = true;
                    endMonitoringDate = (*monitorDates)[j-1];
                }
                j++;
            }
            if (!foundEndPeriod)
            {
                throw ModelException("Accumulator::getEndCurrentMonitoringPeriod related to ", " date " + today.toString() + 
                " failed to find the related end of period");
            }
        }
        else
        {
            throw ModelException("Accumulator::getEndCurrentMonitoringPeriod related to ", " date " + today.toString() + 
                " failed to recognize this date as a monitoring date");
        }
    
        return endMonitoringDate;
    }

    DateTime getStartCurrentMonitoringPeriod(const DateTime& today) const
    {
        DateTimeArraySP monitorDates = getDatesFromCFL(monitoringDatesAndStrikes);
        DateTime startMonitoringDate;
        DateTime payDate;

        bool foundMonitor = false;
        bool foundStartPeriod = false;
        int i = 0;
        int j = 0;

        while (!foundMonitor && i < monitorDates->size())
        {
            if ( today.equals((*monitorDates)[i]))
            {
                foundMonitor = true;
                payDate = (*paymentDates)[i];
            }
            i++;
        }

        if(foundMonitor)
        {
            j = i;
            while (!foundStartPeriod && j >= 0)
            {
                // case where the period considered is the first
                if( j==0)
                {
                    foundStartPeriod = true;
                    startMonitoringDate = (*monitorDates)[0];
                }

                if( (*paymentDates)[j]< payDate)
                {
                    foundStartPeriod = true;
                    startMonitoringDate = (*monitorDates)[j+1];
                }
                
                j--;
            }
            if (!foundStartPeriod)
            {
                throw ModelException("Accumulator::getStartCurrentMonitoringPeriod related to ", " date " + today.toString() + 
                " failed to find the related start of period");
            }
        }
        else
        {
            throw ModelException("Accumulator::getStartCurrentMonitoringPeriod related to ", " date " + today.toString() + 
                " failed to recognize this date as a monitoring date");
        }
    
        return startMonitoringDate;
    }

    // extract the first monitoring date corresponding to this settlement date
    DateTime getFirstMonitoringDate(const DateTime& settleDate) const
    {
        DateTimeArraySP monitorDates = getDatesFromCFL(monitoringDatesAndStrikes);
        DateTime monitoringDate;

        bool foundMonitor = false;
        int i = 0;

        while (!foundMonitor && i < monitorDates->size())
        {
            if ( settleDate.equals((*paymentDates)[i]))
            {
                foundMonitor = true;
                monitoringDate = (*monitorDates)[i];
            }
            i++;
        }

        if (!foundMonitor)
        {
            throw ModelException("Accumulator::getFirstMonitoringDate related to ", " date " + settleDate.toString() + 
            " failed to find the first monitoring date");
        }
    
        return monitoringDate;
    }

    DateTimeArraySP static getDatesFromCFL(const CashFlowArray& cfl)
    {
        DateTimeArraySP dta(new DateTimeArray(cfl.size()));
        for (int i=0; i<dta->size(); i++) { 
            (*dta)[i] = cfl[i].date;
        }
        return dta;
    }

private:
    friend class AccumulatorHelper;
    friend class AccumulatorFDProd;

protected:
    bool hasExtraParticipation() const
    {
        return (!!extraParticipation && extraParticipation->size() !=0 );
    }

    static void load(CClassSP& clazz);
    Accumulator();

    ScheduleSP              monitoringSchedule;
    CashFlowArray           monitoringDatesAndStrikes;
    DateTimeArraySP            monitoringDates;
    DateTimeArraySP            paymentDates;
    DoubleArraySP            histMonSamples;         // level of underlying on historical sample dates
    DoubleArray             participation;
    DoubleArray             barrierLevel;
    DoubleArray                ecoBarrierLevel;        // economic barrier used for reporting only
    DoubleArray             rebate;

    // extra participation available if spot reach the strike
    DoubleArraySP           extraParticipation;  

    bool                    isUpAndOut;
    bool                    isContinuouslyMonitored;
    bool                    isGlobalKO;
    bool                    includeDivs;
    int                        productAccumulated;
    mutable bool            isKO; // $unregistered
    mutable DateTime        knockedOutDate; // $unregistered
    mutable double            KoRebate; // $unregistered
    mutable double            KOlevel; // $unregistered
    mutable bool            stopSameDaySettle;

};
typedef smartPtr<Accumulator> AccumulatorSP;

// helpers
class AccumulatorHelper {
public:
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(Accumulator, clazz);
        SUPERCLASS(Generic1Factor);
        EMPTY_SHELL_METHOD(defaultAccumulator);
        // same as vanilla
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(ISensitiveStrikes);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(LegalTerms::Shift);
        FIELD(monitoringSchedule,        "monitoring Schedule");
        FIELD_MAKE_OPTIONAL(monitoringSchedule); // to be eventually removed
        FIELD(paymentDates, "payment dates: one for each monitoring date");
        FIELD_MAKE_OPTIONAL(paymentDates);
        FIELD(histMonSamples, "level of underlying on hist sample dates");
        FIELD_MAKE_OPTIONAL(histMonSamples);
        FIELD(participation, "Gearing at different strike dates.");
        FIELD(barrierLevel,        "Barrier Level for Knock-Out Condition");
        FIELD(ecoBarrierLevel,        "Reported Barrier");
        FIELD_MAKE_OPTIONAL(ecoBarrierLevel);
        FIELD(rebate,        "Rebate in case of Knock-Out");
        FIELD(isUpAndOut,        "Is product Up and Out ? If not, product is Down and Out");
        FIELD_MAKE_OPTIONAL(isUpAndOut);
        FIELD(isContinuouslyMonitored,   "Is Barrier Continuously monitored ?");
        FIELD_MAKE_OPTIONAL(isContinuouslyMonitored);
        FIELD(isGlobalKO,   "Does a Knock-Out affect all the product ?");
        FIELD_MAKE_OPTIONAL(isGlobalKO);
        FIELD(includeDivs,   "Are intermediary dividends included ?");
        FIELD_MAKE_OPTIONAL(includeDivs);
        FIELD(monitoringDatesAndStrikes, "monitoring CashFlowArray");
        FIELD_MAKE_OPTIONAL(monitoringDatesAndStrikes);
        FIELD(monitoringDates, "monitoring dates");
        FIELD_MAKE_OPTIONAL(monitoringDates);
        FIELD(productAccumulated, "What we accumulated can be a forward, a call or a put");
        FIELD_MAKE_OPTIONAL(productAccumulated);
        FIELD(stopSameDaySettle, "If true, same day settlements are not allowed");
        FIELD_MAKE_OPTIONAL(stopSameDaySettle);
        FIELD(extraParticipation, "Extra participation on a day if spot is above/below strike for downKO/upKO respectively");
        FIELD_MAKE_OPTIONAL(extraParticipation);
    }

    static IObject* defaultAccumulator(){
        return new Accumulator();
    }
};


CClassConstSP const Accumulator::TYPE = CClass::registerClassLoadMethod(
    "Accumulator", typeid(Accumulator), AccumulatorHelper::load);


void Accumulator::GetMarket(const IModel*         model, 
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
    if (premiumSettle.get())
    {
        premiumSettle->getMarket(model, market.get());
    }
}

bool Accumulator::avoidVegaMatrix(const IModel* model)
{
    return CTree1fLV::TYPE->isInstance(model);
}

// constructor: initialize optional ks
Accumulator::Accumulator(): Generic1Factor(TYPE), 
        isUpAndOut(true),
        isContinuouslyMonitored(false),
        isGlobalKO(false),
        includeDivs(false),
        productAccumulated(1),
        isKO(false),
        knockedOutDate(),
        stopSameDaySettle(false){}

void Accumulator::validatePop2Object(){
    static const string method = "Accumulator::validatePop2Object";
    static const string errMsg = "Monitoring Schedule must be either empty or contain the same dates and values as Monitoring Dates and Strikes";
   
    int num = monitoringDatesAndStrikes.size();
    if (num == 0)
    {
        if( !monitoringSchedule || monitoringSchedule->length()==0 )
            throw ModelException(method, "Monitoring schedule must be provided if there is no monitoringDatesAndStrikes.");

        const DateTimeArray& dts = monitoringSchedule->getDates();
        const DoubleArray& vals = monitoringSchedule->getValues();
        monitoringDatesAndStrikes.resize(dts.size());
        for (int i = 0; i < dts.size(); i++)
        {
            monitoringDatesAndStrikes[i].date = dts[i];
            monitoringDatesAndStrikes[i].amount = vals[i];
        }
        num = dts.size();
    }
    else if ( !!monitoringSchedule && monitoringSchedule->length() > 0 )
    {
        const DateTimeArray& dts = monitoringSchedule->getDates();
        const DoubleArray& vals = monitoringSchedule->getValues();

        // if the cashflow array monitoringDatesAndStrikes is used, it must match the monitoringSchedule
        if( dts.size() != num)
            throw ModelException(method, errMsg);
   
        // same size, need to check each element
        bool match = true;
        int i = 0;
        while (match && i < dts.size())
        {
            if (dts[i] != monitoringDatesAndStrikes[i].date || 
                vals[i] != monitoringDatesAndStrikes[i].amount)
            {
                match = false;
            }
            i++;
        }
        if (!match)
        {
            throw ModelException(method, errMsg);
        }
    }

    // check for array dimensions
    if( !paymentDates || paymentDates->size() != num ||
        participation.size() != num ||
        barrierLevel.size() != num ||
        rebate.size() != num )
        throw ModelException(method, "Paymentdate, participation, barrier and rebate must be same size as monitor dates.");

    // may have extra participation
    if( hasExtraParticipation() )
    {
        if( extraParticipation->size() != monitoringDatesAndStrikes.size() )
            throw ModelException(method, "If has extraParticipation, extraParticipation must be same size as Monitoring Schedule.");  
    }

    if( includeDivs && productAccumulated != 1 )
        throw ModelException(method, "IncludeDivs flag only allowed if accumulate fwd contract.");  

    // validation of the monitoringDates field
    if (!!monitoringDates)
    {
        if  (monitoringDates->size() != num)
            throw ModelException(method, "monitoringDates must be the same size as monitoringDatesAndStrikes");
        
        // same size, need to check each element
        for(int i = 0; i < num; i++)
        {
            if (monitoringDatesAndStrikes[i].date != (*monitoringDates)[i])
                throw ModelException(method, "monitoringDates must match the dates in monitoringDatesAndStrikes");
        }
    }
}

void Accumulator::Validate()
{
    static const string method = "Accumulator::Validate";
    // just check the things that aren't/cannot be checked in 
    // validatePop2Object
    if (!asset){
        throw ModelException(method, "Asset is null");
    }
    if (!discount){
        throw ModelException(method, "Discount YC is null");
    }

    // validate against struck until thoroughly tested and agreed methodology
    if (ccyTreatment == Asset::CCY_TREATMENT_STRUCK) {
        throw ModelException(method, "ccy struck type not supported.");
    }

    // validation against forward start
    if ( fwdStarting ) {
        throw ModelException(method, "Forward Starting feature not supported");
    }

    // validation against one coupon and forward start
    if ( fwdStarting && oneContract ) {
        throw ModelException(method, "Cannot be One Contract and Forward Starting");
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

    DateTimeArraySP monitorDates = getDatesFromCFL(monitoringDatesAndStrikes);

    DateTime start = fwdStarting? startDate : valueDate;
    int numPastSampleDates = monitorDates->size() - start.numFutureDates(*monitorDates);

    if (numPastSampleDates > 0)
    {
        // there are some historical dates so check we have sample values 
        if (!histMonSamples.get())
        {
            throw ModelException(method,
                                 "There are " + Format::toString(numPastSampleDates) +
                                 " historical monitoring dates but no historical"
                                 " sample values have been provided.");
        }
        
        if (histMonSamples->size() < numPastSampleDates)
        {
            throw ModelException(method,
                                 "There are " + Format::toString(numPastSampleDates) +
                                 " historical monitoring dates but only "+
                                 Format::toString(histMonSamples->size()) +
                                 " historical sample values have been provided.");
        }
        
        // check that all the historical samples are non-zero 
        int idx;
        for (idx = 0; idx < numPastSampleDates; idx++)
        {
            if (!Maths::isPositive((*histMonSamples)[idx]))
            {
                throw ModelException(method,
                                     "Historical sample " + Format::toString(idx+1)+
                                     " for monitoring date " + (*monitorDates)[idx].toString()+
                                     " must be greater than zero.");
            }
        }
    }

    if (numPastSampleDates ==0)
    {
        /* Create an empty histMonSamples list. This has the benefit
           of discarding a bogus list passed in when none are expected.
           This can cause problems when updating the past samples during theta twk*/
        histMonSamples = DoubleArraySP(new DoubleArray(0));
    }
    else if (numPastSampleDates < histMonSamples->size())
    {
        histMonSamples->resize(numPastSampleDates);
    }

    /* validate that each payment date is not before its associated 
     monitoring point */

    int idx;    

    for ( idx = 0; idx < monitorDates->size(); idx++)
    {
        if ((*monitorDates)[idx].isGreater((*paymentDates)[idx]))
        {
            throw ModelException(method,
                                 "Payment date " + Format::toString(idx+1) + " " + (*paymentDates)[idx].toString()+
                                 " is before it's corresponding monitoring date " +
                                 (*monitorDates)[idx].toString());
        }
    
        // ensure that the payment dates are chronological 
        if (idx > 0 &&
            (*paymentDates)[idx-1].isGreater((*paymentDates)[idx]))
        {
            throw ModelException(method,
                                 "Payment dates are not chronological. The payment "
                                 "date for monitoring date " + Format::toString(idx+1)+ " "+
                                 (*paymentDates)[idx].toString()+
                                 " is before the payment date " + (*paymentDates)[idx-1].toString());
        }
    }

    // validate that there is only one strike per monitoring period
    for (idx = 0; idx < monitorDates->size()-1; idx++)
    {
        if( (*paymentDates)[idx].equals((*paymentDates)[idx+1])
            && (!Maths::isZero(
                                monitoringDatesAndStrikes[idx].amount
                                - monitoringDatesAndStrikes[idx + 1].amount )))
        {
            throw ModelException(method,
                                 "monitoring dates " + Format::toString(idx) + " and "+
                                 Format::toString(idx+1)+ " must have same strike because they are in the same monitoring period ");
        }
        
    }

    // validate that there is a barrier level
    if(barrierLevel.size() != 0)
    {
        if(barrierLevel.size() != monitorDates->size())
        {
            throw ModelException(method,
                                 " the number of barrier levels must be equal to the number of monitoring dates ");
        }

    }
    else
    {
        throw ModelException(method,
                                 " one barrier level must be entered for every monitoring date ");
    }

    // validate default value of product accumulation
    if(productAccumulated != 0)
    {
        if( !((productAccumulated == 1)||
            (productAccumulated == 2)||
            (productAccumulated == 3) ) )
        {
            throw ModelException(method,
                                 " product accumulated must be 1 (meaning Forward), 2 (Call) or 3 (Put) ");
        }
    }
    else
    {
        throw ModelException(method,
                                 " type of product accumulated must be declared ");
    }

    // validate against same day settlements
    if(stopSameDaySettle)
    {
        for (idx = 0; idx < monitorDates->size()-1; idx++)
        {
            if( (*monitorDates)[idx].equals((*paymentDates)[idx],false))
            {
                throw ModelException(method,
                                     "monitoring date " + Format::toString(idx) + " not allowed to settle the same day ");
            }
        }
    }

    // validate that the economic barrier (the one reported) has the same size
    // validate that there is a barrier level
    if(ecoBarrierLevel.size() != 0)
    {
        if(ecoBarrierLevel.size() != barrierLevel.size())
        {
            throw ModelException(method,
                                 " reported economic barrier and real barrier must have the same size ");
        }
    }
    else // populate it in case it is blank
    {
        ecoBarrierLevel.resize(barrierLevel.size());
        int i =0;
        for(i=0;i<ecoBarrierLevel.size();i++)
        {
            ecoBarrierLevel[i] = barrierLevel[i];
        }
    }
}



/** returns the current value date */
DateTime Accumulator::getValueDate() const
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
CashFlowArraySP Accumulator::getKnownCashFlows() const
{
    DateTime settleDate,pastPayDate,monitorDate, pastDate = valueDate;
    DateTimeArraySP monitorDates = getDatesFromCFL(monitoringDatesAndStrikes);
    CashFlowArraySP cfl(new CashFlowArray(0));
    
    // return empty cfl if all monitor dates are in future
    if ( valueDate < monitorDates->front() ) return cfl;
    
    // computing the scaling factor
    double fwdStrt = fwdStarting?asset->fwdValue(startDate):0.0;
    double scalingFactor = InstrumentUtil::scalePremium(oneContract,fwdStarting,notional,fwdStrt,initialSpot);

    // in case product is Knocked Out
    if (isKO) 
    {       
        // adding KO rebate if needed
        pastDate = knockedOutDate; 
        pastPayDate = getPayDate(knockedOutDate);
        settleDate = instSettle->settles(pastPayDate, asset.get()); 
        cfl->push_back(CashFlow(settleDate, scalingFactor*KoRebate));
    }
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
    int iEnd = monitorDates->size()-1, i, j; 
    while( iEnd >= 0 && (*monitorDates)[iEnd] > pastDate ) 
        iEnd--;
    if( iEnd < 0 )    // sanity check, should never happen
        throw ModelException("Accumulator::getKnownCashFlows", "Internal error.");
    
    j = 0;
    for(i=0; i<=iEnd; i++)
    {
        bool isEndPeriod = (i == (monitorDates->size()-1) || (*paymentDates)[i] != (*paymentDates)[i+1] );
        
        // for cash settle, add periodic payment
        // if KO or end of monitor period, include current period accumulated contract.
        // otherwise, only include previous fully historical periods        
        if( !instSettle->isPhysical() &&
            ( isEndPeriod || ( iEnd==i && isKO ) ))
        {
            monitorDate = (*monitorDates)[i];
            settleDate = instSettle->settles((*paymentDates)[i], asset.get());
            double histContracts = getKoHistContracts(monitorDate,true);
            double strike = monitoringDatesAndStrikes[i].amount;
            double prodVal = histContracts*productValue((*histMonSamples)[i], strike);
            if( !Maths::isZero(prodVal) )
                cfl->push_back(CashFlow(settleDate, scalingFactor*prodVal));
        }
        
        // adding dividend payments as they become known
        if( includeDivs && exDivDates->size() > 0 )
        {
            DateTime monitorStart = (*monitorDates)[i];
            DateTime monitorEnd = isEndPeriod?monitorStart:(*monitorDates)[i+1];
            while( j < exDivDates->size() && 
                (*exDivDates)[j].getDate() >= monitorStart.getDate() &&
                (*exDivDates)[j].getDate() < monitorEnd.getDate() &&
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

PhysicalDeliveryArray Accumulator::getPhysicalDelivery() const
{
    DateTime pastDate = valueDate, pastPayDate, monitorDate, settleDate;
    DateTimeArraySP monitorDates = getDatesFromCFL(monitoringDatesAndStrikes);
    PhysicalDeliveryArray pda;
    
    // return empty cfl if all monitor dates are in future
    if ( valueDate < monitorDates->front() ) return pda;
    
    // computing the scaling factor
    double fwdStrt = fwdStarting?asset->fwdValue(startDate):0.0;
    double scalingFactor = InstrumentUtil::scalePremium(oneContract,fwdStarting,notional,fwdStrt,initialSpot);
    
    // find the latest monitor date, can be pastDate
    int iEnd = monitorDates->size()-1, i; 
    while( iEnd >= 0 && (*monitorDates)[iEnd] > pastDate ) iEnd--;
    
    // sanity check, should never happen
    if( iEnd < 0 )
        throw ModelException("Accumulator::getKnownCashFlows", "Internal error.");
    
    for(i=0; i<=iEnd; i++)
    {
        bool isEndPeriod = (i == (monitorDates->size()-1) || (*paymentDates)[i] != (*paymentDates)[i+1] );
        
        // for cash settle, add periodic payment
        // if KO or end of monitor period, include current period accumulated contract.
        // otherwise, only include previous fully historical periods        
        if( isEndPeriod || ( iEnd==i && isKO ) )
        {
            monitorDate = (*monitorDates)[i];
            pastPayDate = (*paymentDates)[i];
            settleDate = instSettle->settles(pastPayDate, asset.get());
            
            double histContracts = getKoHistContracts(monitorDate,true);
            double strike = monitoringDatesAndStrikes[i].amount;
            double spotAtMat = (*histMonSamples)[i];
            
            if( !Maths::isZero(histContracts) &&
                (productAccumulated == 1 ||
                (productAccumulated == 2 && spotAtMat > strike) ||
                (productAccumulated == 3 && spotAtMat < strike) ))
            {
                if( productAccumulated == 3 ) histContracts = - histContracts;
                PhysicalDelivery delivery(scalingFactor*histContracts,strike, monitorDate, settleDate);
                pda.push_back(delivery);
            }    
        }
        
    }
     
    return pda;
}

DateTimeArraySP Accumulator::getPaymentDates() const
{
    DateTimeArraySP payDates(new DateTimeArray(paymentDates->size()));
    
    for(int i=0; i<paymentDates->size(); i++)
        (*payDates)[i] = instSettle->settles((*paymentDates)[i], asset.get());

    return payDates;
}

/** price a dead instrument until settlement - exercised, expired, knocked out etc.
returns true if it is dead (and priced), false if it is not dead */
bool Accumulator::priceDeadInstrument(CControl* control, CResults* results) const
{
    static const string method = "Accumulator::priceDeadInstrument";

    try  
    {
        bool deadInstrument = false;
        DateTime settlementDate = (*paymentDates)[paymentDates->size() - 1];

        if(isKO)
        {   
            DateTime endDate = getPayDate(knockedOutDate);
            if (valueDate >= endDate)
            {
                results->storePrice(0.0, discount->getCcy());
                addOutputRequests(control, results, 0.0, 0.0);
                deadInstrument = true;
            }
        }

        if (valueDate >= settlementDate)
        {   // settled already
            results->storePrice(0.0, discount->getCcy());
            addOutputRequests(control, results, 0.0, 0.0);
            deadInstrument = true;
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
bool Accumulator::sensShift(Theta* shift)
{    
    DateTime newDate = shift->rollDate(valueDate);

    DateTime matDate = (*paymentDates)[paymentDates->size() - 1]; 

//    if (fwdStarting && newDate.isGreaterOrEqual(startDate) &&
//        startDate.isGreaterOrEqual(valueDate))
//    {
//        fwdStarting = false;
//        initialSpot = asset->getThetaSpotOnDate(shift, startDate);
//        monitoringSchedule->scale(initialSpot);
//    }

    /* also need to fill in historical levels for any samples that are in the
           future with respect to the value date but are in the past with respect
           to the new date */

    DateTimeArraySP monitorDates = getDatesFromCFL(monitoringDatesAndStrikes);

    if (!histMonSamples)
    {
        histMonSamples = DoubleArraySP( new DoubleArray(monitorDates->size(), 0.0) );
    }
    else
    {
        histMonSamples->resize(monitorDates->size());
    }


    DoubleArray wgts(monitorDates->size(), 1.0);
    SampleList monDatesValues(*monitorDates, *histMonSamples, wgts);
    monDatesValues.roll(asset.get(),valueDate,newDate,true); 

    monitorDates = DateTimeArraySP( new DateTimeArray(monDatesValues.getDates()) );
    histMonSamples = DoubleArraySP( new DoubleArray(monDatesValues.getValues()) );

    // roll today 
    valueDate = newDate;
    
    return true;
};

/** Satisfy LegalTerms::Shift interface */
bool Accumulator::sensShift(LegalTerms* shift) {
    // Set the barriers for pricing equal to the economic barriers
    // ecoBarrierLevel is the same size as barrierLevel (checked in Validate())
    for(int i=0; i<ecoBarrierLevel.size(); i++)
    {
        barrierLevel[i] = ecoBarrierLevel[i];
    }

    return true;
}

/** when to stop tweaking */
DateTime Accumulator::endDate(const Sensitivity* sensControl) const {
    DateTime matDate = monitoringDatesAndStrikes.back().date;
    DateTime instEnd  = instSettle->settles(matDate, asset.get());
    DateTime assetEnd = asset->settleDate(matDate);
    DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
    return end;
}

void Accumulator::addOutputRequests(Control* control,
                                 Results* results,
                                 const double& fairValue,
                                 const double& indVol) const
{
    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        DateTime matDate = (*paymentDates)[paymentDates->size() - 1];
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
class AccumulatorFDProd: public LatticeProdEDRIns
{
public:
    AccumulatorFDProd( const Accumulator * acc, FDModel * mdl ) :
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
                                                 (*inst->paymentDates)[inst->paymentDates->size() - 1]);

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
        static const string method = "AccumulatorFDProd::Init";
        try {
            int i;

            // all exe dates are critical dates
            DateTimeArray critDates;
            int num = inst->monitoringDatesAndStrikes.size();
            for (i=0; i<num; i++) {          
                critDates.push_back(inst->monitoringDatesAndStrikes[i].date);
                critDates.push_back((*(inst->paymentDates))[i]);
            }

            if( tree1f ) 
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
        static const string method = "AccumulatorFDProd::initProd";
        try {
            int i = 0;
            int j = 0;
            // timeline exists at this point
            int lastStep = model->getLastStep();

            initSlices( numPrices );
            initInsertNode();

            DateTime endDate(inst->valueDate.getDate()+364, inst->valueDate.getTime()); 

            stepIsMonitoringDate.resize(lastStep+1);
            stepIsKODate.resize(lastStep+1);
            stepIsPayDate.resize(lastStep+1);
            stepStrikeLevel.resize(lastStep+1);
            stepParticipation.resize(lastStep+1);
            stepBarrierLevel.resize(lastStep+1);
            stepRebateLevel.resize(lastStep+1);
            stepIsDividendDate.resize(lastStep+1);
            stepDividend.resize(lastStep+1);
            stepExtraParticipation.resize(lastStep+1);

            // past call dates are not taken into account

            DateTimeArraySP monitorDatesSP = inst->getDatesFromCFL(inst->monitoringDatesAndStrikes);
            DateTimeArray dts = *monitorDatesSP;
             int monitorIdx = inst->getIndexOnOrAfterDate(dts, model->getDate(0));

            HolidayConstSP hol = AssetUtil::getHoliday(inst->asset.get());
            
            for (i=0; i < lastStep; i++) 
            {
                const DateTime treeDate = model->getDate(i);

                stepIsMonitoringDate[i] = false;
                stepIsKODate[i] = false;
                stepStrikeLevel[i] = 0.0;
                stepParticipation[i] = 0.0;
                stepExtraParticipation[i] = 0.0;
                stepRebateLevel[i] = 0.0;
                stepIsDividendDate[i] = false;
                stepDividend[i] = 0.0;

                if (inst->isUpAndOut)
                {    
                    stepBarrierLevel[i] = 9999999999999.0;
                }
                else
                {
                    stepBarrierLevel[i] = 0.0;
                }
                
                if (monitorIdx < dts.size() && 
                    treeDate.equals( dts[monitorIdx],false ) )
                {
                    stepBarrierLevel[i] = inst->barrierLevel[monitorIdx];
                    stepRebateLevel[i] = inst->rebate[monitorIdx];
                    stepParticipation[i] = inst->participation[monitorIdx];
                    stepIsKODate[i] = true; /* discrete monitoring uses modified continous monitoring barrier level
                                                 but the algorithmic treatment is the same */
                    stepStrikeLevel[i] = inst->monitoringDatesAndStrikes[monitorIdx].amount;
                    if( inst->hasExtraParticipation() )
                        stepExtraParticipation[i] = (*inst->extraParticipation)[monitorIdx];

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
            while ( (payIdx < inst->paymentDates->size())
                        && ((*(inst->paymentDates))[payIdx] < model->getDate(0)) )
            {
                payIdx ++;
            }

            for (i=0; i<=lastStep; i++) 
            {
                const DateTime treeDate = model->getDate(i);

                stepIsPayDate[i] = false;

                if ( payIdx < inst->paymentDates->size() && 
                    treeDate.equals( (*(inst->paymentDates))[payIdx]) )
                {
                    stepIsPayDate[i] = true; 
                    stepStrikeLevel[i] = inst->monitoringDatesAndStrikes[payIdx].amount;
                    
                    //payIdx ++;
                    int k = payIdx;
                    
                    while ( (payIdx < inst->paymentDates->size())
                        && ((*(inst->paymentDates))[payIdx].equals( (*(inst->paymentDates))[k]) ))
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
                                                         (*inst->paymentDates)[inst->paymentDates->size() - 1]);

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
            DateTime       matDate = (*inst->paymentDates)[inst->paymentDates->size() - 1]; ;
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
                        double size = inst->ecoBarrierLevel.size();
                        DateTimeArraySP monitorDates = inst->getDatesFromCFL(inst->monitoringDatesAndStrikes);
                        
                        // use the barrierWindow routine 
                        DateTime upperDate = BarrierLevel::barrierWindow(inst->valueDate);


                        for (int i=0; i<size; i++)
                        {
                            barlevel = inst->ecoBarrierLevel[i];
                            barDate = (*monitorDates)[i];
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
        DateTime matDate   = inst->monitoringDatesAndStrikes.back().date;
        double   volStrike = inst->monitoringDatesAndStrikes.back().amount;

        CVolRequestConstSP volRequest(
            new LinearStrikeVolRequest(volStrike, getStartDate(), 
                                       matDate, inst->fwdStarting));
        return volRequest;
    }

    /** called before each update() */
    virtual void preCalc(int step)
    {
        static const string method("AccumulatorFDProd::preCalc");
        try {
            if( tree1f )
            {
                int idx = tree1f->getSliceIndex(step);
                const DateTime treeDate = model->getDate(step);

                if ( ( treeDate <= (*inst->paymentDates)[inst->paymentDates->size() - 1])) 
                                                    // && (stepIsKODate[step]) )
                                                  
                {
                    if (treeDate.equals((*inst->paymentDates)[inst->paymentDates->size() - 1]))
                    {
                        inst->KOlevel =  inst->barrierLevel.back();
                    }
                    else
                    {
                        inst->KOlevel =  stepBarrierLevel[step];
                    }

                    // making a correction in case barrier monitoring is discrete
                    if (!inst->isContinuouslyMonitored) {
                        vector<double> vol;
                        // adjust barrier if needed
                        tree1f->GetStepVol(step, vol, &inst->KOlevel, 0, 0); // get vol at barrier
                        Barrier::BarrierAdjustment(vol[0], inst->isUpAndOut, inst->KOlevel);         
                    }

                    if (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION) 
                    { //assuming Lower Barrier is index 0 and Upper uses 1 !!!
                        tree1f->SetInsertNode(idx, 0, inst->KOlevel, 0); // insert barrier level
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
        const TreeSlice & s = payoffIndex->getValue( step );
        int bot, top;
        s.getCalcRange( bot, top );

        const vector< TreeSliceSP > & price = slices;
        int pStart = 0, pEnd = price.size() - 1;

        if (type == FDProduct::BWD_T){
            prod_BWD_T( s,
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
            prod_BWD( s,
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
        static const string method("AccumulatorFDProd::prod_BWD_T");
        try {
            double * s = spot.getValues();
            const vector< double * > & p = getValues( price );

            int j;

            // get strike and maturity date from instrument
            DateTime matDate = inst->monitoringDatesAndStrikes.back().date;
            double matStrike  = inst->monitoringDatesAndStrikes.back().amount;
            const DateTime payDate = (*inst->paymentDates)[inst->paymentDates->size() - 1];    
            double matParticipation = inst->participation.back();

            nextPayDate = payDate;
            double barLevel = inst->KOlevel; //inst->barrierLevel.back();
            double rebate = inst->rebate.back();
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
                    if( (inst->isUpAndOut&& (s[j] >= barLevel)) || (!inst->isUpAndOut && (s[j] <= barLevel)) )
                    {        
                        p[0][j] = rebate;
                    }
                    else
                    {
                        p[0][j] = p[1][j]*matParticipation;
                        if( inst->hasExtraParticipation() )
                        {
                            if( (inst->isUpAndOut && s[j]<matStrike) || (!inst->isUpAndOut && s[j]>matStrike) )
                                p[0][j] += p[1][j] * inst->extraParticipation->back();
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
        static const string method("AccumulatorFDProd::prod_BWD");
        try {
            double * s = spot.getValues();
            const vector< double * > & p = getValues( price );

            int j;
            const DateTime treeDate = model->getDate(step);    

            double barLevel = inst->KOlevel;
            double rebate = stepRebateLevel[step];
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
                    if( (inst->isUpAndOut&& (s[j] >= barLevel)) || (!inst->isUpAndOut && (s[j] <= barLevel)) )
                    {            
                        p[0][j] = (inst->isGlobalKO?0.0:p[3][j]) + rebate*settlementPV;
                    }
                }
            }

            if(stepIsMonitoringDate[step] && step > 0)
            {
                for (j=bot; j<=top; j++)
                {
                    // adding payment of the forward contract maturing today and checking Knock-Out
                    if( (inst->isUpAndOut&& (s[j] >= barLevel)) || (!inst->isUpAndOut && (s[j] <= barLevel)) )
                    {            
                        p[0][j] =  (inst->isGlobalKO?0.0:p[3][j]) + rebate*settlementPV;
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
    const Accumulator* inst;
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

    double KoHistContracts;
    double KOrebate;
    DateTime nextPayDate;
    DateTimeArray        divDates;
    DoubleArray            divAmounts;

    double              knownCFValue;
};

/** create a fd payoff tree - new for all fd/tree state variable interface */
FDProductSP Accumulator::createProduct(FDModel* model) const
{
    return FDProductSP( new AccumulatorFDProd(this, model) );
}

/** returns all strikes on the vol surface to which 
    this instrument is sensitive */
DoubleArraySP Accumulator::getSensitiveStrikes(OutputNameConstSP outputName,
                                               const IModel*      model)
{
    DoubleArraySP sensStrikes = DoubleArraySP(new DoubleArray(0));

    if (avoidVegaMatrix(model)) {
        throw ModelException("Accumulator::getSensitiveStrikes", 
                             "VEGA_MATRIX is not valid for this instrument");
    }

    // get start date for vol interpolation
    DateTime imntStartDate = fwdStarting ? startDate:valueDate;
    
    // get last exercise date in exercise schedule
    DateTime  maturityDate =  monitoringDatesAndStrikes.back().date;

    double   strike       = monitoringDatesAndStrikes.back().amount;


    // create a vol request object to be passed on
    LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(strike, 
                                                                   imntStartDate, 
                                                                   maturityDate,
                                                                   fwdStarting));

    SensitiveStrikeDescriptor sensStrikeDesc;
    sensStrikeDesc.forwardOnly = false;

    asset->getSensitiveStrikes(volRequest.get(), outputName, 
                               sensStrikeDesc, sensStrikes);

    return sensStrikes;
}

extern bool AccumulatorLoad()
{
    return true && Accumulator::TYPE;
}

DRLIB_END_NAMESPACE
