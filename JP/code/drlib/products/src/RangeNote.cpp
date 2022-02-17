
//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RangeNote.cpp
//
//   Description : Range Note Instrument (generic 1 factor)
//
//   Author      : Stephen Hope
//
//   Date        : 5 Sep 2001
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/Vanilla.hpp"
#include "edginc/ErrorHandler.hpp"
#include "edginc/ProtEquity.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/SwapLegIntFace.hpp"
#include "edginc/Tree1fLV.hpp"
#include "edginc/VegaParallel.hpp"

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"

DRLIB_BEGIN_NAMESPACE
class RangeNote;

/** Abstract base class for the different kinds of RangeNote. 
    Note: This is a private class heiarchy and is used only to seperate
    the pricing functionality of each type of Range Note. */
class RNtype{
public:
    friend class RangeNote;
    friend class RangeNoteFDProd;

protected:
    RNtype(const RangeNote* rn);
    virtual ~RNtype(){}

    virtual double priceFutureSample() = 0;
    virtual double pricePastSample() = 0;
    virtual DoubleArray* getSensStrikes(OutputNameConstSP outputName) = 0;

    double singleRangeAsDigitalFactor(const Schedule& rangeSchedule,
                                      bool upsideRange);
    
    double singleRangeAsCallSpreadFactor(const Schedule& rangeSchedule,
                                         double spread,
                                         bool upsideRange,
                                         double& lowStrike,
                                         double& highStrike);

    void getSpreadStrikes(const string& spreadType,
                          double rangeLevel,
                          bool upsideRange,
                          double strikeSpread,
                          double& lowStrike,
                          double& highStrike);

    double priceSpread(const string& spreadType,
                       double rangeLevel,
                       double spread,
                       bool priceAsUpside,
                       double& lowStrike,
                       double& highStrike);
    
    const RangeNote* rn;  /* allows us to see the RangeNote elements
                             even though we are in a different class hierarchy */

private:
    RNtype(const RNtype& rhs);
    RNtype& operator=(const RNtype& rhs);
};

class RNsingleUpside: public RNtype{
public:
    friend class RangeNote;
private:
    RNsingleUpside(const RangeNote* rn):RNtype(rn){};
    virtual ~RNsingleUpside(){};
    RNsingleUpside(const RNsingleUpside& rhs);
    RNsingleUpside& operator=(const RNsingleUpside& rhs);

    virtual double priceFutureSample();
    virtual double pricePastSample();
    virtual DoubleArray* getSensStrikes(OutputNameConstSP outputName);
};

class RNsingleDownside: public RNtype{
public:
    friend class RangeNote;
private:
    RNsingleDownside(const RangeNote* rn):RNtype(rn){};
    virtual ~RNsingleDownside(){};
    RNsingleDownside(const RNsingleDownside& rhs);
    RNsingleDownside& operator=(const RNsingleDownside& rhs);

    virtual double priceFutureSample();
    virtual double pricePastSample();
    virtual DoubleArray* getSensStrikes(OutputNameConstSP outputName);
};

class RNdoubleInside: public RNtype{
public:
    friend class RangeNote;
private:
    RNdoubleInside(const RangeNote* rn):RNtype(rn){};
    virtual ~RNdoubleInside(){};
    RNdoubleInside(const RNdoubleInside& rhs);
    RNdoubleInside& operator=(const RNdoubleInside& rhs);

    virtual double priceFutureSample();
    virtual double pricePastSample();
    virtual DoubleArray* getSensStrikes(OutputNameConstSP outputName);
};

class RNdoubleOutside: public RNtype{
public:
    friend class RangeNote;
private:
    RNdoubleOutside(const RangeNote* rn):RNtype(rn){};
    virtual ~RNdoubleOutside(){};
    RNdoubleOutside(const RNdoubleOutside& rhs);
    RNdoubleOutside& operator=(const RNdoubleOutside& rhs);

    virtual double priceFutureSample();
    virtual double pricePastSample();
    virtual DoubleArray* getSensStrikes(OutputNameConstSP outputName);
};

class RangeNote: public Generic1Factor,
                 virtual public CClosedFormLN::IIntoProduct,
                 public FDModel::IIntoProduct,
                 virtual public LastSensDate{
public:
    static CClassConstSP const TYPE; 

    static const string RANGE_UPSIDE;
    static const string RANGE_DOWNSIDE;
    static const string RANGE_INSIDE;
    static const string RANGE_OUTSIDE;
    static const string SPREAD_NONE;
    static const string SPREAD_INCREASING;
    static const string SPREAD_DECREASING;
    static const string PAY_ON_MONITOR_DATE;
    static const string PAY_ON_SCHED_DATE;
    static const string PAY_ON_SINGLE_DATE;
    static const double MINIMUM_ALPHA;


    virtual void validatePop2Object(){
        static const string method = "RangeNote::validatePop2Object";
        
        // empty
    }

    virtual void Validate(){
        static const string method = "RangeNote::Validate";
        
        try
        {
            // validate the asset specific stuff
            validate();
            
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

            if (!paymentDates.get() && (paymentStyle == PAY_ON_SINGLE_DATE ||
                                    paymentStyle == PAY_ON_SCHED_DATE))
            {
                throw ModelException(method,
                                     "Payment dates must be provided for payment style of ("+
                                     paymentStyle + ")");
            } 
            
            /* ensure that we are pricing for one contract as fixed 
               notional does not make any sense for this model */
            /*if (!oneContract)
            {
                throw ModelException(method,
                                     "Range Notes cannot be priced for fixed notional");
            }*/
            
            ///////////////////////////////////////////////////////////////////////
            //             Validation of monitoring points                       //
            
            if (monitorDates->empty())
            {
                throw ModelException(method,
                                     "You have passed an empty monitoring date list, "
                                     "There isn't a deal to price.");
            }
            
            // check that the first monitoring date is not before the start date
            if ((*monitorDates)[0].isLess(startDate))
            {
                throw ModelException(method, 
                                     "The first monitoring date " + (*monitorDates)[0].toString()+
                                     " is before the start date " + startDate.toString());
            }
            
            DateTime matDate = (*monitorDates)[monitorDates->size()-1];
            /////////////////////////////////////////////////////////////
            //              payRangeType validation                    //
            
            if (!(CString::equalsIgnoreCase(payRangeType, RANGE_UPSIDE)) &&
                !(CString::equalsIgnoreCase(payRangeType, RANGE_DOWNSIDE)) &&
                !(CString::equalsIgnoreCase(payRangeType, RANGE_INSIDE)) &&
                !(CString::equalsIgnoreCase(payRangeType, RANGE_OUTSIDE)))
            {
                throw ModelException(method,
                                     "unknown payRangeType " + payRangeType +
                                     ". Valid types are UPSIDE, DOWNSIDE, "
                                     "INSIDE or OUTSIDE"); 
            }
            
            // check that highRangSchedule has been provided for a double barrier range note
            if ((payRangeType == RANGE_INSIDE) ||
                (payRangeType == RANGE_OUTSIDE))
                
            {
                if (!highRangeSchedule.get())
                {
                    throw ModelException(method,
                                         "highRangeSchedule must be provided for "
                                         "an INSIDE OR OUTSIDE RangeNote instrument");
                }
            }
        
            ////////////////////////////////////////////////////////////////////
            //               Validation of range levels                       //
            
            /* ensure that the rebate and range level schedules are defined from 
               the start date to the end date 
               Set the date from which the ranges should be defined to be the 
               first monitoring date 
               we can take the first element from the monitorDates array as we have 
               already validated that it is chronologically ordered */
            
            DateTime rangeStart = (*monitorDates)[0];
            if (!lowRangeSchedule->coversDateRange(rangeStart, matDate, true))
            {
                throw ModelException(method,
                                     "lowRangeSchedule is not defined from the first "
                                     "monitoring date " + rangeStart.toString() +
                                     " to the last monitoring date " + matDate.toString());
            }
            
            if (payRangeType == RANGE_INSIDE || payRangeType == RANGE_OUTSIDE)
            {
                if (!highRangeSchedule->coversDateRange(rangeStart, matDate, true))
                {
                    throw ModelException(method,
                                         "highRangeSchedule is not defined from the first "
                                         "monitoring date " + rangeStart.toString() +
                                         " to the last monitoring date " + matDate.toString()); 
                }
                
                /* ensure that the high barrier is above the low barrier at every 
                   monitoring point */
                double lowRange, highRange;
                for (int idx = 0; idx < monitorDates->size(); idx++)
                {
                    lowRange = lowRangeSchedule->interpolate((*monitorDates)[idx]);
                    highRange = highRangeSchedule->interpolate((*monitorDates)[idx]);
                    
                    if ((lowRange > highRange) || Maths::equals(lowRange, highRange))
                    {
                        throw ModelException(method,
                                         "lowRange " + Format::toString(lowRange) +
                                             " is not below highRange " + Format::toString(highRange)+
                                             " on monitoring date " + (*monitorDates)[idx].toString());
                    }
                    
                    if (Maths::isNegative(lowRange))
                    {
                        throw ModelException(method,
                                             "lowRange " + Format::toString(lowRange) +
                                             " is negative on monitoring date " + (*monitorDates)[idx].toString());
                    }
                }
            }
            else
            {
                // ensure that the low range is always positive 
                if (!lowRangeSchedule->isNonNegative())
                {
                    throw ModelException(method,
                                         "lowRangeSchedule contains negative points.");
                }
            }
            
            /////////////////////////////////////////////////////////////////////////////////////
            //                    Validation of rebate schedule                                //
            
            if (!rebateSchedule->coversDateRange(rangeStart, matDate, true))
            {
                throw ModelException(method,
                                     "rebateSchedule is not defined from the first "
                                     "monitoring date " + rangeStart.toString() +
                                     " to the last monitoring date " + matDate.toString()); 
            }

            ////////////////////////////////////////////////////////////////////////
            //          Validation of payment style and payment dates             //
            
            if (CString::equalsIgnoreCase(paymentStyle, PAY_ON_SINGLE_DATE))
            {
                // check we have one date in the payment date list 
                if (paymentDates->size() != 1)
                {
                    throw ModelException(method,
                                         "PAY_ON_SINGLE_DATE requested but " + Format::toString(paymentDates->size())+
                                         " pay dates specified");
                }
                
                // check that the payment date is not before the last monitor date 
                if (matDate.isGreater((*paymentDates)[0]))
                {
                    throw ModelException(method,
                                         "The single payment date " + (*paymentDates)[0].toString()+
                                         " is before the last monitoring date "+ matDate.toString());
                }
            }
            else if (CString::equalsIgnoreCase(paymentStyle, PAY_ON_SCHED_DATE))
            {
                if (paymentDates->size() != monitorDates->size())
                {
                    throw ModelException(method,
                                         "PAY_ON_SCHED_DATES therefore " + Format::toString(monitorDates->size())+
                                         " required but " + Format::toString(paymentDates->size())+
                                         " specified.");
                }
            
                /* validate that each payment date is not before its associated 
                   monitoring point */
                for (int idx = 0; idx < monitorDates->size(); idx++)
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
            }
            //ensure that the payment range type is valid 
            else if (CString::equalsIgnoreCase(paymentStyle, PAY_ON_MONITOR_DATE))
            {
                /* In this situation the instrument settlement information is used.
                   Fail if instSettle is not an instance of CashSettlePeriod */
                if (!CashSettlePeriod::TYPE->isInstance(instSettle.get()))
                {
                    throw ModelException(method,
                                         "A paymentStyle of PAY_ON_MONITOR_DATE can only " 
                                         " be combined with an instrument settlement that is"
                                         " a rolling cash settlement (CashSettlePeriod). ");
                }
        }
            else
            {
                throw ModelException(method,
                                     "Unsupported payment style ("+ paymentStyle+
                                     "). Valid payment styles are MONITOR, SCHEDULE OR SINGLE.");  
            }    
            ///////////////////////////////////////////////////////////////////////
            //                    Validate spreadType                            //
            
            if (CString::equalsIgnoreCase(spreadType, SPREAD_INCREASING) ||
                CString::equalsIgnoreCase(spreadType, SPREAD_DECREASING))
            {
                if (Maths::isNegative(lowStrikeSpread) ||
                    lowStrikeSpread < MINIMUM_ALPHA)
                {
                    throw ModelException(method,
                                         "lowStrikeSpread " + Format::toString(lowStrikeSpread)+
                                         " must be greater than DR MINIMUM SPREAD ("+
                                         Format::toString(MINIMUM_ALPHA)+
                                         ")");
                }
            
                /* only check this if we are a RANGE_INSIDE OR RANGE_OUTSIDE */
                if ((payRangeType == RANGE_INSIDE || payRangeType == RANGE_OUTSIDE) &&
                    (Maths::isNegative(highStrikeSpread) ||
                     highStrikeSpread < MINIMUM_ALPHA))
                {
                    throw ModelException(method,
                                         "highStrikeSpread " + Format::toString(highStrikeSpread)+
                                         " must be greater than DR MINIMUM SPREAD ("+
                                         Format::toString(MINIMUM_ALPHA)+
                                         ")");
                }
            }
            else if (!(CString::equalsIgnoreCase(spreadType, SPREAD_NONE)))
            {
                throw ModelException(method,
                                     "Unsupported spread type (" + spreadType+
                                     "). Valid spread types are NONE, INCREASING OR DECREASING.");
            }
        
            if (instSettle->isPhysical())
            {
                throw ModelException(method,
                                     "Physical settlement is not allowed for Range Notes");
            }
        
            // ensure that the rebate schedule is always positive 
            if (!rebateSchedule->isNonNegative())
            {
                throw ModelException(method,
                                     "rebateSchedule contains negative points.");
            }
                        
        }
        catch (exception& e)
        {
            if (rnType) { delete rnType; }
            throw ModelException(e, method);
        }
    }

/** Invoked when Class is 'loaded' */
static void load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(RangeNote, clazz);
    SUPERCLASS(Generic1Factor);
    IMPLEMENTS(CClosedFormLN::IIntoProduct);
    IMPLEMENTS(FDModel::IIntoProduct);
    IMPLEMENTS(LastSensDate);   
    EMPTY_SHELL_METHOD(defaultRangeNote);
    FIELD(payRangeType, "payment range type");
    FIELD(spreadType, "spread type");
    FIELD(paymentStyle, "When payments are made");
    FIELD(lowRangeSchedule, "low barrier schedule");
    FIELD(highRangeSchedule, "high barrier schedule");
    FIELD_MAKE_OPTIONAL(highRangeSchedule);
    FIELD(rebateSchedule, "rebate schedule");
    FIELD(monitorDates, "monitoring dates");
    FIELD(paymentDates, "one or more payment dates");
    FIELD_MAKE_OPTIONAL(paymentDates);
    FIELD(histMonSamples, "level of u/l on hist sample dates");
    FIELD_MAKE_OPTIONAL(histMonSamples);
    FIELD(priceHistAsDigital, "defines historical sample pricing ");
    FIELD(lowStrikeSpread, "low range strike separation");
    FIELD_MAKE_OPTIONAL(lowStrikeSpread);
    FIELD(highStrikeSpread, "high range strike separation");
    FIELD_MAKE_OPTIONAL(highStrikeSpread);
    FIELD(sampleDate, "monitor date");
    FIELD(rebateOnSampleDate, "rebate level on monitor date");
    FIELD(histLevelOnSampDate, "historic level on monitor date");
    FIELD(sampleSettlement, "instrument settlement determined by payment style");
    // extra optional fields for callable case
    FIELD(callSchedule, "Call Schedule");
    FIELD_MAKE_OPTIONAL(callSchedule);
    FIELD(barrierSchedule, "Barrier Schedule for KO");
    FIELD_MAKE_OPTIONAL(barrierSchedule);
    FIELD(optionParticipation, "default(0.0), participation weight for option at maturity");
    FIELD_MAKE_OPTIONAL(optionParticipation);
    FIELD(isCall, " Is final option a call option?");
    FIELD_MAKE_OPTIONAL(isCall);
    FIELD(finalStrike, "Final strike value");
    FIELD_MAKE_OPTIONAL(finalStrike);
    FIELD(floater, "libor funding leg");
    FIELD_MAKE_OPTIONAL(floater);
    FIELD(isUpAndOutKO, "Is KO feature Up and Out ?");
    FIELD_MAKE_OPTIONAL(isUpAndOutKO);
    FIELD(isContMonitoringKO, "Is KO continuously monitored ?");
    FIELD_MAKE_OPTIONAL(isContMonitoringKO);
    FIELD(isCalled, "Is call provision ?");
    FIELD_MAKE_OPTIONAL(isCalled);
    FIELD(calledDate, "Date at which call provision is activated");
    FIELD_MAKE_OPTIONAL(calledDate);
    FIELD(includeNotional, " Is notional paid at maturity ?");
    FIELD_MAKE_OPTIONAL(includeNotional);
    FIELD(KOrebate, " Amount paid in case of Knock-out");
    FIELD_MAKE_OPTIONAL(KOrebate);

    

    // hide from dd interface
    FIELD_MAKE_TRANSIENT(sampleDate); 
    FIELD_MAKE_TRANSIENT(rebateOnSampleDate);
    FIELD_MAKE_TRANSIENT(histLevelOnSampDate);
    FIELD_MAKE_TRANSIENT(sampleSettlement); 
}

static IObject* defaultRangeNote(){
    return new RangeNote();
}
   
private:
    friend class RangeNoteClosedForm;
    friend class RangeNoteFDProd;

    friend class RNtype;     // base class for the different kinds of range notes
    friend class RNsingleUpside;
    friend class RNsingleDownside;
    friend class RNdoubleInside;
    friend class RNdoubleOutside;


    RangeNote():Generic1Factor(TYPE),lowStrikeSpread(0.0), highStrikeSpread(0.0),
                callSchedule(ScheduleSP(   )),
                barrierSchedule(ScheduleSP(   )),optionParticipation(0.0),
                isCall(false),finalStrike(0.0),
                floater(LiborLegSP(   )),isUpAndOutKO(true), isContMonitoringKO(false),
                rebateOnSampleDate(0.0),histLevelOnSampDate(0.0),
                rnType(0),isCalled(false), calledDate(0,0), includeNotional(false),KOrebate(0.0){}

    RangeNote(const RangeNote& rhs);
    RangeNote& operator=(const RangeNote& rhs);

     /* Instantiate the RNtype field with whatever type of RangeNote 'this' type is. 
        Pricing can then be performed as RNtype->price() allowing a clear distinction 
        between which type we are dealing with */
    void instantiateRNtype()const{
        static const string method = "RangeNote::instantiateRNtype";

        if (payRangeType == RANGE_UPSIDE)
        {
            rnType = new RNsingleUpside(this);
        }
        else if (payRangeType == RANGE_DOWNSIDE)
        {
            rnType = new RNsingleDownside(this); 
        }
        else if (payRangeType == RANGE_INSIDE)
        {
            rnType = new RNdoubleInside(this); 
        }
        else if (payRangeType == RANGE_OUTSIDE)
        {
            rnType = new RNdoubleOutside(this); 
        }
        else
        {
            // cant really happen
            throw ModelException(method,
                                 "Unknown payRangeType " + payRangeType +
                                 ". Valid types are UPSIDE, DOWNSIDE, "
                                 "INSIDE or OUTSIDE");
        }
    }

    /** Indicates whether VEGA_MATRIX is sensible for this instrument.
        VEGA_MATRIX is only implemented for FLAT barrier schedules */
    bool avoidVegaMatrix(const IModel* model){

        static const string method = "RangeNote::avoidVegaMatrix";

        bool avoidVegaMatrix = false;
        bool flatHighSchedule = true;
        
        try
        {
            // for lv case, always avoid
            if (CTree1fLV::TYPE->isInstance(model)) {
                avoidVegaMatrix = true;
            } else {                
                bool flatLowSchedule = lowRangeSchedule->isFlat();
            
                // check if its a double barrier
                if (payRangeType == RangeNote::RANGE_INSIDE || payRangeType == RangeNote::RANGE_OUTSIDE)
                {
                    flatHighSchedule = highRangeSchedule->isFlat();
                }
            
                if (!flatLowSchedule || !flatHighSchedule)
                {
                    avoidVegaMatrix = true;
                    ErrorHandler::writeMsg("VEGA_MATRIX only supported for FLAT barrier schedules. "
                                           "VEGA_POINTWISE will be returned instead.");
                }
            }

            return avoidVegaMatrix;
        }
        catch (exception& e)
        {
            if (rnType) { delete rnType; }
            throw ModelException(e, method);
        }
    }

    /** when to stop tweaking */
    DateTime endDate(const Sensitivity* sensControl) const{
        
        /* get the maturity date of the range note - this is equal to the 
           last monitoring date */
        DateTime lastSampleDate = (*monitorDates)[monitorDates->size()-1];

        // get last sens date of underlying 
        DateTime lastUndDate = asset->settleDate(lastSampleDate);

        // also work out what the last payment date would be
        DateTime payDate;
        if (paymentStyle == PAY_ON_SINGLE_DATE)
        {
            payDate = (*paymentDates)[0];
        }
        else if (paymentStyle == PAY_ON_MONITOR_DATE)
        {
            /* settlement delay is applied for this payment type so determine 
               when the trade will settle */
            payDate = instSettle->settles(lastSampleDate,
                                          asset.get());
        }
        else if (paymentStyle == PAY_ON_SCHED_DATE)
        {
            payDate = (*paymentDates)[paymentDates->size()-1];
        }

        DateTime temp = lastSampleDate.isGreater(lastUndDate)? lastSampleDate : lastUndDate;
        DateTime lastSensDate = payDate.isGreater(temp)? payDate : temp;

        return lastSensDate;
    }

    /** Returns all strikes this particular type of RangeNote is sensitive to  */
    DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                      const IModel*       model) {
        static const string method = "RangeNote::getSensitiveStrikes";

        try
        {
            instantiateRNtype();

            DoubleArraySP sensStrikes(new DoubleArray(0));
            SensitiveStrikeDescriptor sensStrikeDesc;

            sensStrikeDesc.forwardOnly = false;
            
            DateTime instStartDate = fwdStarting? startDate : valueDate;
            DateTime maturityDate = (*monitorDates)[monitorDates->size()-1];

            // now get the sensitive strikes specific to the RangeNote
            DoubleArraySP strikeList(rnType->getSensStrikes(outputName));

            // for each strike create a vol request to fire at the asset
            for (int i = 0; i < strikeList->size(); i++)
            {
                LinearStrikeVolRequest volRequest((*strikeList)[i],
                                                  instStartDate,
                                                  maturityDate,
                                                  fwdStarting);
                
                // Generic does the rest
                getSensStrikes(outputName,
                              &volRequest,
                               sensStrikeDesc,
                               sensStrikes);
            }

            delete rnType;

            return sensStrikes;
        }
        catch (exception& e)
        {
            if (rnType) { delete rnType; }      
            throw ModelException(e, method);
        }
    }

    void addOutputRequests(Control* control, 
                           Results* results,
                           double fairValue,
                           double maxPayoff,
                           double pricePCTpayoff)const{
        
        if (control && control->isPricing())
        {
            DateTime matDate = (*monitorDates)[monitorDates->size()-1];

            // Generic1Factor deals with DELAY_PRICE and FWD_AT_MAT if requested
            addRequests(control, 
                        results,
                        fairValue,
                        matDate);

            // MAX_PAYOFF
            InstrumentUtil::recordMaxPayoff(control,
                                            results,
                                            maxPayoff);

            
            // PRICE_PCT_PAYOFF
            InstrumentUtil::recordPricePCTpayoff(control,
                                                 results,
                                                 pricePCTpayoff);

            
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
        }
    }

    /** Call pricing function relevant to particular range note type */
    void price(Control* control, CResults* results)const{
        static const string method = "RangeNote::price";
        
        try
        {
            // determine whether any historical samples exist 
            instantiateRNtype();

            int numFutureDates = valueDate.numFutureDates(*monitorDates);
            
            /* set up anything that only has to be done once per pricing call 
               rather than on every monitoring date */
            DateTime paymentDate;
            
            if (paymentStyle == PAY_ON_SINGLE_DATE)
            {       
                /* set up the payment date variable */
                paymentDate = (*paymentDates)[0];
            }
            
            if (numFutureDates > 0)
            {
                /* set up the settlement object with the correct payment date 
                   if the payment date is the same for every sample */
                if (paymentStyle == PAY_ON_SINGLE_DATE)
                {
                    /* all the sample points on the range note will have the same 
                       payment date so we can create the settlment object here */
                    sampleSettlement = InstrumentSettlementSP(new CashSettleDate((*paymentDates)[0]));
                }
                else if (paymentStyle == PAY_ON_MONITOR_DATE)
                {
                    // The settlement information in the Generic1Factor is valid 
                    sampleSettlement = instSettle;
                }
            }
            
            ///////////////////////// Main Pricing Loop /////////////////////////////////////
            int idx;
            double price = 0.0;
            double samplePrice = 0.0;
            double sampleMaxRebate = 0.0;
            double maxPayoff = 0.0;
            double pricePCTpayoff = 0.0;
            for (idx = 0; idx < monitorDates->size(); idx++)
            {
                //re-initialise the price for this sample 
                samplePrice = 0.0;
                
                // set the sample date
                sampleDate = (*monitorDates)[idx];
                
                // set the historic level on this sample date if any
                if (numFutureDates < monitorDates->size() &&
                    idx < (monitorDates->size()-numFutureDates))
                {
                    histLevelOnSampDate = (*histMonSamples)[idx];
                }
                
                // interpolate the rebate level(s) at this monitor date
                rebateOnSampleDate = rebateSchedule->interpolate(sampleDate);
                
                /* get the payment date 
                   note the payment date for pay on single date has already been set */
                if (paymentStyle == PAY_ON_MONITOR_DATE)
                {
                    // we apply settlement delay in this case only
                    paymentDate = instSettle->settles(sampleDate,
                                                      asset.get());
                }
                else if (paymentStyle == PAY_ON_SCHED_DATE)
                {
                    paymentDate = (*paymentDates)[idx];
                }
                // else PAY_ON_SINGLE_DATE which has already been set
                
                
                
                /* one of the outputs required is the maximum rebate that could
                   be received over the life of the range note. That is the rebate
                   that would be obtained if the sampled level was completely within the 
                   range on every sample date.
                   For historic monitoring dates, the maximum rebate rather than 
                   the rebate indicated by the historic level is captured.
                   The variable sampleMaxRebate holds the level of the maximum
                   rebate that could be received on a given sample date. */
        
                // determine whether the monitoring date is historic 
                if (valueDate.isGreaterOrEqual(sampleDate))
                {
                    // sample is historic - see if we have been paid any monies due 
                    if (valueDate.isGreaterOrEqual(paymentDate))
                    {
                        // payment is also historic so this sample has no value
                        samplePrice = 0.0;
                        sampleMaxRebate = 0.0;
                    }
                    else
                    {
                        // rebate payment has been fixed but not paid 
                        
                        /* get the factor, between zero and one, that indicates 
                           the proportion of the rebate that should be paid given 
                           the relative positions of the historical sample and 
                           range levels */

                        double rebateFactor = rnType->pricePastSample();
                        
                        // The maximum rebate that could have been paid is
                        sampleMaxRebate = rebateOnSampleDate;
                        
                        // no point continuing unless the rebate factor is non-zero
                        if (Maths::isZero(rebateFactor))
                        {
                            samplePrice = 0.0;
                        }
                        else
                        {
                            // discount the rebate from the payment date to the value date
                            double discFactor = discount->pv(valueDate, paymentDate);
                            
                            samplePrice = rebateOnSampleDate * rebateFactor * discFactor;
                        }
                    }
                }
                else
                {
                    /* Payment must be in the future change the
                       settlement object to have the correct rebate 
                       payment date if necessary */
                    sampleSettlement = InstrumentSettlementSP(new CashSettleDate(paymentDate));
                    
                    // price the future sample point
                    
                    samplePrice = rnType->priceFutureSample();
                    
                    // The maximum rebate that could have been paid is
                    sampleMaxRebate = rebateOnSampleDate;
                }
                
                // add the price of this sample to the total price
                price += samplePrice;
                
                maxPayoff += sampleMaxRebate;
                
            }
            
            if (Maths::isZero(price))
            {
                pricePCTpayoff = 0.0;  
            }
            else
            {
                pricePCTpayoff = price/maxPayoff;  
            }
            
            // record the total price in the output 
            
            results->storePrice(price, discount->getCcy());
            
            // add any requested outputs such as FWD_AT_MAT, DELAY_PRICE etc
            addOutputRequests(control,
                              results,
                              price,
                              maxPayoff,
                              pricePCTpayoff);
            
            delete rnType;
        }
        catch (exception& e)
        {
            if (rnType) { delete rnType; }
            throw ModelException(&e, method);
        }
    }
    
    /** Override clone method to copy our extra data over */
    IObject* clone() const{
        // first clone all the registered fields
        IObject*  copy = CObject::clone();
        RangeNote* rn = dynamic_cast<RangeNote*>(copy);
        if (!rn){
            throw ModelException("RangeNote::clone"); // shouldn't happen
        }
        
        rn->rnType = rnType; 
        return copy;
    }

    bool sensShift(Theta* shift){
        // we completely override the Generic1Factor theta shift method
        // and do everything ourselves (cf Extendible.cpp) 
        DateTime newValueDate = shift->rollDate(valueDate);

        double spot = asset->getThetaSpotOnDate(shift,
                                                newValueDate);
        
        /*If fwd start date falls between value date and theta date 
          must multiply the low/high range levels and low/high strike spreads, 
          which are expressed as percentages by the spot to convert to
          absolute values, as will no longer be forward starting.
          Also have to set initial spot to same level. */
        if (fwdStarting &&
            startDate.isGreaterOrEqual(valueDate) && 
            newValueDate.isGreaterOrEqual(startDate))
        {
            // now scale the the range levels and strike spreads 
            lowRangeSchedule->scale(spot);
            lowStrikeSpread *= spot;

            // see if we have a double barrier, if so scale this too
            if (payRangeType == RANGE_INSIDE || payRangeType == RANGE_OUTSIDE)
            {
               highRangeSchedule->scale(spot); 
               highStrikeSpread *= spot;

               if(!!barrierSchedule){
                   barrierSchedule->scale(spot); 
               }
               if(!!callSchedule){
                   callSchedule->scale(spot); 
               }

               finalStrike *= spot;
            }
            // no longer forward starting !
            fwdStarting = false;
        }

        /* also need to fill in historical levels for any samples that are in the
           future with respect to the value date but are in the past with respect
           to the new date */

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
        monDatesValues.roll(asset.get(),valueDate,newValueDate,true); 

        monitorDates = DateTimeArraySP( new DateTimeArray(monDatesValues.getDates()) );
        histMonSamples = DoubleArraySP( new DoubleArray(monDatesValues.getValues()) );


        // taking care of libor leg if it exists
        if(!!floater)
        {
            floater->setFixingforThetaShift( valueDate, 
                                            discount.get(),
                                                newValueDate);
        }
        
        valueDate = newValueDate;
        return true;
    }


    /** Implementation of ClosedFormLN::IntoProduct interface */
    CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;

    /** evaluate range note coupon*/
    double couponToday(const DateTime& today, const double spotToday, bool disc) const; // disc: discount to value date or not
    
    /** find out cash flow amount due at given date */
    double lookUpAmount(const DateTime& today, CashFlowArrayConstSP cfl) const;

    /** find out cash flows due strictly after evalDate */
    double lookUpRemainingAmounts(const DateTime& evalDate, CashFlowArrayConstSP cfl) const;

    /** find out date of end of accrual period corresponding to payment date payDate */
    DateTime lookUpEndAccrualPeriod(const DateTime& payDate) const;

    /** find pay date corresponding to the next monitoring date */
    DateTime getNextPaymentDate(const DateTime& today) const;
    
    /** find payDate of a monitoring date*/
    const DateTime& getPayDate(const DateTime& today) const;

    /** find payDate of a KO monitoring date*/
    const DateTime& getKOPayDate(const DateTime& today) const;

    /* last monitoring date */
    DateTime getLastDate() const
    {
        return (*monitorDates)[monitorDates->size() - 1];
    }

    /** pricing dead instrument */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;
    
    /* to know if product is already KO */
    bool isKnockedOut(const DateTime& today) const;

    /* pv of historical payments prior or equal to today evaluated at evalDate */
    double getRemainingHistPayments(const DateTime& evalDate, const DateTime& today) const;

    /* pv of historical payments Strictly prior to today evaluated at evalDate */
    double getRemainingStrictHistPayments(const DateTime& evalDate, const DateTime& today) const;

    /* recorded historical price */
    double getHistPrice(const DateTime& histDate) const;

    /* check if monitoring date is at the end of a monitoring period */
    bool isEndMonitoringPeriod( const DateTime& testDate) const;

    /* check if testdate is a call payment date */
    bool isCallPayDate (const DateTime& testDate) const;

    /* gets the end of current monitoring period */
    const DateTime& getEndPeriodMonitoringDate (const DateTime& today) const;

    
    /* to know certain cash flows */
    CashFlowArraySP getKnownCashFlows() const;

    /* to know certain payment dates */
    DateTimeArraySP getPaymentDates() const;


private:
    // registered elements

    string               payRangeType;           // inside, outside, upside or downside
    string               spreadType;             // increasing, decreasing or none
    string               paymentStyle;           // payOnMonitorDate, payOnPaySchedDate or payOnSingleDate
    ScheduleSP           lowRangeSchedule;       // defines one level in the range
    ScheduleSP           highRangeSchedule;      // defines second level in the range
    ScheduleSP           rebateSchedule;         // rebate payment if u/l in range on monitor date
    DateTimeArraySP      monitorDates;           // payment is triggered if u/l in range on these dates
    DateTimeArraySP      paymentDates;           // one or more payment dates
    DoubleArraySP        histMonSamples;         // level of u/l on hist sample dates
    bool                 priceHistAsDigital;     // defines historical sample pricing 
    double               lowStrikeSpread;        // low range strike separation for (increasing or decreasing)
    double               highStrikeSpread;       // high range strike separation for (increasing or decreasing)

    // extra optional fields for callable version
    
    // extra fields
    ScheduleSP callSchedule;  //call Schedule with strike level.
    ScheduleSP barrierSchedule; // for ko
    double optionParticipation;
    bool isCall;
    double finalStrike;
    LiborLegSP    floater;
    bool isUpAndOutKO;
    bool isContMonitoringKO;
    
    

    // --- debug ---
    static DoubleArraySP d_prices;
    
    // transients
    mutable DateTime               sampleDate;
    mutable double                 rebateOnSampleDate;
    mutable double                 histLevelOnSampDate;    
    mutable InstrumentSettlementSP sampleSettlement;     // depends on payment style
    mutable RNtype*                rnType; // $unregistered
    mutable DateTime               KnockedOutDate; // useful for priceDeadInstrument $unregistered
    mutable double                   KOlevel; // $unregistered
    mutable double                   fwdAtStart; // $unregistered
    mutable bool                   isCalled; // useful for priceDeadInstrument
    mutable DateTime               calledDate;
    mutable bool                   includeNotional;
    mutable double                   KOrebate;
    mutable bool                   updateLiborValue; // used in PayoffBeforeMat() $unregistered

};

const string RangeNote::RANGE_UPSIDE = "UPSIDE";
const string RangeNote::RANGE_DOWNSIDE = "DOWNSIDE";
const string RangeNote::RANGE_INSIDE = "INSIDE";
const string RangeNote::RANGE_OUTSIDE = "OUTSIDE";
const string RangeNote::SPREAD_NONE = "NONE";
const string RangeNote::SPREAD_INCREASING = "INCREASING";
const string RangeNote::SPREAD_DECREASING = "DECREASING";
const string RangeNote::PAY_ON_MONITOR_DATE = "MONITOR";
const string RangeNote::PAY_ON_SCHED_DATE = "SCHEDULE";
const string RangeNote::PAY_ON_SINGLE_DATE = "SINGLE";
const double RangeNote::MINIMUM_ALPHA = 0.001;

CClassConstSP const RangeNote::TYPE = CClass::registerClassLoadMethod(
    "RangeNote", typeid(RangeNote), RangeNote::load);

double RangeNote::couponToday(const DateTime& today, const double spotToday, bool disc =true) const
{
    
    DateTime paymentDate;
    // set the sample date
    sampleDate = today;
    
    // set the historic level on this sample date if any
    histLevelOnSampDate = spotToday;

    // interpolate the rebate level(s) at this monitor date
    rebateOnSampleDate = rebateSchedule->interpolate(sampleDate);

    /* get the payment date 
    note the payment date for pay on single date has already been set */
    if (paymentStyle == PAY_ON_MONITOR_DATE)
    {
        // we apply settlement delay in this case only
        paymentDate = instSettle->settles(sampleDate,
            asset.get());
    }
    else if (paymentStyle == PAY_ON_SCHED_DATE)
    {
        paymentDate = getPayDate(today);
    }
    // else PAY_ON_SINGLE_DATE which has already been set

    /* get the factor, between zero and one, that indicates 
    the proportion of the rebate that should be paid given 
    the relative positions of the historical sample and 
    range levels */

    double rebateFactor = rnType->pricePastSample();
    
    // no point continuing unless the rebate factor is non-zero
    if (Maths::isZero(rebateFactor))
    {
        return 0.0;
    }
    else
    {
        // discount the rebate from the payment date to the value date
        double discFactor;
        if (disc) {
             discFactor = discount->pv(valueDate, paymentDate);
        }
        else {discFactor = 1.0 ;}
        double numDaysInPd = 1.0; // to do: init # days in each period
        return rebateOnSampleDate * rebateFactor * discFactor / numDaysInPd ;
    }
}

double RangeNote::lookUpAmount(const DateTime& today, CashFlowArrayConstSP cfl) const
{
    double amount = 0.0;
    bool found = false;
    int i = 0;
    while( !found && i < cfl->size())
    {
        if ((*cfl)[i].date.equals(today))
        {
            found = true;
            amount = (*cfl)[i].amount;
        }
        i++;
    }
    
    if (!found)
    {
        throw ModelException("RangeNote::lookUpAmount", "could not find amount");
    }
    
    return amount;
}

// Discounted PV of remaining Libor cash flows happening strictly after evalDate
double RangeNote::lookUpRemainingAmounts(const DateTime& evalDate, CashFlowArrayConstSP cfl) const
{
    double amount,df;
    int i = 0;
    amount = 0.0;

    while (i < cfl->size())
    {
        if ((*cfl)[i].date > evalDate)
        {
            df = discount->pv(evalDate, (*cfl)[i].date);
            amount += df*(*cfl)[i].amount;
        }
        i++;
    }
    
    return amount;
}

/** find out date of end of accrual period corresponding to payment date payDate */
DateTime RangeNote::lookUpEndAccrualPeriod(const DateTime& payDate) const
{
    DateTime endAccrualDate;
    bool found = false;
    int i = 0;
    while( !found && i < floater->RefixDates.size())
    {
        if (floater->PayDates[i].equals(payDate))
        {
            found = true;
            endAccrualDate = floater->AccrualDates[i+1];
        }
        i++;
    }
    
    if (!found)
    {
        throw ModelException("RangeNote::lookUpEndAccrualPeriod", "could not find end of Accrual Period");
    }
    
    return endAccrualDate;
}

bool RangeNote::priceDeadInstrument(CControl* control, CResults* results) const{
    static const string method = "Callable KO Range Note::priceDeadInstrument";
    
    try  {
        
        double df,KOValue,callValue,Value,scalingFactor,fwdStrt;
        
        bool deadInstrument = false;

        bool isKO = isKnockedOut(valueDate);
        
        const DateTime endDate = (*paymentDates)[paymentDates->size() - 1];
        const DateTime matDate = (*monitorDates)[monitorDates->size() - 1];
                                                        
        // computing the scaling factor
        fwdStrt = 0.0;
        if (fwdStarting)
        {
            fwdStrt = asset->fwdValue(startDate);
        }

        scalingFactor = InstrumentUtil::scalePremium(oneContract,fwdStarting,notional,fwdStrt,initialSpot);

        if (isKO && isCalled)
        {
            throw    ModelException("RangeNote::priceDeadInstrument", "product cannot be both Called and KO");
        }

        // case product is Knocked Out
        if (isKO && !!barrierSchedule) 
        {
            // adding remaining historical coupons to be paid
            // here if a monitoring happens on knocked out it is not 
            // monitored, hence the "strict" getRemainingStrictHistPayments()
            // is used instead of the regular getRemainingHistPayments()
            KOValue = getRemainingStrictHistPayments(valueDate,  KnockedOutDate);
            
            // getting KO rebate settlement date
            DateTime rebateSettle = getKOPayDate(KnockedOutDate);

            if ( (KnockedOutDate <= valueDate) && (valueDate < rebateSettle) )
            {
                df = discount->pv(valueDate, rebateSettle);

                // adding remaining Libor payments: these occur between value date and rebateSettle
                if (!!floater)
                { 
                    double liborPayment = 
                        lookUpRemainingAmounts(valueDate, floater->getCashFlowArray(valueDate,discount.get()))
                        - df*lookUpRemainingAmounts(rebateSettle, floater->getCashFlowArray(valueDate,discount.get())); 

                    KOValue -= liborPayment;            
                }
            
                // adding KOrebate
                KOValue += fwdAtStart*KOrebate*df;                
            }

            KOValue *= scalingFactor;
            results->storePrice(KOValue, discount->getCcy());
            deadInstrument = true;  
        }

        // case product is called
        if (isCalled)
        {
            // validate that calledDate is a call and a monitoring date
            if(!!callSchedule)
            {
                int i = 0;
                bool callDateFound = false;
                bool monitorDateFound = false;

                while ( i< callSchedule->length() && !callDateFound )
                {
                    if( calledDate.equals(callSchedule->getDates()[i]) )
                    {
                        callDateFound = true;
                    }
                    i++;
                }

                i = 0;
                while ( i< monitorDates->size() && !monitorDateFound )
                {
                    if( calledDate.equals((*monitorDates)[i]) )
                    {
                        monitorDateFound = true;
                    }
                    i++;
                }

                if ( (!callDateFound || !monitorDateFound) && isCalled)
                {
                    throw ModelException(method,
                                         "Date Product called " + calledDate.toString()+
                                         " must be a call and monitoring date ");
                }  
            }

            callValue = 0.0;
            DateTime settlementDate = getPayDate(calledDate);
            // Now the end of period date is evaluated: after this date there is no more 
            // uncertainty which means priceDeadInstrument can be used
            DateTime endOfPeriodDate = getEndPeriodMonitoringDate(calledDate); 
            
            // priceDeadInstrument will be effective in a called case only 
            // if we are AFTER the end of period
            if ( valueDate > endOfPeriodDate )
            {
                    // adding value of pending historical coupons observed prior or equal
                    // to endOfPeriodDate, seen from valueDate. There are no Libor payments
                    // to add here because the last libor payment occurs at the endOfPeriodDate
                    callValue += getRemainingHistPayments(valueDate, endOfPeriodDate);

                    // check whether to include strike payment or not
                    if (valueDate < settlementDate)
                    {                                            
                        df = discount->pv(valueDate, settlementDate);                                        
                        callValue += df*callSchedule->interpolate(calledDate);
                    }

                    callValue *= scalingFactor;
                    results->storePrice(callValue, discount->getCcy());
                    deadInstrument = true;                                                             
            }
        }

        // case we are in between the last maturity date and the last payment date:
        if ( (!isCalled) && (!isKO) && (valueDate < endDate) && ( valueDate >= matDate) )
        {
            df = discount->pv(valueDate, endDate);
            
            // pv of remaining Range Note Coupons
            Value = getRemainingHistPayments(valueDate, matDate);
            
            // computing pv of remaining Libor payments
            if (!!floater)
            { 
                double liborPayment = lookUpRemainingAmounts(valueDate, floater->getCashFlowArray(valueDate,discount.get())); 
                Value -= liborPayment;            
            }

            // Capped Vanilla Option payoff
            if (!Maths::isZero(optionParticipation))
            {
                Value += df*Maths::min( initialSpot,
                                        optionParticipation
                                        *GetIntrinsic((*histMonSamples)[monitorDates->size() - 1],fwdAtStart*finalStrike,isCall, true/* this allows fwd */)); 
            }

            // Notional at maturity
            if (includeNotional)
            {
                if (fwdStarting)
                {
                        Value += fwdAtStart*df;                
                }
                else 
                {    
                    Value += initialSpot*df;
                }            
            }

            Value *= scalingFactor;

            // multiplying by scaling and storing result
            results->storePrice(Value, discount->getCcy());
            deadInstrument = true;     
        }

        // case date is too late 
        if (valueDate >= endDate) 
        {
            // all over, worth zero
            results->storePrice(0.0, discount->getCcy());
            if (control && control->isPricing()) {
                results->storePrice(0.0, discount->getCcy());
            }                        
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
        
        return deadInstrument;
    }
    catch (exception& e) 
    {
        if (rnType) { delete rnType; }
        throw ModelException(e, method);
    }            
}

DateTimeArraySP RangeNote::getPaymentDates() const
{
    if (paymentStyle == PAY_ON_MONITOR_DATE) {
        int iDate;
        DateTimeArraySP payDates(new DateTimeArray(monitorDates->size()));
        for(iDate=0; iDate<payDates->size();iDate++) {
            (*payDates)[iDate] = instSettle->settles((*monitorDates)[iDate], asset.get());
        }
        return payDates;
    }
    else {
        DateTimeArraySP payDates(new DateTimeArray(paymentDates->size()));
        payDates = paymentDates;
        return payDates;
    }
}

// cashflows known at the current valueDate
CashFlowArraySP RangeNote::getKnownCashFlows() const
{
    double scalingFactor,fwdStrt,Value;

    DateTime endDate;
    if (paymentStyle == PAY_ON_MONITOR_DATE) {
        endDate = instSettle->settles((*monitorDates)[monitorDates->size() - 1], asset.get());
    }
    else {
         endDate = (*paymentDates)[paymentDates->size() - 1];
    }
    
    const DateTime matDate = (*monitorDates)[monitorDates->size() - 1];
    
    DateTime pastDate = valueDate; // pastDate will be useful to extract relevant historical payments 

    CashFlowArraySP cfl(new CashFlowArray(0));
            
    bool isKO = isKnockedOut(valueDate);
                                                            
    // computing the scaling factor
    fwdStrt = 0.0;
    if (fwdStarting)
    {
        fwdStrt = asset->fwdValue(startDate);
    }
    scalingFactor = InstrumentUtil::scalePremium(oneContract,fwdStarting,notional,fwdStrt,initialSpot);

    // validation that product is not both KO and Called
    if (isKO && isCalled)
    {
        throw    ModelException("RangeNote::getKnownCashFlows", "product cannot be both Called and KO");
    }

    // case product is Knocked Out
    if ((!!barrierSchedule) && isKO && !isCalled) 
    {        
        // adding KO rebate if needed
        DateTime rebateSettle = getKOPayDate(KnockedOutDate);
        cfl->push_back(CashFlow(rebateSettle, scalingFactor*fwdAtStart*KOrebate));        
        pastDate = KnockedOutDate; // date up to which one computes historical payments, used below
    }

    // case product is called
    if ( (!!callSchedule) && isCalled && !isKO )
    {
        DateTime settlementDate = getPayDate(calledDate);    
        cfl->push_back(CashFlow(settlementDate, scalingFactor*callSchedule->interpolate(calledDate)));
        
        // Now it is needed to know up to when are historical payments known
        // After calledDate and up to the end of monitoring period, monitorings contine to count
        // hence one needs to check if valueDate is before the end of the monitoring period or not
        if (valueDate >= calledDate)
        {
            DateTime endOfPeriodDate = getEndPeriodMonitoringDate(calledDate);
            if ( valueDate > endOfPeriodDate)
            {
                pastDate = endOfPeriodDate;
            }
            else 
            {
                pastDate = valueDate;
            }
        }        
    }

    // adding all known historical payments monitored before or on pastDate
    int i = 0;
    DateTime currentDate;
    DateTime settleDate;
    double pastSpotValue = 0.0;

    while ( i<monitorDates->size() && (*monitorDates)[i] <= pastDate ) 
    {
        currentDate = (*monitorDates)[i];
        settleDate = getPayDate(currentDate);
        
        pastSpotValue = (*histMonSamples)[i];

        if ( (!!highRangeSchedule.get())
                 && (!!lowRangeSchedule.get())
                 && (pastSpotValue < fwdAtStart*highRangeSchedule->interpolate(currentDate))
                 && (pastSpotValue > fwdAtStart*lowRangeSchedule->interpolate(currentDate)) )
        {    
            cfl->push_back(CashFlow(settleDate, scalingFactor*rebateSchedule->interpolate(currentDate)));
        }        
        i++;
    }

    // remaining Libor payments
    if (!!floater)
    { 
        // looping through to get pending payments
        i =0;
        while ( i<floater->RefixDates.size() && floater->RefixDates[i] <= pastDate ) 
        {
            currentDate = floater->RefixDates[i];
            settleDate = floater->PayDates[i];
            
            double liborPayment = lookUpAmount(settleDate,floater->getCashFlowArray(valueDate,discount.get()));
            cfl->push_back(CashFlow(settleDate, -scalingFactor*liborPayment));
            
            i++;
        }        
    }
    
    // case we are after the last maturity date 
    if ( (!isCalled) && (!isKO) && ( valueDate >= matDate))
    {                    
        // Capped Vanilla Option payoff is known
        if (!Maths::isZero(optionParticipation))
        {
            Value = scalingFactor*Maths::min( initialSpot,
                                    optionParticipation
                                    *GetIntrinsic((*histMonSamples)[monitorDates->size() - 1],fwdAtStart*finalStrike,isCall, true/* this allows fwd */)); 
        
            cfl->push_back(CashFlow(endDate, Value));
        }
    }

    // Notional at maturity: known as a sure payment only if conditions depending
    // on whether the product has been called and /or KO or not
    if (includeNotional)
    {
        bool knownNotionalIncluded = true;
        
        if (!!callSchedule && (!(callSchedule->length() == 0)))
        {
            if (isCalled || (valueDate <= callSchedule->lastDate()) )
            knownNotionalIncluded = false;
        }

        if (!!barrierSchedule && (!(barrierSchedule->length() == 0))) 
        {
            if (isKO || (valueDate <= barrierSchedule->lastDate()) )
            knownNotionalIncluded = false;
        }

        if (knownNotionalIncluded)
        {
            if (fwdStarting)
            {    
                    cfl->push_back(CashFlow(endDate, scalingFactor*fwdAtStart));
            }
            else 
            {    
                    cfl->push_back(CashFlow(endDate, scalingFactor*initialSpot));
            }
        }        
    }

    return cfl;
}

// Discounted value of historical payments recorded prior or equal to today, but not paid yet at evalDate
// evalDate must be strictly before the payment dates for the corresponding payment to be considered.
// When evalDate is exactly  equal to a payment date, by convention payment is not taken into account
double RangeNote::getRemainingHistPayments(const DateTime& evalDate, const DateTime& today) const
{
    int i = 0;
    double df = 0.0;
    double pastValue = 0.0;
    double pastSpotValue = 0.0;

    DateTime currentDate;
    DateTime settleDate;
    DateTime lastPayDate = getNextPaymentDate(today);

    if (evalDate < lastPayDate) 
    {    
        while ( i<monitorDates->size() && (*monitorDates)[i] <= today ) 
        {
            currentDate = (*monitorDates)[i];
            settleDate = getPayDate(currentDate);
            if ( settleDate > evalDate )
            {
                df = discount->pv(evalDate, settleDate);    
                pastSpotValue = (*histMonSamples)[i];

                if ( (!!highRangeSchedule.get())
                         && (!!lowRangeSchedule.get())
                         && (pastSpotValue < fwdAtStart*highRangeSchedule->interpolate(currentDate))
                         && (pastSpotValue > fwdAtStart*lowRangeSchedule->interpolate(currentDate)) )
                {    
                    pastValue += rebateSchedule->interpolate(currentDate)*df;
                }
            }
            i++;
        }    
    }

    return pastValue;
}

// Same function getRemainingHistPayments() except coupon observed today is not taken into account
double RangeNote::getRemainingStrictHistPayments(const DateTime& evalDate, const DateTime& today) const
{
    int i = 0;
    double df = 0.0;
    double pastValue = 0.0;
    double pastSpotValue = 0.0;

    DateTime currentDate;
    DateTime settleDate;
    DateTime lastPayDate = getNextPaymentDate(today);

    if (evalDate < lastPayDate) 
    {    
        // monitor date must be strictly before today
        while ( i<monitorDates->size() && (*monitorDates)[i] < today ) 
        {
            currentDate = (*monitorDates)[i];
            settleDate = getPayDate(currentDate);
            if ( settleDate > evalDate )
            {
                df = discount->pv(evalDate, settleDate);    
                pastSpotValue = (*histMonSamples)[i];

                if ( (!!highRangeSchedule.get())
                         && (!!lowRangeSchedule.get())
                         && (pastSpotValue < fwdAtStart*highRangeSchedule->interpolate(currentDate))
                         && (pastSpotValue > fwdAtStart*lowRangeSchedule->interpolate(currentDate)) )
                {    
                    pastValue += rebateSchedule->interpolate(currentDate)*df;
                }
            }
            i++;
        }    
    }

    return pastValue;
}

// After ensuring that testDate is a monitoring date, the below 
// function checks whether we are at the end of a monitoring period
// defined as: testDate and the monitoring date immediately after
// have different paymentDates. In case testDate is the last monitoring
// date then isEndMonitoringPeriod() returns true
bool RangeNote::isEndMonitoringPeriod(const DateTime& testDate) const
{
    int i = 0;
    bool isEndPeriod = false;    
    bool monitoringDateFound = false;
    bool nextMonitorDateFound = false;
    DateTime nextMonitoringDate;

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

    if (!monitoringDateFound)
    {
        throw ModelException("RangeNote::isEndMonitoringPeriod", "date entered is not a monitoring date");
    }

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

bool RangeNote::isCallPayDate (const DateTime& testDate) const
{
    int i = 0;
    bool isCallPaymentDate = false;
    
    while ( i< callSchedule->length() && !isCallPaymentDate )
    {
        if ( testDate.equals( getPayDate(callSchedule->getDates()[i]) ) )
        {
            isCallPaymentDate    = true;
        }
    }

    return isCallPaymentDate;
}

// starts from a monitoring date and gives the first equal or future end of period 
// monitoring date which is useful in priceDeadInstrument()
const DateTime& RangeNote::getEndPeriodMonitoringDate (const DateTime& today) const
{
    
    if (isEndMonitoringPeriod(today))
    {
        return today;
    }
    else
    {
        int i = 0;
        bool endDateFound = false;
        while ( i<monitorDates->size() && !endDateFound )
        {
            if(today <= (*monitorDates)[i] )
            {
                if (isEndMonitoringPeriod((*monitorDates)[i]) )
                {
                    endDateFound = true;
                }
            }
            i++;
        }

        if (!endDateFound)
        {
            throw ModelException("RangeNote::getEndPeriodMonitoringDate", "did not find end of period date");
        }

        return (*monitorDates)[i-1];
    }    
}

double RangeNote::getHistPrice(const DateTime& histDate) const
{    
    int i = 0;
    bool priceFound = false;

    while ( i<monitorDates->size() && !priceFound )
    {
        if (histDate.equals((*monitorDates)[i], false))
        {
            priceFound = true;
        }
        i++;
    }

    if (priceFound)
    {
        return    (*histMonSamples)[i-1]; 
    }

    else
    {
        throw ModelException("RangeNote::getHistPrice", "date entered is not a historical date");
    }        
}

// To know whether product is already KO before or at date today
// If it is the case, KnockedOutDate will store the KO date
bool RangeNote::isKnockedOut(const DateTime& today) const 
{
    int i,j;
    bool isKO = false;
    double spotValue;
    DateTime monitoringDate;
    DateTime KOdate;

    i = 0;
    fwdAtStart = 1.0;
    // Get spot at start if needed
    if (fwdStarting)
    {
        fwdAtStart = asset->fwdValue(startDate);
    }

    if(!!(barrierSchedule) && (!(barrierSchedule->length() == 0)))
    {
        while( (!isKO) && (i< monitorDates->size()) && ((*(monitorDates))[i] <= today) )
        {
            j = 0;
            spotValue = (*(histMonSamples))[i];
            monitoringDate = (*(monitorDates))[i];

            while (!isKO && j< barrierSchedule->length() )
            {            
                KOdate = barrierSchedule->getDates()[j];
                KOlevel =  fwdAtStart*barrierSchedule->interpolate(KOdate);

                if ( monitoringDate.equals(KOdate) )                                                   
                {
                    if( ( (isUpAndOutKO)&&(spotValue >= KOlevel) ) 
                          || ( (!isUpAndOutKO)&&(spotValue <= KOlevel) ) )
                    {     
                        isKO = true;
                        KnockedOutDate = KOdate;
                    }
                }
                j++;
            }
            i++;
        }
    }

    return isKO;
}



DoubleArraySP RangeNote::d_prices = DoubleArraySP(   );

/** find pay date in PaymentDates */
const DateTime& RangeNote::getPayDate(const DateTime& today) const
{ 
    int i = 0;
    bool found = false;
    
    if (paymentStyle == PAY_ON_SINGLE_DATE) {
        return (*paymentDates)[0];
    }
    else {
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
            if (paymentStyle == PAY_ON_MONITOR_DATE) {
                return (*monitorDates)[i - 1];
            }
            else {
                return (*paymentDates)[i - 1];
            }
        }
        else
        {
            throw ModelException("CallableRN::getPayDate", "date supplied is not in monitoring dates!");
        }
    }
}

/** find KO pay date in PaymentDates. The particular feature here is that in case
a KO is triggered, the payment date will be the payment date corresponding to the previous
monitoring' pay date. In case this value is already in the past, or KO is on the first monitoring date,
KO pay date is by convention set to today  */
const DateTime& RangeNote::getKOPayDate(const DateTime& today) const
{ 
    int i = 0;
    bool found = false;
    
    while (!found &&  i < monitorDates->size() )
    {
        // don't compare the time since multiple KO steps
        // can fall on the same day
        if ((*monitorDates)[i].equals(today, false)) 
        {
            found = true;
        }

        i++;
    }

    if (found)
    {
        if ( (i > 1) && ((*paymentDates)[i - 2] >= today) )
        {
            return (*paymentDates)[i - 2];
        }
        else
        // previous monitoring's pay date is already in the past or does not exist 
        // hence take today
        {
            return (*paymentDates)[i - 1];
        }
    }
    else
    {
        throw ModelException("CallableRN::getKOPayDate", "could not find KO payment date!");
    }
}


/** private class */
class RangeNoteClosedForm: public CClosedFormLN::IProduct{
private:
    const RangeNote* rngNt; // a reference

public:
    RangeNoteClosedForm(const RangeNote* rngNt): rngNt(rngNt){}

    void price(CClosedFormLN* model,
               Control*    control, 
               CResults*   results)const{
        rngNt->price(control, results);
    }
};

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* RangeNote::createProduct(CClosedFormLN* model) const
{
    return new RangeNoteClosedForm(this);
}

RNtype::RNtype(const RangeNote* rn):rn(rn)
{
    // empty
}

double RNtype::singleRangeAsDigitalFactor(const Schedule& rangeSchedule,
                                          bool upsideRange)
{
    double rebateFactor = 0.0;
    
    // determine the range level on the sample date
    double rangeLevel = rangeSchedule.interpolate(rn->sampleDate);

    /* See if the sample is on the range edge 
       as always get paid in this case */
    if (Maths::equals(rn->histLevelOnSampDate, rangeLevel) ||
        (upsideRange && (rn->histLevelOnSampDate > rangeLevel)) ||
        (!upsideRange && (rn->histLevelOnSampDate < rangeLevel)))
    {
        rebateFactor = 1.0;
    }
    return rebateFactor;
}

double RNtype::singleRangeAsCallSpreadFactor(const Schedule& rangeSchedule,
                                             double spread,
                                             bool upsideRange,
                                             double& lowStrike,
                                             double& highStrike)
{
    lowStrike = 0.0;
    highStrike = 0.0;
    
    double rebateFactor = 0.0;
    
    // determine the range level on the sample date

    double rangeLevel = rangeSchedule.interpolate(rn->sampleDate);
    
    /* Get the strikes to be used in the VanillaSpread product. 
       These will be calculated from the rangeLevel on the sample date and 
       the spread parameters */
    getSpreadStrikes(rn->spreadType,
                     rangeLevel,
                     upsideRange,
                     spread,
                     lowStrike,
                     highStrike);

    // see if sample was in payout range 
    if (upsideRange)
    {
        // determine the call spread payout 
        rebateFactor = Maths::max(rn->histLevelOnSampDate - lowStrike, 0.0) -
            Maths::max(rn->histLevelOnSampDate - highStrike, 0.0);
    }
    else
    {
        // determine the put spread payout 
        rebateFactor = Maths::max(highStrike - rn->histLevelOnSampDate, 0.0) -
            Maths::max(lowStrike  - rn->histLevelOnSampDate, 0.0);
    }
    
    // normalise the rebate factor 
    
    rebateFactor /= (highStrike - lowStrike);

    return rebateFactor;
}


double RNsingleUpside::pricePastSample()
{
    static const string method = "RNsingleUpside::pricePastSample";

    double rebateFactor = 0.0;
    double lowStrike = 0.0;
    double highStrike = 0.0;
    try
    {
        if (rn->priceHistAsDigital)
        {
            rebateFactor = singleRangeAsDigitalFactor(*(rn->lowRangeSchedule),
                                                      true);
        }
        else
        {
            rebateFactor = singleRangeAsCallSpreadFactor(*(rn->lowRangeSchedule),
                                                         rn->lowStrikeSpread,
                                                         true,
                                                         lowStrike,
                                                         highStrike);
            
        }
        return rebateFactor;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

double RNsingleDownside::pricePastSample()
{
    static const string method = "RNsingleDownside::pricePastSample";

    double rebateFactor = 0.0;
    double lowStrike = 0.0;
    double highStrike = 0.0;

    try
    {
        if (rn->priceHistAsDigital)
        {
            rebateFactor = singleRangeAsDigitalFactor(*(rn->lowRangeSchedule),
                                                      false);
        }
        else
        {
            rebateFactor = singleRangeAsCallSpreadFactor(*(rn->lowRangeSchedule),
                                                         rn->lowStrikeSpread,
                                                         false,
                                                         lowStrike,
                                                         highStrike);
        
        }
        return rebateFactor;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

double RNdoubleInside::pricePastSample()
{
    static const string method = "RNdoubleInside::pricePastSample";

    double rebateFactor = 0.0;
    try
    {
        if (rn->priceHistAsDigital)
        {
            // see if we are above the low strike 
            
            /* comment for second parameter - for an inside range the low strike
               an upside edge */
            rebateFactor = singleRangeAsDigitalFactor(*(rn->lowRangeSchedule),
                                                      true);
            
            /* no need to check the other barrier as we are definitely 
               outside the range */
            if (!Maths::isZero(rebateFactor))
            {
                rebateFactor = singleRangeAsDigitalFactor(*(rn->highRangeSchedule),
                                                          false); // downside edge
            }
        }
        else
        {
            double lrLowStrike = 0.0;
            double lrHighStrike = 0.0;
            double hrLowStrike = 0.0;
            double hrHighStrike = 0.0;
            
            // see if we are above the low strike
            
            /* comment for second parameter - for an inside range the low strike
               an upside edge */
            double lrRebateFactor = singleRangeAsCallSpreadFactor(*(rn->lowRangeSchedule),
                                                                  rn->lowStrikeSpread,
                                                                  true,
                                                                  lrLowStrike,
                                                                  lrHighStrike);
            
            /*although in theory we can stop if lrRebateFactor is less than 1.0
              continue as we need to validate that the spreads don't overlap after
              spread */
            
         /* note - can price as a long downside rather than a short upside as 
            we don't sum the rebate factors */
            
            double hrRebateFactor = singleRangeAsCallSpreadFactor(*(rn->highRangeSchedule),
                                                                  rn->highStrikeSpread,
                                                                  false,
                                                                  hrLowStrike,
                                                                  hrHighStrike);
            
            if (hrLowStrike < lrHighStrike)
            {
                throw ModelException(method,
                                     "Call spread strikes for the low and high "
                                     "barriers overlap at monitoring date " 
                                     + rn->sampleDate.toString() +
                                     ". High strike of the call spread for the "
                                     " low barrier is " 
                                     + Format::toString(lrHighStrike)+
                                     " and the low strike of the call spread "
                                     "for the high barrier is "
                                     + Format::toString(hrLowStrike)+
                                  ". Please adjust your spread parameters.");
            }
            
            /* now set the rebate factor 
               for an inside range you need to be above the low barrier 
               AND below the low barrier */
            rebateFactor = Maths::min(lrRebateFactor, hrRebateFactor);
        }
        
        return rebateFactor;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

double RNdoubleOutside::pricePastSample()
{
    static const string method = "RNdoubleOutside::pricePastSample";
    
    double rebateFactor = 0.0;
    try
    {
        if (rn->priceHistAsDigital)
        {
         /* see if we are below the low strike 
            comment for second parameter - for an outside range the low 
            strike is a downside edge */
         
            rebateFactor = singleRangeAsDigitalFactor(*(rn->lowRangeSchedule),
                                                      false);
            
            /* no need to check the other barrier as we are definitely 
               outside the range */
            if (!Maths::equals(rebateFactor, 1.0))
            {
                rebateFactor = singleRangeAsDigitalFactor(*(rn->highRangeSchedule),
                                                          true); // downside edge
            }
        }
        else
        {
            double lrLowStrike = 0.0;
            double lrHighStrike = 0.0;
            double hrLowStrike = 0.0;
            double hrHighStrike = 0.0;
            
            /* see if we are below the low strike 
               comment for second parameter - for an outside range the low 
               strike is a downside edge */
            double lrRebateFactor = singleRangeAsCallSpreadFactor(*(rn->lowRangeSchedule),
                                                                  rn->lowStrikeSpread,
                                                                  false,
                                                                  lrLowStrike,
                                                                  lrHighStrike);
            
            /*although in theory we can stop if lrRebateFactor is less than 1.0
              continue as we need to validate that the spreads don't overlap after
              spread */
            double hrRebateFactor = singleRangeAsCallSpreadFactor(*(rn->highRangeSchedule),
                                                                  rn->highStrikeSpread,
                                                                  true,
                                                                  hrLowStrike,
                                                                  hrHighStrike);
            
            if (hrLowStrike < lrHighStrike)
            {
                throw ModelException(method,
                                     "Call spread strikes for the low and high "
                                     "barriers overlap at monitoring date " 
                                     + rn->sampleDate.toString() +
                                     ". High strike of the call spread for the "
                                     " low barrier is " 
                                     + Format::toString(lrHighStrike)+
                                     " and the low strike of the call spread "
                                     "for the high barrier is "
                                     + Format::toString(hrLowStrike)+
                                     ". Please adjust your spread parameters.");
         }
            
            /* now set the rebate factor 
               For an outside range you get paid if you are below
               the low barrier OR above the high barrier */
            
            rebateFactor = Maths::max(lrRebateFactor, hrRebateFactor);
        }
        return rebateFactor;
     }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

/** returns the low and high strikes for a range level depending on the
    spread type of the range note */
void RNtype::getSpreadStrikes(const string& spreadType,
                              double rangeLevel,
                              bool upsideRange,
                              double strikeSpread,
                              double& lowStrike,
                              double& highStrike)
{
    static const string method = "RNtype::getSpreadStrikes";

    /* Set strikes a factor of alpha either side of the range 
       therefore not changing the area under the payout range */
    if (spreadType == RangeNote::SPREAD_NONE)
    {
        lowStrike  = (1.0 - rn->MINIMUM_ALPHA)*rangeLevel;
        highStrike = Maths::isZero(rangeLevel)? (2 * rn->MINIMUM_ALPHA): (1.0 + rn->MINIMUM_ALPHA)*rangeLevel;
    }
    else
    {
        /* The strike depends on whether we want to increase or 
           decrease the payout range and on which side of the strike you 
           get paid 
           
           NB. If you have sold then when spread, increase the payout 
           range and vice versa if long 
           
           upsideRange == TRUE means that you get paid if the underlying 
           is higher than the range level. 
           upsideRange == FALSE means your paid if the underlying is below
           the range level */
        
        if (((spreadType == RangeNote::SPREAD_DECREASING) && upsideRange) ||
            ((spreadType == RangeNote::SPREAD_INCREASING) && !upsideRange))
        {
            lowStrike  = rangeLevel;
            highStrike = rangeLevel + strikeSpread;
        }
        else
        {
            highStrike = rangeLevel;
            lowStrike  = rangeLevel - strikeSpread;
        }
    }

    // ensure that we have not generated any negative strikes
    if (Maths::isNegative(highStrike) || Maths::isNegative(lowStrike))
    {
        string lo(Format::toString(lowStrike));
        string hi(Format::toString(highStrike));
        string m("One of the strikes generated for spread "
                 "purposes is negative! Low strike = "+lo+
                 " High strike = "+hi);
        
        throw ModelException(method, m);
    }
}

/** price as a call spread */
double RNsingleUpside::priceFutureSample()
{
    static const string method = "RNsingleUpside::priceFutureSample";
    
    double lowStrike = 0.0;
    double highStrike = 0.0;

    try
    {
        /* determine the range level on the sample date 
           the spread strikes will be based on this */
        double rangeLevel = rn->lowRangeSchedule->interpolate(rn->sampleDate);
        
        // price the spread for this range
        double spreadPrice = priceSpread(rn->spreadType,
                                         rangeLevel,
                                         rn->lowStrikeSpread,
                                         true,  // price as long call spread
                                         lowStrike,
                                         highStrike);

        
        return spreadPrice;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

/** price as a put spread */
double RNsingleDownside::priceFutureSample()
{
    static const string method = "RNsingleDownside::priceFutureSample";
    
    double lowStrike = 0.0;
    double highStrike = 0.0;

    try
    {
        /* determine the range level on the sample date 
           the spread strikes will be based on this */
        double rangeLevel = rn->lowRangeSchedule->interpolate(rn->sampleDate);
        
        // price the spread for this range
        double spreadPrice = priceSpread(rn->spreadType,
                                         rangeLevel,
                                         rn->lowStrikeSpread,
                                         false,  // price as long put spread
                                         lowStrike,
                                         highStrike);

        
        return spreadPrice;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

/** price as a long call spread (low range) and a short call spread (high range) */
double RNdoubleInside::priceFutureSample()
{
    static const string method = "RNdoubleInside::priceFutureSample";
    double lrLowStrike = 0.0;
    double lrHighStrike = 0.0;
    double hrLowStrike = 0.0;
    double hrHighStrike = 0.0;

    try
    {
        /* Determine the range levels on the sample date 
           the spread strikes will be based on this */
        double lowRangeLevel = rn->lowRangeSchedule->interpolate(rn->sampleDate);
        double highRangeLevel = rn->highRangeSchedule->interpolate(rn->sampleDate);
        
        // price the low range spread as a long call spread
        double lowSpreadPrice = priceSpread(rn->spreadType,
                                            lowRangeLevel,
                                            rn->lowStrikeSpread,
                                            true,  // call spread
                                            lrLowStrike,
                                            lrHighStrike);
       
        
        /* Price the high range spread as a call spread also.
           This will be shorted later. 
           Since we go short this call spread - to get the correct 
           spread profile once the payoff has been inverted we need 
           to change the spread payoff from inside to outside or 
           vice versa */
        string spreadType;
        if (rn->spreadType == RangeNote::SPREAD_NONE)
        {
            // as this is symetric about the x axis don't do anything
            spreadType = rn->spreadType;
        }
        else if (rn->spreadType == RangeNote::SPREAD_DECREASING)
        {
            spreadType = RangeNote::SPREAD_INCREASING;
        }
        else
        {
            // must be RangeNote::INCREASING so reverse this
            spreadType = RangeNote::SPREAD_DECREASING;
        }

        // price the low range spread as a long call spread
        double highSpreadPrice = priceSpread(spreadType,
                                             highRangeLevel,
                                             rn->highStrikeSpread,
                                             true,  // call spread
                                             hrLowStrike,
                                             hrHighStrike);
        
        // check that there is no overlap
        if (hrLowStrike < lrHighStrike)
        {
            throw ModelException(method,
                                 "Call spread strikes for the low and high "
                                 "barriers overlap at monitoring date " 
                                 + rn->sampleDate.toString() +
                                 ". High strike of the call spread for the "
                                 " low barrier is " 
                                 + Format::toString(lrHighStrike)+
                                 " and the low strike of the call spread "
                                 "for the high barrier is "
                                 + Format::toString(hrLowStrike)+
                                 ". Please adjust your spread parameters.");
        }
        
        return (lowSpreadPrice - highSpreadPrice);
        
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

/** price as a long put spread (low range) and a long call spread (high range) */ 
double RNdoubleOutside::priceFutureSample()
{
    static const string method = "RNdoubleOutside::priceFutureSample";
    double lrLowStrike = 0.0;
    double lrHighStrike = 0.0;
    double hrLowStrike = 0.0;
    double hrHighStrike = 0.0;

    try
    {
        /* Determine the range levels on the sample date 
           the spread strikes will be based on this */
        double lowRangeLevel = rn->lowRangeSchedule->interpolate(rn->sampleDate);
        double highRangeLevel = rn->highRangeSchedule->interpolate(rn->sampleDate);
        
        // price the low range spread as a long put spread
        double lowSpreadPrice = priceSpread(rn->spreadType,
                                            lowRangeLevel,
                                            rn->lowStrikeSpread,
                                            false,  // put spread
                                            lrLowStrike,
                                            lrHighStrike);
        
        

        // price the high range spread as a long call spread
        double highSpreadPrice = priceSpread(rn->spreadType,
                                             highRangeLevel,
                                             rn->highStrikeSpread,
                                             true,  // call spread
                                             hrLowStrike,
                                             hrHighStrike);
        
        // check that there is no overlap
        if (hrLowStrike < lrHighStrike)
        {
            throw ModelException(method,
                                 "Call spread strikes for the low and high "
                                 "barriers overlap at monitoring date " 
                                 + rn->sampleDate.toString() +
                                 ". High strike of the call spread for the "
                                 " low barrier is " 
                                 + Format::toString(lrHighStrike)+
                                 " and the low strike of the call spread "
                                 "for the high barrier is "
                                 + Format::toString(hrLowStrike)+
                                 ". Please adjust your spread parameters.");
        }
        
        return (lowSpreadPrice + highSpreadPrice);
        
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

/** Price a vanilla spread for a single range */
double RNtype::priceSpread(const string& spreadType,
                           double rangeLevel,
                           double spread,
                           bool priceAsUpside,
                           double& lowStrike,
                           double& highStrike)
{
    static const string method = "RNtype::priceSpread";

    try
    {
        /* Get the strikes to be used in the vanilla spread 
           These will be calculated from the rangeLevels on the 
           sample date and the spread parameters */
        getSpreadStrikes(spreadType,
                         rangeLevel,
                         priceAsUpside,   // price an upside range
                         spread,
                         lowStrike,
                         highStrike);

        /* For an upside range level we price a call spread,
           otherwise we price a put spread. */
        
        double spreadPrice = CVanilla::priceSpread(rn->valueDate,
                                                   rn->startDate,
                                                   rn->sampleDate,
                                                   priceAsUpside,   // true=call
                                                   rn->fwdStarting,
                                                   true,   // oneContract
                                                   1.0,    // notional
                                                   1.0,    // initialSpot
                                                   lowStrike,
                                                   highStrike,
                                                   rn->sampleSettlement.get(),
                                                   rn->asset.get(),
                                                   rn->discount.get());

        /* caluate the number of options that we need to price
           to get the correct rebate when in the range */
        double numOptions;
        if (rn->spreadType == RangeNote::SPREAD_NONE)
        {
            numOptions = Maths::isZero(rangeLevel)? (rn->rebateOnSampleDate/ (2 * rn->MINIMUM_ALPHA)) 
                :rn->rebateOnSampleDate/(2.0*rn->MINIMUM_ALPHA*rangeLevel);
        }
        else
        {
            numOptions = rn->rebateOnSampleDate/spread;
        }
        
        if (rn->fwdStarting)
        {
            double fwdAtStart = rn->asset->fwdValue(rn->startDate);
            spreadPrice *=  numOptions/fwdAtStart;
        }
        else
        {
            spreadPrice *= numOptions;   
        }

        return spreadPrice;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

DoubleArray* RNsingleUpside::getSensStrikes(OutputNameConstSP outputName)
{
    static const string method = "RNsingleUpside::getSensStrikes";
    // The barrier must be flat to be executing this code
    double barrierLevel = rn->lowRangeSchedule->lastValue();
    double spread = rn->lowStrikeSpread;
    
    DoubleArraySP sensStrikes(new DoubleArray(0)); 

    if (rn->spreadType == RangeNote::SPREAD_INCREASING)
    {
        sensStrikes->push_back(barrierLevel - spread);
        sensStrikes->push_back(barrierLevel);
    }

    else if (rn->spreadType == RangeNote::SPREAD_DECREASING)
    {
        sensStrikes->push_back(barrierLevel);
        sensStrikes->push_back(barrierLevel + spread);
    }
    
    else if (rn->spreadType == RangeNote::SPREAD_NONE)
    {
        /* dont need to account for fwd starting here because alpha
           is always a percentage */
        spread = 2 * barrierLevel * rn->MINIMUM_ALPHA;
        
        sensStrikes->push_back(barrierLevel - (spread / 2));
        sensStrikes->push_back(barrierLevel + (spread / 2));
    }
    else
    {
        // should never get here
        throw ModelException(method,
                             "Unknown spread type "+
                             rn->spreadType);
    }
    
    return sensStrikes.release();
}

DoubleArray* RNsingleDownside::getSensStrikes(OutputNameConstSP outputName)
{
    static const string method = "RNsingleDownside::getSensStrikes";
    // The barrier must be flat to be executing this code
    double barrierLevel = rn->lowRangeSchedule->lastValue();
    double spread = rn->lowStrikeSpread;
    
    DoubleArraySP sensStrikes(new DoubleArray(0)); 

    if (rn->spreadType == RangeNote::SPREAD_DECREASING)
    {
        sensStrikes->push_back(barrierLevel - spread);
        sensStrikes->push_back(barrierLevel);
    }

    else if (rn->spreadType == RangeNote::SPREAD_INCREASING)
    {
        sensStrikes->push_back(barrierLevel);
        sensStrikes->push_back(barrierLevel + spread);
    }

    else if (rn->spreadType == RangeNote::SPREAD_NONE)
    {
        /* dont need to account for fwd starting here because alpha
           is always a percentage */
        spread = 2 * barrierLevel * rn->MINIMUM_ALPHA;

        sensStrikes->push_back(barrierLevel - (spread / 2));
        sensStrikes->push_back(barrierLevel + (spread / 2));
    }
    else
    {
        // should never get here
        throw ModelException(method,
                             "Unknown spread type "+
                             rn->spreadType);
    }

    return sensStrikes.release();
}

DoubleArray* RNdoubleInside::getSensStrikes(OutputNameConstSP outputName)
{
    static const string method = "RNdoubleInside::getSensStrikes";
    // The barrier must be flat to be executing this code
    double lowBarrierLevel = rn->lowRangeSchedule->lastValue();
    double lowSpread = rn->lowStrikeSpread;
    double highBarrierLevel = rn->highRangeSchedule->lastValue();
    double highSpread = rn->highStrikeSpread;
    
    DoubleArraySP sensStrikes(new DoubleArray(0)); 

    if (rn->spreadType == RangeNote::SPREAD_INCREASING)
    {
        sensStrikes->push_back(lowBarrierLevel - lowSpread);
        sensStrikes->push_back(lowBarrierLevel);
        sensStrikes->push_back(highBarrierLevel);
        sensStrikes->push_back(highBarrierLevel + highSpread);
    }

    else if (rn->spreadType == RangeNote::SPREAD_DECREASING)
    {
        sensStrikes->push_back(lowBarrierLevel);
        sensStrikes->push_back(lowBarrierLevel + lowSpread);
        sensStrikes->push_back(highBarrierLevel - highSpread);
        sensStrikes->push_back(highBarrierLevel);
    }

    else if (rn->spreadType == RangeNote::SPREAD_NONE)
    {
        /* dont need to account for fwd starting here because alpha
           is always a percentage */
        lowSpread = 2 * lowBarrierLevel * rn->MINIMUM_ALPHA;
        highSpread = 2 * highBarrierLevel * rn->MINIMUM_ALPHA;

        sensStrikes->push_back(lowBarrierLevel - (lowSpread / 2));
        sensStrikes->push_back(lowBarrierLevel + (lowSpread / 2));
        sensStrikes->push_back(highBarrierLevel - (highSpread / 2));
        sensStrikes->push_back(highBarrierLevel + (highSpread / 2));
    }
    else
    {
        // should never get here
        throw ModelException(method,
                             "Unknown spread type "+
                             rn->spreadType);
    }

    return sensStrikes.release();
}

DoubleArray* RNdoubleOutside::getSensStrikes(OutputNameConstSP outputName)
{
    static const string method = "RNdoubelOutside::getSensStrikes";
    // The barrier must be flat to be executing this code
    double lowBarrierLevel = rn->lowRangeSchedule->lastValue();
    double lowSpread = rn->lowStrikeSpread;
    double highBarrierLevel = rn->highRangeSchedule->lastValue();
    double highSpread = rn->highStrikeSpread;
    
    DoubleArraySP sensStrikes(new DoubleArray(0)); 

    if (rn->spreadType == RangeNote::SPREAD_DECREASING)
    {
        sensStrikes->push_back(lowBarrierLevel - lowSpread);
        sensStrikes->push_back(lowBarrierLevel);
        sensStrikes->push_back(highBarrierLevel);
        sensStrikes->push_back(highBarrierLevel + highSpread);
    }

    else if (rn->spreadType == RangeNote::SPREAD_INCREASING)
    {
        sensStrikes->push_back(lowBarrierLevel);
        sensStrikes->push_back(lowBarrierLevel + lowSpread);
        sensStrikes->push_back(highBarrierLevel - highSpread);
        sensStrikes->push_back(highBarrierLevel);
    }

    else if (rn->spreadType == RangeNote::SPREAD_NONE)
    {
        /* dont need to account for fwd starting here because alpha
           is always a percentage */
        lowSpread = 2 * lowBarrierLevel * rn->MINIMUM_ALPHA;
        highSpread = 2 * highBarrierLevel * rn->MINIMUM_ALPHA;

        sensStrikes->push_back(lowBarrierLevel - (lowSpread / 2));
        sensStrikes->push_back(lowBarrierLevel + (lowSpread / 2));
        sensStrikes->push_back(highBarrierLevel - (highSpread / 2));
        sensStrikes->push_back(highBarrierLevel + (highSpread / 2));
    }
    else
    {
        // should never get here
        throw ModelException(method,
                             "Unknown spread type "+
                             rn->spreadType);
    }

    return sensStrikes.release();
}
        
// Gives the payment date corresponding to the monitoring date that is equal to 
// or the first after today            
DateTime RangeNote::getNextPaymentDate(const DateTime& today) const
{
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
    if (!found)
    {
        throw ModelException("RangeNote::getNextPaymentDate", " date " + today.toString() + 
            " is after all monitoring dates!");
    }
    
    return nextPayDate;
}

/******************************************************************************************************************************/
// product class
class RangeNoteFDProd: public LatticeProdEDRIns
{
public:

    RangeNoteFDProd(const RangeNote* rn, FDModel* model) :
        LatticeProdEDRIns(model, 1, 3), inst(rn) 
    {
        // first: set discount curve
        if( tree1f )
            tree1f->setDiscountCurve( inst->discount.getSP() );

        // second: create spot payoff
        payoffIndex = model->createProduct( IProdCreatorSP( new
            IndexSpecEQ( inst->asset.getName(), inst->asset, inst->ccyTreatment ) ) );

        calledOnValueDate = false;       
    }

    virtual ~RangeNoteFDProd(){}

    /** vanilla, quanto or struck */
    virtual string getCcyTreatment() const
    {
       return inst->ccyTreatment;
    }
     
    /** calculate at barriers for tree */
    virtual void preCalc(int step); 

    /** product payoff method at maturity */
    void prod_BWD_T(
        const TreeSlice & s,
        int step,
        int bot,
        int top,
        int pStart,
        int pEnd,
        const vector< TreeSliceSP > & price);

    // this is payoff boundary condition, for KO, early exercise etc.
    void prod_BWD(
        const TreeSlice & s,
        int step,
        int bot,
        int top,
        int pStart,
        int pEnd,
        const vector< TreeSliceSP > & price);
    
    /** premium scaling */
    virtual double scalePremium(const double& fairValue, 
                                      YieldCurveConstSP disc);

    /** extra output requests */    
    void recordOutput(Control* control, 
                      YieldCurveConstSP disc, 
                      Results* results); 

    void AdjustDeltaShift(CControl* control);

    virtual CVolRequestConstSP GetLNRequest() const;

    virtual void update(int& step, FDProduct::UpdateType type);

    /** ignore start date if not forward starting */
    virtual DateTime getStartDate() const
    {
        return inst->fwdStarting ? inst->startDate : inst->valueDate;
    }
    
    /** initialisation, called ONCE only before initModel() for each new model instance */
    virtual void init(Control*  control) const;

    /** initialising and setting product variables */
    // this is called per pricing call before each pricing 
    virtual void initProd();

private:
    const RangeNote* inst;
    bool            useInsertNode;
    vector<double>  stepCallValue;
    vector<bool>    stepIsMonitoring;        // whether a step is monitoring date
    vector<bool>    stepIsKO;                // or a KO monitoring date
    vector<bool>    stepIsCallable;            // whether a step is callable
    vector<bool>    stepIsLibor;            // or a  Libor Leg payment date
    vector<double>    stepLiborValue;            // what amount to be paid on Libor Leg
    vector<double>  stepIsPayDate;            // whether a step is a payment date
    vector<double>  stepIsCallPayDate;        // or a call payment date
    vector<double>  stepIsKOPayDate;        // or a KO payment date

    double VolDaysPerYear;
    bool calledOnValueDate;
};

/** create a fd payoff product */
FDProductSP RangeNote::createProduct(FDModel* model) const
{
    return FDProductSP( new RangeNoteFDProd(this, model) );
}

// adjust delta
void RangeNoteFDProd::AdjustDeltaShift(CControl* control) {
    // do nothing
}

/** initialise - allow product customisation */
void RangeNoteFDProd::init(CControl* control) const
{
    static const string method = "RangeNoteFDProd::Init";
    try 
    {
        /*if (instPtr->fwdStarting){
                throw ModelException(method," fwdStarting is currently not supported on the Callable Range Note");
        }*/

        if (!(inst->paymentStyle == "SCHEDULE"))
        {
                throw ModelException(method," Payment Style SCHEDULE must be used with Callable Range Note");
        }
        if (!(inst->spreadType == "NONE"))
        {
                throw ModelException(method," spreadType NONE must be used with Callable Range Note");
        }
        if (!(inst->payRangeType == "INSIDE"))
        {
                throw ModelException(method," Range Type INSIDE must be used with Callable Range Note");
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

            if ( inst->fwdStarting && inst->startDate>inst->valueDate) 
                tree1f->controlSameGridFwdStart(inst->ccyTreatment);
        }

        /** customize tree parameters here and set up the tree */
        DateTimeArray segDates;
        segDates.resize(2);         
 
        segDates[0] = getStartDate();
        segDates[1] = inst->getLastDate();

        IntArray density( 1, 1 );

        // add all monitoring dates to critical dates
        DateTimeArray critDates;
        int num = inst->monitorDates->size();
        for (int i=0; i<num; i++) {
            critDates.push_back((*(inst->monitorDates))[i]);
            critDates.push_back((*(inst->paymentDates))[i]);
        }            
        
        // add critical dates
        model->addCritDates( critDates );

        // prepare timeline set up
        model->initSegments( segDates, density );
    }
    catch (exception& e) 
    {
        if (!!inst->rnType) { delete inst->rnType; }
        throw ModelException(e, method);
    }
} 

/** initialise product specific data */
void RangeNoteFDProd::initProd()
{
    static const string method = "RangeNoteFDProd::InitProd";
    try {

        // determine whether any historical samples exist 
        inst->instantiateRNtype();

        int i,idx;
        int lastStep = model->getLastStep();  
  
        initSlices( numPrices );
        initInsertNode();

        inst->fwdAtStart = 1.0;
        // get spot at start if needed
        if (inst->fwdStarting && model->getDate(0)>inst->valueDate)
        {
            inst->fwdAtStart = inst->asset->fwdValue(inst->startDate);
        }

        DateTime endDate(inst->valueDate.getDate()+364, inst->valueDate.getTime()); 
        VolDaysPerYear = 252.0; // model1F->GetTimeMetric()->volDays(instPtr->valueDate, endDate);      
        
        /////////////////////////////////////////////////////////////////////////////
                /* Validations for callableRangeNote : */

        // ensure that the call and barrier schedule are always positive 
        if(!!(inst->callSchedule) && (!(inst->callSchedule->length() == 0)) && !inst->callSchedule->isNonNegative())
        { 
            throw ModelException(method,
                                 "callSchedule contains negative points.");
        }

        if(!!inst->barrierSchedule && !inst->barrierSchedule->isNonNegative())
        {
            throw ModelException(method,
                                 "barrierSchedule contains negative points.");
        }

         
        DateTime callDate ;
        DateTime barrierDate;
        DateTime monitorDate;
        DateTime paymentDate;
        
        // validate that each call date is a monitoring date
        if(!!(inst->callSchedule) && (!(inst->callSchedule->length() == 0)))
        {
            for (idx = 0; idx < inst->callSchedule->length(); idx++)
            {
                callDate = (inst->callSchedule->getDates())[idx];
                
                bool callDateFound = false;
                int i = 0;
                while ( (i < inst->monitorDates->size()) &&( !callDateFound  ) )
                {
                    monitorDate = (*inst->monitorDates)[i];

                    if (monitorDate.equals(callDate)){
                        callDateFound = true;
                    }

                    i++;
                }
                if (!callDateFound)
                {
                    throw ModelException(method,
                                         "Call Date " + Format::toString(idx+1) + " " + callDate.toString()+
                                         " must be a Monitoring Date ");
                }            
            
            }
        }

        // validate that each KO date is a monitoring date
        if(!!(inst->barrierSchedule) && (!(inst->barrierSchedule->length() == 0)))
        {
            int i = 0;
            for ( idx = 0; idx < inst->barrierSchedule->length(); idx++)
            {
                barrierDate = (inst->barrierSchedule->getDates())[idx];
                
                bool barrierDateFound = false;
                
                while ( (i < inst->monitorDates->size()) &&( !barrierDateFound  ) )
                {
                    monitorDate = (*inst->monitorDates)[i];

                    if (monitorDate.equals(barrierDate)){
                        barrierDateFound = true;
                    }
                    i++;
                }
                if (!barrierDateFound)
                {
                    throw ModelException(method,
                                         "Barrier Date " + Format::toString(idx+1) + " " + barrierDate.toString()+
                                         " must be a Monitoring Date ");
                }            
            
            }
        }

        // validate that Libor pay date is a Range Note payment date
        DateTime liborPayDate;
        if(!!inst->floater)
        {
            int i = 0;
            for ( idx = 0; idx < inst->floater->PayDates.size(); idx++)
            {
                liborPayDate = inst->floater->PayDates[idx];
                
                bool liborPayDateFound = false;
                
                while ( (i < inst->paymentDates->size()) &&( !liborPayDateFound  ) )
                {
                    paymentDate = (*inst->paymentDates)[i];

                    if (paymentDate.equals(liborPayDate)){
                        liborPayDateFound = true;
                    }

                    i++;
                }
                if (!liborPayDateFound)
                {
                    throw ModelException(method,
                                         "Libor Pay Date " + Format::toString(idx+1) + " " + liborPayDate.toString()+
                                         " must be a paymentDate ");
                }            
            
            }
        }

        // validate that Fwd Start must be One contract
        if ( inst->fwdStarting && (!inst->oneContract) )
        {
            throw ModelException(method,
                                         "Fwd Start can be used only in the One Contract case ");
        }

        // validate that product is called only if a Call Schedule exists
        if ( (!(inst->callSchedule) || (inst->callSchedule->length() == 0)) && inst->isCalled )
        {    
            throw ModelException(method,
                                         "Product cannot be called : no Call Schedule ");

        }

        //validate that SpotAtStart is not zero or negative
        if ( !Maths::isPositive(inst->initialSpot) )
        {
            throw ModelException(method,
                                         "Initial Spot value must be positive ");
        }

        // Now creating vectors to capture each step's features

        stepIsMonitoring.resize(lastStep + 1);
        stepIsKO.resize(lastStep + 1);
        stepIsCallable.resize(lastStep + 1);
        stepCallValue.resize(lastStep + 1);
        stepIsLibor.resize(lastStep + 1);
        stepLiborValue.resize(lastStep + 1);
        stepIsPayDate.resize(lastStep + 1);
        stepIsCallPayDate.resize(lastStep + 1);
        stepIsKOPayDate.resize(lastStep + 1);

        // initialize updateLiborValue flag
        inst->updateLiborValue = false;
       
        // past monitoring dates are not taken into account 
        int mnIdx = 0;
        while ( (*(inst->monitorDates))[mnIdx] < model->getDate(0) ){
            mnIdx ++;
        }

        for (i=0; i <= lastStep; i++) {
            const DateTime &treeDate = model->getDate(i);

            stepIsMonitoring[i] = false;
            stepIsCallable[i] = false;
            stepIsKO[i] = false;
            stepIsLibor[i] = false;
            stepLiborValue[i] = 0.0;
            stepIsPayDate[i] = false;
            stepIsCallPayDate[i] = false;
            stepIsKOPayDate[i] = false;


            if (mnIdx < inst->monitorDates->size()) 
            {
                if (treeDate.equals((*(inst->monitorDates))[mnIdx]))
                {
                    stepIsMonitoring[i] = true; 
                
                    if(!!(inst->callSchedule))
                    {
                        // see if date is callable (note: all callable dates are monitoring dates)
                        for (idx = 0; idx < inst->callSchedule->length(); idx++)
                        {
                            if (treeDate.equals ( inst->callSchedule->getDates()[idx]))
                            {
                                stepIsCallable[i] = true;
                                stepCallValue[i] = inst->callSchedule->interpolate(treeDate);
                            }                    
                        }
                    }

                    mnIdx++;
                }
            }

            if(!!(inst->barrierSchedule))
            {
                // see if date is KO date(note: all KO dates are monitoring dates)
                // time is not checked, which means that several successive treeDate
                // happening the same day can be KO dates
                for ( idx = 0; idx < inst->barrierSchedule->length(); idx++)
                {
                    if (treeDate.equals(inst->barrierSchedule->getDates()[idx], false))
                    {
                        stepIsKO[i] = true;
                    }
                }
            }
        }

        // past payment dates are not taken into account 
        int payIdx = 0;
        while ( (*(inst->paymentDates))[payIdx] < model->getDate(0) ){
            payIdx ++;
        }

        for (i=0; i <= lastStep; i++) {
            const DateTime &treeDate = model->getDate(i);

            stepIsLibor[i] = false;
            stepLiborValue[i] = 0.0;
            stepIsPayDate[i] = false;
            stepIsCallPayDate[i] = false;
            stepIsKOPayDate[i] = false;


            if ( payIdx < inst->paymentDates->size() && 
                treeDate.equals( (*(inst->paymentDates))[payIdx]) )
            {
                stepIsPayDate[i] = true; 
                
                if(!!(inst->callSchedule))
                {
                    // see if pay date is call pay date 
                    for (idx = 0; idx < inst->callSchedule->length(); idx++)
                     {
                        if (treeDate.equals (inst->getPayDate(inst->callSchedule->getDates()[idx]) ) )
                        {
                            stepIsCallPayDate[i] = true;
                        }                    
                    }
                }

                if(!!(inst->barrierSchedule))
                {
                    // see if pay date is KO pay date 
                    for (idx = 0; idx < inst->barrierSchedule->length(); idx++)
                     {
                        if (treeDate.equals (inst->getPayDate(inst->barrierSchedule->getDates()[idx]) ) )
                        {
                            stepIsKOPayDate[i] = true;
                        }                    
                    }
                }
                

                if(!!inst->floater)
                {
                    // see if date is LiborPayment date
                    for (idx = 0; idx < inst->floater->PayDates.size(); idx++)
                     {
                        if (treeDate.equals (inst->floater->PayDates[idx]) )
                        {
                            stepIsLibor[i] = true;
                            // now determining Libor payment size
                            stepLiborValue[i] = 
                                inst->lookUpAmount(treeDate, 
                                                   inst->floater->getCashFlowArray(
                                                   inst->valueDate,inst->discount.get()));
                        }                    
                    }
                }
                int k = payIdx;
                
                while ( (payIdx < inst->paymentDates->size())
                    && ((*(inst->paymentDates))[payIdx].equals( (*(inst->paymentDates))[k]) ))
                {
                    payIdx++;
                }
            }    
        }  
        useInsertNode = (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION);
        delete inst->rnType;
    }
    catch (exception& e) 
    {
        if (!!inst->rnType) { delete inst->rnType; }
        throw ModelException(e, method);
    }
}

/** calculate barriers and place barriers at inserted node if needed */
void RangeNoteFDProd::preCalc(int step)
{
    // set barrier, make adjustment for once a day monitoring if needed
    // to do: forward start adjustment
 
    int idx = tree1f->getSliceIndex(step);
    const DateTime& treeDate = model->getDate(step);
    
    if( !!(inst->barrierSchedule) && ( !(inst->barrierSchedule->length() == 0) )
                                  && ( treeDate <= inst->barrierSchedule->lastDate())
                                  && ( treeDate >= inst->barrierSchedule->firstDate()) 
                                  && (stepIsKO[step]) )
    {   
        inst->KOlevel =  inst->fwdAtStart*inst->barrierSchedule->interpolate(treeDate);

        // making a correction in case barrier monitoring is discrete
        if (!inst->isContMonitoringKO) {
            vector<double> vol;
            // adjust barrier if needed
            tree1f->GetStepVol(step, vol, &inst->KOlevel, 0, 0); // get vol at barrier
            Barrier::BarrierAdjustment(vol[0], inst->isUpAndOutKO, inst->KOlevel);
        }

        if (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION) 
        { //assuming Lower Barrier is index 0 and Upper uses 1 !!!
            tree1f->SetInsertNode(idx, 0, inst->KOlevel, 0); // inser barrier level
        }

    }
    else
    {
        tree1f->SetInsertNode(idx, 0, -1.0, 0); // if not a barrier date, insert dummy;
    }
}

/** product payoff method at maturity */
void RangeNoteFDProd::prod_BWD_T(const TreeSlice & spot,
                                 int step, 
                                 int bot, 
                                 int top, 
                                 int pStart, 
                                 int pEnd,
                                 const vector< TreeSliceSP > & price) 
{
    static const string method("RangeNoteFDProd::prod_BWD_T");

    double * s = spot.getValues();
    const vector< double * > & p = getValues( price );

    const DateTime& matDate = (*inst->monitorDates)[inst->monitorDates->size() - 1];
    const DateTime& payDate = inst->getPayDate(matDate);
    int j;
    
    // KO check
    // specific pv treatment: KO pay date is previous monitoring's pay date
    double pvKOSettle = inst->discount->pv(matDate, inst->getKOPayDate(matDate));

    vector<bool> KoFlag (top - bot + 1,false);
    // this flag is to know at each spot value if product is knocked-out or not
    
    if(!!(inst->barrierSchedule) && (!(inst->barrierSchedule->length() == 0)))
    {
        if ((inst->barrierSchedule->lastDate() >= matDate)&&(inst->barrierSchedule->firstDate() <= matDate))
        {
            if(stepIsKO[step])  
            {
                double spotValue;

                for (j=bot; j<=top; j++)
                {
                    spotValue = s[j]; 
                    
                    if( ( (inst->isUpAndOutKO)&&(spotValue >= inst->KOlevel) ) || ( (!inst->isUpAndOutKO)&&(spotValue <= inst->KOlevel) ) )
                    {
                        (p[0])[j] = pvKOSettle*inst->KOrebate*inst->fwdAtStart;
                        // if there is a KO, value is equal to the discounted KOrebate plus the 
                        // expected discounted value of the remaining libor payments which is 
                        // captured by (p[2])[j]

                        KoFlag[-bot + j] = true; // KoFlag is to know if we count the monitoring that happens
                        // on this KO date. This flag will be useful for example to have the monitoring date not 
                        // counted when a KO has been triggered. This information will be useful below when adding
                        // coupon range
                    }
                    else
                    {
                        (p[0])[j] = 0.0;
                    }
                }
            }
        }
    }
    
    // settlement is exactly at the payment date
    double settlementPV = inst->discount->pv(matDate,payDate );

    for (j=bot; j<=top; j++)
    {
        if (!KoFlag[-bot + j])
        {
            // vanilla payoff in case matDate is the last monitoring date        
            if (!Maths::isZero(inst->optionParticipation))
            {
                (p[0])[j] = settlementPV
                              *Maths::min( inst->initialSpot,
                                        inst->optionParticipation
                                        *GetIntrinsic(s[j],  inst->fwdAtStart*inst->finalStrike, inst->isCall, true/* this allow fwd */)); 
            
            }
            else 
            {
                (p[0])[j] = 0.0;
            }
        }
        
        (p[1])[j] = 0.0;  
        (p[2])[j] = 0.0; // both vectors are used for optimal exercise strategy if the product is callable

        // range note payoff              
        if ( (!!inst->highRangeSchedule.get())
             && (!!inst->lowRangeSchedule.get())
             && (s[j] < inst->fwdAtStart*inst->highRangeSchedule->interpolate(matDate))
             && (s[j] > inst->fwdAtStart*inst->lowRangeSchedule->interpolate(matDate)) )
        {           
            (p[0])[j] += inst->rebateSchedule->interpolate(matDate)*settlementPV;
        }
        
    }

    // libor leg payment: goes AGAINST the range leg payment
    if (!!inst->floater)
    { 
        if(inst->floater->PayDates.size() != 0)
        {
            if (payDate.equals(inst->floater->PayDates[inst->floater->PayDates.size()-1]))
            {   // adding libor leg payment to (p[0]) 
                    double liborPayment = inst->lookUpAmount(inst->getPayDate(matDate),
                                            inst->floater->getCashFlowArray(
                                                inst->valueDate,inst->discount.get()));
                    for (j=bot; j<=top; j++)
                    {    
                        (p[0])[j] -= liborPayment*settlementPV;
                    }
            }    
        }
    }

    // paying notional at maturity if needed
    if (inst->includeNotional)
    {
        if (inst->fwdStarting)
        {
            for (j=bot; j<=top; j++)
            {    
                (p[0])[j] += inst->fwdAtStart*settlementPV;
            }
        }
        else
        {        
            for (j=bot; j<=top; j++)
            {    
                (p[0])[j] += inst->initialSpot*settlementPV;
            }
        }
    }

    // checking if payment of call provision happens on last payment date    
    if(!!(inst->callSchedule) && (!(inst->callSchedule->length() == 0)))
    {
        if( payDate.equals(inst->getPayDate(inst->callSchedule->lastDate())) )
        { 
            for (j=bot; j<=top; j++)
            {         
                (p[1])[j] = (p[0])[j];                                    
            }
        }
    }
    
    // debug
#if 0
        {
            static bool doneit = false;
            if (!doneit)
            {
                FILE * file = fopen("C:\\tmp\\ payoffatmat.txt","w");
                fprintf(file, "\n ------------------ PAYOFF AT MATURITY ------------ \n");
                for (j=bot; j<=top; j++)
                {
                    fprintf(file, "j = %d \t s[j] = %f \t payoff = %f", j, s[j], (p[0])[j]);
                    fprintf(file, "\n");
                }
                fclose(file);
                doneit = true;
            }
        }
            
#endif

    inst->d_prices = DoubleArraySP(new DoubleArray(top - bot + 1));
    for (int d_i = bot; d_i <= top; d_i++)
    {
        (*inst->d_prices)[d_i-bot] = (p[0])[d_i];
    }
}
      
/** product payoff method at steps earlier than maturity */
void RangeNoteFDProd::prod_BWD(const TreeSlice & spot,
                               int step, 
                               int bot, 
                               int top, 
                               int pStart, 
                               int pEnd,
                               const vector< TreeSliceSP > & price) 
{
    static const string method = "RangeNoteFDProd::prod_BWD";
    try 
    {
        double * s = spot.getValues();
        const vector< double * > & p = getValues( price );

        int j;
        const DateTime& treeDate = model->getDate(step);    
        double df,callValue;       
        double pastValue,strictPastValue,pastSpotValue,spotValue,pvSettle,pvKOSettle;
        bool usePastSpotValue = false;
        bool isCalledToday = false;
        bool todayIsValueDate = treeDate.equals(inst->valueDate);

        vector<bool> KoFlag (top - bot + 1,false);
        // this flag is to know at each spot value if product is knocked-out or not
            
        if (stepIsMonitoring[step])
        {
            pvSettle = inst->discount->pv(treeDate, inst->getPayDate(treeDate));
        }
        
        // KO check
        if(!!(inst->barrierSchedule) && (!(inst->barrierSchedule->length() == 0)))
        {
            if ((inst->barrierSchedule->lastDate() >= treeDate)&&(inst->barrierSchedule->firstDate() <= treeDate))
            {
                if(stepIsKO[step])  
                {
                    // specific pv treatment: KO pay date is previous monitoring's pay date
                    pvKOSettle = inst->discount->pv(treeDate, inst->getKOPayDate(treeDate));

                    if (step == 0)
                    {
                        // use the spot value entered by the user if it is not zero
                        pastSpotValue = inst->getHistPrice(treeDate);
                        if (! Maths::isZero(pastSpotValue)) 
                        { 
                            usePastSpotValue = true;
                            spotValue = pastSpotValue;
                        }
                    }
                                            
                    for (j=bot; j<=top; j++)
                    { 
                        if (!usePastSpotValue)
                        {
                            spotValue = s[j]; 
                        }    
    
                        if( ( (inst->isUpAndOutKO)&&(spotValue >= inst->KOlevel) ) || ( (!inst->isUpAndOutKO)&&(spotValue <= inst->KOlevel) ) )
                        {
                            (p[0])[j] = (p[2])[j] + pvKOSettle*inst->KOrebate*inst->fwdAtStart;
                            // if there is a KO, value is equal to the discounted KOrebate plus the 
                            // expected discounted value of the remaining libor payments which is 
                            // captured by (p[2])[j]

                            KoFlag [-bot + j] = true; // KoFlag is to know if we count the monitoring that happens
                            // on this KO date. This flag will be useful for example to have the monitoring date not 
                            // counted when a KO has been triggered. This information will be useful below when adding
                            // coupon range
                        }
                    }
                }
            }
        }

        //  Call Check
        if(!!(inst->callSchedule) && (!(inst->callSchedule->length() == 0)))
        {
            // First part:looking whether treeDate is the last monitoring date of a given period            
            // If it is the case then tree values are stored in (p[1])[j] to be used for exercise decision
            // (p[1])[j] is the discounted expected value of the RangeNote at the end of 
            // the monitoring period (excluding already scheduled payments that are treated separately).
            if (stepIsMonitoring[step] && inst->isEndMonitoringPeriod(treeDate))
            {        
                for (j=bot; j<=top; j++)
                {         
                    (p[1])[j] = (p[0])[j];                                    
                }    
                
                // and update flag to stop counting Libor payments in (p[2])[j]
                inst->updateLiborValue = false;
            }

            // Second part of the algorithm: if we are at a Call notification Date, (p[1])[j] (the
            // discounted expected value of the Range Note at the next end of period) 
            // is compared to the discounted call value plus discounted expected value of Libor payments
            // to determine optimal exercise
            if(stepIsCallable[step]) 
            {    
                callValue = pvSettle*inst->callSchedule->interpolate(treeDate);

                for (j=bot; j<=top; j++)
                { 
                    if ( (callValue + (p[2])[j]< (p[1])[j]) && !KoFlag [-bot + j] )
                    // The call provision is exercised by the issuer if total called value 
                    // is less than current price and if the product has not been knocked out already
                    {
                        isCalledToday = true;
                        if (todayIsValueDate && j == 0) { calledOnValueDate = true; }
                        (p[0])[j] = callValue + (p[2])[j] + ((p[0])[j] - (p[1])[j]);
                        // when product is called, its value is the discounted callValue plus 
                        // the discounted expected value of libor payments between end of monitoring period
                        // and call payment date ( captured by (p[2])[j]) plus
                        // the discounted expected value of the remaining
                        // payments up to the end of the current monitoring period
                        // which is equal to ((p[0])[j] - (p[1])[j]).
                        // Coupons not paid yet corresponding to previously observed
                        // values are added later in the backpropagation.
                    }                    
                }
            }
        }

        // if a Call PaymentDate or a KO PaymentDate is crossed, and this value is not
        // on an end of monitoring period, updateLiborValue flag is set to true 
        // and (p[2])[j] is initialized to zero
        if ( (stepIsCallPayDate[step] || stepIsKOPayDate[step]) 
              && !(stepIsMonitoring[step] && inst->isEndMonitoringPeriod(treeDate)) )
        {
            inst->updateLiborValue = true;
            for (j=bot; j<=top; j++)
            {
                (p[2])[j] = 0.0;
            }
        }

        // Libor Leg: payment goes AGAINST range note leg
        if (!!inst->floater)
        {
            if(stepIsLibor[step] && (step >0)) 
            // if we are exactly at value date, by convention we assume the payment has already been settled
            // hence the condition step >0
            {    
                for (j=bot; j<=top; j++)
                {
                    (p[0])[j] -= stepLiborValue[step];
                    // If updateLiborValue is swtiched on, (p[2])[j] (which captures the 
                    // discounted expected value of libor streams between call of KO payment date
                    // and end of monitoring period) takes into account the libor payment
                    if ( inst->updateLiborValue)
                    {
                        (p[2])[j] -= stepLiborValue[step];
                    }
                }
            }
        }

        if ((stepIsMonitoring[step]) && (step > 0)) 
        // Adding Range coupon payment. The condition (step > 0) is enforced
        // because a special treatement is performed below for the case step == 0
        {    
            for (j=bot; j<=top; j++)
            {     
                // in case product is not KO, one adds range note value
                if ( !(stepIsKO[step] && KoFlag [-bot + j]) )    
                {
                    if ( (!!inst->highRangeSchedule.get())
                     && (!!inst->lowRangeSchedule.get())
                     && (s[j] < inst->fwdAtStart*inst->highRangeSchedule->interpolate(treeDate))
                     && (s[j] > inst->fwdAtStart*inst->lowRangeSchedule->interpolate(treeDate)) )
                    {                            
                        (p[0])[j] += pvSettle*inst->rebateSchedule->interpolate(treeDate);        
                    }
                }                
            }
        }
        
        // Adding past monitoring value not paid yet, including today's if it is a monitoring value
        // the special case where the product has been called (isCalled true) is treated here
        if (step == 0)
        {
            if (inst->isCalled && !isCalledToday)
            // If this is the case product will be called at the end of the current monitoring period
            // hence the value in the product is equal to the discounted expected value of future coupons
            // plus the call strike payment plus past coupons not paid yet plus remaining libor payments
            // The flag calledToday is to avoid repeating the same operation twice
            {
                // compute the discounted value of the strike price
                df = inst->discount->pv(treeDate, inst->getPayDate(inst->calledDate));        
                callValue = df*inst->callSchedule->interpolate(inst->calledDate);

                for (j=bot; j<=top; j++)
                {     
                    (p[0])[j] = ((p[0])[j]-(p[1])[j]) + callValue + (p[2])[j];            
                }                
            }
            
            // Adding historical coupons not paid yet. Today's coupon
            // is not paid only if there is a KO today and product is not Called
            strictPastValue = inst->getRemainingStrictHistPayments(treeDate,treeDate);
            pastValue = inst->getRemainingHistPayments(treeDate,treeDate);
            
            for (j=bot; j<=top; j++)
            { 
                if (stepIsKO[step] && KoFlag [-bot + j] && !inst->isCalled)
                {
                    (p[0])[j] += strictPastValue;
                }
                else
                {
                    (p[0])[j] += pastValue;
                }
            }            
        }      
    }

    catch (exception& e) {
        throw ModelException(e, method);
    }

    // debug
#if 0
        {
            static bool doneit = false;
            if (!doneit)
            {
                FILE * file = fopen("C:\\tmp\\rn.txt","w");
                fprintf(file, "\n ------------------ STEP IS KO ------------ \n");
                for (int i = 0;i < stepIsKO.size(); i++)
                {
                    DateTime::MonthDayYear mdy = model->getDate(i).toMDY();
                     int ansiDate = (mdy.day + 100 * mdy.month + 10000 * mdy.year);

                    fprintf(file, "i = %d \t stepIsKO = %d \t Date = %d", i, stepIsKO[i]? 1 : 0, 
                        ansiDate);
                    fprintf(file, "\n");
                }
                fclose(file);
                doneit = true;
            }
        }
            
#endif

}

void RangeNoteFDProd::recordOutput(Control* control, 
                                   YieldCurveConstSP disc, 
                                   Results* results)                                    
{
     // get prices at t=0
    // save price
    double price = scalePremium(model->getPrice0( *slices[0] ), disc);
    results->storePrice(price, disc->getCcy());
    
    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        DateTime       matDate = inst->getLastDate(); // to do: verify ok
        double         indVol;
        // calculate indicative vol
        try{
            if ( matDate.isGreater(inst->valueDate) )
            {

                DateTime imntStartDate = inst->fwdStarting? inst->startDate: inst->valueDate;

                // get vol request
                CVolRequestConstSP lnVolRequest = GetLNRequest();

                // interpolate the vol
                CVolProcessedSP  vol(inst->asset->getProcessedVol(lnVolRequest.get()));
                // cast to the type of vol we're expecting
                CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);
                // this should never happen if our get market data has worked properly
                if (!vol){
                    throw ModelException("VanillaTreeFDProd::recordOutput", 
                                         "No Black Scholes Vol");
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
                            indVol,
                            999.000); // to do: what do we need here?


        OutputNameConstSP lastStep(new OutputName("LastStepDebug"));
        results->storeGreek(inst->d_prices, Results::DEBUG_PACKET, lastStep);

        OutputNameConstSP exerciseTodayOut(new OutputName("exerciseToday"));
        results->storeGreek(IObjectSP(CDouble::create(calledOnValueDate? 1.0 : 0.0)), Results::DEBUG_PACKET, exerciseTodayOut);

    }
}

/** returns a vol request for log-normal vol */
CVolRequestConstSP RangeNoteFDProd::GetLNRequest() const
{
    // get strike and maturity date from instrument
    DateTime matDate = inst->getLastDate();

    double volStrike  = inst->finalStrike;

    DateTime imntStartDate = inst->fwdStarting? inst->startDate: inst->valueDate;

    CVolRequestConstSP volRequest(
        new LinearStrikeVolRequest(volStrike, imntStartDate, 
                                   matDate, inst->fwdStarting));
    return volRequest;
}

/** premium scaling */
double RangeNoteFDProd::scalePremium(const double& fairValue, 
                                     YieldCurveConstSP disc)
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

    return (fairValue*scalingFactor*fwdStartDF);
}

void RangeNoteFDProd::update(int& step, 
                             FDProduct::UpdateType type)
{
    // we assume just need one und level for spot here
    const TreeSlice & s = payoffIndex->getValue( step );
    int bot, top;
    s.getCalcRange( bot, top );

    const vector< TreeSliceSP > & price = slices;
    int pStart = 0, pEnd = price.size() - 1;

    if (type == FDProduct::BWD_T)
    {
        prod_BWD_T(s,
                   step,
                   bot,
                   top,
                   pStart, 
                   pEnd,
                   price);

        //insert nodes
        if (tree1f && tree1f->NumOfInsertNode>0)
        {
            prod_BWD_T(*insNodes,
                        step,
                        0,
                        tree1f->NumOfInsertNode-1,
                        pStart, 
                        pEnd,
                       *insPrices);
        }
    }
    else if(type == FDProduct::BWD)
    {
        prod_BWD(s,
                 step,
                 bot,
                 top,
                 pStart, 
                  pEnd,
                 price);

        //insert nodes
        if (tree1f && tree1f->NumOfInsertNode>0)
        {
            prod_BWD(*insNodes,
                      step,
                      0,
                       tree1f->NumOfInsertNode-1,
                      pStart, 
                      pEnd,
                     *insPrices);
        }
    }   
}

/* for class loading */
bool RangeNoteLoad() {
    return (RangeNote::TYPE != 0);
}

DRLIB_END_NAMESPACE






