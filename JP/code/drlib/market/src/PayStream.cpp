/*****************************************************************************
 *
 *    Group       : Equity Derivatives Research
 *
 *    Description : Defines a Payment stream
 *
 *    Author      : Stephen Hope
 *
 *    Date        : 29 May 2001
 *
 *
 *******************************************************************************/

#include "edginc/config.hpp"
#include "edginc/PayStream.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Addin.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/BadDayFollowing.hpp"
#include "edginc/Business252.hpp"
#include "edginc/NonPricingModel.hpp"

DRLIB_BEGIN_NAMESPACE

// wrapper class that can be used in an IMS generic template
class PayStreamLite: public CObject,
                     virtual public ITypeConvert {
public:
    static CClassConstSP const TYPE;
    friend class PayStreamLiteHelper;

    /** create a proper PayStream */
    virtual void convert(IObjectSP&    object,
                         CClassConstSP requiredType) const {
        static const string method = "PayStreamLite::convert";
        try {
            if (requiredType != PayStream::TYPE) {
                throw ModelException(method, 
                                     "Cannot convert a PayStreamLite into "
                                     "object of type "+requiredType->getName());               
            }
            DayCountConventionSP dayCount(DayCountConventionFactory::make(dcc));
            object = IObjectSP(new PayStream(dayCount.get(),
                                             compoundFlat,
                                             floatRate.get(),
                                             *(notionals),
                                             *(refixDates),
                                             *(accrualDates),
                                             *(payDates),
                                             *(fixings),
                                             *(spreads),
                                             *(weights),
                                             *(compounding)));
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    
private:
    PayStreamLite():CObject(TYPE){}
    
    PayStreamLite(const PayStreamLite& rhs);
    PayStreamLite& operator=(const PayStreamLite& rhs);   
   
    DateTimeArraySP refixDates;
    DateTimeArraySP accrualDates;
    DateTimeArraySP payDates;
    DoubleArraySP   fixings;
    DoubleArraySP   spreads;
    DoubleArraySP   weights;
    DoubleArraySP   notionals;
    CBoolArraySP    compounding;
    FloatRateSP     floatRate;
    string          dcc;
    bool            compoundFlat;
};

class PayStreamLiteHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PayStreamLite, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ITypeConvert);
        EMPTY_SHELL_METHOD(defaultPayStreamLite);
        FIELD(refixDates, "refix dates");
        FIELD(accrualDates, "accrual dates");
        FIELD(payDates, "pay dates"); 
        FIELD(fixings, "fixings");
        FIELD(spreads, "spreads");
        FIELD(weights, "weights");
        FIELD(notionals, "notionals");
        FIELD(compounding, "compounding");
        FIELD(floatRate, "float rate");
        FIELD(dcc, "dcc");
        FIELD(compoundFlat, "compound flat");
        FIELD_MAKE_OPTIONAL(floatRate);
    }

    static IObject* defaultPayStreamLite(){
        return new PayStreamLite();
    }
};

CClassConstSP const PayStreamLite::TYPE = CClass::registerClassLoadMethod(
    "PayStreamLite", typeid(PayStreamLite), PayStreamLiteHelper::load);

class PayStreamXLInterface: public CObject,
                            public IXLInterfaceMap{
public:
    static CClassConstSP const TYPE;
    friend class PayStream;  
    friend class PayStreamXLInterfaceHelper;

    /** Map the object */
    virtual IObjectSP map()const{
        static const string method = "PayStreamXLInterface::map";
        try
        {
            PayStreamSP payStream(new PayStream(this->payStreamDCC.get(),
                                                this->compoundFlat,
                                                this->floatRate.get(),
                                                *(this->notionals),
                                                *(this->refixDates),
                                                *(this->accrualDates),
                                                *(this->payDates),
                                                *(this->fixings),
                                                *(this->spreads),
                                                *(this->weights),
                                                *(this->compounding)));
            
            return payStream;
            
        }
        catch (exception& e) 
        {
            throw ModelException(e, method);
        }
    }

    
private:
    PayStreamXLInterface():CObject(TYPE){
        //empty
    }
    
    PayStreamXLInterface(const PayStreamXLInterface &rhs);
    PayStreamXLInterface& operator=(const PayStreamXLInterface& rhs);
    
    /** overrides default */
    void validatePop2Object(){
        static const string method = "PayStreamXLInterface::validatePop2Object";

        try
        { 
            int numItems = payDates->size();
            if (numItems != refixDates->size() ||
                numItems != fixings->size() ||
                numItems != spreads->size() ||
                numItems != weights->size() ||
                numItems != compounding->size() ||
                numItems != notionals->size() ||
                numItems != (accrualDates->size()-1))
            {
                throw ModelException(method,
                                     "The number of pay dates (size = " + Format::toString(payDates->size()) +
                                     "), refix dates (size =" + Format::toString(refixDates->size()) +
                                     "), fixings (size =" + Format::toString(fixings->size()) + 
                                     "), \nspreads (size =" + Format::toString(spreads->size()) + 
                                     "), weights (size =" + Format::toString(weights->size()) + 
                                     "), compounding (size =" + Format::toString(compounding->size()) + 
                                     "), notionals (size =" + Format::toString(notionals->size()) + 
                                     ") must be equal.  \nThe number of accrualDates (size =" + Format::toString(accrualDates->size()) +
                                     ") must be one greater than this. \n");
            }
            
            
            if (numItems)
            {
                if ((*compounding)[compounding->size()-1])
                {
                    throw ModelException(method,
                                         "The last period cannot be compounded.");
                }
            }
        }
        catch (exception& e) 
        {
            throw ModelException(e, method);
        }
    }
    
            

    PayStreamXLInterface(const DayCountConvention* payStreamDCC,
                         bool  compoundFlat,
                         const FloatRate* floatRate,
                         const DoubleArray* notionals,
                         const DateTimeArray* refixDates,
                         const DateTimeArray* accrualDates,
                         const DateTimeArray* payDates,
                         const DoubleArray*   fixings,
                         const DoubleArray*   spreads,
                         const DoubleArray*   weights,
                         const BoolArray*     compounding):
        CObject(TYPE), compoundFlat(compoundFlat){
        
        static const string method = "PayStreamXLInterface::PayStreamXLInterface";
        
        try
        {
            this->payStreamDCC = DayCountConventionSP(copy(payStreamDCC));
            if (floatRate)
            {
                this->floatRate = FloatRateSP(copy(floatRate));
            }
            else
            {
                this->floatRate = FloatRateSP(   );
            }
            this->notionals = DoubleArraySP(copy(notionals));
            this->refixDates = DateTimeArraySP(copy(refixDates));
            this->accrualDates = DateTimeArraySP(copy(accrualDates));
            this->payDates = DateTimeArraySP(copy(payDates));
            this->fixings = DoubleArraySP(copy(fixings));
            this->spreads = DoubleArraySP(copy(spreads));
            this->weights = DoubleArraySP(copy(weights));
            this->compounding = CBoolArraySP(copy(compounding));
            
            this->validatePop2Object();
            
        }
        catch (exception& e) 
        {
            throw ModelException(&e, method);
        }
    }


    DateTimeArraySP      refixDates;
    DateTimeArraySP      accrualDates;
    DateTimeArraySP      payDates;
    DoubleArraySP        fixings;
    DoubleArraySP        spreads;
    DoubleArraySP        weights;
    DoubleArraySP        notionals;
    CBoolArraySP         compounding;
    FloatRateSP          floatRate;
    DayCountConventionSP payStreamDCC;
    bool                 compoundFlat;
};

class PayStreamXLInterfaceHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PayStreamXLInterface, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IXLInterfaceMap);
        EMPTY_SHELL_METHOD(defaultPayStreamXLInterface);
        FIELD(refixDates, "refix dates");
        FIELD(accrualDates, "accrual dates");
        FIELD(payDates, "pay dates"); 
        FIELD(fixings, "fixings");
        FIELD(spreads, "spreads");
        FIELD(weights, "weights");
        FIELD(notionals, "notionals");
        FIELD(compounding, "compounding");
        FIELD(floatRate, "float rate");
        FIELD(payStreamDCC, "pay stream DCC");
        FIELD(compoundFlat, "compound flat");
        FIELD_MAKE_OPTIONAL(floatRate);
    }

    static IObject* defaultPayStreamXLInterface(){
        return new PayStreamXLInterface();
    }
};

CClassConstSP const PayStreamXLInterface::TYPE = CClass::registerClassLoadMethod(
    "PayStreamXLInterface", typeid(PayStreamXLInterface), PayStreamXLInterfaceHelper::load);

typedef smartConstPtr<PayStreamXLInterface> PayStreamXLInterfaceConstSP;
typedef smartPtr<PayStreamXLInterface> PayStreamXLInterfaceSP;



/* for reflection */
PayStream::PayStream(): 
    CObject(TYPE), accruedInterest(0.0)
{
    // empty
}


/** Return last accrue start date */
const DateTime& PayStream::getLastAccrueDate()const
{
    PayStream::PaymentSP lastPayment = (*payments)[payments->size()-1];
    PayStream::Payment::RefixSP lastRefix = (*(lastPayment->refixes))[lastPayment->refixes->size()-1];
            
    return lastRefix->accrueEnd;
}


/** Curtail the PayStream for callable LiborStreams */
void PayStream::curtail(const DateTime& callDatePlusOffset)
{
    static const string method = "PayStream::curtail";

    bool chop = false;
    int pi, ri;


    // loop backwards through the payments to find where the callDate lies
    for (pi = payments->size()-1; pi >= 0; pi--)
    {
        PayStream::PaymentSP currPayment = (*payments)[pi];
        // loop backwards through this payment's refix objects
        for (ri = currPayment->refixes->size()-1; ri >= 0; ri--)
        {
            PayStream::Payment::RefixSP currRefix = (*(currPayment->refixes))[ri];
            if (callDatePlusOffset.isGreater(currRefix->accrueStart) &&
                currRefix->accrueEnd.isGreaterOrEqual(callDatePlusOffset))
            {
                chop = true;
                currPayment->refixes->resize(ri + 1);
                currPayment->payDate = callDatePlusOffset;
                break;
            }
        }
        if (chop)
        {
            break;
        }
    }
    
    if (chop)
    {
        payments->resize(pi + 1);
        
        // now look back through earlier payments to modify the payDate if necessary
        for (pi = payments->size()-1; pi >= 0; pi--)
        {
            PayStream::PaymentSP currPayment = (*payments)[pi];
            if (currPayment->payDate.isGreater(callDatePlusOffset))
            {
                currPayment->payDate = callDatePlusOffset;
            }
        }
    }

}

void PayStream::getMarketData(const IModel* model, const CMarketDataSP market)
{
    floatRate->getMarketData(model, market);
}

/** Constructor with flat list type interface */
PayStream::PayStream(const DayCountConvention* payDayCountConv,
                     bool  compoundFlat,
                     const FloatRate* floatRate,
                     const DoubleArray& notionals,
                     const DateTimeArray& refixDates,
                     const DateTimeArray& accrualDates,
                     const DateTimeArray& payDates,
                     const DoubleArray&   refixLevels,
                     const DoubleArray&   spreads,
                     const DoubleArray&   weights,
                     const BoolArray&     compound)   // TRUE = compound period
    : CObject(TYPE)
{
    static const string method = "PayStream::PayStream";
    int idx = 0;
    accruedInterest = 0.0;

    try
    {
        int numItems = payDates.size();
        if (numItems != refixDates.size() ||
            numItems != refixLevels.size() ||
            numItems != spreads.size() ||
            numItems != weights.size() ||
            numItems != compound.size() ||
            numItems != notionals.size() ||
            numItems != (accrualDates.size()-1))
        {
            throw ModelException(method,
                                 "The number of pay dates (size = " + Format::toString(payDates.size()) +
                                 "), refix dates (size =" + Format::toString(refixDates.size()) +
                                 "), refixLevels (size =" + Format::toString(refixLevels.size()) + 
                                 "), \nspreads (size =" + Format::toString(spreads.size()) + 
                                 "), weights (size =" + Format::toString(weights.size()) + 
                                 "), compound (size =" + Format::toString(compound.size()) + 
                                 "), notionals (size =" + Format::toString(notionals.size()) + 
                                 ") must be equal.  \nThe number of accrualDates (size =" + Format::toString(accrualDates.size()) +
                                 ") must be one greater than this. \n");
        }

        
        if (numItems)
        {
            if (compound[compound.size()-1])
            {
                throw ModelException(method,
                                     "The last period cannot be compounded.");
            }
        }
        
        PaymentArraySP payArray(new PaymentArray);
        
        idx = 0;
        while (idx < numItems)
        {
            double tempNotional = notionals[idx];
            int startIdx = idx;
            do
            {
                if (!(Maths::equals(tempNotional, notionals[idx])))
                {
                    throw ModelException(method,
                                         "notionals must be equal across "
                                         "compounding periods.");
                }
                idx++;
                
            }while (compound[idx-1]);                   
            
            // construct this payment
            PaymentSP payment(new Payment(notionals[startIdx],
                                          refixDates,
                                          accrualDates,
                                          payDates[idx-1],
                                          refixLevels,
                                          spreads,
                                          weights,
                                          startIdx,
                                          idx-1));
            
            payArray->push_back(payment);
        }
        
        // construct the payStream
        this->payments = PaymentArraySP(payArray.release());
        if (floatRate)
        {
            this->floatRate = FloatRateSP(copy(floatRate));
        }
        this->dayCountConv = DayCountConventionSP(copy(payDayCountConv));
        this->compoundFlat = compoundFlat;

        this->validatePop2Object();
    }
    catch (exception& e) 
    {
        throw ModelException(&e, method);
    }
}

PayStream::Payment::Payment(double notional,
                            const DateTimeArray& refixDates,
                            const DateTimeArray& accrualDates,
                            const DateTime& payDate,
                            const DoubleArray& refixLevels,
                            const DoubleArray& spreads,
                            const DoubleArray& weights,
                            int   startIdx,
                            int   endIdx)
    : CObject(TYPE)
{
    static const string method = "PayStream::Payment::Payment"; 
    
    try
    {
        Payment::RefixArraySP refixArray(new Payment::RefixArray);

        // if payment has single accrual period
        if (startIdx == endIdx)
        {
            Payment::RefixSP refix(
                new Payment::Refix(refixDates[startIdx],
                                   accrualDates[startIdx],
                                   accrualDates[startIdx+1],
                                   refixLevels[startIdx],
                                   spreads[startIdx],
                                   weights[startIdx]));
            
            refixArray->push_back(refix);
        }
        else  // compounding within the payment
        {
            for (int i = startIdx; i < endIdx+1; i++)
            {
                Payment::RefixSP refix(
                    new Payment::Refix(refixDates[i],
                                       accrualDates[i],
                                       accrualDates[i+1],
                                       refixLevels[i],
                                       spreads[i],
                                       weights[i]));
                    
                refixArray->push_back(refix);
            }
        }

        // construct payment
        this->payDate = payDate;
        this->refixes = Payment::RefixArraySP(refixArray.release());
        this->notional = notional;
    }
    catch (exception& e) 
    {
        throw ModelException(&e, method);
    }
}

/* for reflection */
PayStream::Payment::Payment(): 
    CObject(TYPE)
{
    // empty
}

/* for reflection */
PayStream::Payment::Refix::Refix(): 
    CObject(TYPE)
{
    // empty
}

PayStream::Payment::Refix::Refix(const DateTime& refixDate,
                                 const DateTime& accrueStart,
                                 const DateTime& accrueEnd,
                                 double refixLevel,
                                 double spread,
                                 double weight)
    : CObject(TYPE), refixDate(refixDate), accrueStart(accrueStart), 
    accrueEnd(accrueEnd), refixLevel(refixLevel), spread(spread), weight(weight)
{
    // empty
}


bool PayStream::load() 
{
    return (PayStream::TYPE && 
            PayStream::Payment::TYPE && 
            PayStream::Payment::Refix::TYPE);
}

/** return TRUE if floatingRate */
bool PayStream::isFloating()const
{
    bool result = false;

    if (floatRate.get())
    {
        result = true;
    }

    return result;
}

/** when to stop tweaking */
DateTime PayStream::lastSensDate()const
{
    if (!payments->empty()) {
        DateTime lastPayDate = (*payments)[payments->size()-1]->payDate;
        DateTime lastSensDate = lastPayDate;

        if (isFloating()) {
            lastSensDate = floatRate->findMatDate(lastPayDate);
        }

        return lastSensDate;
    }
    else {
        // no payments, so no point tweaking
        DateTime lastSensDate;
        return lastSensDate;
    }        
}

/** when do payments occur ? */
DateTimeArraySP PayStream::paymentDates() const {
    DateTimeArraySP paydates(new DateTimeArray(payments->size()));

    for (int i = 0; i < paydates->size(); i++) {
        (*paydates)[i] = (*payments)[i]->payDate;
    }
    return paydates;
}

/** and what are they ? */
CashFlowArrayConstSP PayStream::knownCashflows(
    const DateTime&    today,
    const DoubleArray* notionals,
    bool               isCallable,
    const DateTime&    cutoff,
    double             callFwdRate) {
    if (!isFloating()) {
        if (!isCallable) {
            return cashFlowList(today, notionals);
        }
        return cashFlowList(today, cutoff, callFwdRate, notionals);
    }
    else {
        CashFlowArraySP cfl(new CashFlowArray(0));
        double accInt;
        for (int i = 0; i < payments->size(); i++) {
            if ((*payments)[i]->isKnown(today, true)) {
                if (!isCallable) {
                    cfl->push_back((*payments)[i]->cashFlow(dayCountConv.get(),
                                                            compoundFlat,
                                                            floatRate.get(),
                                                            today,
                                                            &accInt,
                                                            notionals ? &((*notionals)[i]) : 0));
                }
                else {
                    cfl->push_back((*payments)[i]->cashFlow(dayCountConv.get(),
                                                            compoundFlat,
                                                            floatRate.get(),
                                                            today,
                                                            &accInt,
                                                            notionals ? &((*notionals)[i]) : 0,
                                                            cutoff,
                                                            callFwdRate));
                }
            }
        }
        return cfl;
    }
}

/** returns cashflows corresponding to periods whose accrual start date is before cutoff */
CashFlowArrayConstSP PayStream::startedKnownCashflows (const DateTime& today, const DoubleArray* notionals) const {

    CashFlowArraySP cfl(new CashFlowArray(0));
    double accInt;
    for (int i = 0; i < payments->size(); i++) {
        if ((*payments)[i]->isKnownAndStartedAccruing(today, isFloating())) {
                cfl->push_back((*payments)[i]->cashFlow(dayCountConv.get(),
                    compoundFlat,
                    floatRate.get(),
                    today,
                    &accInt,
                    notionals ? &((*notionals)[i]) : 0));
        }
    }

    return cfl;

}

/** Map the object */
IObjectSP PayStream::map()const
{
    static const string method = "PayStream::map";

    try
    {
        // create the required arrays
        DateTimeArraySP refixes(new DateTimeArray(0));
        DateTimeArraySP accruals(new DateTimeArray(0));
        DateTimeArraySP pays(new DateTimeArray(0));
        DoubleArraySP fixes(new DoubleArray(0));
        DoubleArraySP sprds(new DoubleArray(0));
        DoubleArraySP wts(new DoubleArray(0));
        DoubleArraySP ntns(new DoubleArray(0));
        CBoolArraySP cpnd(new BoolArray(0));
        
        int pi, ri;
        bool compoundPeriod = true;

        // loop through the payments 
        for (pi = 0; pi < payments->size(); pi++)
        {
            PayStream::PaymentSP currPayment = (*payments)[pi];
            DateTime currPayDate = currPayment->payDate;
            double   currNotional = currPayment->notional;

            // loop through this payments refixes adding to the above arrays
            for (ri = 0; ri < currPayment->refixes->size(); ri++)
            {
                PayStream::Payment::RefixSP currRefix = (*(currPayment->refixes))[ri];
                
                refixes->push_back(currRefix->refixDate);
                fixes->push_back(currRefix->refixLevel);
                sprds->push_back(currRefix->spread);
                wts->push_back(currRefix->weight);
                pays->push_back(currPayDate);
                ntns->push_back(currNotional);

                compoundPeriod = (ri == currPayment->refixes->size()-1)? false : true; 
                cpnd->push_back(compoundPeriod);


                // only need to add accrueStart for the first refix period of the first payment
                if (ri == 0 && pi == 0)
                {
                    accruals->push_back(currRefix->accrueStart);
                }

                accruals->push_back(currRefix->accrueEnd);
            }
        }
            // now call PayStreamXLInterface constructor
        IObjectSP payStreamXLInterface(
            new PayStreamXLInterface(this->dayCountConv.get(),
                                     this->compoundFlat,
                                     this->floatRate.get()? this->floatRate.get(): 0,
                                     ntns.get(),
                                     refixes.get(),
                                     accruals.get(),
                                     pays.get(),
                                     fixes.get(),
                                     sprds.get(),
                                     wts.get(),
                                     cpnd.get()));
        
        return payStreamXLInterface;
        
    }
    catch(exception& e)
    {
        throw ModelException(&e,
                             method);
    }   
}

/** overrides default */
void PayStream::validatePop2Object()
{
    static const string method = "PayStream::validatePop2Object";

    try
    {  
        // validate the float rate if it's there
        if (floatRate.get())
        {
            floatRate->validatePop2Object();
        }

        bool isFloating = this->isFloating();
        // validate each payment in the list
        for (int idx = 0; idx < payments->size(); idx++)
        {
            bool isAdjust = isFloating ? floatRate->isinArrears() : false;
            (*payments)[idx]->validate(isFloating, isAdjust);

            /* If multiple payments, check that the last accrueEnd of one 
               payment is equal to the first accrue start of the next payment */
            if (idx)
            {
                Payment::RefixArray refixes1;
                Payment::RefixArray refixes2;
                int        sizeOfRefixes1 = 0;

                refixes1 = (*(*payments)[idx-1]->refixes);
                refixes2 = (*(*payments)[idx]->refixes);
                sizeOfRefixes1 = refixes1.size();
                
                if (!(refixes1[sizeOfRefixes1 - 1]->accrueEnd.equals(
                    refixes2[0]->accrueStart)))
                {
                    throw ModelException(method,
                                         "The first accrue start date of a payment " + 
                                         refixes2[0]->accrueStart.toString() +
                                         " must equal the last accrue end date of "
                                         "the previous payment " +
                                         refixes1[sizeOfRefixes1 - 1]->accrueEnd.toString());
                }
            }
        }
    }
    catch(exception& e)
    {
        throw ModelException(&e,
                             method);
    }
}

void PayStream::Payment::validate(bool isFloating, const bool isAdjust)const
{
    static const string method = "PayStream::Payment::validate";
    
    try
    {
        for (int idx = 0; idx < refixes->size(); idx++)
        {
            // validate each Refixing
            (*refixes)[idx]->validate(isFloating, isAdjust);
            if (idx)
            {
                Payment::RefixSP refix1 = (*refixes)[idx-1];
                Payment::RefixSP refix2 = (*refixes)[idx];

                if (!(refix2->accrueStart.equals(refix1->accrueEnd)))
                {
                    throw ModelException(method,
                                         "The accrue start date " +
                                         refix2->accrueStart.toString() +
                                         " must match the previous accrue "
                                         "end date " +
                                         refix1->accrueEnd.toString());
                }
                
                if (refix1->refixDate.isGreater(refix2->refixDate))
                {
                    throw ModelException(method,
                                         "Refix dates must be in order. " 
                                         "Refix date " +
                                         refix2->refixDate.toString() +
                                         " is before refix date " +
                                         refix1->refixDate.toString());
                                         
                }
            }
        }
    }
    catch (exception& e)
    {
        throw ModelException(&e,
                             method);
    }
}


/** overrides default */
void PayStream::Payment::validatePop2Object()
{
    // empty
}

void PayStream::Payment::Refix::validate(bool isFloating, const bool isAdjust)const
{
    static const string method = "Refix::validate";
    
    if (accrueStart.isGreater(accrueEnd))
    {
        throw ModelException(method,
                             "accrue start " +
                             accrueStart.toString() +
                             " is after accrue end " +
                             accrueEnd.toString() +
                             " within the same payment.");
    }

    if (refixDate.isGreater(accrueStart) && isFloating)
    {
        // Can't be in arrears if they are on the same day
        if (!(refixDate.equals(accrueStart, false))
           && !isAdjust)
        {
            throw ModelException(method,
                                 "FloatRate requires to apply convexity/delay adjustment for in arrears! \n" 
                                 " The refix date ("+
                                 refixDate.toString()+
                                 ") is after the accrue start date ("+
                                 accrueStart.toString()+
                                 "). \n Libor streams that need to price in arrears "
                                 " must be priced with FloatRate with convex/delay adjust.");
        }
    }
} 

/** overrides default */
void PayStream::Payment::Refix::validatePop2Object()
{
    // empty
} 

/** Calculate the fwd rate between the callDatePlusOffset and
    the accrueEnd date of the relevant refix period. Note that
    the discount curve is passed since this may be a fixed PayStream
    i.e FloatRate is empty ! */
const double PayStream::callPlusOffsetFwdRate(const DateTime& callDatePlusOffset,
                                              const YieldCurveWrapper& curve)const
{
    static const string method = "PayStream::callFwdRate";
    double callFwdRate = 0.0;
    try
    {
        DateTime lastAccrueEnd = getLastAccrueDate();
		
		if(!curve.isEmpty()) {
			callFwdRate = curve->fwd(callDatePlusOffset,lastAccrueEnd,
										dayCountConv.get(),CompoundBasis::SIMPLE);
		}
		else {
			callFwdRate = floatRate->fwd(callDatePlusOffset,lastAccrueEnd, dayCountConv.get());
		}
        return callFwdRate;
    }                           
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}

/** return the ACCRUED_INTEREST. This is always calculated 
    and stored as a non-registered field in the PayStream object */
double PayStream::getAccruedInterest()const
{
    return accruedInterest;
}

// get the accrued as of 'when'
double PayStream::getAccruedInterest(const DateTime&    today,
                                     const DateTime&    when,
                                     const DoubleArray* notionals) const {
    static const string method("PayStream::getAccruedInterest");
    try {
        int             i;
        double          accrued = 0.0;
        double          accInt  = 0.0;          
        double          historicAccruedInterest = 0.0;
        CashFlowArraySP refixList;

        for (i = 0; i < payments->size(); i++) {
            accInt = 0.0;

            /* if it is a floating rate pay stream then calculate
               the forward rates for dates after today */
            if (floatRate.get()) {
                refixList = CashFlowArraySP((*payments)[i]->getRefixes());
                floatRate->calculateRates(today, *refixList, *paymentDates());
            }
            
            DoubleArray fixings((*payments)[i]->refixes->size());
            
            for (int j = 0; j < (*payments)[i]->refixes->size(); j++) {
                const Payment::Refix* refix = (*(*payments)[i]->refixes)[j].get();

                fixings[j] = floatRate.get() ? (*refixList)[j].amount : 
                    refix->refixLevel;
            }

            CashFlow cf = (*payments)[i]->cashFlow(dayCountConv.get(),
                                                   compoundFlat,
                                                   &fixings.front(),
                                                   when,
                                                   &accInt,
                                                   notionals ? &((*notionals)[i]) : 0);

            PayStream::PaymentSP thisPayment = (*payments)[i];
            if (thisPayment->payDate.isGreater(when) ) {
                /* catch case where 'when' is between last accrueEnd of a payment
                   and the payDate. */
                int numRefixesThisPayment = thisPayment->refixes->size();
                PayStream::Payment::RefixSP lastRefixThisPayment = 
                    (*(thisPayment->refixes))[numRefixesThisPayment-1];
                if (when.isGreater(lastRefixThisPayment->accrueEnd))
                {
                    historicAccruedInterest += cf.amount;
                }
                accrued += accInt;
            }            
        }   
        accrued += historicAccruedInterest;
        return accrued;
    }                    
    catch (exception& e) {
        throw ModelException(e, method);
    }
}



/** Calculate the cash flow list arising from a payment stream.
    The number of notionals must match the number of payments */
CashFlowArrayConstSP PayStream::cashFlowList(const DateTime& today,
                                             const DoubleArray* notionals)
{
    static const string method = "PayStream::cashFlowList";
    double accInt = 0.0;          // ACCRUED_INTEREST  (up to today) 
    double historicAccruedInterest = 0.0;
    
    // This is the non-registered field in the PayStream
    accruedInterest = 0.0;

    if (notionals && (notionals->size() != payments->size()))
    {
        throw ModelException(method,
                             "Number of supplied notionals is different "
                             "to the number of payments.");
    }

    CashFlowArraySP cfl(new CashFlowArray(payments->size()));

    // Just loop through each payment 
    for (int idx = 0; idx < payments->size(); idx++)
    {
        accInt = 0.0;

        (*cfl)[idx] = (*payments)[idx]->cashFlow(dayCountConv.get(),
                                                 compoundFlat,
                                                 floatRate.get(),
                                                 today,
                                                 &accInt,
                                                 notionals ? &((*notionals)[idx]) : 0);

        PayStream::PaymentSP thisPayment = (*payments)[idx];
        if ( thisPayment->payDate.isGreater(today) ) {
            /* catch case where today is between last accrueEnd of a payment
               and the payDate. */
            int numRefixesThisPayment = thisPayment->refixes->size();
            PayStream::Payment::RefixSP lastRefixThisPayment = 
                (*(thisPayment->refixes))[numRefixesThisPayment-1];
            if (today.isGreater(lastRefixThisPayment->accrueEnd))
            {
                historicAccruedInterest += (*cfl)[idx].amount;
            }
            accruedInterest += accInt;
        }
    }

    accruedInterest += historicAccruedInterest;
    
    return cfl;
}

/** Calculate the cash flow list arising from a payment stream.
    The number of notionals must match the number of payments. */
void PayStream::cashFlowList(const DateTime&    today,     // (I) 
                             const DoubleArray* notionals, // (I), optional
                             const DoubleArray& fixings,   // (I)
                             CashFlowArray&     cfl)       // (M) the result
{
    static const string method = "PayStream::cashFlowList";
    double accInt = 0.0;          // ACCRUED_INTEREST  (up to today) 
    double historicAccruedInterest = 0.0;
    
    // This is the non-registered field in the PayStream
    accruedInterest = 0.0;

    if (notionals && (notionals->size() != payments->size()))  {
        throw ModelException(method,
                             "Number of supplied notionals is different "
                             "from the number of payments.");
    }
    if (cfl.size() != payments->size()) { 
        throw ModelException(method,
                             "Length of supplied CashFlowArray is different "
                             "from the number of payments.");
    }

    // Just loop through each payment 
    int fixingsOffset = 0;
    for (int idx = 0; idx < payments->size(); idx++)
    {
        accInt = 0.0;

        cfl[idx] = (*payments)[idx]->cashFlow(dayCountConv.get(),
                                              compoundFlat,
                                              &fixings[fixingsOffset], 
                                              today,
                                              &accInt,
                                              notionals ? &((*notionals)[idx]) : 0);

        PayStream::PaymentSP thisPayment = (*payments)[idx];
        int numRefixesThisPayment = thisPayment->refixes->size();
        fixingsOffset += numRefixesThisPayment;
        if ( thisPayment->payDate.isGreater(today) ) {
            /* catch case where today is between last accrueEnd of a payment
               and the payDate. */
            PayStream::Payment::RefixSP lastRefixThisPayment = 
                (*(thisPayment->refixes))[numRefixesThisPayment-1];
            if (today.isGreater(lastRefixThisPayment->accrueEnd))
            {
                historicAccruedInterest += cfl[idx].amount;
            }
            accruedInterest += accInt;
        }
    }

    accruedInterest += historicAccruedInterest;
}


/** Calculate the cash flow list arising from a payment stream 
    containing payments after 'fromDate' */
CashFlowArrayConstSP PayStream::cashFlowListAfterDate(
    const DateTime& today,        /* (I) value date */
    const DateTime& fromDate)     /* (I) date from */
{
    CashFlowArraySP cfl(new CashFlowArray(0));

    for (int idx = 0; idx < payments->size(); idx++) {
        // Only payments after fromDate
        if ( (*payments)[idx]->payDate.isGreater(fromDate) ) {
            double accInt = 0.0; // not used
            CashFlow cf = (*payments)[idx]->cashFlow(dayCountConv.get(),
                                                     compoundFlat,
                                                     floatRate.get(),
                                                     today,
                                                     &accInt,
                                                     0); // always use notionals from payment
            cfl->push_back(cf);
        }
    }
    return cfl;
}

/** Calculate the cash flow list arising from a payment stream.
    The number of notionals must match the number of payments.
    For use when the parent instrument has a call date */
CashFlowArrayConstSP PayStream::cashFlowList(
    const DateTime& today,               
    const DateTime& callDatePlusOffset,
    double callFwdRate,
    const DoubleArray* notionals)
{
    static const string method = "PayStream::cashFlowList";
    double accInt = 0.0;          // ACCRUED_INTEREST  (up to today) 
    double historicAccruedInterest = 0.0;
    
    // This is the non-registered field in the PayStream
    accruedInterest = 0.0;
    
    if (notionals && (notionals->size() != payments->size()))
    {
        throw ModelException(method,
                             "Number of supplied notionals is different "
                             "to the number of payments.");
    }

    CashFlowArraySP cfl(new CashFlowArray(payments->size()));
    
    // Just loop through each payment 
    for (int idx = 0; idx < payments->size(); idx++)
    {
        accInt = 0.0;
        
        (*cfl)[idx] = (*payments)[idx]->cashFlow(dayCountConv.get(),
                                                 compoundFlat,
                                                 floatRate.get(),
                                                 today,
                                                 &accInt,
                                                 notionals ? &((*notionals)[idx]) : 0,
                                                 callDatePlusOffset,
                                                 callFwdRate);
        
        PayStream::PaymentSP thisPayment = (*payments)[idx];
        if ( thisPayment->payDate.isGreater(today) ) {
            /* catch case where today is between last accrueEnd of a payment
               and the payDate. */
            int numRefixesThisPayment = thisPayment->refixes->size();
            PayStream::Payment::RefixSP lastRefixThisPayment = 
                (*(thisPayment->refixes))[numRefixesThisPayment-1];
            if (today.isGreater(lastRefixThisPayment->accrueEnd))
            {
                historicAccruedInterest += (*cfl)[idx].amount;
            }
            accruedInterest += accInt;
        }
    }
    
    accruedInterest += historicAccruedInterest;
    
    return cfl;
}  

/** Calculate the cash flow list arising from a payment stream.,
    but only up to the cut off payment date passed in. i.e
    generate a cash flow list from a stream within a stream */
CashFlowArrayConstSP PayStream::cashFlowList(const DateTime& fromDate,
                                             const DateTime& lastPayDate)const
{
    static const string method = "PayStream::cashFlowList";
    double accInt = 0.0;  // not used

    CashFlowArraySP cfl(new CashFlowArray(0));
    PayStream::PaymentSP payment;
    int payIdx = 0;
    
    // make sure lastPayDate is actually one of the pay dates
    for (payIdx = 0; payIdx < payments->size(); payIdx++)
    {
        payment = (*payments)[payIdx];
        if (payment->payDate.equals(lastPayDate))
        {
            break;
        }
    }

    if (payIdx < payments->size())
    {
        // Just loop through each payment 
        DateTime payDate;
        CashFlow cf;
        int i = 0;
        do
        {
            payDate = (*payments)[i]->payDate;

            cf = (*payments)[i]->cashFlow(dayCountConv.get(),
                                          compoundFlat,
                                          floatRate.get(),
                                          fromDate,
                                          &accInt,
                                          0);

            cfl->push_back(cf);
            i++;
        } while(!payDate.equals(lastPayDate));
    }
    else
    {
        throw ModelException(method,
                             "payment date "+
                             lastPayDate.toString() +
                             " was not found in the payment stream.");  
    }
    
    return cfl;
}

// return any cashflow due in the next CPND_DATE_WINDOW calendar days
// together with coupon due date and the due date adjusted to next valid business day
AccrualCalendarArraySP PayStream::couponDue(
    const DateTime&      fromDate, 
    bool                 isCallable, 
    const DateTime&      cutoff,
    double               callFwdRate) const 
{
    static const string method("PayStream::couponDue");

    try {
        AccrualCalendarArraySP couponsDue(new AccrualCalendarArray(0));

        // roll forward CPND_DATE_WINDOW calendar days
        DateTime toDate = fromDate.rollDate(AccrualCalendar::CPND_DATE_WINDOW);

        int i=0, j=0;
        double accInt=0.0;          // Not needed
        CashFlowArraySP cfl(new CashFlowArray(0));

        // Get known cashflows falling within CPND_DATE_WINDOW window
        for (i=0; i<payments->size() && (*payments)[i]->payDate < toDate; i++) {

            if ((*payments)[i]->payDate >= fromDate) {
                
                // sufficient check to use fromDate in isKnown and cashflow() when fromDate maybe tomorrow?
                if ((*payments)[i]->isKnown(fromDate, isFloating())) {
                    if (!isCallable) {
                        cfl->push_back((*payments)[i]->cashFlow(dayCountConv.get(),
                                                                compoundFlat,
                                                                floatRate.get(),
                                                                fromDate,
                                                                &accInt,
                                                                0));
                    }
                    else {
                        cfl->push_back((*payments)[i]->cashFlow(dayCountConv.get(),
                                                                compoundFlat,
                                                                floatRate.get(),
                                                                fromDate,
                                                                &accInt,
                                                                0,
                                                                cutoff,
                                                                callFwdRate));
                    }   
                }
            }
        }

        DateTime payDate;
        HolidaySP hols(Holiday::weekendsOnly());
        BadDayFollowing bdf;

        // Push cashflows into results object
        for (j=0; j<cfl->size(); j++) {
            payDate = (*cfl)[j].date;
            AccrualCalendar couponData(payDate, bdf.adjust(payDate, hols.get()), (*cfl)[j].amount);
            couponsDue->push_back(couponData);
        }
        return couponsDue;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// return ACCA_DATE_WINDOW calendar days of accrued interest data starting with fromDate 
AccrualCalendarArraySP PayStream::accrualCalendar(const DateTime& date) const {
    static const string method("PayStream::accrualCalendar");
    try {

        // Accrual calcs are independent of time of day so use SOD
        DateTime fromDate(date.getDate(), DateTime::timeConvert(DateTime::START_OF_DAY));

        AccrualCalendarArraySP accrualData(new AccrualCalendarArray(0));
        //CashFlowArraySP nextCoupons(getCoupons(fromDate));

        DateTime accDate;
        int i, j;
        bool foundNext = false;
        for (i=0; i<AccrualCalendar::ACCA_DATE_WINDOW; i++) {

            foundNext = false;
            accDate = fromDate.rollDate(i);

            // Need the strictly next cashflow date 
            for (j=0; j<payments->size(); j++) {
                if (accDate < (*payments)[j]->payDate) {
                    foundNext = true;
                    break;
                }
            }

            // Populate provided accDate falls between datedDate and final coupon date
            if (foundNext && (accDate <= getLastAccrueDate())) {
                AccrualCalendar accCal(accDate, (*payments)[j]->payDate, getAccruedInterest(fromDate, accDate, 0));
                accrualData->push_back(accCal);
            }
        }

        return accrualData;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Calculate the cash flow arising from a payment.
    If we have a FloatRate then the future rates will be calculated, 
    otherwise the rates inside the payment will be used.
    For the notional, either pass NULL to use the notionals recorded inside
    the payments, otherwise pass a pointer to a double containing the single
    notional to use for that payment. */   
CashFlow PayStream::Payment::cashFlow(const   DayCountConvention* payDayCountConv, 
                                      bool    compoundFlat,                        
                                      const   FloatRate* floatRate,                
                                      const   DateTime& today,
                                      double *accInt,
                                      const double* theNotional)const
{
    CashFlow cf;

    /* The computational method for the two different compounding types is
       quite different. Index rate flat uses an iterative formula. The rate
       variable holds the previous value and this is used to calculate the 
       rate needed for the next period. Specifically
       rate_i = rate_(i-1) * (x_i + 1) + (x_i + s_i)
       here x = yearFrac * refix Level, s = yearFrac * spread.
       Hence we need rate_0 = 0.
       In words, we grow the rate accumulated from the last period by the 
       rate for this period (no spread) and add on the rate for this period
       which includes the spread. (It's not too hard to show that the above
       method is equivalent to the formula in the spec.
       
       For index rate + spread, the formula is just a simple product of
       the rates for each of the periods - this needs 1 to be subtracted
       from the final rate. */

    CashFlowArraySP refixList;

    /* if it is a floating rate pay stream then calculate
       the forward rates for dates after today */
    if (floatRate) {
        refixList = CashFlowArraySP(getRefixes());
        DateTimeArraySP payDates = DateTimeArraySP(new DateTimeArray(refixes->size(),payDate));
        floatRate->calculateRates(today, *refixList, *payDates);
    }

    DoubleArray fixings(refixes->size());

    for (int i = 0; i < refixes->size(); i++) {
        fixings[i] = floatRate ? (*refixList)[i].amount : (*refixes)[i]->refixLevel;
    }

    return cashFlow(payDayCountConv,
                    compoundFlat,
                    &fixings.front(),
                    today,
                    accInt,
                    theNotional);
}

CashFlow PayStream::Payment::cashFlow(
    const DayCountConvention* payDayCountConv,   /* (I) Pay date count convention */
    bool  compoundFlat,                          /* (I) TRUE = Index rate flat,
                                                    FALSE = Index rate + spread */  
    const double* fixings,                       // (I) fix/flt fixings. "double*" allows more flexibility
    const DateTime& today,                       /* (I) value date */                       
    double* accInt,                              /* (M) accruedInterest so far */
    const double* theNotional)const                 /* (I) Use NULL to use notional inside
                                                    payment otherwise this overrides it */
{
    //need to modify if there'll be diff ways to calc the CF with Bus/252
    if(dynamic_cast<const Business252*>(payDayCountConv)){
        return cashFlowBrazil(
            payDayCountConv,   /* (I) Pay date count convention */
            compoundFlat,                          /* (I) TRUE = Index rate flat,
                                                            FALSE = Index rate + spread */  
            fixings,                       // (I) fix/flt fixings. "double*" allows more flexibility
            today,                       /* (I) value date */                       
            accInt,                              /* (M) accruedInterest so far */
            theNotional);                 /* (I) Use NULL to use notional inside*/
    }else{
        return cashFlowStandard(
            payDayCountConv,   /* (I) Pay date count convention */
            compoundFlat,                          /* (I) TRUE = Index rate flat,
                                                            FALSE = Index rate + spread */  
            fixings,                       // (I) fix/flt fixings. "double*" allows more flexibility
            today,                       /* (I) value date */                       
            accInt,                              /* (M) accruedInterest so far */
            theNotional);                 /* (I) Use NULL to use notional inside
                                                    payment otherwise this overrides it */
    }
}


/* it's for Brazil swap
    The day Count convention is business day/252
    Interest amount is calculated different than others, daily compounding
    [(1 + refix Level )*(1+spread)]^(yearFrac)-1
    weights =(* refix->weight) are used as the fraction of CDI agreed in the contract. 
    {[(1 + refix Level )^(1/252)-1]*weight + 1}^{nb Bus Day}
  */
CashFlow PayStream::Payment::cashFlowBrazil(
    const DayCountConvention* payDayCountConv,   /* (I) Pay date count convention */
    bool  compoundFlat,                          /* (I) TRUE = Index rate flat,
                                                    FALSE = Index rate + spread */  
    const double* fixings,                       // (I) fix/flt fixings. "double*" allows more flexibility
    const DateTime& today,                       /* (I) value date */                       
    double* accInt,                              /* (M) accruedInterest so far */
    const double* theNotional)const                 /* (I) Use NULL to use notional inside
                                                    payment otherwise this overrides it */
{
    static const string method("PayStream::Payment::cashFlowBrazil");
    try {
        CashFlow cf;

        /* The computational method for the two different compounding types is
           quite different. Index rate flat uses an iterative formula. The rate
           variable holds the previous value and this is used to calculate the 
           rate needed for the next period. 

           In words, we grow the rate accumulated from the last period by the 
           rate for this period (no spread) and add on the rate for this period
           which includes the spread. (It's not too hard to show that the above
           method is equivalent to the formula in the spec.
       
           For index rate + spread, the formula is just a simple product of
           the rates for each of the periods - this needs 1 to be subtracted
           from the final rate. */

        

        double rate = compoundFlat? 0.0: 1.0;
        double fixRate = 0.0;
        double yearFrac = 0.0;
        double notionalToUse = 0.0;
        double rateAccSoFar = 0.0;
        double yearFracAccSoFar = 0;
        bool   calcAccrueSoFar = false;

        // loop through each refix period
        for (int idx = 0; idx < refixes->size(); idx++)
        {
            RefixSP refix = (*refixes)[idx];

            // get the year fraction for each accrue period 
            yearFrac = payDayCountConv->years(refix->accrueStart,
                                              refix->accrueEnd);

            //if there is no bus day between accrueStart and accrueEnd
            //don't need to accrue rate
            if (yearFrac !=0.0){

                fixRate = fixings[idx] ;

                /* check if today is between this accrueStart and accrueEnd.
                   If so, we need to calculate ACCRUED_INTEREST so far in this 
                   accrual period */
                if (today.isGreater(refix->accrueStart) &&
                    refix->accrueEnd.isGreaterOrEqual(today))
                {
                    yearFracAccSoFar = payDayCountConv->years(refix->accrueStart,
                                                              today);
                    calcAccrueSoFar = true;
                }

                double yearFracToUse;
                double rateAmt;

                if (today.isGreater(refix->accrueStart) && 
                    refix->accrueEnd.isGreaterOrEqual(today)){
                    yearFracToUse = yearFracAccSoFar;
                } else{
                    yearFracToUse = yearFrac;
                }
            
                if (compoundFlat)
                {
                    throw ModelException(method, 
                                         "Compound Flat for Brazil Swap hasn't been implemented yet!");               

                    // rate is zero on first pass through 
                }
                else
                {
                    // rate is 1.0 first time through 
                    double spread = pow( 1.0 + refix->spread, yearFracToUse );
                    //DCC is bus/252
                    double r =  (pow( 1.0 + fixRate, 1.0/252.0 ) - 1.0) * refix->weight + 1.0;
                    r = pow (r, yearFracToUse * 252.0);
                    rateAmt = rate * r * spread;                       
                    
//                    double r = (1.0 + fixRate) * (1.0 + refix->spread);
//                    rateAmt = rate * pow( r, yearFracToUse );                       
                }

                if (today.isGreater(refix->accrueStart) && 
                    refix->accrueEnd.isGreaterOrEqual(today)){
                    rateAccSoFar = rateAmt;
                } 

                rate = rateAmt;
            }
        }

        if (refixes->size())
        {
            if (!compoundFlat)
            {
                /* for compounding using rate + spread above computation 
                   gives us 1 + the rate to use */
                rate -= 1.0;
                if (calcAccrueSoFar)
                {
                    rateAccSoFar -= 1.0;
                }
            }

            if (theNotional)
            {
                notionalToUse = *theNotional;
            }
            else
            {
                notionalToUse = notional;  // within the payment
            }
        
            cf.date = payDate;
            cf.amount = notionalToUse * rate;

            if (calcAccrueSoFar)
            { 
                *accInt = notionalToUse * rateAccSoFar;
            }
        }
    
        return cf;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}



// like the other flavours except fixings are passed in explicitly
CashFlow PayStream::Payment::cashFlowStandard(
    const DayCountConvention* payDayCountConv,   /* (I) Pay date count convention */
    bool  compoundFlat,                          /* (I) TRUE = Index rate flat,
                                                    FALSE = Index rate + spread */  
    const double* fixings,                       // (I) fix/flt fixings. "double*" allows more flexibility
    const DateTime& today,                       /* (I) value date */                       
    double* accInt,                              /* (M) accruedInterest so far */
    const double* theNotional)const                 /* (I) Use NULL to use notional inside
                                                    payment otherwise this overrides it */
{
    static const string method("PayStream::Payment::cashFlowStandard");
    try {
        CashFlow cf;

        /* The computational method for the two different compounding types is
           quite different. Index rate flat uses an iterative formula. The rate
           variable holds the previous value and this is used to calculate the 
           rate needed for the next period. Specifically
           rate_i = rate_(i-1) * (x_i + 1) + (x_i + s_i)
           here x = yearFrac * refix Level, s = yearFrac * spread.
           Hence we need rate_0 = 0.
           In words, we grow the rate accumulated from the last period by the 
           rate for this period (no spread) and add on the rate for this period
           which includes the spread. (It's not too hard to show that the above
           method is equivalent to the formula in the spec.
       
           For index rate + spread, the formula is just a simple product of
           the rates for each of the periods - this needs 1 to be subtracted
           from the final rate. */

        double rate = compoundFlat? 0.0: 1.0;
        double fixRate = 0.0;
        double yearFrac = 0.0;
        double notionalToUse = 0.0;
        double rateAccSoFar = 0.0;
        double yearFracAccSoFar = 0.0;
        int nbDayAccSoFar =0;            
        bool   calcAccrueSoFar = false;

        // loop through each refix period
        for (int idx = 0; idx < refixes->size(); idx++)
        {
            RefixSP refix = (*refixes)[idx];

            fixRate = fixings[idx] * refix->weight;;

            // get the year fraction for each accrue period 
            yearFrac = payDayCountConv->years(refix->accrueStart,
                                              refix->accrueEnd);

            /* check if today is between this accrueStart and accrueEnd.
               If so, we need to calculate ACCRUED_INTEREST so far in this 
               accrual period */
            if (today.isGreater(refix->accrueStart) &&
                refix->accrueEnd.isGreaterOrEqual(today))
            {
                yearFracAccSoFar = payDayCountConv->years(refix->accrueStart,
                                                          today);
                calcAccrueSoFar = true;
            }

        
            if (compoundFlat)
            {
                if (today.isGreater(refix->accrueStart) && 
                    refix->accrueEnd.isGreaterOrEqual(today))
                {
                    rateAccSoFar = rate * (yearFracAccSoFar * fixRate + 1) + 
                        yearFracAccSoFar * (fixRate + refix->spread);
                } 
                // rate is zero on first pass through 
                rate = rate * (yearFrac * fixRate + 1) + 
                    yearFrac * (fixRate + refix->spread);
            }
            else
            {
                if (today.isGreater(refix->accrueStart) && 
                    refix->accrueEnd.isGreaterOrEqual(today))
                {
                    rateAccSoFar = rate * (1.0 + yearFracAccSoFar * 
                                           (fixRate + refix->spread));
                } 

                // rate is 1.0 first time through 
                rate *= 1.0 + yearFrac * (fixRate + refix->spread);
            }
        }

        if (refixes->size())
        {
            if (!compoundFlat)
            {
                /* for compounding using rate + spread above computation 
                   gives us 1 + the rate to use */
                rate -= 1.0;
                if (calcAccrueSoFar)
                {
                    rateAccSoFar -= 1.0;
                }
            }

            if (theNotional)
            {
                notionalToUse = *theNotional;
            }
            else
            {
                notionalToUse = notional;  // within the payment
            }
        
            cf.date = payDate;
            cf.amount = notionalToUse * rate;

            if (calcAccrueSoFar)
            { 
                *accInt = notionalToUse * rateAccSoFar;
            }
        }
    
        return cf;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Calculate the cash flow arising from a payment.
    If call date + offset is within the last refix period 
    then the following behaviour is invoked :
    (A) Accrue over the Notional from accrueStart to accrueEnd as usual.
    (B) Accrue over the Notional from call date + offset to accrueEnd using the fwd rate at call date + offset.
    Subtract (A) from (B) */
CashFlow PayStream::Payment::cashFlow(
    const DayCountConvention* payDayCountConv,
    bool  compoundFlat,                       
    const FloatRate* floatRate,               
    const DateTime& today,                    
    double* accInt,                           
    const double* theNotional,                   
    const DateTime&  callDatePlusOffset,
    const double callFwdRate)const{

    if(dynamic_cast<const Business252*>(payDayCountConv)){
        return cashFlowBrazil(
                    payDayCountConv,   /* (I) Pay date count convention */
                    compoundFlat,      /* (I) TRUE = Index rate flat,
                                              FALSE = Index rate + spread */  
                    floatRate,         /* (I) NULL unless payStream is a float type */
                    today,             /* (I) value date */                       
                    accInt,            /* (M) accruedInterest so far */
                    theNotional,       /* (I) Use NULL to use notional inside
                                              payment otherwise this overrides it */
                    callDatePlusOffset, /* it is callable --> supply the date of call date + offset */ 
                    callFwdRate);       /* fwd rate between callDatePlusOffset and final accrueEnd date */


    }else{
        return cashFlowStandard(
                    payDayCountConv,   /* (I) Pay date count convention */
                    compoundFlat,      /* (I) TRUE = Index rate flat,
                                              FALSE = Index rate + spread */  
                    floatRate,         /* (I) NULL unless payStream is a float type */
                    today,             /* (I) value date */                       
                    accInt,            /* (M) accruedInterest so far */
                    theNotional,       /* (I) Use NULL to use notional inside
                                              payment otherwise this overrides it */
                    callDatePlusOffset, /* it is callable --> supply the date of call date + offset */ 
                    callFwdRate);       /* fwd rate between callDatePlusOffset and final accrueEnd date */
    }
}


/** Calculate the cash flow arising from a payment.
    If call date + offset is within the last refix period 
    then the following behaviour is invoked :
    (A) Accrue over the Notional from accrueStart to accrueEnd as usual.
    (B) Accrue over the Notional from call date + offset to accrueEnd using the fwd rate at call date + offset.
    Subtract (A) from (B) */
CashFlow PayStream::Payment::cashFlowStandard(
    const DayCountConvention* payDayCountConv,
    bool  compoundFlat,                       
    const FloatRate* floatRate,               
    const DateTime& today,                    
    double* accInt,                           
    const double* theNotional,                   
    const DateTime&  callDatePlusOffset,
    const double callFwdRate)const
{
    double rate = compoundFlat? 0.0: 1.0;
    double fixRate = 0.0;
    double yearFrac = 0.0;
    double notionalToUse = 0.0;
    double rateAccSoFar = 0.0;
    double yearFracAccSoFar = 0.0;
    double yearFracTodayToAccEnd = 0.0;
    bool   calcAccrueSoFar = false;
    bool   earlyUnwind = false;
    CashFlowArraySP refixList;

    CashFlow cf;

    /* if it is a floating rate pay stream then calculate
       the forward rates for dates after today */
    if (floatRate)
    {
        refixList = CashFlowArraySP(getRefixes());
        DateTimeArraySP payDates = DateTimeArraySP(new DateTimeArray(refixes->size(),payDate));
        floatRate->calculateRates(today,
                                  *refixList, *payDates);

    }

    earlyUnwind = false;
    // loop through each refix period
    for (int idx = 0; idx < refixes->size(); idx++)
    {
        RefixSP refix = (*refixes)[idx];
        if (callDatePlusOffset.isGreater(refix->accrueStart) &&
            refix->accrueEnd.isGreater(callDatePlusOffset))
        {
            earlyUnwind = true;
            /* note: dont need to do anything special if call date + offset
            is equal to an accrueStart/End */
        }

        if (floatRate)
        {
            fixRate = (*refixList)[idx].amount * refix->weight;
        }
        else
        {
            fixRate = refix->refixLevel * refix->weight;
        }

        // get the year fraction for each accrue period 
        yearFrac = payDayCountConv->years(refix->accrueStart,
                                          refix->accrueEnd);

        /* check if today is between this accrueStart and accrueEnd.
           If so, we need to calculate ACCRUED_INTEREST so far in this 
           accrual period */
        if (today.isGreater(refix->accrueStart) &&
            refix->accrueEnd.isGreaterOrEqual(today))
        {
            yearFracAccSoFar = payDayCountConv->years(refix->accrueStart,
                                                      today);
            yearFracTodayToAccEnd = payDayCountConv->years(today,
                                                           refix->accrueEnd);
            calcAccrueSoFar = true;
        }

        if (compoundFlat)
        {
            if (today.isGreater(refix->accrueStart) && 
                refix->accrueEnd.isGreaterOrEqual(today))
            {
                if (earlyUnwind)
                {
                    // accrue at the same effective rate that is used for callable streams 
                    // i.e we dont use the fwd rate from today to accrue end, we use the callFwdRate
                    rateAccSoFar =  rate * (1 + ((yearFrac * fixRate) - (yearFracTodayToAccEnd * callFwdRate))) +
                        ((yearFrac * (fixRate + refix->spread)) - (yearFracTodayToAccEnd * (callFwdRate + refix->spread)));
                }
                else
                { 
                    rateAccSoFar = rate * (yearFracAccSoFar * fixRate + 1) + 
                        yearFracAccSoFar * (fixRate + refix->spread);
                }
            } 
            if (earlyUnwind)
            {
                double callFwdYearFrac = payDayCountConv->years(callDatePlusOffset,
                                                                refix->accrueEnd);

                rate = rate * (1 + ((yearFrac * fixRate) - (callFwdYearFrac * callFwdRate))) +
                    ((yearFrac * (fixRate + refix->spread)) - (callFwdYearFrac * (callFwdRate + refix->spread)));

            }
            else
            {
                // rate is zero on first pass through 
                rate = rate * (yearFrac * fixRate + 1) + 
                    yearFrac * (fixRate + refix->spread);
            }

        }
        else
        {
            if (today.isGreater(refix->accrueStart) && 
                refix->accrueEnd.isGreaterOrEqual(today))
            {
                if (earlyUnwind)
                {
                    // accrue at the same effective rate that is used for callable streams 
                    // i.e we dont use the fwd rate from today to accrue end, we use the callFwdRate
                    rateAccSoFar = rate * (1.0 + ((yearFrac * (fixRate + refix->spread)) -
                                                  (yearFracTodayToAccEnd * (callFwdRate + refix->spread))));
                }
                else
                {
                    rateAccSoFar = rate * (1.0 + yearFracAccSoFar * 
                                           (fixRate + refix->spread));
                }
            } 

            if (earlyUnwind)
            {
                double callFwdYearFrac = payDayCountConv->years(callDatePlusOffset,
                                                                refix->accrueEnd);

                rate *= 1.0 + ((yearFrac * (fixRate + refix->spread)) -
                               (callFwdYearFrac * (callFwdRate + refix->spread)));
                
            }
            else
            {
                // rate is 1.0 first time through 
                rate *= 1.0 + yearFrac * (fixRate + refix->spread); 
            }
        }
    }

    if (refixes->size())
    {
        if (!compoundFlat)
        {
            /* for compounding using rate + spread above computation 
            gives us 1 + the rate to use */
            rate -= 1.0;
            if (calcAccrueSoFar)
            {
                rateAccSoFar -= 1.0;
            }
        }

        if (theNotional)
        {
            notionalToUse = *theNotional;
        }
        else
        {
            notionalToUse = notional;  // within the payment
        }
        
        cf.date = payDate;
        cf.amount = notionalToUse * rate;

        if (calcAccrueSoFar)
        { 
            *accInt = notionalToUse * rateAccSoFar;
        }
    }
    
    return cf;
}


/** Calculate the cash flow arising from a payment.
    If call date + offset is within the last refix period 
    then the following behaviour is invoked :
    (A) Accrue over the Notional from accrueStart to accrueEnd as usual.
    (B) Accrue over the Notional from call date + offset to accrueEnd using the fwd rate at call date + offset.
    Subtract (A) from (B) */
CashFlow PayStream::Payment::cashFlowBrazil(
    const DayCountConvention* payDayCountConv,
    bool  compoundFlat,                       
    const FloatRate* floatRate,               
    const DateTime& today,                    
    double* accInt,                           
    const double* theNotional,                   
    const DateTime&  callDatePlusOffset,
    const double callFwdRate)const
{
    double rate = compoundFlat? 0.0: 1.0;
    double fixRate = 0.0;
    double yearFrac = 0.0;
    double notionalToUse = 0.0;
    double rateAccSoFar = 0.0;
    double yearFracAccSoFar = 0.0;
    double yearFracTodayToAccEnd = 0.0;
    bool   calcAccrueSoFar = false;
    bool   earlyUnwind = false;
    CashFlowArraySP refixList;

    CashFlow cf;
    
    //to do
    
    throw ModelException("CashFlow PayStream::Payment::cashFlowBrazil",
                         "Callable isn't available for Brazil equity swap.");

    return cf;
}


/** populate the passed cash flow array with the refix dates and
    levels within the payment */
CashFlowArraySP PayStream::Payment::getRefixes()const
{
    CashFlowArraySP refixList(new CashFlowArray(0));

    for (int idx = 0; idx < refixes->size(); idx++)
    {
        CashFlow cf;
        cf.date = (*refixes)[idx]->refixDate;
        cf.amount = (*refixes)[idx]->refixLevel;

        refixList->push_back(cf);
    }
    return refixList;
}           

/** do we have all the data today to calculate this ? */
bool PayStream::Payment::isKnown(const DateTime& today, bool isFloating) const {
    if (!isFloating) {
        return true;  // if it's fixed, must know it
    }

    bool fixed = true;
    for (int i = 0; i < refixes->size() && fixed; i++) {
        fixed = today.isGreaterOrEqual((*refixes)[i]->refixDate);
    }
    return fixed;
}

/** is the payment known and started accruing? */
bool PayStream::Payment::isKnownAndStartedAccruing(const DateTime& today, bool isFloating) const {

    bool refixesAreKnown = false;
    if (!isFloating) {
        refixesAreKnown = true;  // if it's fixed, must know it
    } else {
        bool fixed = true;
        for (int i = 0; i < refixes->size() && fixed; i++) {
            fixed = today.isGreaterOrEqual((*refixes)[i]->accrueStart);
        }
        refixesAreKnown = fixed;
    }

    bool startedAccruing = true;
    for (int i = 0; i < refixes->size() && startedAccruing; i++) {
        startedAccruing = today.isGreaterOrEqual((*refixes)[i]->accrueStart);
    }

    return refixesAreKnown && startedAccruing;
}

/** Returns the floatRate Ccy */
string PayStream::getCcy()const
{
    return floatRate->getCcy();
}

/** Returns the floatRate YC Name */
string PayStream::getYCName() const
{
	return floatRate->getYCName();
}

/* Validation for libor flow. Checks that the payment dateTime is after
   the last refix dateTime for the payment period. */
void PayStream::checkRefixAndPayment()const
{
    static const string method = "PayStream::checkRefixAndPayment";
    
    for (int idx = 0; idx < payments->size(); idx++)
    {
        PaymentSP payment = (*payments)[idx];
        Payment::RefixSP lastRefix = (*(payment->refixes))[payment->refixes->size()-1];
        
        if (payment->payDate.isLess(lastRefix->refixDate))
        {
            throw ModelException(method,
                                 "The payment date time " +
                                 payment->payDate.toString() +
                                 " must be after the last refix date time " +
                                 lastRefix->refixDate.toString());
        }
    }
}

/** Invoked indirectly from THETA shift method. 
    Rolls historical refixes */
void PayStream::Payment::Refix::rollRefix(const DateTime& valueDate,
                                          const DateTime& newValueDate,
                                          const FloatRate& floatRate)
{
    double rate = 0.0;

    /* If there is a refixing between the value date (exclusive) and the 
       next value date (inclusive) , then set the corresponding refixing 
       to today's spot rate */ 
    if ( ( refixDate.isGreater(valueDate) && !refixDate.isGreater(newValueDate)) || 
         ( refixDate.equals(valueDate)     && Maths::isZero(refixLevel))         )
    {
        // get the rate on the valueDate
        rate = floatRate.calculateRate(valueDate);

        refixLevel = rate;
    }
}

   
/** Invoked indirectly by THETA shift method. Rolls the 
    historical refixes for a payment */
void PayStream::Payment::rollRefixes(const DateTime& valueDate,
                                     const DateTime& newValueDate,
                                     const FloatRate& floatRate)
{
    for (int i = 0; i < refixes->size(); i++)
    {
        (*refixes)[i]->rollRefix(valueDate,
                                 newValueDate,
                                 floatRate);
    }
}


/** Invoked by THETA shift method. Roll the historical refixes 
    a day forward */
void PayStream::rollRefixes(const DateTime& valueDate,
                            const DateTime& newValueDate)
{
    // roll refixes for each payment
    for (int i = 0; i < payments->size(); i++)
    {
        (*payments)[i]->rollRefixes(valueDate,
                                    newValueDate,
                                   *floatRate);
    }
}

class RefixHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spread sheet
        REGISTER(PayStream::Payment::Refix, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultRefix);
        FIELD(refixDate, "refixing date");
        FIELD(accrueStart, "accrue start date");
        FIELD(accrueEnd, "accrue end date"); 
        FIELD(refixLevel, "refixing level");
        FIELD(spread, "spread");
        FIELD(weight, "weighting");
    }

    static IObject* defaultRefix(){
        return new PayStream::Payment::Refix();
    }
};

CClassConstSP const PayStream::Payment::Refix::TYPE = CClass::registerClassLoadMethod(
    "Refix", typeid(PayStream::Payment::Refix), RefixHelper::load);


class PaymentHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spread sheet
        REGISTER(PayStream::Payment, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPayment);
        FIELD(payDate, "payment date");
        FIELD(refixes, "array of refixes");
        FIELD(notional, "notional"); 
    }

    static IObject* defaultPayment(){
        return new PayStream::Payment();
    }
};

CClassConstSP const PayStream::Payment::TYPE = CClass::registerClassLoadMethod(
    "Payment", typeid(PayStream::Payment), PaymentHelper::load);


class PayStreamHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PayStream, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IXLInterfaceMap);
        EMPTY_SHELL_METHOD(defaultPayStream);
        FIELD(payments, "array of payments");
        FIELD(floatRate, "floating rate");
        FIELD_MAKE_OPTIONAL(floatRate);
        FIELD(dayCountConv, "payment DCC");
        FIELD(compoundFlat, "TRUE = compoundFlat");
        // transients
        FIELD(accruedInterest, "acruedInterest");
        FIELD_MAKE_TRANSIENT(accruedInterest);
    }

    static IObject* defaultPayStream(){
        return new PayStream();
    }
};

class PayStreamCreateAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // addin parameters
    DayCountConventionSP payDayCountConv;
    bool  compoundFlat;
    FloatRateSP     floatRate;
    DoubleArraySP   notionals;
    DateTimeArraySP refixDates;
    DateTimeArraySP accrualDates;
    DateTimeArraySP payDates;
    DoubleArraySP   refixLevels;
    DoubleArraySP   spreads;
    DoubleArraySP   weights;
    CBoolArraySP    compound;  
    MarketDataConstSP marketData;


    // create a PayStream
    static IObjectSP create(PayStreamCreateAddin *params)
        {
            static const string method = "PayStreamCreateAddin::create";

            // do some validation on the input parameters
            if (params->payDayCountConv.get() == 0 || 
                params->notionals.get() == 0 || params->refixDates.get() == 0 ||
                params->accrualDates.get() == 0 || params->payDates.get() == 0 ||
                params->refixLevels.get() == 0 || params->spreads.get() == 0 ||
                params->weights.get() == 0 || params->compound.get() == 0)
            {
                throw ModelException(method,
                                     "One or more of payDayCountConv, "
                                     "notionals, refixDates, "
                                     "accrualDates, payDates, refixLevels, "
                                     "spreads, weights or compound is NULL.");
            }

	    if (params->marketData.get())
            {
                // May need to fetch holidays for Business/252 DCC
                params->marketData->fetchForNonMarketObjects(params->payDayCountConv, 
                                                             IModelConstSP(new NonPricingModel()),
                                                             "Business252");
            } 
            
            PayStreamSP newPayStream(new PayStream(params->payDayCountConv.get(),
                                                   params->compoundFlat,
                                                   params->floatRate.get(),
                                                   *(params->notionals),
                                                   *(params->refixDates),
                                                   *(params->accrualDates),
                                                   *(params->payDates),
                                                   *(params->refixLevels),
                                                   *(params->spreads),
                                                   *(params->weights),
                                                   *(params->compound)));

            return newPayStream;
        }


    PayStreamCreateAddin(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PayStreamCreateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPayStreamCreateAddin);
        // order of registration effects order of parameters in addin function
        FIELD(payDayCountConv, "day count convention");
        FIELD(compoundFlat, "compoundFlat");
        FIELD(floatRate, "floating rate");
        FIELD_MAKE_OPTIONAL(floatRate);
        FIELD(notionals, "notionals");
        FIELD(refixDates, "refix dates");
        FIELD(accrualDates, "accrual dates");
        FIELD(payDates, "payment dates");
        FIELD(refixLevels, "refixing levels");
        FIELD(spreads, "spreads");
        FIELD(weights, "weights");
        FIELD(compound, "TRUE = compounding period");
        FIELD(marketData, "market data");
        FIELD_MAKE_OPTIONAL(marketData);

        Addin::registerClassObjectMethod("PAY_STREAM",
                                         Addin::UTILITIES,
                                         "Creates a handle to a payment stream ",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)create);
    }
    
    static IObject* defaultPayStreamCreateAddin(){
        return new PayStreamCreateAddin();
    }
 
};

CClassConstSP const PayStreamCreateAddin::TYPE = CClass::registerClassLoadMethod(
    "PayStreamCreateAddin", typeid(PayStreamCreateAddin), load);


CClassConstSP const PayStream::TYPE = CClass::registerClassLoadMethod(
    "PayStream", typeid(PayStream), PayStreamHelper::load);

typedef PayStream::PaymentArray PayStreamPaymentArray; // msvc7 bug
DEFINE_TEMPLATE_TYPE_WITH_NAME("PaymentArray", PayStreamPaymentArray);

typedef PayStream::Payment::RefixArray PayStreamPaymentRefixArray; // msvc7 bug
DEFINE_TEMPLATE_TYPE_WITH_NAME("RefixArray", PayStreamPaymentRefixArray);

DEFINE_TEMPLATE_TYPE(PayStreamArray);


DRLIB_END_NAMESPACE


