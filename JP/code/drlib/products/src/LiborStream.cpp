//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LiborStream.cpp
//
//   Description : LiborStream instrument
//
//   Author      : Stephen Hope
//
//   Date        : 30 May 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/LiborStream.hpp"
#include "edginc/OutputRequestUtil.hpp"

DRLIB_BEGIN_NAMESPACE

/* for reflection */
LiborStream::LiborStream(): 
    CInstrument(TYPE), callDateOffset(0)
{
    ccyTreatment = CAsset::CCY_TREATMENT_NONE;
}

void LiborStream::Validate(){
    static const string method = "LiborStream::Validate";
    try {
        if (!discount) {
            throw ModelException(method, "Discount YC is null");
        }
        
		if (!(ccyTreatment == CAsset::CCY_TREATMENT_NONE  || 
			ccyTreatment == CAsset::CCY_TREATMENT_VANILLA || 
			ccyTreatment == CAsset::CCY_TREATMENT_STRUCK)) {
			throw ModelException(method, "Only None, Vanilla or Struck are currently supported as currency treatment");
		}

        payStream->checkRefixAndPayment();

        // Callable feature checks

        // need to hardcode this incase hol inside discount is NoHolidays
        HolidaySP weekendsOnly(Holiday::weekendsOnly());
        if (isCallable) {
            if (callDateOffset < 0)
            {
                throw ModelException(method,
                                     "call date offset cannot be negative !");

            }
            // if call date is in past move to today or the next business day
            if (valueDate.isGreater(callDate))
            {
                callDate = valueDate;
            }
            // incase its a weekend
            while(!weekendsOnly->isBusinessDay(callDate))
            {
                callDate = callDate.rollDate(1);
            }
            // we should adjust for the yield curve hols here ??
            // i.e. Monday might be a holiday according to YC Hol

            // No sense in SLS being callable if call date + offset
            // would give a date after the very last accrue end date
            DateTime lastAccrueDate = payStream->getLastAccrueDate();
            DateTime callableCutOffDate = lastAccrueDate.rollDate(-callDateOffset);
            
            if (callDate.isGreaterOrEqual(callableCutOffDate))
            {
                isCallable = false;
            }
            
            if (isCallable)
            {
                DateTime callDatePlusOffset = callDate.rollDate(callDateOffset);
                // incase its a weekend
                while (!weekendsOnly->isBusinessDay(callDatePlusOffset))
                {
                    callDatePlusOffset = callDatePlusOffset.rollDate(1);
                }

                payStream->curtail(callDatePlusOffset);
            }
        }

        // if there are fees, these must end on or before the paystream
        if (fees.get()) {
            if (isCallable) {
                // until we know what we're doing
                throw ModelException(method, "callable libor stream with fees "
                                     "not yet supported");
            }

            CashFlow::ensureDatesIncreasing(*fees, "fees", false);
            int numfee = fees->size();
            if (numfee > 0) {
                const DateTime& lastFee = (*fees)[numfee-1].date;
                if (lastFee > endDate(0)) {
                    throw ModelException(method, "last fee (" + 
                                         lastFee.toString() + 
                                         ") is after stream maturity (" + 
                                         endDate(0).toString() + ")");
                }
            }
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

// copy market data relevant to the instrument
void LiborStream::GetMarket(const IModel* model, const CMarketDataSP market){
	static const string method = "SimpleCashFlowStream::GetMarket";
    market->GetReferenceDate(valueDate);
    discount.getData(model, market);
	string discountName = discount->getName();
	string underlyingName = "no und";

	if (!underlyingCurve.isEmpty()) {
		underlyingCurve.getData(model, market);		
		underlyingName = underlyingCurve->getName();
	}

    /* The FloatRate is not represented as an asset inside the instrument
       but we may still have to get the yield curve from the market data cache */
    if (payStream->isFloating()) {
		payStream->getMarketData(model, market);
		string floatName = payStream->getYCName();
		if (underlyingName != "no und" && underlyingName != floatName) {
		throw ModelException(method, "When floating, underlyingCurve " + underlyingName
			+ " has to be empty or equal to floatRate " + floatName + " ");
		}
		underlyingName = floatName;
    }
	
	if (underlyingName != "no und" && underlyingName != discountName) {
		string fxName = market->getFXName(underlyingName, discountName);
		MarketObjectSP fxObj = market->GetData(fxName, FXAsset::TYPE);
		fxObj->getMarket(model, market.get());
		fxAsset = FXAssetSP::dynamicCast(fxObj);
		ccyTreatment = CAsset::CCY_TREATMENT_STRUCK;  // need to remove this line when ccyTreatment will possibly be quanto
	}
}

void LiborStream::price(Control* control,
                        CResults* results) const
{
    static const string method = "LiborStream::price";
    
    try
    {
        DateTime             callDatePlusOffset;
        double               callFwdRate = 0.0;
        CashFlowArrayConstSP cfl;
        if (isCallable)
        {
            callDatePlusOffset = callDate.rollDate(callDateOffset);
            // in case its a weekend
            HolidaySP weekendsOnly(Holiday::weekendsOnly());
            while(!weekendsOnly->isBusinessDay(callDatePlusOffset))
            {
                callDatePlusOffset = callDatePlusOffset.rollDate(1);
            }
            /* Calculate the fwd rate between the call date + offset
               and the accrueEnd date of the relevant refix period.*/
			if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK) {
				callFwdRate = payStream->callPlusOffsetFwdRate(callDatePlusOffset, underlyingCurve);
			}
			else {
				callFwdRate = payStream->callPlusOffsetFwdRate(callDatePlusOffset, discount);
			}

            // calculate the cash flow list from the payment stream
            cfl = payStream->cashFlowList(valueDate,
                                          callDatePlusOffset,
                                          callFwdRate,
                                          0);
        }
        else
        {
            // calculate the cash flow list from the payment stream
            cfl = payStream->cashFlowList(valueDate, 0);
        }
        
        CashFlowArraySP opcf(new CashFlowArray()); // for output only
        double value = 0.0;
        for (int idx = 0; idx < cfl->size(); idx++)
        {
            // exclude cash flows in the past 
            DateTime cfDate = (*cfl)[idx].date;
            if (cfDate.isGreater(valueDate))
            {
                // pv the cash flow
                value += (*cfl)[idx].amount * discount->pv(cfDate)
					* ((ccyTreatment == CAsset::CCY_TREATMENT_STRUCK) ? fxAsset->fwdValue(cfDate) : 1.0);
                if (control && control->isPricing())
                {
                    opcf->push_back(CashFlow(cfDate, (*cfl)[idx].amount)); 
                }
            }
        }

        // fold in any fees
        if (fees.get()) {
            for (int i = 0; i < fees->size(); i++) {
                const DateTime& feedate = (*fees)[i].date;
                // unclear what happens for callable here
                // do we pull the fee after call date to call date 
                // or drop it
                // validate against callable for now
#if 0
                if (feedate > valueDate && 
                    (!isCallable || 
                     (isCallable && feedate <= callDatePlusOffset))) {
                    value += (*fees)[i].amount * discount->pv(feedate)
						* ((ccyTreatment == CAsset::CCY_TREATMENT_STRUCK) ? fxAsset->fwdValue(feedate) : 1.0);
                }
#endif
                if (feedate > valueDate) {
                    value += (*fees)[i].amount*discount->pv(feedate)
						* ((ccyTreatment == CAsset::CCY_TREATMENT_STRUCK) ? fxAsset->fwdValue(feedate) : 1.0);
                }
            }
        }

        results->storePrice(value, discount->getCcy());

        if (control && control->isPricing()) {
            requests(control,results,isCallable,callDatePlusOffset,callFwdRate);
            if (opcf->size() > 0) {   
                const string name = "LIBOR";
                OutputNameConstSP irOutput(new OutputName(name));
                results->storeGreek(opcf, Results::DEBUG_PACKET, irOutput);
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Shifts the object using given shift. */
bool LiborStream::sensShift(Theta* shift)
{
    DateTime newValueDate = shift->rollDate(valueDate);

    if (payStream->isFloating()) {
        payStream->rollRefixes(valueDate, newValueDate);
    }

    valueDate = newValueDate;

    return true;
}

/** what's today ? */
DateTime LiborStream::getValueDate() const 
{
    return valueDate;
}

/** when to stop tweaking */
DateTime LiborStream::endDate(const Sensitivity* sensControl) const 
{
    return payStream->lastSensDate();
}

CreditSupportSP LiborStream::createCreditSupport(CMarketDataSP market){
    return CreditSupportSP(new LiborStreamCreditSupport(this, market));
}

/** Support for ITaxableInstBasic */
const DateTime LiborStream::getFinalPaymentDate() const {
    if (isCallable) {
        throw ModelException("LiborStream::getFinalPaymentDate",
                             "Tax is not supported for callable LiborStreams");
    }
    CashFlowArrayConstSP cfl = payStream->cashFlowList(valueDate, 0);
    return cfl->back().date;
}

/** Support for ITaxableInstWithCoupons */
CashFlowArrayConstSP LiborStream::getCoupons() const {
    // it's already done for us ...
    CashFlowArrayConstSP cfl = payStream->cashFlowList(valueDate, 0);
    return cfl;
}

void LiborStream::requests(Control*        control, 
                           CResults*       results,
                           bool            isCallable,
                           const DateTime& cutoff,
                           double          callFwdRate) const {
    static const string method = "LiborStream::requests";
    try {
        OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTimeArraySP dates(payStream->paymentDates());

            if (fees.get()) {
                for (int i = 0; i < fees->size(); i++) {
                    dates->push_back((*fees)[i].date);
                }
            }
            OutputRequestUtil::recordPaymentDates(control,results,dates.get()); 
       }
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            CashFlowArrayConstSP payCFL(payStream->knownCashflows(valueDate,
                                                                  0,
                                                                  isCallable,
                                                                  cutoff,
                                                                  callFwdRate));
            // handle fees
            CashFlowArraySP feeCFL(new CashFlowArray(0));
            if (fees.get()) {
                // what happens for callable ??
                for (int i = 0; i < fees->size(); i++) {
                    feeCFL->push_back((*fees)[i]);
                }
            }

            CashFlowArraySP payCFLnc(new CashFlowArray(0));
            for (int i = 0; i < payCFL->size(); i++) {
                payCFLnc->push_back((*payCFL)[i]);           
            }

            CashFlowArraySP cfl(CashFlow::merge(payCFLnc, feeCFL));

            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    payCFL.get()); 
        }
        
        if (control->requestsOutput(OutputRequest::ACCRUED_INTEREST, request)) {
            results->storeRequestResult(request, payStream->getAccruedInterest());
        }

        // COUPON_DUE 
        request = control->requestsOutput(OutputRequest::COUPON_DUE);
        if (request) {
            AccrualCalendarArraySP coupons = payStream->couponDue(valueDate,
                                                                  isCallable,
                                                                  cutoff,
                                                                  callFwdRate);
            if (coupons->empty()) {
                results->storeNotApplicable(request); 
            } else {
                results->storeRequestResult(request, coupons); 
            }
        }    

        // ACCRUAL_CALENDAR 
        request = control->requestsOutput(OutputRequest::ACCRUAL_CALENDAR);
        if (request) {
            AccrualCalendarArraySP accData = payStream->accrualCalendar(valueDate);
            if (accData->empty()) {
                results->storeNotApplicable(request); 
            } else {
                results->storeRequestResult(request, accData); 
            } 
        }    
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}



/** private class */
class LiborStreamClosedForm: public ClosedForm::IProduct{
private:
    const LiborStream* lib; // a reference

public:
    LiborStreamClosedForm(const LiborStream* lib): lib(lib){}

    void price(ClosedForm* model,
               Control*    control, 
               CResults*   results) const{
        lib->price(control, results);
    }
};

/** Implementation of ClosedForm::IntoProduct interface */
ClosedForm::IProduct* LiborStream::createProduct(
    ClosedForm* model) const{
    return new LiborStreamClosedForm(this);
}


/** Returns the name of the instrument's discount currency */
string LiborStream::discountYieldCurveName() const {
    return discount.getName();
}


class LiborStreamHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(LiborStream, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(ClosedForm::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(Theta::IShift);
        EMPTY_SHELL_METHOD(defaultLiborStream);
        FIELD(valueDate, "value date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(isCallable, "TRUE = callable");
        FIELD(callDate, "call date");
        FIELD_MAKE_OPTIONAL(callDate);
        FIELD(callDateOffset, "call date offset");
        FIELD_MAKE_OPTIONAL(callDateOffset);
        FIELD(discount, "identifies discount curve");
        FIELD(payStream, "payment stream");
        FIELD(fees, "fees");
        FIELD_MAKE_OPTIONAL(fees);
		FIELD(ccyTreatment, "ccyTreatment");
        FIELD_MAKE_OPTIONAL(ccyTreatment);					
		FIELD(underlyingCurve, "underlyingCurve");	
        FIELD_MAKE_OPTIONAL(underlyingCurve);
		FIELD(fxAsset, "fx");
		FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(fxAsset);
    }

    static IObject* defaultLiborStream(){
        return new LiborStream();
    }
};

CClassConstSP const LiborStream::TYPE = CClass::registerClassLoadMethod(
    "LiborStream", typeid(LiborStream), LiborStreamHelper::load);
bool  LiborStreamLoad() {
    return (LiborStream::TYPE != 0);
   }

   

DRLIB_END_NAMESPACE
