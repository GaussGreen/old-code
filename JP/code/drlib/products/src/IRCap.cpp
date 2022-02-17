//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRCap.cpp
//
//   Description : Interest rate cap/floor/straddle
//
//   Author      : Andrew J Swain
//
//   Date        : 18 January 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GenericSimpleIR.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Maths.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/Black.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/VolProcessedBSIR.hpp"
#include "edginc/SwapMaturityVolRequest.hpp"
#include "edginc/SensitiveIRVolPoints.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/ClosedFormIRLN.hpp"

DRLIB_BEGIN_NAMESPACE

class IRCap: public GenericSimpleIR,
             public virtual CClosedFormLN::IIntoProduct,
             public virtual ClosedFormIRLN::IIntoProduct,
             public virtual ISensitiveIRVolPoints,
             public virtual LastSensDate {
public:
    static CClassConstSP const TYPE; 

    virtual void validatePop2Object(){
        static const string method = "IRCap::validatePop2Object";
    }

    virtual void Validate() {
        static const string method = "IRCap::Validate";
        try {
            int i;

            if (!vol.get()) {
                throw ModelException(method, "no IR vol supplied");
            }

            if (Maths::isNegative(capLevel)) {
                throw ModelException(method,
                                     "negative cap/floor level (" + 
                                     Format::toString(capLevel) + ")");
            }

            if (Maths::isNegative(capNotional)) {
                throw ModelException(method,
                                     "negative cap/floor notional (" + 
                                     Format::toString(capNotional) + ")");
            }
 
            if (accrualDates->size() != resetDates->size()+1) {
                throw ModelException(method,
                                     "must have 1 more accrual dates "
                                     "than reset dates");
            }

            if (refixLevels.size() != resetDates->size()) {
                throw ModelException(method,
                                     "must have same number of refix levels"
                                     " as reset dates");
            }

            // check that the resetDates and accrualDates are in ascending order
            for (i = 1; i < resetDates->size(); i++) {
                if ((*resetDates)[i-1].isGreater((*resetDates)[i])) {
                    throw ModelException(method,
                                         "reset dates are not in ascending "
                                         "order (" + 
                                         (*resetDates)[i-1].toString() + "), ("+
                                         (*resetDates)[i].toString() + ")");                    
                }                  
            } 
            for (i = 1; i < accrualDates->size(); i++) {
                if ((*accrualDates)[i-1].isGreater((*accrualDates)[i])) {
                    throw ModelException(method,
                                         "accrual dates are not in ascending "
                                         "order (" + 
                                         (*accrualDates)[i-1].toString() + "), ("+
                                         (*accrualDates)[i].toString() + ")"); 
                }                   
            } 
 
            // check that the accrual dates are >= to the reset dates 
            for (i = 0; i < resetDates->size(); i++) {
                if ((*accrualDates)[i].isLess((*resetDates)[i])) {
                    throw ModelException(method,
                                         "accrual date (" + 
                                         (*accrualDates)[i].toString() + 
                                         ") is before reset date ("+
                                         (*resetDates)[i].toString() + ")"); 
                }                   
            } 
 
            // the refix levels should all be zero unless 
            // it is a valid historic refix
            for (i = 0; i < resetDates->size(); i++) {
                if (Maths::isNegative(refixLevels[i])) {
                    throw ModelException(method,
                                         "refix level at (" + 
                                         (*resetDates)[i].toString() + 
                                         ") is negative (" + 
                                         Format::toString(refixLevels[i])+")");
                }
                else if (Maths::isPositive(refixLevels[i]) && 
                         valueDate.isLess((*resetDates)[i])) {
                    throw ModelException(method,
                                         "non-zero refix levels (" + 
                                         Format::toString(refixLevels[i]) +
                                         ") at (" + 
                                         (*resetDates)[i].toString() + 
                                         ") should only be supplied for "
                                         "historic fixings");
                }
                else if (Maths::isZero(refixLevels[i]) &&
                         valueDate.isGreaterOrEqual((*resetDates)[i])) {
                    throw ModelException(method,
                                         "refix levels must be supplied for "
                                         "reset dates in the past or on the "
                                         "value date - check fixing at (" + 
                                         (*resetDates)[i].toString() + ")");
                }
            }
        }    
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(IRCap, clazz);
        SUPERCLASS(GenericSimpleIR);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(ClosedFormIRLN::IIntoProduct);
        IMPLEMENTS(ISensitiveIRVolPoints);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultIRCap);
        FIELD(resetDates, "resetDates");
        FIELD(accrualDates, "accrualDates");
        FIELD(refixLevels, "refixLevels");
        FIELD(dayCountConv, "dayCountConv");
        FIELD(capType, "capType");
        FIELD(capLevel, "capLevel");
        FIELD(capNotional, "capNotional");
    }

    static IObject* defaultIRCap(){
        return new IRCap();
    }
    
private:
    static const string CAP;
    static const string FLOOR;
    static const string STRADDLE;

    friend class IRCapClosedForm;

    IRCap():GenericSimpleIR(TYPE) {}; 
    IRCap(const IRCap& rhs);
    IRCap& operator=(const IRCap& rhs);

    int capletFrequency(const DateTime& d1, const DateTime& d2) const {
        // capletFreq = (int)(((accEnd - accStart)/365)*12 + 0.5)
        // This equation ensures that the capletFreq is always an exact
        // number of months
        
        int capletStart = d1.getDate();
        int capletEnd   = d2.getDate();
        int capletFreq  = int((capletEnd-capletStart)/365.0*12.0+0.5);
        return capletFreq;
    }

    void requests(Control* control, CResults* results) const {
        static const string method = "CashFlowStream::requests";
        try {
            OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request) {
                DateTimeArraySP dates(new DateTimeArray(resetDates->size()));
                for (int i = 0; i < resetDates->size(); i++) {
                    (*dates)[i] = (*accrualDates)[i+1];
                }
                OutputRequestUtil::recordPaymentDates(control,results,dates.get()); 
            }
        
            request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request) {
                int           i;
                bool          isCall;
                CashFlowArray payments;

                isCall = !CString::equalsIgnoreCase(capType, FLOOR);

                for (i = 0; i < resetDates->size(); i++) {
                    // if  historic fixing can be used
                    if (valueDate.isGreaterOrEqual((*resetDates)[i])) {
                        double premium;
                        if (isCall) {
                            premium = Maths::max(refixLevels[i] - capLevel, 0.0);
                        }
                        else {
                            premium = Maths::max(capLevel - refixLevels[i], 0.0);
                        }
                        payments.push_back(CashFlow((*accrualDates)[i+1], premium));

                        if (CString::equalsIgnoreCase(capType, STRADDLE)) {
                            // add on the put
                            premium = Maths::max(capLevel - refixLevels[i], 0.0);
                            payments.push_back(CashFlow((*accrualDates)[i+1], premium));
                        }                            
                    }
                }
                
                OutputRequestUtil::recordKnownCashflows(control,
                                                        results,
                                                        discount->getCcy(),
                                                        &payments); 
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    void price(Control* control, CResults* results)const{
        static const string method = "IRCap::price";
        
        try  {
            double value = 0.0;
            double premium = 0.0;
            double pv;
            int    i;
            int    j;
            bool   isCall;
            int    loops;

            DayCountConventionSP dcc(DayCountConventionFactory::make(dayCountConv));

            DoubleArray fixings(refixLevels.size());

            if (CString::equalsIgnoreCase(capType, CAP)) {
                isCall = true;
                loops  = 1;
            }
            else if (CString::equalsIgnoreCase(capType, FLOOR)) {
                isCall = false;
                loops  = 1;
            }
            else if (CString::equalsIgnoreCase(capType, STRADDLE)) {
                isCall = true;
                loops  = 2;
            }
            else {
                throw ModelException(method,
                                     "invalid cap type (" + capType + ")");
            }

            // find the caplet frequency from adjacent accrualDates - any two
            // will do.
            int capletFreq  = capletFrequency((*accrualDates)[0], 
                                              (*accrualDates)[1]);

            MaturityPeriodSP tenor(new MaturityPeriod(capletFreq, "M"));

            for (i = 0; i < fixings.size(); i++) {
                fixings[i] = refixLevels[i];
                if (Maths::isZero(fixings[i])) {
                    fixings[i] = discount->fwd((*accrualDates)[i],
                                               (*accrualDates)[i+1],
                                               dcc.get(),
                                               CompoundBasis::SIMPLE);
                }
            }

            for (i = 0; i < loops; i++) {
                for (j = 0; j < resetDates->size(); j++) {
                    // if caplet is active and reset date in the future
                    if (valueDate.isLess((*resetDates)[j])) {
                        pv = discount->pv((*accrualDates)[j+1]);

                        CVolRequestSP volRequest(new SwapMaturityVolRequest(tenor.get()));
                        
                        CVolProcessedSP volCurve(vol->getProcessedVol(volRequest.get(), 0));

                        CVolProcessedBS& volBS = dynamic_cast<CVolProcessedBS&>(*volCurve);
                        
                        double variance = volBS.CalcVar(valueDate, (*resetDates)[j]);
                        
                        premium = Black::price(isCall,
                                               fixings[j],
                                               capLevel,
                                               pv,
                                               variance);
                    }
                    // if caplet active but historic fixing can be used
                    else if (valueDate.isGreaterOrEqual((*resetDates)[j]) &&
                             valueDate.isLess((*accrualDates)[j+1])) {
                        pv = discount->pv((*accrualDates)[j+1]);

                        if (isCall) {
                            premium = Maths::max(fixings[j] - capLevel, 0.0);
                        }
                        else {
                            premium = Maths::max(capLevel - fixings[j], 0.0);
                        }
                        premium *= pv;
                    }
                    else {
                                // the caplet is inactive
                        continue;
                    }

                    premium *= capNotional * dcc->years((*accrualDates)[j],
                                                        (*accrualDates)[j+1]);

                    value += premium;
                }
                if (CString::equalsIgnoreCase(capType, STRADDLE)) {
                    isCall = false;
                }
            }

            results->storePrice(value, discount->getCcy());       

            if (control && control->isPricing() ) {
                requests(control, results);
            }    
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
    
    DateTime endDate(const Sensitivity* sensControl) const {
        return (*accrualDates)[accrualDates->size()-1];
    }

    virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
        OutputNameConstSP outputName,
        const IModel*      model) const {
        static const string method = "IRCap::getSensitiveIRVolPoints";
        try {
            int capletFreq = capletFrequency((*accrualDates)[0],
                                             (*accrualDates)[1]);

            MaturityPeriod tenor(capletFreq, "M");
            SwapMaturityVolRequest volRequest(&tenor);
            return VolProcessedBSIR::sensitiveIRVolPoints(vol.get(), 
                                                          discount.get(),
                                                          &volRequest, 
                                                          *resetDates);

#if 0
            IRGridPointArraySP points(new IRGridPointArray(0));

            int capletFreq = capletFrequency((*accrualDates)[0],
                                             (*accrualDates)[1]);

            MaturityPeriodSP tenor(new MaturityPeriod(capletFreq, "M"));
            CVolRequestSP    volRequest(new SwapMaturityVolRequest(tenor.get()));
            
            vol->sensitiveIRVolPoints(volRequest.get(),
                                      outputName,
                                      *resetDates,
                                      points);

            return points;
#endif
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

/** Implementation of ClosedFormLN::IntoProduct interface */
    CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;
    ClosedFormIRLN::IProduct* createProduct(ClosedFormIRLN* model) const;

private:
    DateTimeArraySP resetDates;        // interest rate reset dates 
    DateTimeArraySP accrualDates;      // accrue start/pay dates 
    DoubleArray     refixLevels;       // refix levels 
    string          dayCountConv;      // day count convention 
    string          capType;           // C = cap, F = floor, S = straddle 
    //(long cap + long floor)
    double          capLevel;          // cap/floor strike level 
    double          capNotional;       // notional 
};

const string IRCap::CAP      = "C";
const string IRCap::FLOOR    = "F";
const string IRCap::STRADDLE = "S";

CClassConstSP const IRCap::TYPE = CClass::registerClassLoadMethod(
    "IRCap", typeid(IRCap), IRCap::load);


/** private class */
class IRCapClosedForm: public virtual CClosedFormLN::IProduct,
                       public virtual ClosedFormIRLN::IProduct {
private:
    const IRCap* cap; // a reference

public:
    IRCapClosedForm(const IRCap* cap): cap(cap){}

    void price(CClosedFormLN* model,
               Control*       control, 
               CResults*      results)const{
        cap->price(control, results);
    }

    void price(ClosedFormIRLN* model,
               Control*        control, 
               CResults*       results) const{
        cap->price(control, results);
    }
};

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* IRCap::createProduct(CClosedFormLN* model) const
{
    return new IRCapClosedForm(this);
}
ClosedFormIRLN::IProduct* IRCap::createProduct(ClosedFormIRLN* model) const
{
    return new IRCapClosedForm(this);
}

// for class loading 
bool IRCapLoad() {
    return (IRCap::TYPE != 0);
}

DRLIB_END_NAMESPACE
