//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CashFlow.cpp
//
//   Description : Well, it's a cashflow
//
//   Author      : Andrew J Swain
//
//   Date        : 7 February 2001
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#define QLIB_CASHFLOW_CPP
#include "edginc/CashFlow.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Maths.hpp"
#include <sstream>

DRLIB_BEGIN_NAMESPACE

INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<CashFlowArray _COMMA_ CashFlowArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<CashFlowList>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<CashFlowCluster>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<CashFlowCluster>);

CashFlow::CashFlow(const DateTime& date, double amount): 
    CObject(TYPE), date(date), amount(amount) {
    // done
}

CashFlow::CashFlow() : CObject(TYPE), amount(0.0) {}

CashFlow::~CashFlow() {}

//// Returns true if dates are equal and if amounts are equal using
//// Maths::isZero(. - .). Method added to support instantiating array template
bool CashFlow::operator==(const CashFlow& rhs) const{
    return (date == rhs.date && Maths::isZero(amount - rhs.amount));
}

class CashFlowHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CashFlow, clazz);
        SUPERCLASS(CObject);
        clazz->enableCloneOptimisations();
        EMPTY_SHELL_METHOD(defaultCashFlow);
        FIELD(date, "date");
        FIELD(amount, "amount");
    }
    static IObject* defaultCashFlow(){
        return new CashFlow();
    }
};

/** Validates that dates within the cashflow are in increasing
    order. If failIfEmpty is true, then throws exception if cfArray is
    empty. The description string is used in the exception messages
    should describe what the cashflow represents */
void CashFlow::ensureDatesIncreasing(const CashFlowArray& cfArray,
                                     const string&        description,
                                     bool                 failIfEmpty){
    static const string routine("CashFlow::ensureDatesIncreasing");
    if (failIfEmpty && cfArray.empty()){
        throw ModelException(routine, description + " (cashflow array) is "
                             "empty");
    }
    for (int i = 1; i < cfArray.size(); i++){
        if (cfArray[i-1].date.isGreater(cfArray[i].date)){
            throw ModelException(routine, "Dates in " + description + 
                                 " (array of cashflows) are not increasing: [" + 
                                 Format::toString(i) + "] = " + cfArray[i-1].date.toString() + 
                                 " is after [" + Format::toString(i+1) + "] = " +
                                 cfArray[i].date.toString() + "]");
        }
    }
}

/** Aggragates cashflows that occur at the same date and time. Modifies
    existing cashflow (so could end up with shorter array) */
/* TODO: this function has O(n^2) worst case time complexity (erase from the middle of a vector is costly)
   rewrite it similar to remove() algorithm to have guaranteed complexity.
   For now fix the begin()+1 error, by guaranteeing that begin()+1 is <= end() */
void CashFlow::aggregate(CashFlowArray& cfArray){
    if (cfArray.size() < 2)
        return;
    for (vector<CashFlow>::iterator iter = cfArray.begin() + 1; 
         iter < cfArray.end(); /* ++ in loop body */){
        if (iter->date.equals((iter-1)->date)){
            (iter-1)->amount += iter->amount;
            iter = cfArray.erase(iter);
        } else {
            ++iter;
        }
    }
}

CClassConstSP const CashFlow::TYPE = CClass::registerClassLoadMethod(
    "CashFlow", typeid(CashFlow), CashFlowHelper::load);


/** helper functions for sorting */
bool CashFlow::lessThenForDates(const CashFlow& x, const CashFlow& y)
{
    return (x.date < y.date);
}
bool CashFlow::lessThenForAmounts(const CashFlow& x, const CashFlow& y)
{
    return (x.amount < y.amount);
}

/** Simple Linear Interpolation routine */
double CashFlow::interpolate(const ExpiryArray& exp,
                             const DoubleArray& vals,
                             const DateTime& today,
                             const DateTime& date,
                             const bool extendFlat)
{
    static const string method("CashFlow::interpolate");
    try
    {
        int sizeExp = exp.size();
        int sizeVal = vals.size();
        if (sizeExp != sizeVal)
        {
            throw ModelException(method, "Expiry Array and Value Array must be of equal size");
        }
        
        if (sizeExp < 1) {
            throw ModelException(method, "Expiry Array is empty");
        }
        
        int dt = date.getDate();
        double rate;

        if (sizeExp == 1){ 
            rate = vals[0];
        }else if (exp[0]->toDate(today).getDate() >= dt) {
            // Do *flat* extrapolation only when going backwards.
            rate = vals[0];
        } else  if (extendFlat && exp[sizeExp - 1]->toDate(today).getDate() <= dt) {
            // extrapolate flat off end of zero curve
            rate = vals[sizeExp - 1];
        } else {
            int  lo;
            int  hi;
            
            // Do a binary search to find the lower and upper bounds
            int mid;
            lo = 0;
            hi = sizeExp -1;
            
            while ((hi - lo) > 1) {
                mid = (hi + lo) >> 1;  // compute a mid point
                if (dt >= exp[mid]->toDate(today).getDate()) {
                    lo = mid;
                }
                else {
                    hi = mid;
                }
            }

            if (exp[lo]->toDate(today).getDate() == dt) {
                rate = vals[lo];
            }
            else if (exp[hi]->toDate(today).getDate() == dt) {
                rate = vals[hi];
            }
            else {
                // mimic the year fraction method of Actual365F
                double hi_lo = (exp[hi]->toDate(today).daysDiff(exp[lo]->toDate(today))/365.0);
                double dt_lo = (date.daysDiff(exp[lo]->toDate(today))/365.0);
                
                rate = vals[lo]+((vals[hi]-vals[lo])/hi_lo)*dt_lo;
            }
        }
        return rate;

    }
    catch (exception &e) 
    {
        throw ModelException(e, method);
    }
}

/** Interpolate a rate on a date, boolean gives option of extending 
    flat off the curve or linearly interpolating */
double CashFlow::interpolate(const CashFlowArray& cf,
                             const DateTime & date,
                             const bool extendFlat)
{
    static const string method("CashFlow::interpolate");
    try
    {
        int size = cf.size();
        if (size < 1) {
            throw ModelException(method, "no points in CashFlow");
        }
        int    dt = date.getDate();
        double rate;
        
        if (size == 1){ // can size be different from getLength() ? Yes
            rate = cf[0].amount;
        } else if (cf[0].date.getDate() >= dt) {
            // Do *flat* extrapolation only when going backwards.
            rate = cf[0].amount;
        } else  if (extendFlat && cf[size-1].date.getDate() <= dt) {
            // extrapolate flat off end of zero curve
            rate = cf[size-1].amount;
        } else {
            int  lo;
            int  hi;

            // Do a binary search to find the lower and upper bounds
            int mid;
            lo = 0;
            hi = size -1;
            
            while ((hi - lo) > 1) {
                mid = (hi + lo) >> 1;  // compute a mid point
                if (dt >= cf[mid].date.getDate()) {
                    lo = mid;
                }
                else {
                    hi = mid;
                }
            }
            
            if (cf[lo].date.getDate() == dt) {
                rate = cf[lo].amount;
            }
            else if (cf[hi].date.getDate() == dt) {
                rate = cf[hi].amount;
            }
            else {
                // mimic the year fraction method of Actual365F
                double hi_lo = (cf[hi].date.daysDiff(cf[lo].date)/365.0);
                double dt_lo = (date.daysDiff(cf[lo].date)/365.0);
                
                rate = cf[lo].amount+((cf[hi].amount-cf[lo].amount)/hi_lo)*dt_lo;
            }
        }
        return rate;
    }
    catch (exception &e) 
    {
        throw ModelException(e, method);
    }
}

/** This one requires an exact match on 'date' */
double CashFlow::interpolate(const CashFlowArray& cf,
                             const DateTime& date)
{
    static const string method("CashFlow::interpolate");
    try
    {
        int size = cf.size();
        if (size < 1) {
            throw ModelException(method, "no points in CashFlowArray");
        }
        double rate;
        
        if (size == 1){ // can size be different from getLength() ? Yes
            if (cf[0].date != date) {
                throw ModelException(method, "No value at date " + date.toString());
            }
            rate = cf[0].amount;
        } else if (date < cf[0].date ||
                   date > cf.back().date) {
            throw ModelException(method, "Date " + date.toString() + 
                                 "is out of bounds: " + cf[0].date.toString() +
                                 " to " + cf.back().date.toString());
        } else {
            int  lo;
            int  hi;

            // Do a binary search to find the lower and upper bounds
            int mid;
            lo = 0;
            hi = size -1;
            
            while ((hi - lo) > 1) {
                mid = (hi + lo) >> 1;  // compute a mid point
                if (date >= cf[mid].date) {
                    lo = mid;
                }
                else {
                    hi = mid;
                }
            }
            
            if (cf[lo].date == date) {
                rate = cf[lo].amount;
            }
            else if (cf[hi].date == date) {
                rate = cf[hi].amount;
            }
            else {
                throw ModelException(method, "No entry found for " + date.toString());
            }
        }
        return rate;
    }
    catch (exception &e) 
    {
        throw ModelException(e, method);
    }
}


/** Returns the dates inside the CashFlowArray */
DateTimeArray CashFlow::dates(const CashFlowArray& cf){
    DateTimeArray dates(cf.size());
    for (int i = 0; i < dates.size(); i++){
        dates[i] = cf[i].date;
    }
    return dates;
}

/** Returns the amounts inside the CashFlowArray */
CDoubleArraySP CashFlow::amounts(const CashFlowArray& cf) {
    CDoubleArraySP amounts(new CDoubleArray(cf.size()));
    for (int i = 0; i < amounts->size(); i++){
        (*amounts)[i] = cf[i].amount;
    }
    return amounts;
}

/** Constructs a CashFlowArray from dates and amounts arrays */
CashFlowArraySP CashFlow::createCashFlowArray(
    const DateTimeArray& dates,
    const CDoubleArray& amounts)
{
    if (dates.size() != amounts.size()) {
        throw ModelException(
            "CashFlow::createCashFlowArray",
             "Attempting to create a cash flow array with "
             + Format::toString(dates.size())
             + " date(s) and "
             + Format::toString(amounts.size())
             + " amount(s).");        
    }
    
    CashFlowArraySP result(new CashFlowArray(dates.size()));
    for(int i=0;i<dates.size();i++) {
        (*result)[i] = CashFlow(dates[i], amounts[i]);
    }
    
    return result;
}

/** Cooperate with strings */
/* Format is "DateTime::toString(),Amount" */
CashFlowArray CashFlow::fromStringArray(const StringArray& asStrings) {
    static const string method("CashFlow::fromStringArray");
    CashFlowArray cfa(0);
    for(int i=0; i<asStrings.size(); i++) {
        const string& str = asStrings[i];
        string dateString;
        string timeString("EOD"); // default - this piece is optional
        double amt;
        const char* pos = str.c_str();
        int paramId = 0;
        while (*pos) {
            char   c;
            // skip whitespace
            while ((c = *pos) == ' ' || c == '\t'){
                pos++;
            }
            if (c == ','){
                // skip over our chosen separator and continue
                pos++;
                // count params
                paramId++;
                if (paramId>1) {
                    throw ModelException(method, "Too many commas in " + str);
                }
            } else {
                // next contiguous block is the text we want
                const char* varBegin = pos;
                do{
                    pos++; /* Get another character. */
                    c = *pos;
                } while (c != '\0' && c != ',' && c != ' ' && c != '\t');
                const char* varEnd = pos;
                switch (paramId) {
                case 0:
                    if (dateString.empty()) {
                        // first time through get the date
                        dateString = string(varBegin, varEnd - varBegin);
                    } else {
                        // if there's a next time, it'll be the time
                        timeString = string(varBegin, varEnd - varBegin);
                    }
                    break;
                case 1:
                {
                    string amts(varBegin, varEnd - varBegin);
                    char*   endPos;
                    amt = strtod(amts.c_str(), &endPos);
                    if (*endPos != '\0') {
                        throw ModelException(method, 
                                             "Badly formed number : str");
                    }
                }
                break;
                default:
                    throw ModelException(method, 
                                         "Badly formed CashFlow string representation : " + str);
                }
            } 
        }
        if (paramId!=1) {
            throw ModelException(method, 
                                 "Badly formed CashFlow string representation : " + str);
        }
        DateTime aDate = DateTime(dateString, timeString);
        cfa.push_back(CashFlow(aDate, amt));
    }
    // XXX ?? should we ensureIncreasing?
    return cfa;
}
StringArray CashFlow::toStringArray(const CashFlowArray& cfa) {
    StringArray my(0);
    for(int i=0; i<cfa.size(); i++) {
        string dt(cfa[i].date.toString());
        ostringstream ost;
        ost << cfa[i].amount;
        string amt(ost.str()); 
        string cfstr(dt + "," + amt);
        my.push_back(cfstr);
    }
    return my;
}

/** Accumulates the cashflows in the CashFlowArraySP passed in, so that 
    the cashflows in the returned CashFlowArraySP are:
    cf_out[i].amount = SUM_j{cf_in[j].amount} for all j <= i */
CashFlowArraySP CashFlow::accumulateCashFlows(CashFlowArraySP inp) {
    CashFlowArraySP cumulativeCashFlows(inp.clone());
    int numCashFlows = cumulativeCashFlows->size();
    double cumulativeAmount = 0.0;
    for (int i=0; i < numCashFlows; ++i) {
        // Add cumulativeAmount to (*cumulativeCashFlows)[i].amount and store 
        // that as the new cumulativeAmount
        cumulativeAmount = 
            ((*cumulativeCashFlows)[i].amount += cumulativeAmount);
    }
    return cumulativeCashFlows;
}



/** specialisations of arrayObjectCast */
/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<CashFlow>::toIObject(
    const CashFlow& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<CashFlow>::toIObject(CashFlow& value){
    return CashFlowSP::attachToRef(&value);
}

/** Turns the IObjectSP into a CashFlow */
const CashFlow& arrayObjectCast<CashFlow>::fromIObject(IObjectSP& value){
    CashFlow *cfPtr = dynamic_cast<CashFlow *>(value.get());
    if (!cfPtr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " a CashFlow");
    }
    return *cfPtr;
}

// explicit clone for cashflow arrays - for performance
IObject* arrayClone<CashFlow>::clone(const CArray* arrayToClone){
    const CashFlowArray& theArray = 
        static_cast<const CashFlowArray&>(*arrayToClone);
    return new CashFlowArray(theArray);
}

//template<> CClassConstSP const CashFlowArray::TYPE = CClass::registerClassLoadMethod(
//    "CashFlowArray", typeid(CashFlowArray), load);
DEFINE_TEMPLATE_TYPE(CashFlowArray);

/** specialisations of arrayObjectCast for CashFlowCluster */
/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<CashFlowArray>::toIObject(
    const CashFlowArray& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<CashFlowArray>::toIObject(CashFlowArray& value){
    return IObjectSP::attachToRef(&value);
}

/** Turns the IObjectSP into a CashFlowArray */
const CashFlowArray& arrayObjectCast<CashFlowArray>::fromIObject(
    IObjectSP& value){
    CashFlowArray *dtPtr = dynamic_cast<CashFlowArray *>(value.get());
    if (!dtPtr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " a CashFlowArray");
    }
    return *dtPtr;
}

// explicit clone for arrays of cashflow arrays - for performance
IObject* arrayClone<CashFlowArray>::clone(const CArray* arrayToClone){
    const CashFlowCluster& theCluster = 
        static_cast<const CashFlowCluster&>(*arrayToClone);
    return new CashFlowCluster(theCluster);
}

//template<> CClassConstSP const CashFlowCluster::TYPE = 
//CClass::registerClassLoadMethod(
//    "CashFlowArrayArray", typeid(CashFlowCluster), load);
DEFINE_TEMPLATE_TYPE_WITH_NAME("CashFlowArrayArray", CashFlowCluster);

class CashFlow::CFAddin: public CObject{
    static CClassConstSP const TYPE;

    /** addin takes two parameters  */
    DateTimeArraySP     dates;
    DoubleArraySP       amounts;

    /** the 'addin function' - builds array of correct type */
    static IObjectSP createCFArray(CFAddin* params){
        if (params->dates->size() != params->amounts->size()){
            throw ModelException("CFAddin::createCFArray", "Inconsistent "
                                 "number of dates (" + 
                                 Format::toString(params->dates->size()) + 
                                 ") and amounts (" + 
                                 Format::toString(params->amounts->size())
                                 + ")");
        }
        CashFlowArraySP cfArray(new CashFlowArray(params->dates->size()));
        for (int i = 0; i < params->dates->size(); i++){
            (*cfArray)[i].date = (*params->dates)[i];
            (*cfArray)[i].amount = (*params->amounts)[i];
        }
        return cfArray;
    }
    /** for reflection */
    CFAddin():  CObject(TYPE){}

    static IObject* defaultCFAddin(){
        return new CFAddin();
    }
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CFAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCFAddin);
        FIELD(dates, "Dates");
        FIELD(amounts, "Amounts");

        Addin::registerClassObjectMethod("CASHFLOW_ARRAY",
                                         Addin::UTILITIES,
                                         "Creates a CashFlowArray",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)createCFArray);
                                         
    }
};
    
CClassConstSP const CashFlow::CFAddin::TYPE = CClass::registerClassLoadMethod(
    "CashFlow::CFAddin", typeid(CFAddin), load);

class CashFlow::GetCFAddin: public CObject{
    static CClassConstSP const TYPE;

    /** addin takes one parameter */
    CashFlowArraySP cashFlows;

    /** the 'addin function' - builds arrays of dates and amounts */
    static IObjectSP getCashFlows(GetCFAddin* params){
        
        ObjectArraySP output(new ObjectArray(2)); // 2 columns
        
        int cfSize = params->cashFlows->size();
        DateTimeArraySP dates(new DateTimeArray(cfSize));
        DoubleArraySP amounts(new DoubleArray(cfSize));
        for (int i = 0; i < cfSize; i++){
            (*dates)[i] = (*params->cashFlows)[i].date;
            (*amounts)[i] = (*params->cashFlows)[i].amount;
        }
        (*output)[0] = dates;
        (*output)[1] = amounts;
        
        return output;
    }
    
    /** for reflection */
    GetCFAddin():  CObject(TYPE){}

    static IObject* defaultGetCFAddin(){
        return new GetCFAddin();
    }
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(GetCFAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetCFAddin);
        FIELD(cashFlows, "Cash Flow array");

        Addin::registerClassObjectMethod("GET_CASHFLOWS",
                                         Addin::UTILITIES,
                                         "Creates a date and an amount array from a CashFlowArray",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)getCashFlows);
    }
};
    
CClassConstSP const CashFlow::GetCFAddin::TYPE = CClass::registerClassLoadMethod(
    "CashFlow::GetCFAddin", typeid(GetCFAddin), load);

CashFlowArraySP CashFlow::merge(const CashFlowArrayConstSP x, 
                                const CashFlowArrayConstSP y) {

    static const string method = "CashFlow::merge";
    try {
        CashFlowArraySP merged(new CashFlowArray(0));

        if (!x) {
            if (!y) {
                return merged; // 0 length array
            }
            else {
                return CashFlowArraySP(y.clone());
            }
        }
        if (!y) {
            return CashFlowArraySP(x.clone());
        }
        
        ensureDatesIncreasing(*x, "1st cashflows", false); 
        ensureDatesIncreasing(*y, "2nd cashflows", false); 

        int ix = 0; // iterate over x
        int iy = 0; // iterate over y

        while (ix < x->size() && iy < y->size()) {
            const DateTime& xDate = (*x)[ix].date;
            const DateTime& yDate = (*y)[iy].date;

            if (xDate < yDate) {
                merged->push_back((*x)[ix++]);
            }
            else if (xDate > yDate) {
                merged->push_back((*y)[iy++]);
            }
            else {
                merged->push_back((*x)[ix]);
                (*merged)[merged->size()-1].amount = (*x)[ix++].amount+
                                                     (*y)[iy++].amount;
            }
        }
                
        while (ix < x->size()) {
            merged->push_back((*x)[ix++]);
        }

        while (iy < y->size()) {
            merged->push_back((*y)[iy++]);
        }
                        
        return merged;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


CashFlowList::CashFlowList(const CashFlowArray* cfl): 
    CObject(TYPE), cfl(copy(cfl)) {}

/** scale by factor x */
void CashFlowList::scale(double x) {
    for (int i = 0; i < cfl->size(); i++) {
        (*cfl)[i].amount *= x;
    }
}

/** add CashFlowList to this result (Implementation of
    CombinableResult) */
void CashFlowList::add(const CombinableResult& x, double scaleFactor) {
    // gcc bug: force to IObject before dynamic cast
    const CashFlowList& toAdd= dynamic_cast<const CashFlowList&>(static_cast<const IObject&>(x));

    // handle adding arrays of different lengths & dates
    CashFlowListSP rhs(copy(&toAdd));
    rhs->scale(scaleFactor);

    cfl = CashFlow::merge(cfl, rhs->cfl);
}

CashFlowList::CashFlowList() : CObject(TYPE) {}

class CashFlowListHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CashFlowList, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(CombinableResult);
        EMPTY_SHELL_METHOD(defaultCashFlowList);
        FIELD(cfl, "cashflows");
    }
    static IObject* defaultCashFlowList(){
        return new CashFlowList();
    }
};

CClassConstSP const CashFlowList::TYPE = CClass::registerClassLoadMethod(
    "CashFlowList", typeid(CashFlowList), CashFlowListHelper::load);


DRLIB_END_NAMESPACE
