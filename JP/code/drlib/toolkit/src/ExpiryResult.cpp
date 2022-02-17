//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ExpiryResult.cpp
//
//   Description : Captures result of a result at a given expiry, ie stores
//                 an expiry and a double
//
//   Author      : Mark A Robson
//
//   Date        : 5 March 2001
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#define QLIB_EXPIRYRESULT_CPP
#include "edginc/ExpiryResult.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

/** Captures result of a result at a given expiry, ie stores an expiry
    and a double */

ExpiryResult::~ExpiryResult(){}

/* for simplicity allow array of ExpiryResults to be an array of
   structures rather than an array of pointers - dictates public
   default constructor */
ExpiryResult::ExpiryResult(): CObject(TYPE), result(0){}

ExpiryResult::ExpiryResult(const ExpiryConstSP& expiry, double result): 
    CObject(TYPE), expiry(copy(expiry.get())), result(result) {}


/** returns the expiry associated with this result */
ExpiryConstSP ExpiryResult::getExpiry() const{
    return expiry;
}

/** returns the double detailing this result */
double ExpiryResult::getResult() const{
    return result;
}

//// Returns true if expiries are equal and if "result"s are equal
//// using Maths::isZero(. - .). Method added to support
//// instantiating array template
bool ExpiryResult::operator==(const ExpiryResult& rhs) const{
    return (expiry->equals(rhs.expiry.get()) &&
            Maths::isZero(result-rhs.result));
}

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void ExpiryResult::outputWrite(const string& linePrefix,
                               const string& prefix,
                               ostream&      stream) const{
    ios::fmtflags oldFlags = stream.flags();     // save settings
    long               oldPrecision = stream.precision(); // save settings
    if (result < 1 && result > -1){
        // 8 dp
        stream.flags(ios::fixed);
    } else {
        // 8 sf
        stream.unsetf(ios::fixed);
    }        
    stream.precision(8);

    stream << linePrefix << prefix << "_" << expiry->toString() 
           << ": " << result << endl;
    stream.flags(oldFlags); // restore settings
    stream.precision(oldPrecision);
}

/** scale by factor x */
void ExpiryResult::scale(double x) {
    result *= x;
}

/** add ExpiryResult object to this result (Implementation of
    CombinableResult) */
void ExpiryResult::add(const CombinableResult& x, double scaleFactor) {
    // gcc bug: force to IObject before dynamic cast
    const ExpiryResult& ex = dynamic_cast<const ExpiryResult&>(static_cast<const IObject&> (x));
    if (!expiry->equals(ex.expiry.get())) {
        // shouldn't happen really
        throw ModelException("ExpiryResult::add", 
                             "can't add " + ex.expiry->toString() +
                             " to " + expiry->toString());
    }
    result += scaleFactor * ex.result;           
}

/** adds the supplied value to this object */
void ExpiryResult::addToResult(double value){
    result += value;
}

/** Constructs an ExpiryResultArray from dates and amounts arrays */
ExpiryResultArraySP ExpiryResult::createExpiryResultArray(
    const ExpiryArray& expiries,
    const CDoubleArray& amounts)
{
    if (expiries.size() != amounts.size()) {
        throw ModelException(
            "ExpiryResult::createExpiryResultArray",
             "Attempting to create a cash flow array with "
             + Format::toString(expiries.size())
             + " date(s) and "
             + Format::toString(amounts.size())
             + " amount(s).");        
    }
    
    ExpiryResultArraySP result(new ExpiryResultArray(expiries.size()));
    for(int i=0;i<expiries.size();i++) {
        (*result)[i] = ExpiryResult(expiries[i], amounts[i]);
    }
    
    return result;
}


/** Returns index of the expiry in supplied ExpiryResultArray.
    The supplied base date is used in the comparison so that it is 
    possible to search for eg MaturityPeriods in ExpiryResultArrays 
    were the expiries are BenchmarkDates */
int ExpiryResult::search(const DateTime& base,
                         const ExpiryConstSP expiry,
                         const ExpiryResultArraySP expiryResults)
{
    static const string routine = "ExpiryResult::search(base...)";

    const CClassConstSP expiryClass = expiry->getClass();
    const DateTime& expiryDate = expiry->toDate(base);
    ExpirySP expiryIdx;

    for (int idx = 0; idx < expiryResults->size(); idx++){
        expiryIdx = (*expiryResults)[idx].expiry; // for ease
        if (expiryClass == expiryIdx->getClass()) {
            if (expiry->equals(expiryIdx.get())) {
                return idx;
            }
        }
        else {
            if (expiryDate.equals(expiryIdx->toDate(base))) {
                return idx;
            }
        }
    }
    throw ModelException(routine, 
                         "Expiry ("+expiry->toString()+") not found");   
}

/** Reports if the supplied ExpiryArray is a subset of the expiries in the
    ExpiryResultArray - NB: not the other way around!
    The supplied base date is used in the comparison so that it is 
    possible to compare eg MaturityPeriods and ExpiryResultArrays 
    were the expiries are BenchmarkDates */
int ExpiryResult::isSubset(const DateTime& base,
                           const ExpiryArrayConstSP part,
                           const ExpiryResultArraySP full)
{
    int pSize = part->size();
    int fSize = full->size();
    int f = 0;
    for (int p=0; p<pSize; p++) {
        ExpirySP partExp = (*part)[p];
        // Search all elements of part in full sequentially (ie, they must be
        // in the same order - typically, in increasing order)
        for (/*f*/; f<fSize; f++) {
            ExpirySP fullExp = (*full)[f].expiry;        

            if (partExp->getClass() == fullExp->getClass()) {
                if (partExp->equals(fullExp.get())) {
                    break; // Get out of the for loop in f
                }
                // else skip expiries in full which are not present in part
            }
            else {
                const DateTime& fDate = fullExp->toDate(base);
                const DateTime& pDate = partExp->toDate(base);
                if (fDate.equals(pDate)) {
                    break; // Get out of the for loop in f
                }
                // else skip expiries in full which are not present in part
            }
        }
        if (f == fSize) {
            // It means we did not find part[p] in full, ie, part is not a subset
            return false;
        }
    }
    return true;
}


/** specialisations of arrayObjectCast */
/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<ExpiryResult>::toIObject(
    const ExpiryResult& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<ExpiryResult>::toIObject(ExpiryResult& value){
    return IObjectSP::attachToRef(&value);
}

/** Turns the IObjectSP into a ExpiryResult */
const ExpiryResult& arrayObjectCast<ExpiryResult>::fromIObject(
    IObjectSP& value){
    ExpiryResult *ptr = dynamic_cast<ExpiryResult *>(value.get());
    if (!ptr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " an ExpiryResult");
    }
    return *ptr;
}

class ExpiryResultHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ExpiryResult, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(CombinableResult);
        EMPTY_SHELL_METHOD(defaultExpiryResult);
        FIELD(expiry, "Indicates which point was tweaked ");
        FIELD(result, "The result of the expiry tweak");
    }
    static IObject* defaultExpiryResult(){
        return new ExpiryResult();
    }
};

CClassConstSP const ExpiryResult::TYPE = CClass::registerClassLoadMethod(
    "ExpiryResult", typeid(ExpiryResult), ExpiryResultHelper::load);

//template<> CClassConstSP const ExpiryResultArray::TYPE = CClass::registerClassLoadMethod(
//    "ExpiryResultArray", typeid(ExpiryResultArray), load);
DEFINE_TEMPLATE_TYPE(ExpiryResultArray);

DRLIB_END_NAMESPACE
