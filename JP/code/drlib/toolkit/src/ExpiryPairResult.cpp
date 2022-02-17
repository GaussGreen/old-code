/**
 * @file ExpiryPairResult.cpp
 */

#include "edginc/config.hpp"
#define QLIB_EXPIRYPAIRRESULT_CPP
#include "edginc/ExpiryPairResult.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

/** Captures result of a result at a given expiry pair, ie stores an expiry pair
    and a double */

ExpiryPairResult::~ExpiryPairResult(){}

/* for simplicity allow array of ExpiryPairResults to be an array of
   structures rather than an array of pointers - dictates public
   default constructor */
ExpiryPairResult::ExpiryPairResult(): CObject(TYPE), result(0){}

ExpiryPairResult::ExpiryPairResult(const ExpiryPairConstSP& expiryPair, double result): 
    CObject(TYPE), expiryPair(copy(expiryPair.get())), result(result) {}

    
/** returns the expiry associated with this result */
ExpiryPairConstSP ExpiryPairResult::getExpiryPair() const{
    return expiryPair;
}

/** returns the double detailing this result */
double ExpiryPairResult::getResult() const{
    return result;
}

//// Returns true if expiries are equal and if "result"s are equal
//// using Maths::isZero(. - .). Method added to support
//// instantiating array template
bool ExpiryPairResult::operator==(const ExpiryPairResult& rhs) const{
    return (expiryPair->equals(rhs.expiryPair.get()) &&
            Maths::isZero(result-rhs.result));
}

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void ExpiryPairResult::outputWrite(const string& linePrefix,
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

    stream << linePrefix << prefix << "_" << expiryPair->toString() 
           << ": " << result << endl;
    stream.flags(oldFlags); // restore settings
    stream.precision(oldPrecision);
}

/** scale by factor x */
void ExpiryPairResult::scale(double x) {
    result *= x;
}

/** add ExpiryPairResult object to this result (Implementation of
    CombinableResult) */
void ExpiryPairResult::add(const CombinableResult& x, double scaleFactor) {
    // gcc bug: force to IObject before dynamic cast
    const ExpiryPairResult& ex = dynamic_cast<const ExpiryPairResult&>(static_cast<const IObject&> (x));
    if (!expiryPair->equals(ex.expiryPair.get())) {
        // shouldn't happen really
        throw ModelException("ExpiryPairResult::add", 
                             "can't add " + ex.expiryPair->toString() +
                             " to " + expiryPair->toString());
    }
    result += scaleFactor * ex.result;           
}

/** adds the supplied value to this object */
void ExpiryPairResult::addToResult(double value){
    result += value;
}

/** Constructs an ExpiryPairResultArray from dates and amounts arrays */
ExpiryPairResultArraySP ExpiryPairResult::createExpiryPairResultArray(
    const ExpiryPairArray& expiries,
    const CDoubleArray& amounts)
{
    if (expiries.size() != amounts.size()) {
        throw ModelException(
            "ExpiryPairResult::createExpiryPairResultArray",
             "Attempting to create a cash flow array with "
             + Format::toString(expiries.size())
             + " date(s) and "
             + Format::toString(amounts.size())
             + " amount(s).");        
    }
    
    ExpiryPairResultArraySP result(new ExpiryPairResultArray(expiries.size()));
    for(int i=0;i<expiries.size();i++) {
        (*result)[i] = ExpiryPairResult(expiries[i], amounts[i]);
    }
    
    return result;
}


/** Returns index of the expiry in supplied ExpiryPairResultArray.
    The supplied base date is used in the comparison so that it is 
    possible to search for eg MaturityPeriods in ExpiryPairResultArrays 
    were the expiries are BenchmarkDates */
int ExpiryPairResult::search(const DateTime& base,
                         const ExpiryPairConstSP expiryPair,
                         const ExpiryPairResultArraySP expiryResults)
{
    static const string routine = "ExpiryPairResult::search(base...)";

    const CClassConstSP expiryClass = expiryPair->getClass();
    const DateTime& expiryDate = expiryPair->optExpiry->toDate(base);
    ExpiryPairSP expiryIdx;

    for (int idx = 0; idx < expiryResults->size(); idx++){
        expiryIdx = (*expiryResults)[idx].expiryPair; // for ease
        if (expiryClass == expiryIdx->getClass()) {
            if (expiryPair->equals(expiryIdx.get())) {
                return idx;
            }
        }
        else {
            if (expiryDate.equals(expiryIdx->optExpiry->toDate(base))) {
                return idx;
            }
        }
    }
    throw ModelException(routine, 
                         "Expiry ("+expiryPair->toString()+") not found");   
}

/** Reports if the supplied ExpiryPairArray is a subset of the expiries in the
    ExpiryPairResultArray - NB: not the other way around!
    The supplied base date is used in the comparison so that it is 
    possible to compare eg MaturityPeriods and ExpiryPairResultArrays 
    were the expiries are BenchmarkDates */
int ExpiryPairResult::isSubset(const DateTime& base,
                           const ExpiryPairArrayConstSP part,
                           const ExpiryPairResultArraySP full)
{
    int pSize = part->size();
    int fSize = full->size();
    int f = 0;
    for (int p=0; p<pSize; p++) {
        ExpiryPairSP partExp = (*part)[p];
        // Search all elements of part in full sequentially (ie, they must be
        // in the same order - typically, in increasing order)
        for (/*f*/; f<fSize; f++) {
            ExpiryPairSP fullExp = (*full)[f].expiryPair;        

            if (partExp->getClass() == fullExp->getClass()) {
                if (partExp->equals(fullExp.get())) {
                    break; // Get out of the for loop in f
                }
                // else skip expiries in full which are not present in part
            }
            else {
                const DateTime& fDate = fullExp->optExpiry->toDate(base);
                const DateTime& pDate = partExp->optExpiry->toDate(base);
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
IObjectConstSP arrayObjectCast<ExpiryPairResult>::toIObject(
    const ExpiryPairResult& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<ExpiryPairResult>::toIObject(ExpiryPairResult& value){
    return IObjectSP::attachToRef(&value);
}

/** Turns the IObjectSP into a ExpiryPairResult */
const ExpiryPairResult& arrayObjectCast<ExpiryPairResult>::fromIObject(
    IObjectSP& value){
    ExpiryPairResult *ptr = dynamic_cast<ExpiryPairResult *>(value.get());
    if (!ptr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " an ExpiryPairResult");
    }
    return *ptr;
}

class ExpiryPairResultHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ExpiryPairResult, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(CombinableResult);
        EMPTY_SHELL_METHOD(defaultExpiryPairResult);
        FIELD(expiryPair, "Indicates which pair point was tweaked ");
        FIELD(result, "The result of the expiry tweak");
    }
    static IObject* defaultExpiryPairResult(){
        return new ExpiryPairResult();
    }
};

CClassConstSP const ExpiryPairResult::TYPE = CClass::registerClassLoadMethod(
    "ExpiryPairResult", typeid(ExpiryPairResult), ExpiryPairResultHelper::load);

DEFINE_TEMPLATE_TYPE(ExpiryPairResultArray);

DRLIB_END_NAMESPACE
