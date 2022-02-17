//----------------------------------------------------------------------------
//
//   Group       : GCCT Derivatives Research
//
//   Filename    : ExpiryResultss.cpp
//
//   Description : Captures double results of a result at a given expiry, ie stores
//                 an expiry and a double array. Strickly based on ExpiryResultss
//
//   Author      : Sean Chen
//
//   Date        : 17 May 2005
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/ExpiryResults.hpp"


DRLIB_BEGIN_NAMESPACE

/** Captures results of a result at a given expiry, ie stores an expiry
    and a double array */

ExpiryResults::~ExpiryResults(){}

/* for simplicity allow array of ExpiryResultss to be an array of
   structures rather than an array of pointers - dictates public
   default constructor */
ExpiryResults::ExpiryResults(): CObject(TYPE){}

ExpiryResults::ExpiryResults(const ExpiryConstSP& expiry, const CDoubleArray& result): 
    CObject(TYPE), expiry(copy(expiry.get())), result(result) {}


/** returns the expiry associated with this result */
ExpiryConstSP ExpiryResults::getExpiry() const
{
    return expiry;
}

/** returns the double detailing this result */
CDoubleArray ExpiryResults::getResult() const
{
    return result;
}

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void ExpiryResults::outputWrite(const string& linePrefix,
                               const string& prefix,
                               ostream&      stream) const{
    ios::fmtflags oldFlags = stream.flags();     // save settings
    long               oldPrecision = stream.precision(); // save settings
    if (result[0] < 1 && result[0] > -1){
        // 8 dp
        stream.flags(ios::fixed);
    } else {
        // 8 sf
        stream.unsetf(ios::fixed);
    }        
    stream.precision(8);

    stream << linePrefix << prefix << "_" << expiry->toString();
	for (int i=0; i<result.size(); i++)
           stream << ": " << result[i];
	stream << endl;
    stream.flags(oldFlags); // restore settings
    stream.precision(oldPrecision);
}

/** scale by factor x */
void ExpiryResults::scale(double x) {
    for (int i=0; i<result.size(); i++)
		result[i] *= x;
}

/** add ExpiryResults object to this result (Implementation of
    CombinableResult) */
void ExpiryResults::add(const CombinableResult& x, double scaleFactor) {
    // gcc bug: force to IObject before dynamic cast
    const ExpiryResults& ex = dynamic_cast<const ExpiryResults&>((IObject&)x);
    if (!expiry->equals(ex.expiry.get())) {
        // shouldn't happen really
        throw ModelException("ExpiryResults::add", 
                             "can't add " + ex.expiry->toString() +
                             " to " + expiry->toString());
    }

	if (result.size() != ex.result.size()) 
	{
        // shouldn't happen really
        throw ModelException("ExpiryResults::add", 
                             "can't add " + ex.expiry->toString() +
                             " to " + expiry->toString() + " element size different");
    }
	for (int i=0; i<result.size(); i++)
		result[i] += scaleFactor * (ex.result)[i];           
}

/** adds the supplied value to this object */
void ExpiryResults::addToResult(const CDoubleArray& value)
{
	if (result.size() != value.size()) 
	{
        // shouldn't happen really
        throw ModelException("ExpiryResults::addToResult", 
                             "can't add value to result: element size different");
    }

	for (int i=0; i<result.size(); i++)
		result[i] += value[i];           
}

/** specialisations of arrayObjectCast */
/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<ExpiryResults>::toIObject(
    const ExpiryResults& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<ExpiryResults>::toIObject(ExpiryResults& value){
    return IObjectSP::attachToRef(&value);
}

/** Turns the IObjectSP into a ExpiryResults */
ExpiryResults arrayObjectCast<ExpiryResults>::fromIObject(IObjectSP& value){
    ExpiryResults *ptr = dynamic_cast<ExpiryResults *>(value.get());
    if (!ptr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " an ExpiryResults");
    }
    return *ptr;
}

class ExpiryResultsHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ExpiryResults, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(CombinableResult);
        EMPTY_SHELL_METHOD(defaultExpiryResults);
        FIELD(expiry, "Indicates which point was tweaked ");
        FIELD(result, "The result of the expiry tweak");
    }
    static IObject* defaultExpiryResults(){
        return new ExpiryResults();
    }
};

CClassConstSP const ExpiryResults::TYPE = CClass::registerClassLoadMethod(
    "ExpiryResults", typeid(ExpiryResults), ExpiryResultsHelper::load);

DEFINE_TEMPLATE_TYPE(ExpiryResultsArray);

DRLIB_END_NAMESPACE
