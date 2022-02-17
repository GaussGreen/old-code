//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MaturityTimePeriod.cpp
//
//   Description : Defines floating expiries used to define yield curve & 
//                 vol surface points e.g. 1M, 5Y with a fixed time of day
//
//   Author      : Andrew J Swain
//
//   Date        : 15 February 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include "edginc/Writer.hpp"


DRLIB_BEGIN_NAMESPACE


MaturityTimePeriod::MaturityTimePeriod(const string& period, int time) : 
    Expiry(TYPE), period(new MaturityPeriod(period)), time(time){}

/** Explicit implementation for performance */
IObject* MaturityTimePeriod::clone() const{
    int& count = getRefCount();
    if (count == 0){
        return new MaturityTimePeriod(*this);
    } else {
        return  const_cast<IObject*>((const IObject*)this);
    }
}

MaturityTimePeriod::MaturityTimePeriod(const MaturityTimePeriod& matTime): 
    Expiry(TYPE), period(matTime.period.clone()), time(matTime.time){}


/** overrides default */
void MaturityTimePeriod::validatePop2Object(){
    // currently empty - may want to validate the time
}

/** Returns true if given expiry matches this. Overridden for performance */
bool MaturityTimePeriod::equalTo(const IObject* expiry) const{
    if (this == expiry){
        return true;
    }
    if (!expiry || expiry->getClass() != TYPE){
        return false;
    }
    const MaturityTimePeriod* matPeriod =
        STATIC_CAST(MaturityTimePeriod, expiry);
    return (time == matPeriod->time && 
            period->equals(matPeriod->period.get()));
}

/** Returns true if given expiry matches this */
bool MaturityTimePeriod::equals(const Expiry* expiry) const{
    return equalTo(expiry);
}

/** returns a hashcode for the object based upon the period and time.
    Overridden for performance */
int MaturityTimePeriod::hashCode() const{
    return (((size_t) TYPE) ^ period->hashCode() ^ time);
}

MaturityTimePeriod::~MaturityTimePeriod() {
    // empty
}

string MaturityTimePeriod::toString() const {
    return (period->toString() + " " + DateTime::timeFormat(time));
}

DateTime MaturityTimePeriod::toDate(const DateTime& aDate) const {
    static const string method = "MaturityTimePeriod::toDate";
    try {
        DateTime       date;
        int            dt;
        
        date = period->toDate(aDate);
        
        dt = date.getDate();

        return (DateTime(dt, time));
    }
    catch (exception& e) {
        throw ModelException(e, method, "Failed");
    } 
}

/** write object out to writer */
void MaturityTimePeriod::write(const string& tag, Writer* writer) const {
    try {
        IObjectConstSP obj(writer->objectStart(tag, "", this, false));
        if (obj.get()){
            writer->write(toString());
        }
        writer->objectEnd(tag, this);
    }
    catch (exception& e){
        throw ModelException(e, "MaturityTimePeriod::write");
    }
}

/** populate an empty object from reader */
void MaturityTimePeriod::import(Reader::Node* elem, Reader* reader) {
    static const string method = "MaturityTimePeriod::import";
    try {
        string    formatted = elem->value();

        // find the space that separates period & time
        string::size_type i = formatted.rfind(' ');

        if (i == string::npos) {
            throw ModelException(method, formatted + 
                                 " not in correct 'period time' format");
        }
        // import has privileged access
        const_cast<MaturityPeriodConstSP&>(period) = MaturityPeriodConstSP(
            new MaturityPeriod(formatted.substr(0, i)));
        const_cast<int&>(time) = DateTime::timeConvert(formatted.substr(i+1));
     }
    catch (exception& e) {
        throw ModelException(e, method, "Failed");
    }
} 

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void MaturityTimePeriod::outputWrite(
    const string& linePrefix,
    const string& prefix,
    ostream&      stream) const{
    stream << linePrefix << prefix << ": " << toString() << endl;
}

/** Returns the maturity portion of the MaturityTimePeriod */
string  MaturityTimePeriod::getMaturity() const{
    return period->toString();
}

/** Returns the time portion of the MaturityTimePeriod */
int  MaturityTimePeriod::getTime() const{
    return time;
}

/** Returns the MaturityPeriod portion of the MaturityTimePeriod */
MaturityPeriodConstSP MaturityTimePeriod::getMaturityPeriod() const{
    return period;
}

/* for reflection */
MaturityTimePeriod::MaturityTimePeriod():Expiry(TYPE), time(0){}

// we never put a public/private interface on MaturityTimePeriod, yet 
// we don't want to expose the int that represents the time. Too late to 
// change, so for DR Interface purposes provide a proxy which acts as the 
// interface we want to expose
class MaturityTimePeriodProxy: public CObject {
public:
    static CClassConstSP const TYPE;
private:
    string period;     
    string time;

    static IObjectSP toMTP(const IObjectConstSP& obj) {
        const MaturityTimePeriodProxy* proxy = 
            dynamic_cast<const MaturityTimePeriodProxy*>(obj.get());
        if (!proxy) {
            throw ModelException("MaturityTimePeriodProxy::toMTP",
                                 "object is not a MaturityTimePeriodProxy");
        }
    
        return IObjectSP(new MaturityTimePeriod(proxy->period, 
                                                DateTime::timeConvert(proxy->time)));
    }

    static IObjectSP fromMTP(const IObjectConstSP& obj) {
        const MaturityTimePeriod* mtp = 
            dynamic_cast<const MaturityTimePeriod*>(obj.get());
        if (!mtp) {
            throw ModelException("MaturityTimePeriodProxy::fromMTP",
                                 "object is not a MaturityTimePeriod");
        }
    
        return IObjectSP(new MaturityTimePeriodProxy(mtp->getMaturity(),
                                                     DateTime::timeFormat(mtp->getTime())));
    }

    MaturityTimePeriodProxy(const string& period, const string& time) : CObject(TYPE),
        period(period), time(time) {}

    /** for reflection */
    MaturityTimePeriodProxy(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); 
        clazz->setDescription("An expiry with a fixed time of day");
        REGISTER(MaturityTimePeriodProxy, clazz);
        SUPERCLASS(CObject);
        clazz->enableCloneOptimisations();
        EMPTY_SHELL_METHOD(defaultMaturityTimePeriodProxy);
        registerObjectProxy(MaturityTimePeriod::TYPE,
                            MaturityTimePeriodProxy::TYPE,
                            fromMTP,
                            toMTP);
        FIELD(period, "period");
        FIELD(time, "time (either HH-MM-SS format or BEX, XBS, SOD, EOD)");
    }

    static IObject* defaultMaturityTimePeriodProxy(){
        return new MaturityTimePeriodProxy();
    }    
};

CClassConstSP const MaturityTimePeriodProxy::TYPE = 
CClass::registerClassLoadMethod("MaturityTimePeriodProxy",
                                typeid(MaturityTimePeriodProxy), load);

class MaturityTimePeriodHelper{
public:
/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDRIProxyType(MaturityTimePeriodProxy::TYPE); // for dri
        REGISTER(MaturityTimePeriod, clazz);
        SUPERCLASS(Expiry);
        EMPTY_SHELL_METHOD(defaultMaturityTimePeriod);
        FIELD(period, "period eg 1W");
        FIELD(time, "time");
    }
    
    static IObject* defaultMaturityTimePeriod(){
        return new MaturityTimePeriod();
    }
};

CClassConstSP const MaturityTimePeriod::TYPE = CClass::registerClassLoadMethod(
    "MaturityTimePeriod", typeid(MaturityTimePeriod), 
    MaturityTimePeriodHelper::load);

DRLIB_END_NAMESPACE
