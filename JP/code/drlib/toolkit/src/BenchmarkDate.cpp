//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BenchmarkDate.cpp
//
//   Description : Defines fixed date expiries used to define yield curve & 
//                 vol surface points
//
//   Author      : Andrew J Swain
//
//   Date        : 29 January 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_BENCHMARKDATE_CPP
#include "edginc/BenchmarkDate.hpp"
#include "edginc/Writer.hpp"

DRLIB_BEGIN_NAMESPACE
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<BenchmarkDate>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<BenchmarkDate>);

BenchmarkDate::BenchmarkDate(const DateTime& date) : Expiry(TYPE), date(date) {
    // empty
}

BenchmarkDate::BenchmarkDate(const BenchmarkDate& bm): 
    Expiry(TYPE), date(bm.date) {}

/** Explicit implementation for performance */
IObject* BenchmarkDate::clone() const{
    int& count = getRefCount();
    if (count == 0){
        return new BenchmarkDate(*this);
    } else {
        return  const_cast<IObject*>((const IObject*)this);
    }
}

BenchmarkDate::~BenchmarkDate() {
    // empty
}


string BenchmarkDate::toString() const {
    return date.toString();
}

DateTime BenchmarkDate::toDate(const DateTime& aDate) const {
    return date;
}

/** write object out in XML format */
void BenchmarkDate::write(const string& tag, Writer* writer) const {
    try {
        IObjectConstSP obj(writer->objectStart(tag, "", this, false));
        if (obj.get()){
            writer->write(toString());
        }
        writer->objectEnd(tag, this);
    }
    catch (exception& e){
        throw ModelException(&e, "BenchmarkDate::write");
    }
}

/** Returns true if given object matches this. Overridden for performance */
bool BenchmarkDate::equalTo(const IObject* expiry) const{
    if (this == expiry){
        return true;
    }
    if (!expiry || expiry->getClass() != TYPE){
        return false;
    }
    const BenchmarkDate* bMark = STATIC_CAST(BenchmarkDate, expiry);
    return (date.equals(bMark->date));
}


/** Returns true if given expiry matches this */
bool BenchmarkDate::equals(const Expiry* expiry) const{
    return equalTo(expiry);
}

/** returns a hashcode for the object based upon the date and time.
    Overridden for performance */
int BenchmarkDate::hashCode() const{
    return (((size_t) TYPE) ^ date.hashCode());
}

/** specific to BenchmarkDate - return this expiry as an absolute date */
DateTime BenchmarkDate::toDate() const{
    return date;
}

/** populate an empty object from XML description */
void BenchmarkDate::import(Reader::Node* elem, Reader* reader) {
    static const string method = "BenchmarkDate::import";
    try {
        // xmlImport has privileged access
        const_cast<DateTime&>(date).import(elem, reader);
    }
    catch (exception& e) {
        throw ModelException(e, method, "Failed");
    }
} 

/* for reflection */
BenchmarkDate::BenchmarkDate():Expiry(TYPE){}

// we never put a public/private interface on BenchmarkDate, yet 
// we don't want to expose the insides. Too late to change, so for 
// DR Interface purposes provide a proxy which acts as the 
// interface we want to expose
class BenchmarkDateProxy: public CObject {
public:
    static CClassConstSP const TYPE;
private:
    string date;     
    string time;

    static IObjectSP toBMD(const IObjectConstSP& obj) {
        const BenchmarkDateProxy* proxy = 
            dynamic_cast<const BenchmarkDateProxy*>(obj.get());
        if (!proxy) {
            throw ModelException("BenchmarkDateProxy::toBMD",
                                 "object is not a BenchmarkDateProxy");
        }
    
        DateTime dt(proxy->date, proxy->time);

        return IObjectSP(new BenchmarkDate(dt));
    }

    static IObjectSP fromBMD(const IObjectConstSP& obj) {
        const BenchmarkDate* bmd = 
            dynamic_cast<const BenchmarkDate*>(obj.get());
        if (!bmd) {
            throw ModelException("BenchmarkDateProxy::fromBMD",
                                 "object is not a BenchmarkDate");
        }
    
        DateTime dt = bmd->toDate();
        return IObjectSP(new BenchmarkDateProxy(DateTime::dateFormat(dt.getDate()),
                                                DateTime::timeFormat(dt.getTime())));
    }

    BenchmarkDateProxy(const string& date, const string& time) : CObject(TYPE),
        date(date), time(time) {}

    /** for reflection */
    BenchmarkDateProxy(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); 
        REGISTER(BenchmarkDateProxy, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBenchmarkDateProxy);
        registerObjectProxy(BenchmarkDate::TYPE,
                            BenchmarkDateProxy::TYPE,
                            fromBMD,
                            toBMD);
        FIELD(date, "date (in DD-MMM-YYYY format)");
        FIELD(time, "time (either HH-MM-SS format or BEX, XBS, SOD, EOD)");
    }

    static IObject* defaultBenchmarkDateProxy(){
        return new BenchmarkDateProxy();
    }    
};

CClassConstSP const BenchmarkDateProxy::TYPE = 
CClass::registerClassLoadMethod("BenchmarkDateProxy",
                                typeid(BenchmarkDateProxy), load);

class BenchmarkDateHelper{
public:
/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDRIProxyType(BenchmarkDateProxy::TYPE); // use proxy for dri
        REGISTER(BenchmarkDate, clazz);
        SUPERCLASS(Expiry);
        clazz->enableCloneOptimisations();
        EMPTY_SHELL_METHOD(defaultBenchmarkDate);
        FIELD(date, "benchmark date");
    }
    
    static IObject* defaultBenchmarkDate(){
        return new BenchmarkDate();
    }
};

CClassConstSP const BenchmarkDate::TYPE = CClass::registerClassLoadMethod(
    "BenchmarkDate", typeid(BenchmarkDate), BenchmarkDateHelper::load);


DRLIB_END_NAMESPACE
