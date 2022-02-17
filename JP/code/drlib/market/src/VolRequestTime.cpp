//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolRequestTime.cpp
//
//   Description : Vol request for dealing with [trading] time
//
//   Author      : Mark A Robson
//
//   Date        : 1 Nov 2004
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#define QLIB_VOLREQUESTTIME_CPP
#include "edginc/VolRequestTime.hpp"

DRLIB_BEGIN_NAMESPACE

VolRequestTime::~VolRequestTime(){}

VolRequestTime::VolRequestTime(): CVolRequest(TYPE){}

class VolRequestTime::Processed: public CObject,
                      public virtual IVolProcessed{
public:
    static CClassConstSP const TYPE;

    Processed(const string&     name, 
              TimeMetricConstSP metric): 
        CObject(TYPE), name(name), metric(metric){}

    /** identifies the market data name of the volatility */
    virtual string getName() const{
        return name;
    }

    /** calculates the trading time between two dates */
    virtual double calcTradingTime(const DateTime &date1, 
                                   const DateTime &date2) const{
        return metric->yearFrac(date1, date2);
    }

    /** retieve time measure for the vol */
    virtual TimeMetricConstSP GetTimeMetric()const{
        return metric;
    }

    virtual ~Processed(){}
private:
    static void load(CClassSP& clazz){
        REGISTER(Processed, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IVolProcessed);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(name, "Vol's name");
        FIELD(metric, "How to handle trading time");
    }

    static IObject* defaultConstructor(){
        return new VolRequestTime();
    }

    /// fields
    string            name;
    TimeMetricConstSP metric;
};

CClassConstSP const VolRequestTime::Processed::TYPE = 
CClass::registerClassLoadMethod("VolRequestTime::Processed",
                                typeid(VolRequestTime::Processed), load);

/** Creates a VolProcessed which supports only those methods that 
    VolProcessed does (ie suitable for dealing with [trading] time */
IVolProcessed* VolRequestTime::createVolProcessed(const string&     name, 
                                                  TimeMetricConstSP metric){
    return new Processed(name, metric);
}

/** Creates a VolProcessed which supports only those methods that 
    VolProcessed does (ie suitable for dealing with [trading] time).
    Here, simple DateTime::yearFrac is used */
IVolProcessed* VolRequestTime::createVolProcessed(const string&     name){
    HolidaySP hols(Holiday::noHolidays());
    return new Processed(name, 
                         TimeMetricConstSP(new TimeMetric(1.0, hols.get())));
}


void VolRequestTime::load(CClassSP& clazz){
    REGISTER(VolRequestTime, clazz);
    SUPERCLASS(CVolRequest);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

IObject* VolRequestTime::defaultConstructor(){
    return new VolRequestTime();
}

CClassConstSP const VolRequestTime::TYPE = 
CClass::registerClassLoadMethod("VolRequestTime", typeid(VolRequestTime), load);

DRLIB_END_NAMESPACE

