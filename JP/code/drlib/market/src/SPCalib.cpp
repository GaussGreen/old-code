//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPCalib.cpp
//
//   Description : smile params for basis spread
//
//   Date        : 
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SPCalib.hpp"
#include "edginc/Format.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/MRSpotVolRequest.hpp"

DRLIB_BEGIN_NAMESPACE


/// TYPE'd but not storing any reflection information
CClassConstSP const SPCalib::TYPE = 
CClass::registerClassLoadMethod(
    "SPCalib", typeid(SPCalib), load);

CClassConstSP const SPCalib::Processed::TYPE = 
CClass::registerClassLoadMethod(
    "SPCalib::Processed", typeid(SPCalib::Processed), load);

void SPCalib::Processed::spotVol(
            const DateTime&      initialStartDate,
            const DateTimeArray& subsequentDates,
            DoubleArray&         spotvol) const
{
    spotvol = DoubleArray(subsequentDates.size(), vol->flatSpotVol);
}
/** calculates the trading time between two dates */
double SPCalib::Processed::calcTradingTime(const DateTime &date1, 
                                const DateTime &date2) const
{
    return date1.yearFrac(date2);
}

/** retrieve time measure for the vol - no time metric here */
TimeMetricConstSP SPCalib::Processed::GetTimeMetric() const
{
    HolidaySP noHols(Holiday::noHolidays());
    return TimeMetricConstSP(new TimeMetric(1.0, noHols.get()));
}

SPCalib::SPCalib() :
    MarketObject(TYPE),
    qSpread(0.0),
    flatSpotVol(0.0),
    meanRev(0.0),
    name()
{}

bool SPCalib::recognizedVolRequest(const CVolRequest* volRequest)
{
    return  volRequest->getClass() == SPCalib::SPCalibRequest::TYPE ||
            volRequest->getClass() == MRSpotVolRequest::TYPE;
}


CVolProcessed* SPCalib::getProcessedVol(const CVolRequest* volRequest) const
{
    if (recognizedVolRequest(volRequest))
    {
        return new Processed(this);
    }
    throw ModelException("SPCalib::getProcessedVol",
                            "Request of type "+
                            volRequest->getClass()->getName()+
                            " not supported");
}

CClassConstSP const SPCalib::SPCalibRequest::TYPE = 
CClass::registerClassLoadMethod("SPCalib::SPCalibRequest", typeid(SPCalib::SPCalibRequest), load);

void SPCalib::SPCalibRequest::load(CClassSP& clazz) 
{
    // class deliberately private - no need for clients to know it exists
    REGISTER(SPCalib::SPCalibRequest, clazz);
    SUPERCLASS(CVolRequest);
    EMPTY_SHELL_METHOD(defaultConstructor);
}






void SPCalib::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(SPCalib, clazz);
    SUPERCLASS(MarketObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(name, "Name for this vol");
    FIELD(flatSpotVol, "Spot vol");
    FIELD(meanRev, "Mean reversion");
    FIELD(qSpread, "Smile parameter: 0.0 - normal, 1.0 - lognormal");
}


// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(SPCalibWrapper);


void SPCalib::Processed::load(CClassSP& clazz)
{
    // class deliberately private - no need for clients to know it exists
//    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(SPCalib::Processed, clazz);
    SUPERCLASS(MRSpotVolProcessed);
    EMPTY_SHELL_METHOD(defaultConstructor);
}


/* external symbol to allow class to be forced to be linked in */
bool SPCalibLinkIn(){
    return SPCalib::TYPE != NULL && SPCalib::Processed::TYPE != NULL;
}



DRLIB_END_NAMESPACE
