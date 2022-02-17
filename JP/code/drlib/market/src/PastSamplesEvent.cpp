
#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/PastSamplesEvent.hpp"
#include "edginc/EventResults.hpp"

DRLIB_BEGIN_NAMESPACE

/* for reflection */
PastSamplesEvent::PastSamplesEvent(): Event(TYPE) {}

// Invoked when Class is 'loaded'
void PastSamplesEvent::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(PastSamplesEvent, clazz);
    SUPERCLASS(Event);
    IMPLEMENTS(IEvent);
    EMPTY_SHELL_METHOD(defaultPastSamplesEvent);
    FIELD(observationType, "Observation Type e.g. Open, Close, High, Low etc");
    FIELD(source, "The source from which the sample is taken");
    FIELD(samplingDates, "The dates sampling was attempted (prior to ISDA adjustment)");
    FIELD(observationDates, "The dates the sample was actually observed");
    FIELD(overridden, "Were samples overridden at the instrument level");
    FIELD(levels, "The levels of the asset which were observed");
    FIELD(context, "Some more contexr for the fixing");
    FIELD(overrideInts, "Transient field to get round iterator issue");
    FIELD_MAKE_TRANSIENT(overrideInts);
}

PastSamplesEvent::PastSamplesEvent(const DateTime& eDate, const string source, 
                                   const string obsType, const FixingType* context) :
                        Event(TYPE, eDate), observationType(obsType), source(source),
                        context(copy(context)) {
     samplingDates = DateTimeArraySP(new DateTimeArray(0));
     observationDates = DateTimeArraySP(new DateTimeArray(0));
     levels = DoubleArraySP(new DoubleArray(0));
     overrideInts = IntArraySP(new IntArray(0));
}

IObject* PastSamplesEvent::defaultPastSamplesEvent(){
    return new PastSamplesEvent();
}

void PastSamplesEvent::buildBoolArray() {
    overridden = BoolArraySP(new BoolArray(overrideInts->size()));
    for (int i = 0; i < overrideInts->size(); ++i) {
        (*overridden)[i] = ((*overrideInts)[i] == 0) ? false : true;
    }
    // don't need it anymore so we might as well clear it
    overrideInts->clear();
}

void PastSamplesEvent::addSample(const DateTime& sampleDate, 
                                 const DateTime& obsDate, 
                                 double level, bool isOverride) {
    try {
        // let's see if we've already got it/ work out where it goes
        vector<DateTime>::iterator sampleIter(samplingDates->begin());
        vector<DateTime>::iterator obsIter(observationDates->begin());
        vector<double>::iterator levelIter(levels->begin());
        vector<int>::iterator overrideIter(overrideInts->begin());

        int overrideInt = isOverride ? 1 : 0;
        for( ; sampleIter != samplingDates->end()
                        && (*sampleIter) <= sampleDate; ) {
            if ((*sampleIter) == sampleDate && (*obsIter) == obsDate) {
                if ((*levelIter) != level) {
                    throw ModelException("Inconsistent levels found for "
                                            "asset for sample date (" +
                                            sampleDate.toString() +
                                            ") and observationDate (" +
                                            obsDate.toString() +
                                            "). Two levels found - " +
                                            Format::toString((*levelIter)) +
                                            " and " + Format::toString(level));
                }
                if ((*overrideIter) != overrideInt) {
                    double override = isOverride ? level : (*levelIter);
                    double central = isOverride ? (*levelIter) : level;
                    throw ModelException("Inconsistent levels found for "
                                            "asset for sample date (" +
                                            sampleDate.toString() +
                                            ") and observationDate (" +
                                            obsDate.toString() +
                                            "). Found centralised level " +
                                            Format::toString(central) +
                                            " and an override level " + 
                                            Format::toString(override));
                }
                return; //it's already here - our job is done
            } else {
                ++sampleIter;
                ++obsIter;
                ++levelIter;
                ++overrideIter;
            }
        }
        // if we got here it must be a new sample so let's stick it in
        sampleIter = samplingDates->insert(sampleIter, sampleDate);
        obsIter = observationDates->insert(obsIter, obsDate);
        levelIter = levels->insert(levelIter, level);
        overrideIter = overrideInts->insert(overrideIter, overrideInt);
    } catch (exception& e) {
        throw ModelException(e, "PastSamplesEvent::addSample");
    }
}    

CClassConstSP PastSamplesEvent::getEventInterface() const {
    return IEventHandler::TYPE;
}

void PastSamplesEvent::getEvents(IObjectConstSP       obj,
                          IModel*              model,
                          const DateTime&      eDate,
                          EventResults*        events) const {
    const IEventHandler& eventObj = dynamic_cast<const IEventHandler&>(*obj);
    eventObj.getEvents(this, model, eDate, events);
}

// if we get risk and legal events we need to know whether a risk
// event will just be a duplicated legal event (and hence removable)
bool PastSamplesEvent::isValidRiskEvent(const DateTime& today) const {
    // all past samples events can only be legal
    return false;
}

CClassConstSP const PastSamplesEvent::TYPE = CClass::registerClassLoadMethod(
    "PastSamplesEvent", typeid(PastSamplesEvent), load);

CClassConstSP const PastSamplesEvent::IEventHandler::TYPE = CClass::registerInterfaceLoadMethod(
    "PastSamplesEvent::IEventHandler", typeid(PastSamplesEvent::IEventHandler), 0);

DEFINE_TEMPLATE_TYPE(PastSamplesEventArray);

/*****************************************************************************/

CClassConstSP const PastSamplesCollector::TYPE = CClass::registerClassLoadMethod(
    "PastSamplesCollector", typeid(PastSamplesCollector), load);

// Invoked when Class is 'loaded'
void PastSamplesCollector::load(CClassSP& clazz){
//    clazz->setPublic(); // commented out to be private
    REGISTER(PastSamplesCollector, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultPastSamplesCollector);
    FIELD(samples, "The collected samples");
    FIELD_MAKE_TRANSIENT(samples);
    FIELD(numResultSets, "Number of different sample sets");
    FIELD_MAKE_TRANSIENT(numResultSets);
    FIELD(today, "The collected samples");
    FIELD_MAKE_TRANSIENT(today);
}

/* for reflection */
PastSamplesCollector::PastSamplesCollector(): CObject(TYPE),
            numResultSets(0) {
    samples = PastSamplesEventArraySP(new PastSamplesEventArray(0));
}

PastSamplesCollector::PastSamplesCollector(const DateTime& today): 
            CObject(TYPE), numResultSets(0), today(today) {
    samples = PastSamplesEventArraySP(new PastSamplesEventArray(0));
}
IObject* PastSamplesCollector::defaultPastSamplesCollector(){
    return new PastSamplesCollector();
}

// the main method
// adds a new past sample event to the list.
// it will put it in the correct object by asset/src/obsType
// and chekc for duplication
void PastSamplesCollector::addSample(const DateTime&          sampleDate,
                                    const DateTime&          obsDate,
                                    const ObservationSource* source,
                                    const ObservationType*   obsType,
                                    double                   level,
                                    const FixingType*        fixType,
                                    bool                     isOverride) {
    try {
        // use primary source and to string on the type/source
        const string sourceName = source->getPrimarySource();
        const string obsName = obsType->toString();

        // first see if we've already got a result set for this 
        // underlying/source and observationType
        int idx = 0;
        bool exists = false;
        for (idx = 0; idx < numResultSets; ++idx) {
            if (sourceName == ((*samples)[idx])->source &&
                    obsName == ((*samples)[idx])->observationType &&
                    fixType->equals(*(((*samples)[idx])->context))) {
                exists = true;
                break;
            }
        }
        if (!exists) {
            // make a new object for this new underlying/source/obsType
            PastSamplesEventSP sampleObj = 
                PastSamplesEventSP(new PastSamplesEvent(today, sourceName, obsName, fixType));
            samples->append(sampleObj);
            idx = numResultSets++;
        }
        // and add the sample to the list - checking for dodgy duplicates
        ((*samples)[idx])->addSample(sampleDate, obsDate, level, isOverride);
    } catch (exception& e) {
        throw ModelException(e, "PastSamplesCollector::addSample");
    }
}    

// put the finished event results in the results structure
void PastSamplesCollector::finaliseResults(EventResults* events) {
    for (int i = 0; i < numResultSets; ++i) {
        (*samples)[i]->buildBoolArray();
        events->addEvent((*samples)[i].get());
    }
}

DRLIB_END_NAMESPACE
