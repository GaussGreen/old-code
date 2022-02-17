//----------------------------------------------------------------------------
//
//   Filename    : ObservationType.cpp
//
//   Description : Classes for observation types - for use with AssetHistory
//
//   Author      : Ian Stares   
//
//   Date        : January 19, 2006
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/ObservationType.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

// base class
bool ObservationType::equals(const ObservationType& other) const {
    return getClass() == other.getClass();
}

ObservationTypeSP ObservationType::make(const string& suffix) {
    static const string routine = "ObservationType::make";
    try {
        // Create class by reflection
        string className = "Observation" + suffix;
        CClassConstSP clazz = CClass::forName(className);
        IObjectSP obj(clazz->newInstance());
        ObservationTypeSP obsType = ObservationTypeSP::dynamicCast(obj);
        return obsType;
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

ObservationType::ObservationType(CClassConstSP clazz) : CObject(clazz) {}

void ObservationType::load(CClassSP& clazz) {
    REGISTER(ObservationType, clazz);
    SUPERCLASS(CObject);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

CClassConstSP const ObservationType::TYPE = CClass::registerClassLoadMethod(
    "ObservationType", typeid(ObservationType), ObservationType::load);

DEFINE_TEMPLATE_TYPE(ObservationTypeArray);

/*****************************************************************************/

void ObservationOpen::load(CClassSP& clazz){
    REGISTER(ObservationOpen, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultOpenObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationOpen::ObservationOpen(): ObservationType(TYPE) {}

string ObservationOpen::toString() const {
    return "Open";
}

string ObservationOpen::indicativeTime() const {
    return DateTime::START_OF_DAY;
}

IObject* ObservationOpen::defaultOpenObsType() {
    return new ObservationOpen();
}

CClassConstSP const ObservationOpen::TYPE = CClass::registerClassLoadMethod(
    "ObservationOpen", typeid(ObservationOpen), ObservationOpen::load);
   
/*****************************************************************************/

void ObservationClose::load(CClassSP& clazz){
    REGISTER(ObservationClose, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultCloseObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationClose::ObservationClose(): ObservationType(TYPE) {}

string ObservationClose::toString() const {
    return "Close";
}

string ObservationClose::indicativeTime() const {
    return DateTime::END_OF_DAY;
}

IObject* ObservationClose::defaultCloseObsType() {
    return new ObservationClose();
}

CClassConstSP const ObservationClose::TYPE = CClass::registerClassLoadMethod(
    "ObservationClose", typeid(ObservationClose), ObservationClose::load);
   
/*****************************************************************************/

void ObservationHigh::load(CClassSP& clazz){
    REGISTER(ObservationHigh, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultHighObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationHigh::ObservationHigh(): ObservationType(TYPE) {}

string ObservationHigh::toString() const {
    return "High";
}

string ObservationHigh::indicativeTime() const {
    return DateTime::END_OF_DAY;
}

IObject* ObservationHigh::defaultHighObsType() {
    return new ObservationHigh();
}

CClassConstSP const ObservationHigh::TYPE = CClass::registerClassLoadMethod(
    "ObservationHigh", typeid(ObservationHigh), ObservationHigh::load);
   
/*****************************************************************************/

void ObservationLow::load(CClassSP& clazz){
    REGISTER(ObservationLow, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultLowObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationLow::ObservationLow(): ObservationType(TYPE) {}

string ObservationLow::toString() const {
    return "Low";
}

string ObservationLow::indicativeTime() const {
    return DateTime::END_OF_DAY;
}

IObject* ObservationLow::defaultLowObsType() {
    return new ObservationLow();
}

CClassConstSP const ObservationLow::TYPE = CClass::registerClassLoadMethod(
    "ObservationLow", typeid(ObservationLow), ObservationLow::load);
   
/*****************************************************************************/

void ObservationMid::load(CClassSP& clazz){
    REGISTER(ObservationMid, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultMidObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationMid::ObservationMid(): ObservationType(TYPE) {}

string ObservationMid::toString() const {
    return "Mid";
}

string ObservationMid::indicativeTime() const {
    return DateTime::END_OF_DAY;
}

IObject* ObservationMid::defaultMidObsType() {
    return new ObservationMid();
}

CClassConstSP const ObservationMid::TYPE = CClass::registerClassLoadMethod(
    "ObservationMid", typeid(ObservationMid), ObservationMid::load);
   
/*****************************************************************************/

void ObservationExact::load(CClassSP& clazz){
    REGISTER(ObservationExact, clazz);
    SUPERCLASS(ObservationType);
    FIELD(time, "The time of day at which observation is made");
    FIELD_MAKE_OPTIONAL(time); // (temporarily for backwards compatibility with rates code)
    FIELD(id, "Name of the type to distinguish between ones with same time");
    FIELD_MAKE_OPTIONAL(id); // (temporarily for backwards compatibility with rates code)
    EMPTY_SHELL_METHOD(defaultExactObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

bool ObservationExact::equals(const ObservationType& other) const {
    bool equal = false;
    if (getClass() == other.getClass()) {
        const ObservationExact& myOther = dynamic_cast<const ObservationExact&>(other);
        equal = (id == myOther.id && time == myOther.time);
    }
    return equal;
}

ObservationExact::ObservationExact(): ObservationType(TYPE) {}

string ObservationExact::toString() const {
    return id.empty() ? "Exact" : id;
}

// temporary for rates code which has some empty ObservationExact objects
bool ObservationExact::hasTimeAndId() const {
    return (!id.empty() && ! time.empty());
}

string ObservationExact::indicativeTime() const {
    if (!time.empty()) {
        return time;
    }
    throw ModelException("ObservationExact::indicativeTime",
                        "No time available");;
}

IObject* ObservationExact::defaultExactObsType() {
    return new ObservationExact();
}

CClassConstSP const ObservationExact::TYPE = CClass::registerClassLoadMethod(
    "ObservationExact", typeid(ObservationExact), ObservationExact::load);
   
/*****************************************************************************/

void ObservationVWAP::load(CClassSP& clazz){
    REGISTER(ObservationVWAP, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultVWAPObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationVWAP::ObservationVWAP(): ObservationType(TYPE) {}

string ObservationVWAP::toString() const {
    return "VWAP";
}

string ObservationVWAP::indicativeTime() const {
    return DateTime::END_OF_DAY;
}

IObject* ObservationVWAP::defaultVWAPObsType() {
    return new ObservationVWAP();
}

CClassConstSP const ObservationVWAP::TYPE = CClass::registerClassLoadMethod(
    "ObservationVWAP", typeid(ObservationVWAP), ObservationVWAP::load);
   
/*****************************************************************************/

void ObservationOSPOption::load(CClassSP& clazz){
    REGISTER(ObservationOSPOption, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultOSPOptionObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationOSPOption::ObservationOSPOption(): ObservationType(TYPE) {}

string ObservationOSPOption::toString() const {
    return "OSP Option";
}

string ObservationOSPOption::indicativeTime() const {
    return DateTime::END_OF_DAY;
}

IObject* ObservationOSPOption::defaultOSPOptionObsType() {
    return new ObservationOSPOption();
}

CClassConstSP const ObservationOSPOption::TYPE = CClass::registerClassLoadMethod(
    "ObservationOSPOption", typeid(ObservationOSPOption), ObservationOSPOption::load);
   
/*****************************************************************************/

void ObservationOSPOptionOpen::load(CClassSP& clazz){
    REGISTER(ObservationOSPOptionOpen, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultOSPOptionOpenObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationOSPOptionOpen::ObservationOSPOptionOpen(): ObservationType(TYPE) {}

string ObservationOSPOptionOpen::toString() const {
    return "OSP Option Open";
}

string ObservationOSPOptionOpen::indicativeTime() const {
    return DateTime::START_OF_DAY;
}

IObject* ObservationOSPOptionOpen::defaultOSPOptionOpenObsType() {
    return new ObservationOSPOptionOpen();
}

CClassConstSP const ObservationOSPOptionOpen::TYPE = CClass::registerClassLoadMethod(
    "ObservationOSPOptionOpen", typeid(ObservationOSPOptionOpen), ObservationOSPOptionOpen::load);
/*****************************************************************************/

void ObservationOSPFutureOpen::load(CClassSP& clazz){
    REGISTER(ObservationOSPFutureOpen, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultOSPFutureOpenObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationOSPFutureOpen::ObservationOSPFutureOpen(): ObservationType(TYPE) {}

string ObservationOSPFutureOpen::toString() const {
    return "OSP Future Open";
}

string ObservationOSPFutureOpen::indicativeTime() const {
    return DateTime::START_OF_DAY;
}

IObject* ObservationOSPFutureOpen::defaultOSPFutureOpenObsType() {
    return new ObservationOSPFutureOpen();
}

CClassConstSP const ObservationOSPFutureOpen::TYPE = CClass::registerClassLoadMethod(
    "ObservationOSPFutureOpen", typeid(ObservationOSPFutureOpen), ObservationOSPFutureOpen::load);
   
/*****************************************************************************/

void ObservationFuture::load(CClassSP& clazz){
    REGISTER(ObservationFuture, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultFutureObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationFuture::ObservationFuture(): ObservationType(TYPE) {}

string ObservationFuture::toString() const {
    return "Future";
}

string ObservationFuture::indicativeTime() const {
    return DateTime::END_OF_DAY;
}

IObject* ObservationFuture::defaultFutureObsType() {
    return new ObservationFuture();
}

CClassConstSP const ObservationFuture::TYPE = CClass::registerClassLoadMethod(
    "ObservationFuture", typeid(ObservationFuture), ObservationFuture::load);
   
/*****************************************************************************/

void ObservationAverage::load(CClassSP& clazz){
    REGISTER(ObservationAverage, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultAverageObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationAverage::ObservationAverage(): ObservationType(TYPE) {}

string ObservationAverage::toString() const {
    return "Average";
}

string ObservationAverage::indicativeTime() const {
    return DateTime::END_OF_DAY;
}

IObject* ObservationAverage::defaultAverageObsType() {
    return new ObservationAverage();
}

CClassConstSP const ObservationAverage::TYPE = CClass::registerClassLoadMethod(
    "ObservationAverage", typeid(ObservationAverage), ObservationAverage::load);
   
/*****************************************************************************/

void ObservationOther::load(CClassSP& clazz){
    REGISTER(ObservationOther, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultOtherObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationOther::ObservationOther(): ObservationType(TYPE) {}

string ObservationOther::toString() const {
    return "Other";
}

string ObservationOther::indicativeTime() const {
    return DateTime::END_OF_DAY;
}

IObject* ObservationOther::defaultOtherObsType() {
    return new ObservationOther();
}

CClassConstSP const ObservationOther::TYPE = CClass::registerClassLoadMethod(
    "ObservationOther", typeid(ObservationOther), ObservationOther::load);
   
/*****************************************************************************/

void ObservationNotUsed::load(CClassSP& clazz){
    REGISTER(ObservationNotUsed, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultNotUsedObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationNotUsed::ObservationNotUsed(): ObservationType(TYPE) {}

string ObservationNotUsed::toString() const {
    return "Not Used";
}

string ObservationNotUsed::indicativeTime() const {
    return DateTime::END_OF_DAY;
// XXX We should really complain about being asked for the time of a 'NotUsed' type,
//     but this has been disabled to make 6.6.0.0 work with the observation types
//     supplied for ObservationBuilderDaily, which is getting NotUsed. 
//    throw ModelException("ObservationNotUsed::indicativeTime",
//                         "Internal error - method not implemented");;
}

IObject* ObservationNotUsed::defaultNotUsedObsType() {
    return new ObservationNotUsed();
}

CClassConstSP const ObservationNotUsed::TYPE = CClass::registerClassLoadMethod(
    "ObservationNotUsed", typeid(ObservationNotUsed), ObservationNotUsed::load);
   
/*****************************************************************************/

void ObservationMIBOpen::load(CClassSP& clazz){
    REGISTER(ObservationMIBOpen, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultMIBOpenObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationMIBOpen::ObservationMIBOpen(): ObservationType(TYPE) {}

string ObservationMIBOpen::toString() const {
    return "MIB Open";
}

string ObservationMIBOpen::indicativeTime() const {
    return DateTime::START_OF_DAY;
}

IObject* ObservationMIBOpen::defaultMIBOpenObsType() {
    return new ObservationMIBOpen();
}

CClassConstSP const ObservationMIBOpen::TYPE = CClass::registerClassLoadMethod(
    "ObservationMIBOpen", typeid(ObservationMIBOpen), ObservationMIBOpen::load);
   
/*****************************************************************************/

void ObservationNikkei::load(CClassSP& clazz){
    REGISTER(ObservationNikkei, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultNikkeiObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationNikkei::ObservationNikkei(): ObservationType(TYPE) {}

string ObservationNikkei::toString() const {
    return "Nikkei";
}

string ObservationNikkei::indicativeTime() const {
    return DateTime::START_OF_DAY;
}

IObject* ObservationNikkei::defaultNikkeiObsType() {
    return new ObservationNikkei();
}

CClassConstSP const ObservationNikkei::TYPE = CClass::registerClassLoadMethod(
    "ObservationNikkei", typeid(ObservationNikkei), ObservationNikkei::load);
 /*****************************************************************************/

void ObservationNikkei300::load(CClassSP& clazz){
    REGISTER(ObservationNikkei300, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultNikkei300ObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationNikkei300::ObservationNikkei300(): ObservationType(TYPE) {}

string ObservationNikkei300::toString() const {
    return "Nikkei 300";
}

string ObservationNikkei300::indicativeTime() const {
    return DateTime::START_OF_DAY;
}

IObject* ObservationNikkei300::defaultNikkei300ObsType() {
    return new ObservationNikkei300();
}

CClassConstSP const ObservationNikkei300::TYPE = CClass::registerClassLoadMethod(
    "ObservationNikkei300", typeid(ObservationNikkei300), ObservationNikkei300::load);
 /*****************************************************************************/

void ObservationTime::load(CClassSP& clazz){
    REGISTER(ObservationTime, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultTimeObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationTime::ObservationTime(): ObservationType(TYPE) {}

string ObservationTime::toString() const {
    return "Time";
}

string ObservationTime::indicativeTime() const {
    return DateTime::END_OF_DAY;
}

IObject* ObservationTime::defaultTimeObsType() {
    return new ObservationTime();
}

CClassConstSP const ObservationTime::TYPE = CClass::registerClassLoadMethod(
    "ObservationTime", typeid(ObservationTime), ObservationTime::load);
 /*****************************************************************************/

void ObservationOptionExpiry::load(CClassSP& clazz){
    REGISTER(ObservationOptionExpiry, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultOptionExpiryObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationOptionExpiry::ObservationOptionExpiry(): ObservationType(TYPE) {}

string ObservationOptionExpiry::toString() const {
    return "Option Expiry";
}

string ObservationOptionExpiry::indicativeTime() const {
    return DateTime::END_OF_DAY;
}

IObject* ObservationOptionExpiry::defaultOptionExpiryObsType() {
    return new ObservationOptionExpiry();
}

CClassConstSP const ObservationOptionExpiry::TYPE = CClass::registerClassLoadMethod(
    "ObservationOptionExpiry", typeid(ObservationOptionExpiry), ObservationOptionExpiry::load);
 /*****************************************************************************/

void ObservationPrecio::load(CClassSP& clazz){
    REGISTER(ObservationPrecio, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultPrecioObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationPrecio::ObservationPrecio(): ObservationType(TYPE) {}

string ObservationPrecio::toString() const {
    return "Precio";
}

string ObservationPrecio::indicativeTime() const {
    return DateTime::END_OF_DAY;
}

IObject* ObservationPrecio::defaultPrecioObsType() {
    return new ObservationPrecio();
}

CClassConstSP const ObservationPrecio::TYPE = CClass::registerClassLoadMethod(
    "ObservationPrecio", typeid(ObservationPrecio), ObservationPrecio::load);
 /*****************************************************************************/

void ObservationPremi::load(CClassSP& clazz){
    REGISTER(ObservationPremi, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultPremiObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationPremi::ObservationPremi(): ObservationType(TYPE) {}

string ObservationPremi::toString() const {
    return "Premi";
}

string ObservationPremi::indicativeTime() const {
    return DateTime::END_OF_DAY;
}

IObject* ObservationPremi::defaultPremiObsType() {
    return new ObservationPremi();
}

CClassConstSP const ObservationPremi::TYPE = CClass::registerClassLoadMethod(
    "ObservationPremi", typeid(ObservationPremi), ObservationPremi::load);
 /*****************************************************************************/

void ObservationPrezzo::load(CClassSP& clazz){
    REGISTER(ObservationPrezzo, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultPrezzoObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationPrezzo::ObservationPrezzo(): ObservationType(TYPE) {}

string ObservationPrezzo::toString() const {
    return "Prezzo Ufficiale";
}

string ObservationPrezzo::indicativeTime() const {
    return DateTime::END_OF_DAY;
}

IObject* ObservationPrezzo::defaultPrezzoObsType() {
    return new ObservationPrezzo();
}

CClassConstSP const ObservationPrezzo::TYPE = CClass::registerClassLoadMethod(
    "ObservationPrezzo", typeid(ObservationPrezzo), ObservationPrezzo::load);
 /*****************************************************************************/

void ObservationRifer::load(CClassSP& clazz){
    REGISTER(ObservationRifer, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultRiferObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationRifer::ObservationRifer(): ObservationType(TYPE) {}

string ObservationRifer::toString() const {
    return "Riferimento";
}

string ObservationRifer::indicativeTime() const {
    return DateTime::END_OF_DAY;
}

IObject* ObservationRifer::defaultRiferObsType() {
    return new ObservationRifer();
}

CClassConstSP const ObservationRifer::TYPE = CClass::registerClassLoadMethod(
    "ObservationRifer", typeid(ObservationRifer), ObservationRifer::load);
 /*****************************************************************************/

void ObservationSQSP500::load(CClassSP& clazz){
    REGISTER(ObservationSQSP500, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultSQSP500ObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationSQSP500::ObservationSQSP500(): ObservationType(TYPE) {}

string ObservationSQSP500::toString() const {
    return "SQ S&P500";
}

string ObservationSQSP500::indicativeTime() const {
    return DateTime::START_OF_DAY;
}

IObject* ObservationSQSP500::defaultSQSP500ObsType() {
    return new ObservationSQSP500();
}

CClassConstSP const ObservationSQSP500::TYPE = CClass::registerClassLoadMethod(
    "ObservationSQSP500", typeid(ObservationSQSP500), ObservationSQSP500::load);
 /*****************************************************************************/

void ObservationTopix::load(CClassSP& clazz){
    REGISTER(ObservationTopix, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultTopixObsType);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

ObservationTopix::ObservationTopix(): ObservationType(TYPE) {}

string ObservationTopix::toString() const {
    return "Topix";
}

string ObservationTopix::indicativeTime() const {
    return DateTime::START_OF_DAY;
}

IObject* ObservationTopix::defaultTopixObsType() {
    return new ObservationTopix();
}

CClassConstSP const ObservationTopix::TYPE = CClass::registerClassLoadMethod(
    "ObservationTopix", typeid(ObservationTopix), ObservationTopix::load);
 /*****************************************************************************/

void StruckObservationType::load(CClassSP& clazz){
    REGISTER(StruckObservationType, clazz);
    SUPERCLASS(ObservationType);
    EMPTY_SHELL_METHOD(defaultStruckObsType);
    FIELD(primaryObsType, "The obs type for the plain asset");
    FIELD(fxObsType, "The obs type for the fx asset");
    // NOTE THIS ONE IS CURRENTLY EXPOSED PUBLICLY
}

StruckObservationType::StruckObservationType(): ObservationType(TYPE) {}

StruckObservationType::StruckObservationType(ObservationType* primary,
                                             ObservationType* fx): 
                        ObservationType(TYPE) {
    primaryObsType = ObservationTypeSP(primary);
    fxObsType = ObservationTypeSP(fx);
}

string StruckObservationType::toString() const {
    return primaryObsType->toString();
}

string StruckObservationType::indicativeTime() const {
    // Use the obs rule of the underlying, not the fx
    return primaryObsType->indicativeTime();
}

IObject* StruckObservationType::defaultStruckObsType() {
    return new StruckObservationType();
}

ObservationType* StruckObservationType::getPrimaryObsType() const {
    return primaryObsType.get();
}

ObservationType* StruckObservationType::getFXObsType() const {
    return fxObsType.get();
}

CClassConstSP const StruckObservationType::TYPE = CClass::registerClassLoadMethod(
    "StruckObservationType", typeid(StruckObservationType), StruckObservationType::load);
   
/*****************************************************************************/
DRLIB_END_NAMESPACE




