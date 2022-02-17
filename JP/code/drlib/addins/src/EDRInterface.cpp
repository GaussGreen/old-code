//////////////////////////////////////////////////////////////
// Equity Derivatives Research: Models Library Interface
//
//
//////////////////////////////////////////////////////////////

#include "edginc/config.hpp"
#include "edginc/EDRInterface.h"
#include "edginc/DataDictionary.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Malloc.hpp"
#include "edginc/Library.hpp"
#include "edginc/AddinLib.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CompositeInstrument.hpp"
#include "edginc/OutputFile.hpp"
#include "edginc/DRWrapper.hpp"
#include "edginc/Modifier.hpp"
#include "edginc/Version.hpp"
#include "edginc/XMLReader.hpp"
#include "edginc/PrivateObject.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE
/** We want to work with smart pointers but need to hide these from
    the interface - so wrap then in a structure and return pointers to
    that */
class ObjectWrapper{
public:
    // constructor from smart pointer
    ObjectWrapper(const IObjectSP& object): object(object){
        if (!object){
            throw ModelException("ObjectWrapper", "NULL pointer supplied");
        }
    }

    // constructor from pointer
    ObjectWrapper(IObject* object): object(object){
        if (!object){
            throw ModelException("ObjectWrapper", "NULL pointer supplied");
        }
    }

    // occassionally useful to allow null - may be better off banning this
    // though
    ObjectWrapper() {}

    // single field
    IObjectSP object;

};

class DataDictIterator: public CObject{
public:
    static CClassConstSP const TYPE;
    // fields
    CDataDictionary::ObjectHashTable::iterator    iterator; // $unregistered
    CDataDictionary::ObjectHashTable::iterator    iteratorEnd; // $unregistered
    ObjectWrapper                                 objWrapper; // $unregistered

    // constructor
    DataDictIterator(const CDataDictionarySP& dd):
        CObject(TYPE), 
        iterator(dd->iterationBegin()),
        iteratorEnd(dd->iterationEnd()) {}

    /** Invoked when this class is 'loaded' */
    static void load(CClassSP& classToLoad){
        // empty for now
    }
};

CClassConstSP const DataDictIterator::TYPE = CClass::registerClassLoadMethod(
    "DataDictIterator", typeid(DataDictIterator), load);
typedef smartPtr<DataDictIterator> DataDictIteratorSP;

class ResultsIterator: public CObject{
public:
    static CClassConstSP const TYPE;
    // fields
    ResultsSP             results; // $unregistered
    vector<const string*> packets;     // $unregistered
    vector<pair<OutputNameConstSP, IObjectConstSP> > packetContents;  // $unregistered
    unsigned int          packetIdx; // $unregistered
    unsigned int          contentsIdx; // $unregistered

    ObjectWrapper  nameWrapper; // $unregistered
    ObjectWrapper  objWrapper; // $unregistered

    // constructor
    ResultsIterator(const ResultsSP& results): CObject(TYPE), 
        results(results),packetIdx(0), contentsIdx(0){

        packets = results->listPyramidPackets();
        if (!packets.empty()){
            packetContents = results->listPacketResults(*packets[0]);
            if (packetContents.empty()){
                increment();
            }
        }
    }

    // increment
    void increment(void){
        contentsIdx++;
        if (contentsIdx == packetContents.size()){
            do {
                packetIdx++;
                contentsIdx = 0;
                if (packetIdx < packets.size()){
                    packetContents = 
                        results->listPacketResults(*packets[packetIdx]);
                }
            } while (packetIdx < packets.size() && packetContents.size() == 0);
        }
    }
    /** Invoked when this class is 'loaded' */
    static void load(CClassSP& classToLoad){
        // empty for now
    }
};

CClassConstSP const ResultsIterator::TYPE = CClass::registerClassLoadMethod(
    "ResultsIterator", typeid(ResultsIterator), load);
typedef smartPtr<ResultsIterator> ResultsIteratorSP;

DRLIB_END_NAMESPACE

USING_DRLIB_NAMESPACE

static char* stringCopy(const char* in) {
    char *out = 0;

    if ((out = (char*)malloc(sizeof(char) * strlen(in) + 1))) {
        strcpy(out, in);
    }
    return out;
}

// build a composite instrument
static CompositeInstrumentSP makeComposite(
    int                  lengthOfArrays,
    EDRObject*           models,
    EDRObject*           instruments,
    EDRObject*           controls,
    double*              multipliers,
    double*              weights,
    EDRObject            scenario,
    EDRMarketDataCache   marketData) {
    const static char routine[] = "makeComposite";

    CompositeInstrumentSP composite;

    try {
        int i;
        int size = lengthOfArrays;
        if (!models){
            throw ModelException(routine, "NULL models");
        }
        if (!instruments){
            throw ModelException(routine, "NULL instruments");
        }
        if (!controls){
            throw ModelException(routine, "NULL controls");
        }
        if (!weights){
            throw ModelException(routine, "NULL weights");
        }
        if (!marketData){
            throw ModelException(routine, "NULL marketData");
        }
        if (size <= 0){
            throw ModelException(routine, "array lengths must be > 0");
        }

        ObjectWrapper* scenWrapper = (ObjectWrapper*)scenario;
        ObjectWrapper* dataWrapper = (ObjectWrapper*)marketData;

        ScenarioSP scenarioSP;
        if (scenWrapper) {
            scenarioSP = ScenarioSP::dynamicCast(scenWrapper->object);
        }

        CMarketDataSP dataSP(CMarketDataSP::dynamicCast(dataWrapper->object));
 
        ObjectArraySP modelSP(new ObjectArray(size));
        ObjectArraySP instSP(new ObjectArray(size));
        ObjectArraySP ctrlSP(new ObjectArray(size));
        DoubleArraySP multsSP(new DoubleArray(size));
        DoubleArraySP weightSP(new DoubleArray(size));

        for (i = 0; i < size; i++) {
            ObjectWrapper* modelWrapper = (ObjectWrapper*)models[i];
            ObjectWrapper* instWrapper = (ObjectWrapper*)instruments[i];
            ObjectWrapper* ctrlWrapper = (ObjectWrapper*)controls[i];
            (*modelSP)[i]  = modelWrapper->object;
            (*instSP)[i]   = instWrapper->object;
            (*ctrlSP)[i]   = ctrlWrapper->object;
            (*multsSP)[i]    = multipliers[i];
            (*weightSP)[i] = weights[i];
        }
  
        composite = CompositeInstrumentSP(new CompositeInstrument(scenarioSP,
                                                                  modelSP,
                                                                  instSP,
                                                                  ctrlSP,
                                                                  multsSP,
                                                                  weightSP,
                                                                  dataSP));
    } 
    catch (exception& e) {
        throw ModelException(e, routine);
    }
    return composite;
}

// how do we handle exceptions ?
// we catch them and call EASErrorHandler - this stores the exception
// inside a static "error". We also set another static "handled" which
// indicates if EAS have dealt with the error. We also increment yet
// another static "muppet" which shows how may exceptions occured since
// EAS bothered to handle it.
// We use this draconian approach to force EAS to handle errors rather
// than blithely carry on - no functions in this file will work unless
// "ok2run" returns "true" which indicates that all exceptions have been
// handled

// initial error condition only cleared when library initialised 
static ModelException error = 
ModelException("Library not initalised yet (Call EdrStartup)");
static bool handled = false;
static int  muppet  = 0;
static bool turnedOn = false; // true when library initialised
static bool turnedOff = false; // true when library shut down

static void EASErrorHandler(exception& e, const char* message) {

    error   = ModelException(e, message);
    handled = false;  // error occured, not dealt with yet
}

// has any error been handled, so we're good to go
static bool ok2run(const char* method) {
    if (!handled) {
        try {
            muppet++;
            string msg = "can't run " + string(method) + 
                " as an error has not been handled";
            error.addMsg(msg);

            msg = "the library has been called " + Format::toString(muppet) + 
                " times since the last error was ignored";

            error.addMsg(msg);
        }
        catch (exception& ) {
            // do nowt
        }        
    }
    return handled;
}

/** Returns true if there are no errors outstanding */
DLL_EXPORT EDRBool EdrErrorCleared(){
    return handled;
}

/** return stack trace as a dynamic C style string - 
    use free() to cleanup */
DLL_EXPORT char* EdrErrorStackTrace() {
    handled = true;  // tried to deal with error
    muppet = 0;
    try {
        return error.stackTrace();
    }
    catch (exception& ) {
        return 0;
    }
}
/** VITAL to call this before DOING ANYTHING AT ALL - return TRUE/FALSE */
DLL_EXPORT EDRBool EdrStartup(void)
{
    static const char routine[] = "EdrStartup";
    
    try {
        Library::startup();
        // ensure all symbols are linked in
        AddinLib::linkInAllLibs();

        if (turnedOff && !handled) {
            // foil attempt to circumvent error handing
            throw ModelException(routine, 
                                 "DR library has been shut down with "
                                 "unhandled errors");
        }

        handled = true;
        turnedOn = true;
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}

DLL_EXPORT void EdrShutdown(void)
{
    try {
        if (turnedOn){
            Library::shutdown();
            turnedOff = true;
        }
    } catch (exception& ) {
        // and do what exactly ?
    }
}

/** return version as a dynamic C style string - 
    use free() to cleanup */
DLL_EXPORT char* EdrVersion() {
    const static char routine[] = "EdrVersion";

    if (!ok2run(routine)) {
        return 0;
    }
    try {
        return (stringCopy(CVersion::DRLibVersion().c_str()));

    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return 0;
    }
}


/************* Creation of Atomic Objects  *********************/

/** Creates a wrapper to an integer. Returns NULL on failure */
DLL_EXPORT EDRObject EdrIntNew(int value){
    static const char routine[] = "EdrIntNew";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* obj;
    try{
        obj = new ObjectWrapper(CInt::create(value)); 
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}

/** Creates a wrapper to a double. Returns NULL on failure */
DLL_EXPORT EDRObject EdrDoubleNew(double value){
    static const char routine[] = "EdrDoubleNew";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* obj;
    try{
        obj = new ObjectWrapper(CDouble::create(value)); 
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}

/** Creates a wrapper to a string. Returns NULL on failure */
DLL_EXPORT EDRObject EdrStringNew(const char* value){
    static const char routine[] = "EdrStringNew";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* obj;
    try{
        if (!value){
            throw ModelException(routine, "NULL string provided");
        }
        obj = new ObjectWrapper(CString::create(value)); 
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}

/** Creates a wrapper to a bool. value must be 0 or 1. Returns NULL on
    failure */
DLL_EXPORT EDRObject EdrBoolNew(int value){
    const static char routine[] = "EdrBoolNew";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* obj;
    try{
        bool boolValue;
        if (value == 0){
            boolValue = false;
        }else if (value == 1){
            boolValue = true;
        } else {
            throw ModelException(routine, Format::toString(value)+" passed for"
                                 " boolean value");
        }

        obj = new ObjectWrapper(CBool::create(boolValue)); 
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}

/** Creates a wrapper to a double matrix. The matrix is created
    with all elements set to 0. Returns NULL on failure */
DLL_EXPORT EDRObject EdrDoubleMatrixEmptyNew(int             numCols,
                                             int             numRows){
    const static char routine[] = "EdrDoubleMatrixEmptyNew";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* obj;
    try{
        obj = new ObjectWrapper(new CDoubleMatrix(numCols, numRows)); 
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}

/** Sets a value in a double matrix created with via
    EdrDoubleMatrixEmptyNew */
DLL_EXPORT EDRBool EdrDoubleMatrixSet(EDRObject       matrix,
                                      int             col,
                                      int             row,
                                      double          value){
    const static char routine[] = "EdrDoubleMatrixSet";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }

    try {
        if (!matrix){
            throw ModelException(routine, "NULL object supplied");
        }
        ObjectWrapper* wrapper = (ObjectWrapper*)matrix;
        CDoubleMatrixSP matrixSP(
            CDoubleMatrixSP::dynamicCast(wrapper->object));
        (*matrixSP)[col][row] = value;
    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}


/** Creates a datetime object. Date is in dd-mmm-yyyy format.
    Time is hh:mm:ss or "SOD" or "EOD". Returns NULL on failure */
DLL_EXPORT EDRObject EdrDateTimeNew(const char* date, const char* time){
    const static char routine[] = "EdrDateTimeNew";
    ObjectWrapper* obj;

    if (!ok2run(routine)) {
        return 0;
    }

    try{
        obj = new ObjectWrapper(new DateTime(date, time)); 
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}

/** Creates an expiry which corresponds to a fixed date. Date is
    in dd-mmm-yyyy format.  Time is hh:mm:ss or "SOD" or "EOD".
    Returns NULL on failure */
DLL_EXPORT EDRObject EdrBenchmarkDateNew(const char* date, const char* time){
    const static char routine[] = "EdrBenchmarkDateNew";
    ObjectWrapper* obj;

    if (!ok2run(routine)) {
        return 0;
    }

    try{
        obj = new ObjectWrapper(new BenchmarkDate(DateTime(date, time))); 
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}

/** Creates a wrapper to an enum which is specifed by the type of the
    enum and a string (null terminated) representing the value to
    use. (No references are kept to either strings) . Returns NULL on
    failure.  A list of possible values can be obtained via
    EdrEnumValues. */
DLL_EXPORT EDRObject EdrEnumNew(const char* enumType,
                                const char* value){
    static const char routine[] = "EdrEnumNew";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* obj;
    try{
        if (!enumType){
            throw ModelException(routine, "NULL enumType provided");
        }
        if (!value){
            throw ModelException(routine, "NULL value provided");
        }
        // turn enumType into CClass
        CClassConstSP clazz = CClass::forName(string(enumType));
        obj = new ObjectWrapper(Enum::create(clazz, string(value))); 
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}

/** Creates an expiry (e.g. "6M"). Returns NULL on failure */
DLL_EXPORT EDRObject EdrMaturityPeriodNew(const char* period){
    const static char routine[] = "EdrMaturityPeriodNew";
    ObjectWrapper* obj;

    if (!ok2run(routine)) {
        return 0;
    }

    try{
        obj = new ObjectWrapper(new MaturityPeriod(period)); 
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}

/** Creates an expiry (e.g. "6M") with a fixed time 
    (either hh-mm-ss or "SOD" or "EOD"). Returns NULL on failure */
DLL_EXPORT EDRObject EdrMaturityTimePeriodNew(const char* period, 
                                              const char* time)
{
    const static char routine[] = "EdrMaturityTimePeriodNew";
    ObjectWrapper* obj;

    if (!ok2run(routine)) {
        return 0;
    }

    try{
        int tm = DateTime::timeConvert(time);
        obj = new ObjectWrapper(new MaturityTimePeriod(period, tm)); 
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}

/************ Decomposing Atomic Objects *******************/

DLL_EXPORT EDRBool EdrDoubleGet(const EDRObject object,
                                double*         d)
{
    const static char routine[] = "EdrDoubleGet";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }

    *d = 0.0;

    try {
        if (!object){
            throw ModelException(routine, "NULL object supplied");
        }
        ObjectWrapper* wrapper = (ObjectWrapper*)object;
        CDoubleSP dblSP(CDoubleSP::dynamicCast(wrapper->object));

        *d = dblSP->doubleValue();
    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}

/** returns dynamic string - use free */
DLL_EXPORT char* EdrExpiryGet(const EDRObject object)
{
    const static char routine[] = "EdrExpiryGet";

    if (!ok2run(routine)) {
        return 0;
    }

    try {
        if (!object){
            throw ModelException(routine, "NULL object supplied");
        }
        ObjectWrapper* wrapper = (ObjectWrapper*)object;
        ExpirySP expirySP(ExpirySP::dynamicCast(wrapper->object));

        return (stringCopy(expirySP->toString().c_str()));

    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return 0;
    }
}

DLL_EXPORT char* EdrDateTimeGet(const EDRObject object)
{
    const static char routine[] = "EdrDateTimeGet";

    if (!ok2run(routine)) {
        return 0;
    }

    try {
        if (!object){
            throw ModelException(routine, "NULL object supplied");
        }
        ObjectWrapper* wrapper = (ObjectWrapper*)object;
        DateTimeSP dtSP(DateTimeSP::dynamicCast(wrapper->object));

        return (stringCopy(dtSP->toString().c_str()));

    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return 0;
    }
}

DLL_EXPORT char* EdrStringGet(const EDRObject object)
{
    const static char routine[] = "EdrStringGet";

    if (!ok2run(routine)) {
        return 0;
    }

    try {
        if (!object){
            throw ModelException(routine, "NULL object supplied");
        }
        ObjectWrapper* wrapper = (ObjectWrapper*)object;
        CStringSP sSP(CStringSP::dynamicCast(wrapper->object));

        return (stringCopy(sSP->stringValue().c_str()));

    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return 0;
    }
}

DLL_EXPORT EDRBool EdrBoolGet(const EDRObject object,
                              EDRBool*        b)
{
    const static char routine[] = "EdrBoolGet";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }

    *b = EDR_FALSE;

    try {
        if (!object){
            throw ModelException(routine, "NULL object supplied");
        }
        ObjectWrapper* wrapper = (ObjectWrapper*)object;
        CBoolSP bSP(CBoolSP::dynamicCast(wrapper->object));

        *b = bSP->boolValue();
    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}

EDRBool EdrIntGet(const EDRObject object,
                  int*            i)
{
    const static char routine[] = "EdrIntGet";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }

    *i = 0;

    try {
        if (!object){
            throw ModelException(routine, "NULL object supplied");
        }
        ObjectWrapper* wrapper = (ObjectWrapper*)object;
        CIntSP iSP(CIntSP::dynamicCast(wrapper->object));

        *i = iSP->intValue();
    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}

/** Returns the number of rows and columns of the matrix. The
    supplied object must be a wrapped double matrix (see
    EdrDoubleMatrixEmptyNew). Returns EDR_FALSE on failure*/
DLL_EXPORT EDRBool EdrDoubleMatrixGetSize(const EDRObject object,   /* (I) */
                                          int*            numCols,  /* (O) */
                                          int*            numRows){ /* (O) */
    const static char routine[] = "EdrDoubleMatrixGetSize";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }

    try {
        if (!object){
            throw ModelException(routine, "NULL object supplied");
        }
        ObjectWrapper* wrapper = (ObjectWrapper*)object;
        CDoubleMatrixSP
            matrixSP(CDoubleMatrixSP::dynamicCast(wrapper->object));
        *numCols = matrixSP->numCols();
        *numRows = matrixSP->numRows();
    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}

/** Returns an element of a double matrix together. The supplied
    object must be a wrapped double matrix (see
    EdrDoubleMatrixNew). Returns EDR_FALSE on failure */
DLL_EXPORT EDRBool EdrDoubleMatrixGet(const EDRObject object,  /* (I) */
                                      int             col,     /* (I) */
                                      int             row,     /* (I) */
                                      double*         value){  /* (O) */
    const static char routine[] = "EdrDoubleMatrixGet";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }

    try {
        if (!object){
            throw ModelException(routine, "NULL object supplied");
        }
        ObjectWrapper* wrapper = (ObjectWrapper*)object;
        CDoubleMatrixSP
            matrixSP(CDoubleMatrixSP::dynamicCast(wrapper->object));
        *value = (*matrixSP)[col][row];
    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}

/** returns a const char* style string (do NOT use free) representing
    the value of the supplied enum which must be a wrapped enum
    (see EdrEnumNew). Returns NULL on failure */
DLL_EXPORT const char* EdrEnumGet(const EDRObject object)
{
    const static char routine[] = "EdrEnumGet";

    if (!ok2run(routine)) {
        return 0;
    }

    try {
        if (!object){
            throw ModelException(routine, "NULL object supplied");
        }
        ObjectWrapper* wrapper = (ObjectWrapper*)object;
        EnumSP enumSP(EnumSP::dynamicCast(wrapper->object));

        return (enumSP->enumValueAsString().c_str());

    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return 0;
    }
}

/************* Creation of Array Objects  *********************/

/** Creates an empty array of object whose components are of
    componentType. The length of the array created is determined by
    numElements. Returns NULL on failure. */
DLL_EXPORT EDRObject EdrArrayNew(const char* componentType,
                                 int         numElements){
    const static char routine[] = "EdrArrayNew";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* obj = 0;
    try{
        if (!componentType){
            throw ModelException(routine, "NULL componentType");
        }
        // look up class of componentType
        CClassConstSP componentClass= CClass::forName(string(componentType));
        if (!Modifier::isPublic(componentClass->getModifiers())){
            throw ModelException(routine, string(componentType) + " is not a "
                                 "public class - instantiation disallowed");
        }
        IArraySP array(componentClass->
                       newArrayInstanceByComponent(numElements));
        obj = new ObjectWrapper(array);
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}

/** Sets the supplied object in the supplied array at the
    specified index location (index must lie in [0,numElements-1]).
    Returns EDR_FALSE on failure */
DLL_EXPORT EDRBool EdrArraySet(EDRObject array,    /* (M) */
                               int       index,    /* (I) */
                               EDRObject object){  /* (I) */
    const static char routine[] = "EdrArraySet";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }
    try{
        if (!array){
            throw ModelException(routine, "NULL array");
        }
        if (!object){
            throw ModelException(routine, "NULL object");
        }

        ObjectWrapper* arrayWrapper = (ObjectWrapper*)array;
        ObjectWrapper* objWrapper   = (ObjectWrapper*)object;
        IArraySP array(IArraySP::dynamicCast(arrayWrapper->object));

        IObjectSP objToSet(objWrapper->object);
        //// convert any public objects to private objects etc if needed
        CObject::checkType(objToSet, array->getClass()->getComponentType());
        array->set(index, objToSet);
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}

/** Appends the supplied object to the end of supplied array. The
    length of the array is increased by 1. (aka push_back). 
    Returns EDR_FALSE on failure*/
DLL_EXPORT EDRBool EdrArrayAppend(EDRObject array,
                                  EDRObject object){
    const static char routine[] = "EdrArrayAppend";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }
    try{
        if (!array){
            throw ModelException(routine, "NULL array");
        }
        if (!object){
            throw ModelException(routine, "NULL object");
        }
        ObjectWrapper* arrayWrapper = (ObjectWrapper*)array;
        ObjectWrapper* objWrapper   = (ObjectWrapper*)object;
        IArraySP array(IArraySP::dynamicCast(arrayWrapper->object));
        IObjectSP objToAppend(objWrapper->object);
        //// convert any public objects to private objects etc if needed
        CObject::checkType(objToAppend, array->getClass()->getComponentType());
        array->append(objToAppend);
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}


/** how long is an array ? */
DLL_EXPORT EDRBool EdrArrayLength(const EDRObject object,
                                  int*            length)
{
    const static char routine[] = "EdrArrayLength";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }
    *length = 0;
    try{
        if (!object){
            throw ModelException(routine, "NULL object");
        }

        ObjectWrapper* objWrapper = (ObjectWrapper*)object;
        
        IArraySP array(IArraySP::dynamicCast(objWrapper->object));

        *length = array->getLength();
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}

/** get the i'th item from an array (starts at 0) */
DLL_EXPORT EDRObject EdrArrayItem(EDRObject object,
                                  int       index)
{
    const static char routine[] = "EdrArrayItem";

    if (!ok2run(routine)) {
        return 0;
    }
    ObjectWrapper* obj = 0;
    try{
        if (!object){
            throw ModelException(routine, "NULL object");
        }

        ObjectWrapper* objWrapper = (ObjectWrapper*)object;
        
        IArraySP array(IArraySP::dynamicCast(objWrapper->object));
        IObjectSP theElem(array->get(index));
        // make sure we don't give out private objects
        if (IPrivateObject::TYPE->isInstance(theElem)){
            theElem = CObject::convertToPublicRep(theElem);
        }
        obj = new ObjectWrapper(theElem); 
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}

/************* DataDictionary Support *********************/

/** Creates a data dictionary corresponding to the object of the
    supplied type. Returns NULL on failure */
DLL_EXPORT EDRDataDict EdrDataDictNew(const char* type){
    const static char routine[] = "EdrDataDictNew";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper*     obj = 0;
    try{
        if (!type){
            throw ModelException(routine, "NULL type supplied");
        }
        obj = new ObjectWrapper(CDataDictionary::create(type)); 
    } catch (exception& e){
        EASErrorHandler(e, "EdrDataDictNew");
        return 0;
    }
    return (EDRDataDict) obj;
}

/** Adds the specified object to the supplied data dictionary. Returns
    TRUE/FALSE. */
DLL_EXPORT EDRBool EdrDataDictAdd(EDRDataDict  dataDict,
                                  const char*  fieldName,
                                  EDRObject    objectToAdd)
{
    const static char routine[] = "EdrDataDictAdd";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }

    try{
        if (!dataDict){
            throw ModelException(routine, "NULL dataDict");
        }
        if (!fieldName){
            throw ModelException(routine, "fieldName is NULL");
        }

        ObjectWrapper* ddWrapper  = (ObjectWrapper*)dataDict;
        // allow NULL objectToAdd to flow down to underlying method - allows
        // better error messages. ie objWrapper may be null here
        ObjectWrapper* objWrapper = (ObjectWrapper*)objectToAdd;
        
        // have to support dictionaries & DR Wrappers via same interface
        if (CDataDictionary::TYPE->isInstance(ddWrapper->object)){
            CDataDictionarySP dd(CDataDictionarySP::dynamicCast(
                ddWrapper->object));
            dd->put(fieldName, objWrapper? objWrapper->object: IObjectSP(   ));
        }
        else {
            DRWrapperSP dd(DRWrapperSP::dynamicCast(ddWrapper->object));
            dd->put(fieldName, objWrapper? objWrapper->object: IObjectSP(   ));
        }

    } catch (exception& e){
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}

/** Converts supplied data dictionary into its object equivalent */
DLL_EXPORT EDRObject EdrDataDictToObject(EDRDataDict dataDict)
{
    const static char routine[] = "EdrDataDictToObject";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* obj = 0;
    try{
        if (!dataDict){
            throw ModelException(routine, "NULL dataDict");
        }
        ObjectWrapper* ddWrapper = (ObjectWrapper*)dataDict;
        if (!CDataDictionary::TYPE->isInstance(ddWrapper->object)){
            throw ModelException(routine,
                                 "Can't convert DRWrapper to an object");
        }
        CDataDictionarySP dd(CDataDictionarySP::dynamicCast(
            ddWrapper->object));
        obj = new ObjectWrapper(dd->pop2Object());
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}

/** same as above but no deep copies */
DLL_EXPORT EDRObject EdrDataDictToObjectShallow(EDRDataDict dataDict)
{
    const static char routine[] = "EdrDataDictToObjectShallow";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* obj = 0;
    try{
        if (!dataDict){
            throw ModelException(routine, "NULL dataDict");
        }
        ObjectWrapper* ddWrapper = (ObjectWrapper*)dataDict;
        if (!CDataDictionary::TYPE->isInstance(ddWrapper->object)){
            throw ModelException(routine,
                                 "Can't convert DRWrapper to an object");
        }
        CDataDictionarySP dd(CDataDictionarySP::dynamicCast(
            ddWrapper->object));
        dd->copyOnPop2Obj(false);
        obj = new ObjectWrapper(dd->pop2Object());
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}

/** and the other way round */
DLL_EXPORT EDRDataDict EdrObjectToDataDict(EDRObject object)
{
    const static char routine[] = "EdrObjectToDataDict";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* obj = 0;
    try{
        if (!object){
            throw ModelException(routine, "NULL object");
        }
        ObjectWrapper* wrapper = (ObjectWrapper*)object;
        obj = new ObjectWrapper(CDataDictionary::pop2DataDict(wrapper->object));
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRDataDict) obj;
}

/** Get an object out of a data dictionary */
DLL_EXPORT EDRObject EdrDataDictGet(const EDRDataDict  dataDict,
                                    const char*        fieldName)
{
    const static char routine[] = "EdrDataDictGet";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* obj = 0;

    try{
        if (!dataDict){
            throw ModelException(routine, "NULL dataDict");
        }
        if (!fieldName){
            throw ModelException(routine, "fieldName is NULL");
        }

        ObjectWrapper* ddWrapper = (ObjectWrapper*)dataDict;
        
        if (!CDataDictionary::TYPE->isInstance(ddWrapper->object)){
            throw ModelException(routine, "Not supported for DRWrappers");
        }
        CDataDictionarySP dd(CDataDictionarySP::dynamicCast(
            ddWrapper->object));
        obj = new ObjectWrapper(dd->get(fieldName));


    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}


/** Creates an EDRDataDictIterator object (which must eventually be freed)
    This object can then be used to explore the contents of a data
    dictonary */
DLL_EXPORT EDRDataDictIterator EdrDataDictIteratorGet(EDRDataDict  dataDict){
    const static char routine[] = "EdrDataDictIteratorGet";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* obj = 0;

    try{
        if (!dataDict){
            throw ModelException(routine, "NULL dataDict");
        }
        ObjectWrapper* ddWrapper = (ObjectWrapper*)dataDict;
        
        if (!CDataDictionary::TYPE->isInstance(ddWrapper->object)){
            throw ModelException(routine, "Iteration not supported for "
                                 "DR wrappers");
        }
        CDataDictionarySP dd(CDataDictionarySP::dynamicCast(
            ddWrapper->object));
        obj = new ObjectWrapper(new DataDictIterator(dd));
        
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}

/** Returns next <field, object> pair in iteration - field and
    object will be null when there are no more elements in the
    iteration. Neither the field nor the object should be freed */
DLL_EXPORT EDRBool EdrDataDictIteratorNext(
    EDRDataDictIterator iterator, /* (M) */
    const char**        field,    /* (O) */
    EDRObject*          object){  /* (O) */
    const static char routine[] = "EdrDataDictIteratorNext";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }
    try{
        if (!iterator){
            throw ModelException(routine, "NULL iterator");
        }
        ObjectWrapper* iterWrapper = (ObjectWrapper*)iterator;
        DataDictIteratorSP iter(
            DataDictIteratorSP::dynamicCast(iterWrapper->object));
        if (iter->iterator == iter->iteratorEnd){
            *field = 0;
            *object = 0;
        } else {
            *field  = iter->iterator->first->getName().c_str();
            iter->objWrapper = ObjectWrapper(iter->iterator->second);
            *object = (EDRObject)&iter->objWrapper;
            ++iter->iterator;
        }
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}
        

/** Returns a non dynamic string giving the runtime class of the
    object */
DLL_EXPORT const char* EdrObjectGetType(EDRObject object){
    const static char routine[] = "EdrObjectGetType";
 
    if (!ok2run(routine)) {
        return 0;
    }

    const char* typeName = 0;
    try{
        if (!object){
            throw ModelException(routine, "NULL object");
        }
        ObjectWrapper* objWrapper = (ObjectWrapper*)object;
        typeName = objWrapper->object->getClass()->getName().c_str();
    }catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return typeName;
}

/** Returns a null terminated array of const char * containing the
    types within the library - each element in the array gives the
    name of a type. Use free to release the memory used by the
    array. Do not free the elements of the array */
DLL_EXPORT const char** EdrTypeList(void){
    const static char routine[] = "EdrTypeList";

    if (!ok2run(routine)) {
        return 0;
    }

    const char** types = 0;
    try {
        const CClassVec& allClasses = CClass::allClasses();
        types = NEW_ARRAY(const char*, allClasses.size()+1);
        int j = 0;
        for (unsigned int i = 0; i < allClasses.size(); i++){
            // hide private types from EAS
            if (!Modifier::isPrivate(allClasses[i]->getModifiers())){
                types[j] = allClasses[i]->getName().c_str();
                j++;
            }
        }
    } catch (exception& e) {
        FREE(types);
        EASErrorHandler(e, routine);
        return 0;
    }
    return types;
}
    
/** List all types that can be used to build (via the data dictionary
    route) an object of type type or an object derived from type type */
DLL_EXPORT const char** EdrTypeListConstructorTypes(const char* type){
    const static char routine[] = "EdrTypeListConstructorTypes";

    if (!ok2run(routine)) {
        return 0;
    }

    const char** types = 0;
    try {
        if (!type){
            throw ModelException(routine, "NULL type supplied");
        }
        CClassConstSP clazz = CClass::forName(string(type));
        CClassVec  allClasses = clazz->listAllConstructorClasses();
        types = NEW_ARRAY(const char*, allClasses.size()+1);
        for (unsigned int i = 0; i < allClasses.size(); i++){
            types[i] = allClasses[i]->getName().c_str();
        }
    } catch (exception& e) {
        FREE(types);
        EASErrorHandler(e, routine);
        return 0;
    }
    return types;
}
    
/** Returns, via the modifiers parameter, the modifiers for this
    class or interface, encoded in an integer. The values
    ABSTRACT, INTERFACE, PROTECTED and PUBLIC should be used to
    decode the integer using bitwise or.  Returns TRUE/FALSE for
    success/failure */
DLL_EXPORT EDRBool EdrTypeGetModifiers(const char* type,       /* (I) */
                                       int*        modifiers){ /* (O) */
    const static char routine[] = "EdrTypeGetModifiers";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }

    *modifiers = 0;

    try {
        if (!type){
            throw ModelException(routine, "NULL type supplied");
        }
        CClassConstSP clazz = CClass::forName(string(type));
        *modifiers = clazz->getModifiers();
    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}
    
/** find out how many items are in a type description */
DLL_EXPORT EDRBool EdrTypeNumItems(
    const char* type,
    int*        numComponents)     /* (O) */
{
    const static char routine[] = "EdrTypeNumItems";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }

    *numComponents = 0;

    try {
        if (!type){
            throw ModelException(routine, "NULL type supplied");
        }
        CClassConstSP clazz = CClass::forName(string(type));
        CFieldArray fields(Addin::getDataClassFields(clazz));
        *numComponents = fields.size();
    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}

/** return the field name of the i'th item in a type description 
    do NOT free the output of this function */
DLL_EXPORT const char* EdrTypeFieldName(
    const char* type,
    int         index)
{
    const static char routine[] = "EdrTypeFieldName";

    if (!ok2run(routine)) {
        return 0;
    }

    try {
        if (!type){
            throw ModelException(routine, "NULL type supplied");
        }
        CClassConstSP clazz = CClass::forName(string(type));
        CFieldArray fields(Addin::getDataClassFields(clazz));
        
        if (index < 0 || index >= (int)fields.size()) {
            throw ModelException(routine, "index out of range");
        }

        return (fields[index]->getName().c_str());

    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return 0;
    }
}

/** return the type name of the i'th item in a type description 
    do NOT free the output of this function */
DLL_EXPORT const char* EdrTypeFieldType(
    const char* type,
    int         index)
{
    const static char routine[] = "EdrTypeFieldType";

    if (!ok2run(routine)) {
        return 0;
    }

    try {
        if (!type){
            throw ModelException(routine, "NULL type supplied");
        }
        CClassConstSP clazz = CClass::forName(string(type));
        CFieldArray fields(Addin::getDataClassFields(clazz));
        
        if (index < 0 || index >= (int)fields.size()) {
            throw ModelException(routine, "index out of range");
        }

        return (fields[index]->getType()->getName().c_str());

    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return 0;
    }
}

/** is a field in a type optional ? */
DLL_EXPORT EDRBool EdrTypeFieldIsOptional(
    const char* type,
    int         index,
    EDRBool*    optional)     /* (O) */
{
    const static char routine[] = "EdrTypeFieldIsOptional";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }

    *optional = EDR_FALSE;

    try {
        if (!type){
            throw ModelException(routine, "NULL type supplied");
        }
        CClassConstSP clazz = CClass::forName(string(type));
        CFieldArray fields(Addin::getDataClassFields(clazz));
        
        if (index < 0 || index >= (int)fields.size()) {
            throw ModelException(routine, "index out of range");
        }

        *optional = fields[index]->isOptional()? EDR_TRUE: EDR_FALSE;

    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}


/** return the description of the i'th item in a type description 
    do NOT free the output of this function */
DLL_EXPORT const char* EdrTypeFieldDescription(
    const char* type,
    int         index)
{
    const static char routine[] = "EdrTypeFieldDescription";

    if (!ok2run(routine)) {
        return 0;
    }

    try {
        if (!type){
            throw ModelException(routine, "NULL type supplied");
        }
        CClassConstSP clazz = CClass::forName(string(type));
        CFieldArray fields(Addin::getDataClassFields(clazz));
        
        if (index < 0 || index >= (int)fields.size()) {
            throw ModelException(routine, "index out of range");
        }

        return (fields[index]->getDescription().c_str());

    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return 0;
    }
}

/** is a type an array ? */
DLL_EXPORT EDRBool EdrTypeIsArray(
    const char* type,
    EDRBool*    isArray)     /* (O) */
{
    const static char routine[] = "EdrTypeIsArray";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }

    *isArray = EDR_FALSE;

    try {
        if (!type){
            throw ModelException(routine, "NULL type supplied");
        }
        CClassConstSP clazz = CClass::forName(string(type));

        *isArray = clazz->isArray()? EDR_TRUE: EDR_FALSE;

    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}

/** what is the type of element in an array ? 
    do NOT free the output of this function */
DLL_EXPORT const char* EdrArrayElemType(const char* type)
{
    const static char routine[] = "EdrArrayElemType";

    if (!ok2run(routine)) {
        return 0;
    }

    try {
        if (!type){
            throw ModelException(routine, "NULL type supplied");
        }
        CClassConstSP clazz = CClass::forName(string(type));
        CClassConstSP arrayClazz = clazz->getComponentType();

        return (arrayClazz->getName().c_str());

    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return 0;
    }
}

/** Does a type represent an enum? The return value indicates
    SUCCESS/FAILURE */
DLL_EXPORT EDRBool EdrTypeIsEnum(
    const char* type,
    EDRBool*    isEnum)     /* (O) */
{
    const static char routine[] = "EdrTypeIsEnum";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }

    *isEnum = EDR_FALSE;

    try {
        if (!type){
            throw ModelException(routine, "NULL type supplied");
        }
        CClassConstSP clazz = CClass::forName(string(type));

        *isEnum = clazz->isEnum()? EDR_TRUE: EDR_FALSE;

    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}

/** Returns a null terminated array of const char * containing either
    the possible values for the specified enum type (getDescriptions =
    FALSE) or the descriptions for each of these values
    (getDescriptions = TRUE). When getDescriptions = FALSE, each
    element in the array gives a possible value with which an enum of
    the specified type can be built. 
    For either value of getDescriptions, use free to release the memory
    used by the array but do not free the elements of the array.  Returns
    NULL on failure. */
DLL_EXPORT const char** EdrEnumValues(
    const char* enumType,   /* (I) */
    EDRBool     getDescriptions)     /* (I) */
{
    const static char routine[] = "EdrEnumValues";

    if (!ok2run(routine)) {
        return 0;
    }

    const char** values = 0;
    try {
        if (!enumType){
            throw ModelException(routine, "NULL type supplied");
        }
        CClassConstSP clazz = CClass::forName(string(enumType));
        const vector<CClass::EnumValue>* enumValues = clazz->getEnumValues();
        if (!enumValues){
            throw ModelException(routine, "Type '"+string(enumType)+
                                 " is not an enum");
        }
        values = NEW_ARRAY(const char*, enumValues->size()+1);
        for (unsigned int i = 0; i < enumValues->size(); i++){
            values[i] = getDescriptions?
                (*enumValues)[i].comment.c_str():
                (*enumValues)[i].valueAsString.c_str();
        }
    } catch (exception& e) {
        FREE(values);
        EASErrorHandler(e, routine);
        return 0;
    }
    return values;
}

/** Determines if an object whose type is given by typeOfObject is
    derived from typeOfClass */
DLL_EXPORT EDRBool EdrTypeIsAssignableFrom(const char* typeOfClass,
                                           const char* typeOfObject,
                                           EDRBool*    assignable){
    const static char routine[] = "EdrTypeIsAssignableFrom";
    if (!ok2run(routine)) {
        return EDR_FALSE;
    }
    try {
        if (!typeOfClass){
            throw ModelException(routine, "NULL typeOfClass supplied");
        }
        if (!typeOfObject){
            throw ModelException(routine, "NULL typeOfObject supplied");
        }
        CClassConstSP typeClazz = CClass::forName(string(typeOfClass));
        CClassConstSP objClazz = CClass::forName(string(typeOfObject));
        *assignable = typeClazz->isAssignableFrom(objClazz)? 
            EDR_TRUE: EDR_FALSE;
        return EDR_TRUE;
    } catch (exception& e) {
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
}
    

/************* Market Data Cache support ******************/

/** Creates a new cache of market data using the supplied string as
    its id */
DLL_EXPORT EDRMarketDataCache EdrMarketDataCacheNew(EDRObject   valueDate){
    const static char routine[] = "EdrMarketDataCacheNew";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper*     obj = 0;
    try{
        if (!valueDate){
            throw ModelException(routine, "NULL valueDate");
        }
        ObjectWrapper* date = (ObjectWrapper*) valueDate;
        DateTimeSP     dateSP(DateTimeSP::dynamicCast(date->object));
        obj = new ObjectWrapper(new MarketData(*dateSP));
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRMarketDataCache) obj;
}


DLL_EXPORT EDRBool EdrMarketDataCacheAdd(EDRMarketDataCache cache,
                                         EDRObject          objectToAdd){
    const static char routine[] = "EdrMarketDataCacheAdd";

    if (!ok2run(routine)) {
        return 0;
    }

    try{
        if (!cache){
            throw ModelException(routine, "NULL cache");
        }
        if (!objectToAdd){
            throw ModelException(routine, "NULL objectToAdd");
        }
        ObjectWrapper* cacheWrapper = (ObjectWrapper*)cache;
        ObjectWrapper* objWrapper = (ObjectWrapper*)objectToAdd;
        CMarketDataSP cacheSP(CMarketDataSP::dynamicCast(
            cacheWrapper->object));

        cacheSP->put("", objWrapper->object);
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}

/** Creates a copy of a market data cache. Performance wise: the
     contents of the cache are not deep copied, ie references are
     taken to the data within the cache */
DLL_EXPORT EDRMarketDataCache EdrMarketDataCacheClone(
    EDRMarketDataCache cache){
    const static char routine[] = "EdrMarketDataCacheClone";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper*     obj = 0;
    try{
        if (!cache){
            throw ModelException(routine, "NULL cache");
        }
        ObjectWrapper* cacheWrapper = (ObjectWrapper*)cache;
        CMarketDataSP cacheSP(CMarketDataSP::dynamicCast(
            cacheWrapper->object));
        CMarketDataSP copy(cacheSP.clone());
        obj = new ObjectWrapper(copy);

    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRMarketDataCache) obj;
}


/************* DR Wrapper support *********************/
/** Creates a DR Wrapper corresponding to the object of the
    supplied type. Returns NULL on failure */
DLL_EXPORT EDRDataDict EdrDRWrapperNew(const char* type) {
    const static char routine[] = "EdrDRWrapperNew";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* obj = 0;
    try {
        if (!type) {
            throw ModelException(routine, "NULL type supplied");
        }
        obj = new ObjectWrapper(DRWrapper::create(type)); 
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRDataDict) obj;
}

/** Write inputs to file */
DLL_EXPORT EDRBool EdrWriteInputsToFile(const char*          filename,
                                        int                  size,
                                        EDRObject*           models,
                                        EDRObject*           instruments,
                                        EDRObject*           controls,
                                        double*              multipliers,    
                                        double*              weights,
                                        EDRObject            scenario,
                                        EDRMarketDataCache   cache) {
    const static char routine[] = "EdrWriteInputsToFile";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }

    try {       
        CompositeInstrumentSP composite = makeComposite(size,
                                                        models,
                                                        instruments,
                                                        controls,
                                                        multipliers,
                                                        weights,
                                                        scenario,
                                                        cache);

        XMLWriter xml(filename);

        composite->write("COMPOSITE", &xml);

    } catch (exception& e){
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}

/** Writes the supplied object to file. */
DLL_EXPORT EDRBool EdrWriteResultsToFile(const char*         filename,
                                         EDRObject            results){
    const static char routine[] = "EdrWriteResultsToFile";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }

    try {       
        ObjectWrapper* objWrapper = (ObjectWrapper*)results;
        OutputFile   output(filename);
        output.write(objWrapper->object.get());
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}

/** Writes any object to file */
DLL_EXPORT EDRBool EdrWriteObjectToFile(const char* filename,
                                        EDRObject   object){
    const static char routine[] = "EdrWriteObjectToFile";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }

    try {   
        if (!filename){
            throw ModelException(routine, "filename is NULL");
        }
        if (!object){
            throw ModelException(routine, "object is NULL");
        }

        ObjectWrapper* objWrapper = (ObjectWrapper*)object;

        XMLWriter xml(filename);
        objWrapper->object->write("OBJECT", &xml);

    } catch (exception& e){
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}

/** returns object written to file using EdrWriteObjectToFile. 
    Also works for file created using EdrWriteInputsToFile - in this
    case the object returned will be of type "CompositeInstrument". Use
    EDR_TYPE_COMPONENTS("CompositeInstrument") in Excel to see
    what fields this object contains */
DLL_EXPORT EDRObject EdrReadObjectFromFile(const char* filename){
    const static char routine[] = "EdrReadObjectFromFile";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* results = 0;
    try {
        XMLReader reader(filename, true);
        IObjectSP o(reader.read());
        results = new ObjectWrapper(o);
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject)results;
}


/** read wrapper results back from file. Same outputs as EdrMain */
DLL_EXPORT EDRObject EdrDRWrapperRead(const char*  filename){        /* (I) */
    const static char routine[] = "EdrDRWrapperRead";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* results = 0;
    try {
        IObjectSP o(OutputFile::xmlReadResults(filename));
        results = new ObjectWrapper(o);
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject)results;
}

/** return wrapper errors from file as a dynamic C style string - 
    use free() to cleanup */
DLL_EXPORT char* EdrDRWrapperError(const char* filename) {  /* (I) */
    const static char routine[] = "EdrDRWrapperError";

    if (!ok2run(routine)) {
        return 0;
    }
    char * errors = 0;
    try {
        string err(OutputFile::xmlReadErrors(filename));
        errors = stringCopy(err.c_str());
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return errors;
}

/** Parses the supplied string and build the object that the xml describes.
    The xml must conform to the EDR representation */
DLL_EXPORT EDRObject EdrXMLRead(const char* xmlString){
    const static char routine[] = "EdrXMLRead";
    // usual checks
    if (!ok2run(routine)) {
        return 0;
    }
    ObjectWrapper* object = 0;
    try {
        // check inputs
        if (!xmlString){
            throw ModelException(routine, "xmlString is NULL");
        }
        // convert char* into string
        string buffer(xmlString);
        // build XML input stream
        XMLReader reader(buffer, false);
        // parse
        IObjectSP obj(reader.read());
        // create wrapper object
        object = new ObjectWrapper(obj);
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject)object;
}

/** Takes the supplied object and writes it, in xml format, to a string. The
    returned string is dynamic and must be freed using free() */
DLL_EXPORT char* EdrXMLWrite(const EDRObject object){
    const static char routine[] = "EdrXMLWrite";
    // usual checks
    if (!ok2run(routine)) {
        return 0;
    }
    const ObjectWrapper* wrapper = (const ObjectWrapper*)object;
    char*          output  = 0;
    try {
        // check inputs
        if (!wrapper || !wrapper->object){
            throw ModelException(routine, "NULL object supplied");
        }
        // create xml output stream
        string buffer;
        XMLWriter stream(buffer, false);
        // serialise
        wrapper->object->write("EDRObject", &stream);
        // get hold of string
        output = stringCopy(buffer.c_str());
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return output;
}


/************* Pricing/Sensitivities/Scenario support ******************/

/** Main Pricing Call */
DLL_EXPORT EDRObject EdrMain(int                  size,
                             EDRObject*           models,
                             EDRObject*           instruments,
                             EDRObject*           controls,
                             double*              multipliers,    
                             double*              weights,
                             EDRObject            scenario,
                             EDRMarketDataCache   marketData)
{
    const static char routine[] = "EdrMain";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* results = 0;
    try {
        
        CompositeInstrumentSP composite = makeComposite(size,
                                                        models,
                                                        instruments,
                                                        controls,
                                                        multipliers,
                                                        weights,
                                                        scenario,
                                                        marketData);

        results = new ObjectWrapper(composite->run());

    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject)results;
}

/** Generic function for all other actions (e.g. event handling).
    Create object of appropriate type, pass it in, and get results
    back
*/
DLL_EXPORT EDRObject EdrAction(EDRObject input)  /* (I) */
{
    const static char method[] = "EdrAction";

    if (!ok2run(method)) {
        return 0;
    }

    ObjectWrapper* output = 0;

    try {
        if (!input){
            throw ModelException(method, "NULL input");
        }

        ObjectWrapper* inputWrapper = (ObjectWrapper*)input;

        ClientRunnable* runner = dynamic_cast<ClientRunnable*>(inputWrapper->object.get());

        if (!runner){
            throw ModelException(method, "Input is of incorrect type ("+
                                 inputWrapper->object->getClass()->getName()+")"
                                 " Object must be derived from ClientRunnable");
        }

        output = new ObjectWrapper(runner->run());

    } catch (exception& e){
        EASErrorHandler(e, method);
        return 0;
    }
    return (EDRObject)output;
}


/************* Results support *********************/

/** Creates an EDRResultsIterator object (which must eventually be freed)
    This object can then be used to explore the contents of the results
    object */
DLL_EXPORT EDRResultsIterator EdrResultsIteratorGet(EDRObject  results){
    const static char routine[] = "EdrResultsIteratorGet";

    if (!ok2run(routine)) {
        return 0;
    }

    ObjectWrapper* obj = 0;

    try{
        if (!results){
            throw ModelException(routine, "NULL results");
        }
        ObjectWrapper* wrapper = (ObjectWrapper*)results;
        
        ResultsSP resultsSP(ResultsSP::dynamicCast(wrapper->object));
        ResultsIterator* iter = new ResultsIterator(resultsSP);
        obj = new ObjectWrapper(iter);
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return (EDRObject) obj;
}

/** Returns next <packet name, identifier, object> triple in
    iteration - packet, name and object will be null when there
    are no more elements in the iteration. Neither the packet nor
    the name nor the object should be freed */
DLL_EXPORT EDRBool EdrResultsIteratorNext(
    EDRResultsIterator   iterator, /* (M) */
    const char**         packet,   /* (O) */
    EDROutputName*       name,     /* (O) */
    EDRObject*           object){  /* (O) */
    const static char routine[] = "EdrResultsIteratorNext";

    if (!ok2run(routine)) {
        return EDR_FALSE;
    }

    try{
        if (!iterator){
            throw ModelException(routine, "NULL iterator");
        }
        ObjectWrapper* iterWrapper = (ObjectWrapper*)iterator;
        ResultsIteratorSP iter(
            ResultsIteratorSP::dynamicCast(iterWrapper->object));
        if (iter->packetIdx == iter->packets.size()){
            *packet = 0;
            *name = 0;
            *object = 0;
        } else {
            *packet  = iter->packets[iter->packetIdx]->c_str();
            OutputNameConstSP theName = 
                iter->packetContents[iter->contentsIdx].first;
            iter->nameWrapper =ObjectWrapper(OutputNameSP::constCast(theName));
            IObjectConstSP theObj = 
                iter->packetContents[iter->contentsIdx].second;
            iter->objWrapper = ObjectWrapper(IObjectSP::constCast(theObj));
            *name = (EDROutputName) &iter->nameWrapper;
            *object = (EDRObject) &iter->objWrapper;
            iter->increment();
        }
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}
 

/** Returns the number of strings which make up an output name */
DLL_EXPORT EDRBool EdrOutputNameIDCount(
    const EDROutputName outputName,  /* (I) */
    int*                numStrings){ /* (O) */
    const static char routine[] = "EdrOutputNameIDCount";
    if (!ok2run(routine)) {
        return EDR_FALSE;
    }
    try{
        if (!outputName){
            ModelException(routine, "NULL outputName");
        }
        ObjectWrapper* nameWrapper = (ObjectWrapper*)outputName;
        OutputNameSP nameSP(OutputNameSP::dynamicCast(nameWrapper->object));
        *numStrings = nameSP->idCount();
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return EDR_FALSE;
    }
    return EDR_TRUE;
}



/** Returns the component string, identified by index, which
    makes up the supplied EDROutputName */
DLL_EXPORT const char* EdrOutputNameIdGet(const EDROutputName outputName,  /* (I) */
                                          int                 index){      /* (O) */
    const static char routine[] = "EdrOutputNameIdGet";
    const char*  component = 0;
    if (!ok2run(routine)) {
        return 0;
    }
    try{
        if (!outputName){
            ModelException(routine, "NULL outputName");
        }
        ObjectWrapper* nameWrapper = (ObjectWrapper*)outputName;
        OutputNameSP nameSP(OutputNameSP::dynamicCast(nameWrapper->object));
        const string& strComponent = nameSP->idGet(index);
        component = strComponent.c_str();
    } catch (exception& e){
        EASErrorHandler(e, routine);
        return 0;
    }
    return component;
}
        
/************* Memory Support *********************/
    
/** Frees the memory contained within an EDRObject. */
DLL_EXPORT void EdrObjectFree(EDRObject object){
    ObjectWrapper* objWrapper = (ObjectWrapper*)object;
    delete objWrapper; // delete allows NULL's 
}

/** Frees the memory of a C string pointer (char*) returned by EDR.
    this is provided to avoid the need of explict call to free() by client which may 
    give rise problems due to lib versions. */
DLL_EXPORT void EdrCStrFree(char* str){
    if (str != 0)
        free (str);
}

/*** Classes to test interface (complements EDRTester.cpp) ***/   
class EdrMainTester: public CObject {
public:
    static CClassConstSP const TYPE;       
    static IObjectSP run(EdrMainTester* params);
    
private:
    ScenarioSP    scenario;
    IModelSP      model;
    CInstrumentSP inst;
    CControlSP    ctrl;
    CMarketDataSP mkt;
    
    EdrMainTester();
    friend class EdrMainTesterHelper;
};

EdrMainTester::EdrMainTester() : CObject(TYPE) {}

IObjectSP EdrMainTester::run(EdrMainTester* params) {
    static const string method = "EdrMainTester::run";
    try {
        bool ok = false;

        CStringArraySP output(new CStringArray(0));

        EDRMarketDataCache cache1 = 0;
        EDRMarketDataCache cache2 = 0;
        
        ObjectWrapper* edrScenario = 0;
        ObjectWrapper* edrModel    = new ObjectWrapper(params->model);
        ObjectWrapper* edrInst     = new ObjectWrapper(params->inst);
        ObjectWrapper* edrCtrl     = new ObjectWrapper(params->ctrl);

        if (params->scenario.get()) {
            edrScenario = new ObjectWrapper(params->scenario.get());
        }
        
        EDRObject     edrResults = 0;
        double        multiplier = 1.0;
        double        weight     = 1.0;
        EDROutputName name = 0;
        EDRObject     object = 0;
        EDRObject     today = 0;

        EDRResultsIterator iter = 0;

        DateTime valueDate =  params->mkt->GetReferenceDate();

        Library::shutdown();
        if (!EdrStartup()) {
            goto done;
        }

        if (!(today = (EDRObject)new ObjectWrapper(valueDate.clone()))) {
            goto done;
        }

        cache1 = (EDRMarketDataCache)(new ObjectWrapper(params->mkt));

        // test building a cache, adding to it, then copying it 
        if (!(cache2 = EdrMarketDataCacheNew(today))) {
            goto done;
        }

        // add something to it

        // copy input cache
        EdrObjectFree(cache2);

        if (!(cache2 = EdrMarketDataCacheClone(cache1))) {
            goto done;
        }


        edrResults = EdrMain(1,
                             (EDRObject *)&edrModel,
                             (EDRObject *)&edrInst,
                             (EDRObject *)&edrCtrl,
                             &multiplier,
                             &weight,
                             (EDRObject)edrScenario,
                             cache2);

        if (!edrResults) {
            goto done;
        }

        if (!(iter = EdrResultsIteratorGet(edrResults))) {
            goto done;
        }

        const char*    packet;
        EDRBool        iterStatus;

        while ((iterStatus = EdrResultsIteratorNext(iter, &packet,
                                                    &name,
                                                    &object)) &&
               packet) {
            int numStrings;
            if (!(EdrOutputNameIDCount(name, &numStrings))) {
                goto done;
            }

            string label = packet;
            label += "_";
            for (int j = 0; j < numStrings; j++) {
                const char* nm = EdrOutputNameIdGet(name, j);
                label += nm;
            }
            output->push_back(label);
        }


        ok = true;
    done:
        EdrObjectFree(edrResults);
        EdrObjectFree((EDRObject)edrScenario);
        EdrObjectFree((EDRObject)edrModel);
        EdrObjectFree((EDRObject)edrInst);
        EdrObjectFree((EDRObject)edrCtrl);
        EdrObjectFree(cache1);
        EdrObjectFree(cache2);
        EdrObjectFree(iter);
        EdrObjectFree(today);

        if (!EdrErrorCleared()) {
            char* error = EdrErrorStackTrace();
            ModelException e(method, error);
            free(error);
            throw e;
        }

        EdrShutdown();
        Library::startup();
        return output;
    }
    catch (exception& e) {
        EdrShutdown();
        Library::startup();
        throw ModelException(e, method);
    }
}

class  EdrMainTesterHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(EdrMainTester, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultEdrMainTester);

        FIELD(scenario, "scenario");
        FIELD_MAKE_OPTIONAL(scenario);
        FIELD(model, "model");
        FIELD(inst, "instrument");
        FIELD(ctrl, "control");
        FIELD(mkt, "market");

        Addin::registerClassObjectMethod("IFACE_EDR_MAIN_TESTER",
                                         Addin::UTILITIES,
                                         "tests EDR main methods",
                                         EdrMainTester::TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)EdrMainTester::run);


    }

    static IObject* defaultEdrMainTester(){
        return new EdrMainTester();
    }
};

CClassConstSP const EdrMainTester::TYPE = 
CClass::registerClassLoadMethod(
    "EdrMainTester", typeid(EdrMainTester), EdrMainTesterHelper::load);
  

class EdrDataDictTester: public CObject {
public:
    static CClassConstSP const TYPE;       
    static IObjectSP run(EdrDataDictTester* params);
    
private:
    IObjectSP object;

    EdrDataDictTester();
    friend class EdrDataDictTesterHelper;
};

EdrDataDictTester::EdrDataDictTester() : CObject(TYPE) {}

IObjectSP EdrDataDictTester::run(EdrDataDictTester* params) {
    static const string method = "EdrDataDictTester::run";
    try {
        bool ok = false;

        IObjectSP      output(0);
        ObjectWrapper* edrobject = new ObjectWrapper(params->object);
        EDRDataDict    dd1 = 0;
        EDRDataDict    dd2 = 0;
        const char*    typekey = 0;
        const char*    field;
        EDRObject      cmpt = 0;
        EDRObject      out = 0;
        EDRBool        iterStatus;

        EDRDataDictIterator iter = 0;

        Library::shutdown();
        if (!EdrStartup()) {
            goto done;
        }

        if (!(dd1 = EdrObjectToDataDict((EDRObject)edrobject))) {
            goto done;
        }

        if (!(typekey = EdrObjectGetType((EDRObject)edrobject))) {
            goto done;
        }

        if (!(dd2 = EdrDataDictNew(typekey))) {
            goto done;
        }

        if (!(iter = EdrDataDictIteratorGet(dd1))) {
            goto done;
        }

        while ((iterStatus = EdrDataDictIteratorNext(iter, &field, &cmpt)) &&
               field){
            cmpt = EdrDataDictGet(dd1, field); // test get()
            
            if (!(EdrDataDictAdd(dd2, field, cmpt))) {
                goto done;
            }
            EdrObjectFree(cmpt);
        }

        if (!(out = EdrDataDictToObject(dd2))) {
            goto done;
        }

        output = IObjectSP(((ObjectWrapper *)out)->object);

        ok = true;
    done:
        EdrObjectFree((EDRObject)edrobject);
        EdrObjectFree(dd1);
        EdrObjectFree(dd2);
        EdrObjectFree(iter);
        EdrObjectFree(out);

        if (!EdrErrorCleared()) {
            char* error = EdrErrorStackTrace();
            ModelException e(method, error);
            free(error);
            throw e;
        }

        EdrShutdown();
        Library::startup();
        return output;
    }
    catch (exception& e) {
        EdrShutdown();
        Library::startup();
        throw ModelException(e, method);
    }
}

class  EdrDataDictTesterHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(EdrDataDictTester, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultEdrDataDictTester);

        FIELD(object, "object");

        Addin::registerClassObjectMethod("IFACE_EDR_DATADICT_TESTER",
                                         Addin::UTILITIES,
                                         "tests EDR datadict methods",
                                         EdrDataDictTester::TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)EdrDataDictTester::run);


    }

    static IObject* defaultEdrDataDictTester(){
        return new EdrDataDictTester();
    }
};

CClassConstSP const EdrDataDictTester::TYPE = 
CClass::registerClassLoadMethod(
    "EdrDataDictTester", 
    typeid(EdrDataDictTester), 
    EdrDataDictTesterHelper::load);
  
class EdrActionTester: public CObject {
public:
    static CClassConstSP const TYPE;       
    static IObjectSP run(EdrActionTester* params);
    
private:
    IObjectSP object;

    EdrActionTester();
    friend class EdrActionTesterHelper;
};

EdrActionTester::EdrActionTester() : CObject(TYPE) {}

IObjectSP EdrActionTester::run(EdrActionTester* params) {
    static const string method = "EdrActionTester::run";
    try {
        bool ok = false;
        IObjectSP& object = params->object;
        // check that action object is public otherwise EAS can't create it
        if (!Modifier::isPublic(object->getClass()->getModifiers())){
            throw ModelException(method, "Action ("+
                                 object->getClass()->getName()+") is not "
                                 "public");
        }

        IObjectSP      output(0);
        EDRObject      out = 0;
        ObjectWrapper* input = new ObjectWrapper(object);

        Library::shutdown();
        if (!EdrStartup()) {
            goto done;
        }

        if (!(out = EdrAction((EDRObject)input))) {
            goto done;
        }

        output = IObjectSP(((ObjectWrapper *)out)->object);

        ok = true;
    done:
        EdrObjectFree((EDRObject)input);
        EdrObjectFree(out);

        if (!EdrErrorCleared()) {
            char* error = EdrErrorStackTrace();
            ModelException e(method, error);
            free(error);
            throw e;
        }

        EdrShutdown();
        Library::startup();
        return output;
    }
    catch (exception& e) {
        EdrShutdown();
        Library::startup();
        throw ModelException(e, method);
    }
}

class  EdrActionTesterHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(EdrActionTester, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultEdrActionTester);

        FIELD(object, "object");

        Addin::registerClassObjectMethod("IFACE_EDR_ACTION_TESTER",
                                         Addin::UTILITIES,
                                         "tests EDR Action methods",
                                         EdrActionTester::TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)EdrActionTester::run);


    }

    static IObject* defaultEdrActionTester(){
        return new EdrActionTester();
    }
};

CClassConstSP const EdrActionTester::TYPE = 
CClass::registerClassLoadMethod(
    "EdrActionTester", 
    typeid(EdrActionTester), 
    EdrActionTesterHelper::load);
  

