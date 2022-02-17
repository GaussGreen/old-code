#include "edginc/config.hpp"
#define QLIB_RESULTS_CPP
#if defined(_MSC_VER)
// disable warning truncated decorated names information
#pragma warning(disable : 4503)
#endif

#include "edginc/Results.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Writer.hpp"
#include "edginc/DerivativeAsset.hpp"
#include "edginc/CommandLineParams.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/Untweakable.hpp"
#include ext_hash_map
DRLIB_BEGIN_NAMESPACE

// for DR Interface purposes
class ResultsProxy: public CObject {
public:
    static CClassConstSP const TYPE;

private:
    StringArray     packet;
    OutputNameArray identifier;
    ObjectArray     value;

    typedef smartPtr<ResultsProxy> ResultsProxySP;

    static IObjectSP toResults(const IObjectConstSP& obj) {
        const ResultsProxy* proxy = dynamic_cast<const ResultsProxy*>(obj.get());
        if (!proxy) {
            throw ModelException("ResultsProxy::toResults",
                                 "object is not a ResultsProxy");
        }
    
        if (proxy->packet.size() != proxy->identifier.size() ||
            proxy->packet.size() != proxy->value.size()) {
            throw ModelException("ResultsProxy::toResults", 
                                 "internal array size mismatch");
        }

        ResultsSP results(new Results);
        for (int i = 0; i < proxy->packet.size(); i++) {
            results->storeGreek(proxy->value[i], 
                                proxy->packet[i], 
                                proxy->identifier[i]);
        }
        return results;
    }

    static IObjectSP fromResults(const IObjectConstSP& obj) {
        const Results* results = dynamic_cast<const Results*>(obj.get());
        if (!results) {
            throw ModelException("ResultsProxy::fromResults",
                                 "object is not a Results");
        }
    
        ResultsProxySP proxy(new ResultsProxy());

        vector<const string*> packets(results->listPyramidPackets());
        for (int i = 0; i < (int)packets.size(); i++) {
            OutputNameArraySP names(results->packetContents(*packets[i]));
            for (int j = 0; j < names->size(); j++) {
                IObjectConstSP value(results->retrieveGreek(*packets[i],
                                                            (*names)[j]));

                proxy->packet.push_back(*packets[i]);
                proxy->identifier.push_back((*names)[j]);
                proxy->value.push_back(IObjectSP::constCast(value));
            }
        }

        return proxy;
    }

    ResultsProxy() : CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); 
        REGISTER(ResultsProxy, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultResultsProxy);
        registerObjectProxy(Results::TYPE,
                            ResultsProxy::TYPE,
                            fromResults,
                            toResults);
        FIELD(packet, "packet");
        FIELD(identifier, "identifier");
        FIELD(value, "value");
        FIELD_MAKE_OPTIONAL(packet);
        FIELD_MAKE_OPTIONAL(identifier);
        FIELD_MAKE_OPTIONAL(value);
    }

    static IObject* defaultResultsProxy(){
        return new ResultsProxy();
    }  
};

CClassConstSP const ResultsProxy::TYPE = CClass::registerClassLoadMethod(
    "ResultsProxy", typeid(ResultsProxy), load);

/** Currently only used for holding the price - use anywhere else would need
    to be looked at carefully to ensure scaling/addition of results still
    worked */
class ResultWithCcy: public CObject,
                     virtual public CombinableResult {
public:
    static CClassConstSP const TYPE;
    double     result;
    string     resultCcy;
    
    ResultWithCcy(double result, const string& resultCcy):
        CObject(TYPE), result(result), resultCcy(resultCcy){}
    
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const{
        CDoubleSP value(CDouble::create(result));
        value->outputWrite(linePrefix, prefix+"("+resultCcy+")", stream);
    }
    /** scale by factor x */
    virtual void scale(double x){
        result *= x;
    }

    /** add an object (scaled by scaleFactor) to this result.  */
    virtual void add(const CombinableResult& x, double scaleFactor){

          // however: gcc-4.0 has problem with it at runtime
          const ResultWithCcy& resultWithCcy = 
	  // gcc bug: force to IObject before dynamic cast
	  dynamic_cast<const ResultWithCcy&>(static_cast<const IObject&>(x));
          result += scaleFactor * resultWithCcy.result;           
    }

private:
    ResultWithCcy():CObject(TYPE){}

    static IObject* defaultResultWithCcy(){
        return new ResultWithCcy();
    }
    static void load(CClassSP& clazz){
        clazz->setPublic(); 
        REGISTER(ResultWithCcy, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(CombinableResult);
        EMPTY_SHELL_METHOD(defaultResultWithCcy);
        FIELD(result, "The result's value");
        FIELD(resultCcy, "Currency of the result");
    }
};
CClassConstSP const ResultWithCcy::TYPE = CClass::registerClassLoadMethod(
    "ResultWithCcy", typeid(ResultWithCcy), load);

#if 0
// more general version of the above - comes in 2 flavours - combinable and
// non-combinable depending on the abilities of the object it contains
class ResultObjectWithCcy: public CObject {
public:
    static CClassConstSP const TYPE;
    IObjectSP  result;
    string     resultCcy;
    
    ResultObjectWithCcy(IObjectSP result, const string& resultCcy):
        CObject(TYPE), result(result), resultCcy(resultCcy){}

    ResultObjectWithCcy(CClassConstSP clazz, 
                        IObjectSP     result, 
                        const string& resultCcy):
        CObject(clazz), result(result), resultCcy(resultCcy){}
    
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const{
        result->outputWrite(linePrefix, prefix+"("+resultCcy+")", stream);
    }

protected:
    ResultObjectWithCcy(CClassConstSP clazz): CObject(clazz){}

private:
    ResultObjectWithCcy():CObject(TYPE){}

    static IObject* defaultResultObjectWithCcy(){
        return new ResultObjectWithCcy();
    }
    static void load(CClassSP& clazz){
        REGISTER(ResultObjectWithCcy, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultResultObjectWithCcy);
        FIELD(result, "The result");
        FIELD(resultCcy, "Currency of the result");
    }
};
CClassConstSP const ResultObjectWithCcy::TYPE = CClass::registerClassLoadMethod(
    "ResultObjectWithCcy", typeid(ResultObjectWithCcy), load);

class CombiningResultObjectWithCcy: public ResultObjectWithCcy,
                                    virtual public CombinableResult {
public:
    static CClassConstSP const TYPE;
    
    CombiningResultObjectWithCcy(IObjectSP result, const string& resultCcy):
        ResultObjectWithCcy(TYPE, result, resultCcy){}
    
    /** scale by factor x */
    virtual void scale(double x){
        CombinableResult& combiner = dynamic_cast<CombinableResult&>(*result);
        combiner.scale(x);
    }

    /** add an object (scaled by scaleFactor) to this result.  */
    virtual void add(const CombinableResult& x, double scaleFactor){
        CombinableResult& combiner = dynamic_cast<CombinableResult&>(*result);
        combiner.add(x, scaleFactor);
    }

private:
    CombiningResultObjectWithCcy():ResultObjectWithCcy(TYPE){}

    static IObject* defaultCombiningResultObjectWithCcy(){
        return new CombiningResultObjectWithCcy();
    }
    static void load(CClassSP& clazz){
        REGISTER(CombiningResultObjectWithCcy, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(CombinableResult);
        EMPTY_SHELL_METHOD(defaultCombiningResultObjectWithCcy);
        FIELD(result, "The result");
        FIELD(resultCcy, "Currency of the result");
    }
};
CClassConstSP const CombiningResultObjectWithCcy::TYPE = 
CClass::registerClassLoadMethod(
    "CombiningResultObjectWithCcy", typeid(CombiningResultObjectWithCcy), load);

#endif

/* Ideally we'd define PacketData as:
typedef hash_map<OutputNameConstSP, IObjectSP, OutputName::HashUtil,
                 OutputName::HashUtil> PacketData;
but this is causing problems on solaris when we put PacketData into 
PacketHash because the symbols become too long. So fix is to wrap template
into explicit class */
class PacketData: public hash_map<OutputNameConstSP, IObjectSP, 
                  OutputName::HashUtil, OutputName::HashUtil> {
public:
    typedef hash_map<OutputNameConstSP, IObjectSP, 
                     OutputName::HashUtil,
                     OutputName::HashUtil>::const_iterator const_iterator;
    typedef hash_map<OutputNameConstSP, IObjectSP, 
        OutputName::HashUtil, OutputName::HashUtil>::iterator iterator;
};
typedef refCountPtr<PacketData> PacketDataSP;

struct MyStringHash {
    size_t operator()(const string& str) const {
        return (hash_string(str.c_str()));
    }
};

typedef hash_map<string, PacketDataSP, MyStringHash> PacketHash;

/** Note this class is used to hold the hash table and hence needs to be
    registered with the type information */
class ResultsHelper: public CObject{
    ResultWithCcy*  fairValue; // not registered $unregistered
    /*CDoubleSP   fairValue;*/   /* handy reference to the value (aka price)
                                - note: just a cache, true value in packets */
public:
    static CClassConstSP const TYPE;
    static OutputNameSP  nameForFairValue; // initialised on load
    PacketHash  packets; // $unregistered
    
    /** Invoked when ResultsHelper class is 'loaded' */
    static void loadResultsHelperClass(CClassSP& clazz){
        REGISTER(ResultsHelper, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultResultsHelper);
        nameForFairValue = OutputNameSP(new OutputName(Results::VALUE));
        Results::emptyName        = OutputNameSP(new OutputName(""));
    }

    static IObject* defaultResultsHelper(){
        return new ResultsHelper();
    }
    /** Invoked when Results class is 'loaded' */
    static void loadResultsClass(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDRIProxyType(ResultsProxy::TYPE); // use proxy for dri
        REGISTER(Results, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultResults);
        FIELD(data, "Hashtable containing results");
    }

    static IObject* defaultResults(){
        return new Results();
    }
    ResultsHelper():CObject(TYPE), fairValue(0){}

    ~ResultsHelper(){}

    //// updates the cached fair value field
    void updateCachedFairValue(){
        try{
            PacketHash::const_iterator iter1 = 
                packets.find(Results::INSTRUMENT_PACKET);
            if (!(iter1 == packets.end())){
                PacketData::const_iterator iter2 = 
                    iter1->second->find(nameForFairValue);
                if (iter2 != iter1->second->end()){
                    ResultWithCcy& resultWithCcy = 
                        dynamic_cast<ResultWithCcy&>(*iter2->second);
                    setFairValue(&resultWithCcy);
                }
            }
        } catch (exception& e){
            throw ModelException(e, "ResultsHelper::updateCachedFairValue");
        }
    }

    // must override this method as default will fail to copy hash table */
    virtual IObject* clone() const{
        ResultsHelper* copy = new ResultsHelper();
        try{
            for (PacketHash::const_iterator iter1 = packets.begin();
                 !(iter1 == packets.end()); ++iter1){
                const PacketDataSP& packetToCopy = iter1->second;
                PacketDataSP packetData(new PacketData());
                copy->packets[iter1->first] = packetData;
                for (PacketData::const_iterator iter2 = packetToCopy->begin();
                     !(iter2 == packetToCopy->end()); ++iter2){
                    OutputNameSP nameCopy(iter2->first.clone());
                    IObjectSP    objCopy(iter2->second.clone());
                    (*packetData)[nameCopy] = objCopy;
                }
            }
            // also need to 'copy' cached reference to price (if there)
            copy->updateCachedFairValue();
        } catch (exception& e){
            delete copy;
            throw ModelException(e, "ResultsHelper::clone");
        }
        return copy;
    }
            
    // retrieves packet from results - do not free. If quiet then may return
    // null if packet doesn't exist
    PacketData* retrievePacket(const string& packet, bool quiet){
        PacketHash::const_iterator iter = packets.find(packet);
        if (iter == packets.end()){
            if (quiet){
                return 0;
            } else {
                throw ModelException("Results::retrievePacket",
                                     "No results under "+packet);
            }
        }
        return iter->second.get();
    }

    //// scales elements of packet by scaleFactor
    static void scalePacket(PacketData*       packet,       // modify this by
                            double            scaleFactor){ // scaling by this
        // scale combinable results in packet 
        for (PacketData::iterator iter = packet->begin();
             !(iter == packet->end()); ++iter) {
            scaleSingleResult(iter->second, scaleFactor);
        }
    }

    static void scaleSingleResult(const IObjectSP& result,
                                  double           scaleFactor){
        if (result->getClass()->isArray()){
            IArray& arrayToScale = dynamic_cast<IArray&>(*result);
            for (int i = 0; i < arrayToScale.getLength(); i++){
                scaleSingleResult(arrayToScale.get(i), scaleFactor);
            }
        } else if (CombinableResult::TYPE->isInstance(result)) {
            CombinableResult& cr = 
                dynamic_cast<CombinableResult&>(*result);
            cr.scale(scaleFactor);
        }
    }

    //// adds elements of packetToAdd (scaled by scaleFactor) to packet
    static void addPacket(const string&     packetName,   // eg DELTA
                          PacketData*       packet,       // modify this by
                          const PacketData* packetToAdd,  // adding this
                          double            scaleFactor){ // scaled by this
        // add combinable results in packet 
        // need to take special care in case either of the packets has the
        // magic "not applicable" in them
        bool done = false;
        if (packetToAdd->size() == 1){
            PacketData::const_iterator iter = packetToAdd->begin();
            if (NotApplicable::TYPE->isInstance(iter->second)){
                done = !packet->empty(); /* done if there is already a
                                `           result in packet */
            }
        }
        if (!done && packet->size() == 1){
            PacketData::iterator iter = packet->begin();
            if (NotApplicable::TYPE->isInstance(iter->second)){
                packet->erase(iter); // remove it
            }
        } 
        if (!done){
            for (PacketData::const_iterator iter = packetToAdd->begin();
                 !(iter == packetToAdd->end()); ++iter) {
                addSingleResult(packetName, iter->first, packet, 
                                iter->second, scaleFactor);
            }
        }
    }

    //// merge elements of packetToAdd to packet
    static void mergePacket(const string&     packetName,   // eg SPOT_PRICE
                            PacketData*       packet,       // modify this by
                            const PacketData* packetToAdd){ // merging this
        // merge results in packet 
        // need to take special care in case either of the packets has the
        // magic "not applicable" in them
        bool done = false;
        if (packetToAdd->size() == 1){
            PacketData::const_iterator iter = packetToAdd->begin();
            if (NotApplicable::TYPE->isInstance(iter->second)){
                done = !packet->empty(); /* done if there is already a
                                           result in packet */
            }
        } else if (packet->size() == 1){
            PacketData::iterator iter = packet->begin();
            if (NotApplicable::TYPE->isInstance(iter->second)){
                packet->erase(iter); // remove it
            }
        } 
        if (!done){
            for (PacketData::const_iterator iter = packetToAdd->begin();
                 !(iter == packetToAdd->end()); ++iter) {
                const OutputNameConstSP& name = iter->first;
                // does this exist in packet?
                if (packet->find(name) == packet->end()){
                    // doesn't exist in packet so add it
                    try{
                        (*packet)[name] = IObjectSP(iter->second->clone());
                    } catch (exception& e){
                        // record what went wrong
                        (*packet)[name] = IObjectSP(new Untweakable(e));
                    }
                }
            }
        }
    }

    static void addArrayResult(
        const string&            packetName,
        const OutputNameConstSP& outputName,
        PacketData*              packet,       // modify this by
        const IArray&            arrayToAdd,   // adding this
        double                   scaleFactor){ // scaled by this
        static const string method("Results::addArrayResult");
        CClassConstSP  arrayClass = arrayToAdd.getClass();
        if (CombinableResult::TYPE->isAssignableFrom(
            arrayClass->getComponentType())){
            // see if it exists in packet
            PacketData::iterator iter2 = packet->find(outputName);
            if (iter2 == packet->end()){
                // just clone and scale
                IObjectSP arrayCopy(arrayToAdd.clone());
                scaleSingleResult(arrayCopy, scaleFactor);
                (*packet)[outputName] = arrayCopy;
            } else {
                IObjectSP& object = iter2->second;
                CClassConstSP objType = object->getClass();
                if (arrayClass == objType){
                    // if types match exactly then add
                    IArray& arrayDest = dynamic_cast<IArray&>(*object);
                    // merge arrays - assume order matches up
                    if (arrayDest.getLength() != arrayToAdd.getLength()){
                        string m("Could not combine "+ packetName+" for "+
                                 outputName->toString()+" as arrays are "
                                 "of different length");
                        throw ModelException(method, m);
                    }
                    for (int i = 0; i < arrayToAdd.getLength(); i++){
                        IObjectConstSP objToAdd(arrayToAdd.get(i));
                        IObjectSP      objDest(arrayDest.get(i));
                        
                        CombinableResult& crDest = 
                        dynamic_cast<CombinableResult&>(*objDest);

                        const CombinableResult& crToAdd = 
                            dynamic_cast<const CombinableResult&>(*objToAdd);
                        crDest.add(crToAdd, scaleFactor);

                        //Some arrays (e.g. CDoubleArray) contain the 
                        //raw data in the vector rather than smart pointers.
                        //Hence get() returned a copy of the data.
                        arrayDest.set(i, objDest);
                    }
                } else if (objType->isArray() ||
                           !CombinableMixedResult::TYPE->isInstance(object)){
                    // giving up
                    throw ModelException(method, "Can not combine "+
                                         arrayToAdd.getClass()->getName()+
                                         " with "+ objType->getName());
                } else {
                    // finally use CombinableMixedResult
                    CombinableMixedResult& cmr = 
                        dynamic_cast<CombinableMixedResult&>(*object);
                    (*packet)[outputName] = IObjectSP(
                        cmr.addResult(arrayToAdd, scaleFactor));
                }
            }
        }
               
    }   
    static void addSingleResult(
        const string&            packetName,
        const OutputNameConstSP& outputName,
        PacketData*              packet,       // modify this by
        const IObjectConstSP&    resultToAdd,  // adding this
        double                   scaleFactor){ // scaled by this
                                
        static const string routine("Results::addSingleResult");
        if (resultToAdd->getClass()->isArray()){
            const IArray& arrayToAdd = 
                dynamic_cast<const IArray&>(*resultToAdd);
            addArrayResult(packetName, outputName, packet,
                           arrayToAdd, scaleFactor);
        } else if (CombinableResult::TYPE->isInstance(resultToAdd)) {
            const CombinableResult& crToAdd = 
                dynamic_cast<const CombinableResult&>(*resultToAdd);
            // see if it exists in packet
            PacketData::iterator iter2 = packet->find(outputName);
            try{
                IObjectSP obj = iter2 == packet->end()? 
                    IObjectSP(): iter2->second;
                if (!obj){
                    CombinableResult* crCopy = copy(&crToAdd);
                    crCopy->scale(scaleFactor);
                    (*packet)[outputName] = IObjectSP(crCopy);
                } 
                // now see how to add the results together - either same
                // type or need one of them to be a CombinableMixedResult
                else if (crToAdd.getClass() == obj->getClass()){
                    CombinableResult& cr =
                        dynamic_cast<CombinableResult&>(*obj);
                    cr.add(crToAdd, scaleFactor);
                } else if (CombinableMixedResult::TYPE->isInstance(obj)){
                    CombinableMixedResult& cmr = 
                        dynamic_cast<CombinableMixedResult&>(*obj);
                    IObjectSP newResult(cmr.addResult(crToAdd, scaleFactor));
                    (*packet)[outputName] = newResult;
                }else if (CombinableMixedResult::TYPE->isInstance(&crToAdd)){
                    // tricky - first copy it
                    IObjectSP cmrObj(crToAdd.clone());
                    CombinableMixedResult& cmr = 
                        dynamic_cast<CombinableMixedResult&>(*cmrObj);
                    // then scale it
                    cmr.scale(scaleFactor);
                    // then add original to it (scaling by 1.0)
                    IObjectSP newResult(cmr.addResult(*obj, 1.0));
                    // finally store
                    (*packet)[outputName] = newResult;
                } else {
                    string m("Couldn't add results for "+packetName+
                             " for "+ outputName->toString()+ 
                             " because type "+
                             iter2->second->getClass()->getName()+" is"
                             " not combinable");
                    throw ModelException(routine, m);
                }
            } catch (exception& e){
                // record what went wrong
                (*packet)[outputName] = IObjectSP(new Untweakable(e));
            }
        }
    }
    /** Is the fair value stored. This is useful for mtm results (when there
        are no theo results) */
    bool fairValueExists() const {
        return fairValue != 0;
    }

    /** Just returns a reference to the double SP which holds a reference
        to the price, but checks that it is not null */
    ResultWithCcy* getFairValue(){
        if (!fairValue){
            throw ModelException("Results::getFairValue", "No price stored");
        }
        return fairValue;
    }
       
    void setFairValue(ResultWithCcy* fairValue){
        this->fairValue = fairValue;
    }
};
 
OutputNameSP  ResultsHelper::nameForFairValue; // initialised in load
CClassConstSP const ResultsHelper::TYPE = CClass::registerClassLoadMethod(
    "ResultsHelper", typeid(ResultsHelper),
    ResultsHelper::loadResultsHelperClass);

CClassConstSP const Results::TYPE = CClass::registerClassLoadMethod(
    "Results", typeid(Results), ResultsHelper::loadResultsClass);

DEFINE_TEMPLATE_TYPE_WITH_NAME("ResultsArray", CResultsArray);

const string Results::INSTRUMENT_PACKET      = "Instrument";
const string Results::FWD_AT_MAT_PACKET      = "FWD_AT_MAT";
const string Results::DEBUG_PACKET           = "Debug";
const string Results::DEBUG_PACKETS_PREFIX   = "DEBUG_";
const string Results::VALUE                  = "VALUE";
const string Results::ESW_LEG_PRICE          = "ESW_LEG_PRICE";
const string Results::KNOWN_CASHFLOWS_PACKET = "KNOWN_CASHFLOWS";
const string Results::CLEAN_DEFAULT_SPREAD_CURVE_PACKET = "CLEAN_DEFAULT_SPREAD_CURVE";
const string Results::SHIFT_SIZE_POSTFIX     = "_SHIFT_SIZE";
const string Results::SPI_PACKET             = "SPI";
const string Results::TRANCHE_CONTINGENT_LEG_PACKET = "ContingentLeg";
const string Results::TRANCHE_FEE_LEG_PACKET = "FeeLeg";

OutputNameSP Results::emptyName;

/* Results are stored in 'packets'. A packet is identified by a string (
   eg delta, rho etc). All packets are stored in a hash table. Each packet
   is itself a hash table. Here the key for the hash table is of type
   OutputNameConstSP. The values are of type IObjectSP */

Results::Results(): CObject(TYPE){
    data = new ResultsHelper();
}

Results::~Results(){
    delete data;
}

void Results::storeResult(const string&            packet,
                          const OutputNameConstSP& resultName,
                          const IObjectSP&         result){
    PacketHash::const_iterator iter = data->packets.find(packet);
    PacketDataSP packetData;
    if (iter == data->packets.end()){
        packetData = PacketDataSP(new PacketData());
        data->packets[packet] = packetData;
    } else {
        packetData = iter->second;
    }
    (*packetData)[resultName] = result;
}

/** stores double together with ccy ISO code */
void Results::storeRequestResult(OutputRequest*           request,
                                 double                   result,
                                 const string&            ccyISOCode,
                                 const OutputNameConstSP& outputName){
    IObjectSP wrappedResult(new ResultWithCcy(result, ccyISOCode));
    storeRequestResult(request, wrappedResult, outputName);
}

#if 0
/** stores object together with ccy ISO code */
void Results::storeRequestResult(OutputRequest*           request,
                                 const IObjectSP&         result,
                                 const string&            ccyISOCode,
                                 const OutputNameConstSP& outputName) {
    IObjectSP wrappedResult;
    if (dynamic_cast<const CombinableResult*>(result.get())) {
        wrappedResult=IObjectSP(new CombiningResultObjectWithCcy(result, 
                                                                 ccyISOCode));
    }
    else {
        wrappedResult=IObjectSP(new ResultObjectWithCcy(result, ccyISOCode));
    }        
    storeRequestResult(request, wrappedResult, outputName);
}
#endif

void Results::storeRequestResult(OutputRequest*           request,
                                 const IObjectSP&         resultObj,
                                 const OutputNameConstSP& outputName)
{
    string packetName = request->getPacketName();
    try {
        storeResult(packetName, outputName, resultObj);
        // flag successful calculation of requested output
        request->setHasFinished(true);
    }
    catch (exception& e) {
        throw ModelException(e, "Results::storeRequestResult");
    }
}

// lazy version of the above
void Results::storeRequestResult(OutputRequest* request,
                                 double         result,
                                 const string&  outputName) {
    IObjectSP    obj(CDouble::create(result));
    OutputNameSP name(new OutputName(outputName));

    storeRequestResult(request, obj, name);
}

/** for requests that don't need further qualification e.g. IND_VOL */
void Results::storeRequestResult(
    OutputRequest*    request,
    const IObjectSP&  resultObj)
{
    try {
        const string&  packetName = request->getPacketName();
        OutputNameSP outputName(new OutputName(request->getRequestName()));
        storeResult(packetName, outputName, resultObj);
        // flag successful calculation of requested output
        request->setHasFinished(true);
    }
    catch (exception& e) {
        throw ModelException(e, "Results::storeRequestResult");
    }
}

/** for requests that don't need further qualification e.g. IND_VOL */
void Results::storeRequestResult(
    OutputRequest* request,
    double         result)
{
    try {
        IObjectSP resultObj(CDouble::create(result));
        storeRequestResult(request, resultObj);
    }
    catch (exception& e) {
        throw ModelException(e, "Results::storeRequestResult");
    }
}

/** Retrieves a generic output request result. Note returns reference to result */
IObjectConstSP Results::retrieveRequestResult(
              const string& outputRequestName) const{
    try {
        OutputRequestSP outReq(new OutputRequest(outputRequestName));
        OutputNameSP outputName(new OutputName(outReq->getRequestName()));
        return retrieveGreek(outReq->getPacketName(),
                             outputName);
    }
    catch (exception& e) {
        throw ModelException(e, "Results::retrieveRequestResult");
    }

}
//// helper method. Do we have results for when DerivativeAsset's use mtm
bool Results::mtmResultsExist() const{
    OutputNameConstSP mtmName(
        new OutputName(DerivativeAsset::ASSET_MTM_ID));
    return (exists(Results::INSTRUMENT_PACKET, mtmName));
}

//// helper method. Get results for when DerivativeAsset's use mtm
Results* Results::getMtmResults(){
    OutputNameConstSP mtmName(
        new OutputName(DerivativeAsset::ASSET_MTM_ID));
    IObjectSP mtmResultsObj =
        retrieveResult(Results::INSTRUMENT_PACKET, mtmName);
    Results& mtmResults = dynamic_cast<Results&>(*mtmResultsObj);
    return &mtmResults;
}
//// helper method. Get results for when DerivativeAsset's use mtm
const Results* Results::getMtmResults() const{
    Results* nonConstRes = const_cast<Results*>(this);
    return nonConstRes->getMtmResults();
}

/** scale all the Additive results (as determined by the control)
    inside a Results object by supplied factor. The singleInstStatistics
    indicates whether we're trying to caclulate the mean value of
    a result by pricing the same instrument multiple times */
void Results::scale(
    const  CControlSP&   control,
    double               scaleFactor,
    bool                 singleInstStatistics){
    static const string method = "Results::scale";
    try {
        // existence of fair values implies existence of theo results
        if (data->fairValueExists()){
            // first scale the price - this works as fairValue is a reference
            data->getFairValue()->scale(scaleFactor);
        
            // next loop over the requested sensitivities
            SensitivityArrayConstSP sens = control->getSens();
            for (int i = 0; i < sens->size(); i++) {
                // ask each sens to scale them
                (*sens)[i]->scaleResult(this, scaleFactor);
            }
            // then loop over OutputRequests
            OutputRequestArrayConstSP requests = control->getOutputRequests();
            for (int j = 0; j < requests->size(); j++) {
                // ask each OutputRequests to scale them
                (*requests)[j]->scaleResult(this, scaleFactor, 
                                            singleInstStatistics);
            }
            
        }
        // then recurse if we have additional set of results
        if (mtmResultsExist()){
            Results* mtmResults = getMtmResults();
            mtmResults->scale(control, scaleFactor, singleInstStatistics);
        }
    } catch (exception& e) {
        throw ModelException(e, method, "Failed to scale results");
    }
}

/** scale all the CombinableResults inside a Results by a factor */
void Results::scalePostProcess(const  CControlSP&   control,
                               double               scaleFactor) {
    static const string method = "Results::scalePostProcess";
    try {
        // first scale the price - this works as fairValue is a reference
        if (data->fairValueExists()){
            data->getFairValue()->scale(scaleFactor);
        }
        // next loop over the requested sensitivities
        SensitivityArrayConstSP sens = control->getSens();
        for (int i = 0; i < sens->size(); i++) {
            // ask each sens to scale them
            (*sens)[i]->scaleResult(this, scaleFactor);
        }
        // then loop over OutputRequests
        OutputRequestArrayConstSP requests = control->getOutputRequests();
        for (int j = 0; j < requests->size(); j++) {
            // ask each OutputRequests to scale them
            (*requests)[j]->scalePostProcess(this, scaleFactor);
        }

    } catch (exception& e) {
        throw ModelException(e, method, "Failed to scale results");
    }
}


/** scale all the CombinableResult results in the given packet
    by supplied factor */
void Results::scale(const string&       packetName,
                    double              scaleFactor){
    PacketData* packet = data->retrievePacket(packetName, false);
    ResultsHelper::scalePacket(packet, scaleFactor);
}


/** scale the result with given sensName in the given packet
    by the supplied factor */
void Results::scale(const string&       packetName,
                    const string&       sensName,
                    double              scaleFactor){
    OutputNameSP outputName(new OutputName(sensName));
    IObjectSP    result(retrieveResult(packetName, outputName));
    ResultsHelper::scaleSingleResult(result, scaleFactor);
}
    

/** Modify 'this' Results set by adding all results in resultsToAdd
    as indicated by control. Assumes both sets of results exists. Use
    for adding theo or mtm */
void Results::addSingleResultsSet(
    const Results*     resultsToAdd,
    const CControlSP&  control,  // used to obtain resultsToAdd
    double             scaleFactor,
    bool               sameInstrument){ /* are the 2 results for
                                           the same instrument */
    try{
        // add prices 
        if (!data->fairValueExists() && resultsToAdd->data->fairValueExists()){
            storePrice(0.0, resultsToAdd->getCcyName());
        }
        if (resultsToAdd->data->fairValueExists()){
            data->getFairValue()->add(*resultsToAdd->data->getFairValue(),
                                      scaleFactor);
        }
        
        // next loop over the requested sensitivities
        SensitivityArrayConstSP sens = control->getSens();
        for (int i = 0; i < sens->size(); i++) {
            // ask each sens to add them
            (*sens)[i]->addResult(this, resultsToAdd, scaleFactor);
        }
        // then loop over OutputRequests
        OutputRequestArrayConstSP requests =
            control->getOutputRequests();
        for (int j = 0; j < requests->size(); j++) {
            // ask each OutputRequests to scale them
            (*requests)[j]->addResult(this, resultsToAdd, 
                                      scaleFactor, sameInstrument);
        }
    } catch (exception& e){
        throw ModelException(e, "addSingleResultsSet");
    }
}

/** Modify 'this' Results set by adding all results in resultsToAdd
    as indicated by control */
void Results::add(
    const Results*     resultsToAdd,
    const CControlSP&  control,  // used to obtain resultsToAdd
    double             scaleFactor,
    bool               sameInstrument){ /* are the 2 results for
                                           the same instrument */
   static const string method = "Results::add";
    try {
        bool mtmResultsToAddExist = resultsToAdd->mtmResultsExist();
        bool mtmResultsDestExist = mtmResultsExist();
        bool theoResultsToAddExist = resultsToAdd->data->fairValueExists();
        bool theoResultsDestExist = data->fairValueExists();
        // create duplicate copies of anything we need first before we do
        // any adding. Avoids double counting.
        if (!mtmResultsDestExist && mtmResultsToAddExist){
            OutputNameConstSP mtmName(
                new OutputName(DerivativeAsset::ASSET_MTM_ID));
            IObjectSP clonedRes(clone());
            storeGreek(clonedRes, Results::INSTRUMENT_PACKET, mtmName);
        }
        if (!theoResultsDestExist && theoResultsToAddExist && 
            mtmResultsDestExist){
            const Results* source = getMtmResults();
            // set empty fair value
            storePrice(0.0, source->getCcyName());
            addSingleResultsSet(source, control, 
                                scaleFactor, sameInstrument);
        }
        // if separate results don't exist just add theo results */
        if (mtmResultsToAddExist || mtmResultsDestExist){
            const Results* mtmResultsToAdd = mtmResultsToAddExist? 
                resultsToAdd->getMtmResults(): resultsToAdd;
            Results* mtmResults = getMtmResults();
            mtmResults->addSingleResultsSet(mtmResultsToAdd, control, 
                                            scaleFactor, sameInstrument);
        }
        if (theoResultsToAddExist || theoResultsDestExist){
            const Results* theoSource = theoResultsToAddExist?
                resultsToAdd: resultsToAdd->getMtmResults();
            addSingleResultsSet(theoSource, control, 
                                scaleFactor, sameInstrument);
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Modify 'this' Results set by adding all CombinableResult
    results in the packet of given name to this Results */
void Results::add(const string&      packetName,
                  const Results*     resultsToAdd,
                  double             scaleFactor){
    try{
        // look up the two packets and fire at the ResultsHelper
        PacketData* packet = data->retrievePacket(packetName, true);
        if (!packet){
            // results don't exist in running total, so create empty packet
            packet = new PacketData();
            data->packets[packetName] = PacketDataSP(packet);
        }
        const PacketData* packetToAdd =
            resultsToAdd->data->retrievePacket(packetName, false);
        ResultsHelper::addPacket(packetName, packet, packetToAdd, scaleFactor);
    } catch (exception& e){
        throw ModelException(e, "Results::add",
                             "Couldn't combine results for "+ packetName);
    }
}

/** Modify 'this' Results set by merging all objects in the packet of
    given name with those in this Results object */
void Results::merge(const string&      packetName,
                    const Results*     resultsToAdd){
    try{
        // look up the two packets and fire at the ResultsHelper
        PacketData* packet = data->retrievePacket(packetName, true);
        if (!packet){
            // results don't exist in running total, so create empty packet
            packet = new PacketData();
            data->packets[packetName] = PacketDataSP(packet);
        }
        const PacketData* packetToAdd =
            resultsToAdd->data->retrievePacket(packetName, false);
        if (packetToAdd) { 
            ResultsHelper::mergePacket(packetName, packet, packetToAdd);
        }
        
    } catch (exception& e){
        throw ModelException(e, "Results::merge",
                             "Couldn't merge results for "+ packetName);
    }
}

/** Merge all packets from resultsToAdd into this results object */
void Results::merge(const Results* resultsToAdd) { 
    const PacketHash& ph = resultsToAdd->data->packets;
    for (PacketHash::const_iterator it = ph.begin();it != ph.end();it++) { 
        const string& packetName = it->first;
        merge(packetName, resultsToAdd);
    }
}

/** Modify this Results by overwriting all results contained in
    resultsToCopyFrom */
void Results::overwrite(const Results* resultsToCopyFrom){
    for (PacketHash::const_iterator iter = 
             resultsToCopyFrom->data->packets.begin();
         !(iter == resultsToCopyFrom->data->packets.end()); ++iter){
        const string&  packetName = iter->first;
        const PacketDataSP& packet = iter->second;
        PacketDataSP& packetToOverwrite = data->packets[packetName];
        if (!packetToOverwrite){
            packetToOverwrite = PacketDataSP(new PacketData());
        }
        for (PacketData::const_iterator iter = packet->begin();
             !(iter == packet->end()); ++iter){
            (*packetToOverwrite)[iter->first] = iter->second;
        }
    }
    data->updateCachedFairValue();
}
        
//// Combines results from pricing an array of instruments
ResultsSP Results::combineResults(
    const CControlArray&   ctrls,
    const DoubleArray&     weights,
    CResultsArray&         results){
    try{
        // 1. create initially empty Results
        ResultsSP total(new Results());
        // 2. Then fill up with initial fair value
        //total->storePrice(0.0, results[0]->getCcyName());
        // 3. Next, loop over results and add them in
        for (int i = 0; i < results.size(); i++){
            total->add(results[i].get(), ctrls[i], weights[i], false);
        }
        return total;
    }
    catch (exception& e) {
        throw ModelException(e, "Results::combineResults");
    }
}  

        

/** Modify 'this' Results set by adding the result with given
    sensName in the given packet by the supplied factor */
void Results::add(const string&      packetName,
                  const string&      sensName,
                  const Results*     resultsSource,
                  double             scaleFactor){
    // loop up the result to add
    OutputNameSP outputName(new OutputName(sensName));
    IObjectSP    resultToAdd(resultsSource->retrieveResult(packetName,
                                                           outputName));
    // look up which packet to add it to
    PacketData* packet = data->retrievePacket(packetName, false);
    // and do the biz
    ResultsHelper::addSingleResult(packetName, outputName, packet,
                                   resultToAdd, scaleFactor);
}

/** Removes the specified packet if it exists */
void Results::removePacket(const string&       packetName){
    PacketHash::iterator iter1 = data->packets.find(packetName);
    if (iter1 != data->packets.end()){
        // packet exists - so get rid of it
        data->packets.erase(iter1);
    }
}
 
/** Store not applicable for the given packet name and outputName,
    use null for outputName if the whole packet is not applicable */
void Results::storeNotApplicable(const string&       packetName,
                                 const OutputNameSP& outputName,
                                 const IObjectSP&    na){
    try {
        if (!outputName){
            PacketHash::iterator iter1 = data->packets.find(packetName);
            if (iter1 != data->packets.end()){
                // packet exists - get rid of it
                data->packets.erase(iter1);
            }
        }
        storeResult(packetName, !outputName? emptyName: outputName, na);
    }
    catch (exception& e) {
        throw ModelException(e, "Results::storeNotApplicable");
    }
}

/** Store not applicable for the given packet name */
void Results::storeNotApplicable(const string&  packetName){
    storeNotApplicable(packetName, OutputNameSP(), // whole packet
                       IObjectSP(new NotApplicable()));
}

/** If a greek doesn't make sense for a product (e.g. delta for cashflows)
    store a NotApplicable object as the result */
void Results::storeNotApplicable(const Sensitivity* sens) {
    NotApplicable na; // generate default cause
    storeNotApplicable(sens, na.getCause());
}

//// As above but can specify cause
void Results::storeNotApplicable(const Sensitivity*     sens,
                                 NotApplicable::Cause   cause){
    const string& sensName = sens->getSensOutputName();
    const string& packetName = sens->getPacketName();
    OutputNameSP outputName;
    if (packetName != sensName){
        outputName = OutputNameSP(new OutputName(sensName));
    }
    storeNotApplicable(packetName, outputName, 
                       IObjectSP(new NotApplicable(cause)));
}
    
void Results::storeNotApplicable(const OutputRequest*    request, 
                                 const IObjectSP&        na){
    try {
        const string& packetName = request->getPacketName();
        const string& requestName = request->getRequestName();
        OutputNameSP name;
        // that's a bit heavy handed - if the request has it's own
        // packet, we (ie. it does not belong into the instrument
        // packet, we only want to store the packet name. Other-
        // wise, we need the packet name and the request name
        if (packetName != requestName) {
            name = OutputNameSP(new OutputName(requestName));
        }
        storeNotApplicable(packetName, name, na);
    }
    catch (exception& e) {
        throw ModelException(e, "Results::storeNotApplicable");
    }
}

/** If an output request doesn't make sense for a product 
    (e.g. effective strike for a vanilla) store a NotApplicable 
    object as the result */
void Results::storeNotApplicable(const OutputRequest*    request) {
    IObjectSP na(new NotApplicable());
    storeNotApplicable(request, na);
}

//// As above but can specify cause
void Results::storeNotApplicable(const OutputRequest*    request,
                                 NotApplicable::Cause    cause){
    IObjectSP na(new NotApplicable(cause));
    storeNotApplicable(request, na);
}

/** Returns true if an a result is stored for the specified request
    and it is of type NotApplicable */
const NotApplicable* Results::isNotApplicable(
    const OutputRequest* request) const{
    const string& requestName = request->getRequestName();
    const string& packetName = request->getPacketName();
    OutputNameSP outputName;
    if (requestName == packetName){
        outputName = emptyName;
    } else {
        outputName = OutputNameSP (new OutputName(requestName));
    }
    IObjectConstSP obj;
    if (exists(packetName, outputName)){
        obj = retrieveResult(packetName, outputName);
        if (NotApplicable::TYPE->isInstance(obj)){
            return (&dynamic_cast<const NotApplicable&>(*obj));
        }
    }
    return 0;
}
 
// has a greek been stored as NotApplicable?
bool Results::isNotApplicable(SensControl* sens) const {
    const string& sensName = sens->getSensOutputName();
    const string& packetName = sens->getPacketName();
    OutputNameSP outputName;
    if (packetName != sensName){
        outputName = OutputNameSP(new OutputName(sensName));
    }  
    else {
        outputName =  emptyName;
    }
    if (exists(packetName, outputName)){
        IObjectConstSP obj = retrieveResult(packetName, outputName);    

        return (NotApplicable::TYPE->isInstance(obj));
    }    
    return false;
}
    
IObjectSP Results::retrieveResult(const string&              packet,
                                  const OutputNameConstSP&   resultName) const{
    PacketData* thePacket =  data->retrievePacket(packet, false);
    PacketData::const_iterator iter2 = thePacket->find(resultName);
    if (iter2 == thePacket->end()){
        throw ModelException("Results::retrieveResult",
                             "No result found for "+resultName->toString()+
                             " under "+packet);
    }
    return iter2->second;
}

/** Returns the name of currency in which the value is in */
const string& Results::getCcyName() const{
    if (mtmResultsExist()){
        const Results* mtm = getMtmResults();
        return mtm->getCcyName();
    }
    return data->getFairValue()->resultCcy;
}


/** Stores the instrument value */
void Results::storePrice(double         value,
                         const string&  ccyName){
    smartPtr<ResultWithCcy> wrappedFV(new ResultWithCcy(value, ccyName));
    // store in correct place
    storeResult(INSTRUMENT_PACKET, ResultsHelper::nameForFairValue, wrappedFV);
    // save our reference  in local cache
    data->setFairValue(wrappedFV.get());
}

/** Stores a generic greek */ 
void Results::storeGreek(IObjectSP                result,
                         const string&            packetName,
                         const OutputNameConstSP& outputName){
    storeResult(packetName, outputName, result);
}

/** Removes a generic greek if it exists. The packet is also removed
    if it becomes empty */ 
void Results::removeGreek(const string&            packetName,
                          const OutputNameConstSP& outputName){
    PacketHash::iterator iter1 = data->packets.find(packetName);
    if (iter1 != data->packets.end()){
        // packet exists. See if result exists
        PacketData::iterator iter2 = iter1->second->find(outputName);
        if (iter2 != iter1->second->end()){
            iter1->second->erase(iter2);
            if (iter1->second->empty()){
                data->packets.erase(iter1);
            }
        }
    }
}

void Results::storeCcylessResult(IObjectSP                result,
                                 const string&            packetName,
                                 const OutputNameConstSP& outputName){
    storeResult(packetName, outputName, result);
}

/** Stores a generic greek. Equivalent to  
    storeGreek(result,
    sens->getSensOutputName(), 
    sens->getMarketDataName()) */
void Results::storeGreek(IObjectSP      result,
                         SensControl*   sens){  // (I) identifies greek
    storeResult(sens->getSensOutputName(), sens->getMarketDataName(), result);
}
    

/** Stores a 'scalar' greek */
void Results::storeScalarGreek(double         result,
                               SensControl*   sens){  // (I) identifies greek
    CDoubleSP doubleVal(CDouble::create(result));
    storeResult(sens->getSensOutputName(), sens->getMarketDataName(), 
                doubleVal);
}

/** Stores a 'scalar' greek */ 
void Results::storeScalarGreek(double                   result,
                               const string&            packetName,
                               const OutputNameConstSP& outputName){
    CDoubleSP doubleVal(CDouble::create(result));
    storeResult(packetName, outputName, doubleVal);
}


/** Retrieves the instrument price */
double Results::retrievePrice() const{
    return data->getFairValue()->result;
}

/** Whether a price has been set using storePrice() */
bool Results::priceExists() const{
    return data->fairValueExists();
}

/** checks whether results hold a valid results for a given greek */
bool Results::isValidScalarGreek(
    const string&            packetName,
    const OutputNameConstSP& outputName) const // (I) identifies greek
{
    IObjectSP result = retrieveResult(
                            packetName,
                            outputName);

    return (CDouble::TYPE->isInstance(result.get()));
}

/** Retrieves all scalar results (incl price, calc_time etc.) with packet names, output names.
    Returns num of results retrieved. */
int Results::retrieveAllScalarGreeks(StringArray&         packetNames,
                                     OutputNameArray&     outNames,
                                     DoubleArray&         results) const
{
    packetNames.clear();
    outNames.clear();
    results.clear();
    for (PacketHash::const_iterator iter1 = data->packets.begin();
         !(iter1 == data->packets.end()); ++iter1)
    {
        for (PacketData::const_iterator iter2 = iter1->second->begin();
             !(iter2 == iter1->second->end()); ++iter2){
            packetNames.push_back(iter1->first);
            OutputNameSP hack(OutputNameSP::constCast(iter2->first));
            outNames.push_back(OutputNameSP(copy(hack.get())));
            try{
                if (iter1->first == INSTRUMENT_PACKET && 
                    iter2->first->equals(VALUE))
                    results.push_back(data->getFairValue()->result);
                else
                {
                    CDoubleSP res(CDoubleSP::dynamicCast(iter2->second));
                    results.push_back(res->doubleValue());
                }
            }
            catch(ModelException&){
                // continue if the result is not a scalar
                packetNames.erase(packetNames.end()-1);
                outNames.erase(outNames.end()-1);
            }
        }
     }
     return outNames.size();
}

/** Retrieves a 'scalar' greek */
double Results::retrieveScalarGreek(
    CSensControl*  sens) const // (I) identifies greek
{
    CDoubleSP result(CDoubleSP::dynamicCast(
        retrieveResult(sens->getSensOutputName(), sens->getMarketDataName())));
    return result->doubleValue();
}

/** Retrieves a 'scalar' greek */
double Results::retrieveScalarGreek(
    const string&            packetName,
    const OutputNameConstSP& outputName) const
{
    CDoubleSP result(CDoubleSP::dynamicCast(
        retrieveResult(packetName, outputName)));
    return result->doubleValue();
}

/** Retrieves a generic greek. Note returns reference to result */
IObjectConstSP Results::retrieveGreek(
    const string&            packetName,
    const OutputNameConstSP& outputName) const{
    return retrieveResult(packetName, outputName);
}

/** Returns true if given result exists. Equivalent to 
    exists(sens->getSensOutputName(), sens->getMarketDataName()) */
bool Results::exists(SensControl*   sens) const{
    return exists(sens->getSensOutputName(), sens->getMarketDataName());
}


/** Returns true if given result exists.  */
bool Results::exists(const string&            packetName,
                     const OutputNameConstSP& outputName) const{
    PacketHash::const_iterator iter1 = data->packets.find(packetName);
    if (iter1 == data->packets.end()){
        return false;
    }
    PacketData::const_iterator iter2 = iter1->second->find(outputName);
    if ( iter2 == iter1->second->end() ) {
        // check if NotApplicable object exists
        iter2 = iter1->second->find(emptyName);
        if ( iter2 != iter1->second->end() &&
             NotApplicable::TYPE->isInstance(iter2->second)) {
            return true;
        } else {
            return false;
        }
    }
    else {
        // output exists
        return true;
    }
}
/** Returns true if a result exists for the specified output
    request */
bool Results::exists(const OutputRequest* outputRequest) const{
    if (outputRequest->getRequestName() == outputRequest->getPacketName()){
        return packetExists(outputRequest->getPacketName());
    }
    OutputNameSP outputName(new OutputName(outputRequest->getRequestName()));
    return exists(outputRequest->getPacketName(), outputName);
}

/** Returns true if given packet exists.  */
bool Results::packetExists(const string&  packetName) const{
    PacketHash::const_iterator iter1 = data->packets.find(packetName);
    return (iter1 != data->packets.end());
}

/** Returns the list of packets within the Results */
StringArraySP Results::listAllPackets() const {
    StringArraySP packetNames(new StringArray());
    for (PacketHash::const_iterator iter = data->packets.begin();
         !(iter ==  data->packets.end()); ++iter){
        packetNames->push_back(iter->first);
    }
    return packetNames;
}



/** Returns a list of output names within a packet */
OutputNameArraySP Results::packetContents(const string& packetName) const {
    static const string method = "Results::packetContents";
    OutputNameArraySP names(new OutputNameArray(0));
    try {
        PacketHash::const_iterator iter1 = data->packets.find(packetName);
        if (iter1 == data->packets.end()){
            throw ModelException(method, "No results under "+packetName);
        }

        for (PacketData::const_iterator iter2 = iter1->second->begin();
             !(iter2 == iter1->second->end()); ++iter2){
            OutputNameSP hack(OutputNameSP::constCast(iter2->first));
            //names->push_back(OutputNameSP(copy(hack.get())));
            names->push_back(hack);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
    return names;
}


/** Returns the list of packets within the Results for use by
    Pyramid (hides certain packets eg DEBUG) */
vector<const string*> Results::listPyramidPackets() const{
    vector<const string*> packetNames;
    for (PacketHash::const_iterator iter = data->packets.begin();
         !(iter ==  data->packets.end()); ++iter){
        // hide DEBUG_PACKETs and any shift size packets
        if (iter->first != DEBUG_PACKET && 
            !CString::equalsIgnoreCase(iter->first, DEBUG_PACKETS_PREFIX,
                                       DEBUG_PACKETS_PREFIX.size())){
            unsigned int len = SHIFT_SIZE_POSTFIX.size();
            if (iter->first.size() < len ||
                iter->first.substr(len) != SHIFT_SIZE_POSTFIX){
                packetNames.push_back(&iter->first);
            }
        }
    }
    return packetNames;
}
           
/** Returns an array of all the results in the given packet. */
vector<pair<OutputNameConstSP, IObjectConstSP> > Results::listPacketResults(
    const string&  packetName) const{
    PacketHash::const_iterator iter1 = data->packets.find(packetName);
    if (iter1 == data->packets.end()){
        throw ModelException("Results::listPacketResults",
                             "No results under "+packetName);
    }
    vector<pair<OutputNameConstSP, IObjectConstSP> > results;
    for (PacketData::const_iterator iter2 = iter1->second->begin();
         !(iter2 == iter1->second->end()); ++iter2){
        IObjectConstSP objectToReturn = iter2->second;
        results.push_back(pair<OutputNameConstSP, 
                          IObjectConstSP>(iter2->first, objectToReturn));
    }
    return results;
}

/** Helper class for xml read/write for hash table of 
    <OutputNameConstSP, IObjectSP> */
class ResultsEntry: public CObject{
public:
    static CClassConstSP const TYPE;
    OutputNameConstSP  key;
    IObjectSP          value;
    static const string TAG;
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ResultsEntry, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultResultsEntry);
        FIELD(key, "key");
        FIELD(value, "value");
    }

    static IObject* defaultResultsEntry(){
        return new ResultsEntry();
    }
    // for reflection
    ResultsEntry(): CObject(TYPE){}

    ResultsEntry(const OutputNameConstSP& key, 
                 const IObjectSP&         value):
        CObject(TYPE), key(key), value(value){}
};
const string ResultsEntry::TAG = "ResultsEntry";

CClassConstSP const ResultsEntry::TYPE =
CClass::registerClassLoadMethod("ResultsEntry", typeid(ResultsEntry), load);

/** write object out to writer */
void Results::write(const string& tag, Writer* writer) const {
    try {
        IObjectConstSP obj(writer->objectStart(tag, "", this, true));
        if (obj.get()){
            // if we have implemented IPrivateObject then we need to act
            // on obj here which will be of type IPublicObject
            for (PacketHash::const_iterator iter = data->packets.begin();
                 !(iter == data->packets.end()); ++iter){
                /* if PacketData was typed we'd just do
                   iter->second->write(iter->first, xml);
                */
                const string&  packetName = iter->first;
                const PacketDataSP& packet = iter->second;
                writer->objectStart(packetName, "", 0, true);
                for (PacketData::const_iterator iter = packet->begin();
                     !(iter == packet->end()); ++iter){
                    ResultsEntry entry(iter->first, iter->second);
                    // write out entry
                    entry.write(ResultsEntry::TAG, writer);
                }
                writer->objectEnd(packetName, this);
            }
        }
        writer->objectEnd(tag, this);
    }
    catch (exception& e){
        throw ModelException(e, "Results::write");
    }
}

/** populate an empty object from XML description */
void Results::import(Reader::Node* elem, Reader* reader){
    const static string routine = "Results::import";
    try{
        Reader::NodeListSP packets(elem->children());
        for (unsigned int i = 0; i < packets->size(); i++){
            Reader::Node* packetNode = (*packets)[i].get();
            string packetName = packetNode->name();
            // now read contents of packet
            Reader::NodeListSP nodes(packetNode->children());
            for (unsigned int j = 0; j < nodes->size(); j++){
                Reader::Node* entry = (*nodes)[j].get();
                // read in entry
                IObjectSP entryObj(reader->read(entry));
                ResultsEntry& resultsEntry = 
                    dynamic_cast<ResultsEntry&>(*entryObj);
                // at last can add entry
                storeResult(packetName, resultsEntry.key, resultsEntry.value);
            }
        }
        // finally update cached fair price value
        data->updateCachedFairValue();
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}
 
/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void Results::outputWrite(const string& linePrefix,
                          const string& prefix, ostream& stream) const{
    static IObjectSP debugToReport(
        CommandLineParams::getParameterArgs(CommandLineParams::DebugPackets));
    string debugName;
    if (debugToReport.get()){
        debugName = dynamic_cast<CString&>(*debugToReport).stringValue();
    }
    for (PacketHash::const_iterator iter = data->packets.begin();
         !(iter == data->packets.end()); ++iter)
    {
        if ((!debugName.empty() && 
             CString::equalsIgnoreCase(iter->first, debugName,
                                       debugName.size())) ||
            (iter->first != DEBUG_PACKET &&
             !CString::equalsIgnoreCase(iter->first, DEBUG_PACKETS_PREFIX,
                                        DEBUG_PACKETS_PREFIX.size()))){
            const string&  packetNamePrefix = prefix == ""?
                iter->first: prefix + "_" + iter->first;
            const PacketDataSP& packet = iter->second;

            for (PacketData::const_iterator iter = packet->begin();
                 !(iter == packet->end()); ++iter)
            {
                const string&  resultName = iter->first->toString();
                iter->second->outputWrite(linePrefix, packetNamePrefix + 
                                          (!resultName.empty() ? "_" : "") + 
                                          resultName, stream);
            }
        }
    }
}

/** Addin for building handles to arrays of specific types of objects */
class ResultsAddin: public CObject,
                    public virtual ClientRunnable{
    static CClassConstSP const TYPE;

    /** addin takes two parameters - the name of the result to get
        and the result object */
    ResultsSP   results; // holds the results
    string      id1;
    StringArray id2;
    StringArray id3;

    IObjectSP run(){
        return getResult(this);
    }

    /** the 'addin function' - builds array of correct type */
    static IObjectSP getResult(ResultsAddin* params){
        static const string routine = "ResultsAddin::getResult";
        try{
            int size = 3;
            if (params->id3.empty()){
                size--;
                if (params->id2.empty()){
                    size--;
                    if (params->id1.empty()){
                        size--;
                    }
                }
            }
            if (size == 0){
                // list all packets
                return params->results->listAllPackets();
            }
            if (size == 1) {
                if (params->id1 == Results::VALUE){
                    double value = params->results->retrievePrice();
                    return IObjectSP(CDouble::create(value));
                }
                // list all the identifiers in the packet
                return params->results->packetContents(params->id1);
            } 

            // can get multiple results in one go
            if (!params->id3.empty()) {
                if (params->id3.size() != params->id2.size()) {
                    throw ModelException(routine, "lists of names are "
                                         "of different lengths");
                }
            }
            if (params->id2.size() == 1) {
                // just 1 result
                OutputNameSP name = OutputNameSP(size == 3?
                                                 new OutputName(params->id2[0],
                                                                params->id3[0]):
                                                 new OutputName(params->id2[0]));
                IObjectConstSP result = 
                    params->results->retrieveGreek(params->id1, name);
                // return a copy of the results to the addin
                return IObjectSP(result.clone());
            }
            else {
                ObjectArraySP result(new ObjectArray(params->id2.size()));
                for (int i = 0; i < params->id2.size(); i++) {
                    OutputNameSP name = OutputNameSP(size == 3?
                                                     new OutputName(params->id2[i],
                                                                    params->id3[i]):
                                                     new OutputName(params->id2[i]));                    

                    IObjectConstSP greek = 
                        params->results->retrieveGreek(params->id1, name);
                    (*result)[i] = IObjectSP(greek.clone());
                }
                
                return result;
            }
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    ResultsAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(ResultsAddin, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultResultsAddin);
        FIELD(results, "The results from pricing");
        FIELD(id1, "Packet name or VALUE");
        FIELD_MAKE_OPTIONAL(id1);
        FIELD(id2, "First identifier");
        FIELD_MAKE_OPTIONAL(id2);
        FIELD(id3, "Second identifier");
        FIELD_MAKE_OPTIONAL(id3);
        Addin::registerInstanceObjectMethod("GET_RESULT",
                                            Addin::RISK,
                                            "Returns a specific result and "
                                            "unpacks it",
                                            TYPE,
                                            false,
                                            Addin::expandMulti,
                                            (Addin::ObjMethod*)getResult);
        Addin::registerInstanceObjectMethod("GET_RESULT_OBJECT",
                                            Addin::RISK,
                                            "Returns a specific result as a"
                                            " handle",
                                            TYPE,
                                            true,
                                            Addin::returnHandle,
                                            (Addin::ObjMethod*)getResult);
    }

    static IObject* defaultResultsAddin(){
        return new ResultsAddin();
    }
    
};

CClassConstSP const ResultsAddin::TYPE = CClass::registerClassLoadMethod(
    "ResultsAddin", typeid(ResultsAddin), load);



DRLIB_END_NAMESPACE
