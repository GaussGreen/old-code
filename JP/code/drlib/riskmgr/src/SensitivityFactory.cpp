//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SensitivityFactory.cpp
//
//   Description : Class for managing instantiation of Sensitivities
//
//   Author      : Andrew J Swain
//
//   Date        : 2 November 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/Control.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/Addin.hpp"
#include "edginc/MegaShuffleInterface.hpp"

#include ext_hash_map
#include <algorithm>
DRLIB_BEGIN_NAMESPACE

SensitivityFactory::ICreation::ICreation(){}
SensitivityFactory::ICreation::~ICreation(){}
SensitivityFactory::IDefault::IDefault(){}
SensitivityFactory::IDefault::~IDefault(){}
SensitivityFactory::IScalar::IScalar(){}
SensitivityFactory::IScalar::~IScalar(){}

class Lookup {
public:
    Lookup() : factory(), example(0), shiftInterface(0) {}

    Lookup(SensitivityFactory::ICreation* factory,
           Sensitivity*                   example,
           CClassConstSP                  shiftInterface) :
        factory(factory), example(example), shiftInterface(shiftInterface) {
        // empty
    }

    refCountPtr<SensitivityFactory::ICreation> factory;
    SensitivitySP                              example;
    CClassConstSP                              shiftInterface;
};

typedef hash_map<string, Lookup, Hashtable::StringHash> FactoryHashtable;

/** Support ability to build sens controls based upon their name */
static FactoryHashtable factories;

/** Register a sensitivity with the factory. Method takes ownership of
    ICreation* */
void SensitivityFactory::addSens(
    const string& name,            // name for this sensitivity eg DELTA
    ICreation*    factory,         // capable of building Sensitivity
    Sensitivity*  example,         // any old instance of the sensitivity
    CClassConstSP shiftInterface)  /* the associate shift interface for
                                      this tweak eg Delta::Shift */
{
    const static string method("SensitivityFactory::addSens");
    if (!factory){
        throw ModelException(method, "NULL factory");
    }
    if (!example){
        throw ModelException(method, "NULL example");
    }
    FactoryHashtable::const_iterator iter = factories.find(name);
    if (iter != factories.end()){
        throw ModelException(method, name+" already registered!");
    }

    Lookup lookup(factory, example, shiftInterface);
    factories[name] = lookup;
}

Sensitivity* SensitivityFactory::defaultSensitivity(const string& name){
    static const string method("SensitivityFactory::defaultSensitivity");
    FactoryHashtable::const_iterator iter = factories.find(name);
    if (iter == factories.end()){
        throw ModelException(method, "Unknown sensitivity: "+name);
    }

    Lookup     lookup = iter->second;
    ICreation* factory = lookup.factory.get();
    IDefault*  defaultFactory = dynamic_cast<IDefault*>(factory);
    if (!defaultFactory){
        throw ModelException(method, "Default construction for "+name +
                             "not supported");
    }
    return defaultFactory->createDefault();
}

Sensitivity* SensitivityFactory::scalarSensitivity(const string& name,
                                                   double        shiftSize){
    static const string method("SensitivityFactory::scalarSensitivity");
    FactoryHashtable::const_iterator iter = factories.find(name);
    if (iter == factories.end()){
        throw ModelException(method, "Unknown sensitivity: "+name);
    }

    Lookup     lookup = iter->second;
    ICreation* factory = lookup.factory.get();
    IScalar*   scalarFactory = dynamic_cast<IScalar*>(factory);
    if (!scalarFactory){
        throw ModelException("Scalar shift not supported for "+name);
    }
    return scalarFactory->createScalar(shiftSize);
}

// addin functions
class SensListAddin: public CObject {
    static CClassConstSP const TYPE;

    /** the 'addin function' - list all sensitivities */
    static IObjectSP listSens(SensListAddin* params){
        static const string method = "SensListAddin::listSens";
        try {
            CStringArraySP sens(new CStringArray(0));
            FactoryHashtable::const_iterator iter = factories.begin();
            while (iter != factories.end()) {
                sens->push_back(iter->first);
                ++iter;
            }
            sort(sens->begin(), sens->end());
            return sens;
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** for reflection */
    SensListAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(SensListAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSensListAddin);

        Addin::registerClassObjectMethod("LIST_SENS",
                                         Addin::UTILITIES,
                                         "lists all sensitivities",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)listSens);
    }

    static IObject* defaultSensListAddin(){
        return new SensListAddin();
    }
};

CClassConstSP const SensListAddin::TYPE = CClass::registerClassLoadMethod(
    "SensListAddin", typeid(SensListAddin), load);

class SensNameAddin: public CObject {
    static CClassConstSP const TYPE;

    string name;

    // what output packet does this live in ?
    static string packetName(SensNameAddin* params){
        static const string method = "SensNameAddin::packetName";
        try {
            FactoryHashtable::const_iterator iter = factories.find(params->name);
            if (iter == factories.end()) {
                throw ModelException(method,
                                     "Unknown sensitivity: "+params->name);
            }

            Lookup lookup = iter->second;
            return (lookup.example->getPacketName());
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // what class represents this ?
    static string className(SensNameAddin* params){
        static const string method = "SensNameAddin::className";
        try {
            FactoryHashtable::const_iterator iter = factories.find(params->name);
            if (iter == factories.end()) {
                throw ModelException(method,
                                     "Unknown sensitivity: "+params->name);
            }

            Lookup lookup = iter->second;
            return (lookup.example->getClass()->getName());
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // what's the output name for this class ?
    static string outputName(SensNameAddin* params){
        static const string method = "SensNameAddin::outputName";
        try {
            string oname;
            bool   found = false;
            FactoryHashtable::const_iterator iter = factories.begin();
            while (iter != factories.end() && !found) {
                Lookup lookup = iter->second;
                if (params->name ==lookup.example->getClass()->getName()) {
                    found = true;
                    oname = iter->first;
                }
                ++iter;
            }

            if (!found) {
                throw ModelException(method,
                                     "sensitivity " + params->name +
                                     " doesn't exist");
            }
            return oname;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // is this sensitivity defaultable ?
    static bool hasDefault(SensNameAddin* params){
        static const string method = "SensNameAddin::hasDefault";
        try {
            FactoryHashtable::const_iterator iter = factories.find(params->name);
            if (iter == factories.end()){
                throw ModelException(method,
                                     "Unknown sensitivity: "+params->name);
            }

            Lookup     lookup = iter->second;
            SensitivityFactory::ICreation* factory = lookup.factory.get();
            SensitivityFactory::IDefault*  defaultFactory =
                dynamic_cast<SensitivityFactory::IDefault*>(factory);

            return (defaultFactory ? true : false);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // return default shift
    static IObjectSP sensDefault(SensNameAddin* params){
        static const string method = "SensNameAddin::sensDefault";
        try {
            return IObjectSP(SensitivityFactory::defaultSensitivity(params->name));
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }


    /** for reflection */
    SensNameAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(SensNameAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSensNameAddin);
        FIELD(name, "name");
        Addin::registerClassStringMethod("SENS_PACKET_NAME",
                                         Addin::UTILITIES,
                                         "gives packet name for a sensitivity",
                                         TYPE,
                                         (Addin::StringMethod*)packetName);

        Addin::registerClassStringMethod("SENS_CLASS_NAME",
                                         Addin::UTILITIES,
                                         "gives class name for a sensitivity",
                                         TYPE,
                                         (Addin::StringMethod*)className);

        Addin::registerClassStringMethod("SENS_OUTPUT_NAME",
                                         Addin::UTILITIES,
                                         "gives output name for a sensitivity",
                                         TYPE,
                                         (Addin::StringMethod*)outputName);

        Addin::registerClassBoolMethod("SENS_IS_DEFAULTABLE",
                                       Addin::UTILITIES,
                                       "is a sensitivity defaultable",
                                       TYPE,
                                       (Addin::BoolMethod*)hasDefault);

        Addin::registerClassObjectMethod("SENS_DEFAULT",
                                         Addin::UTILITIES,
                                         "builds default sensitivity shift",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)sensDefault);
    }

    static IObject* defaultSensNameAddin(){
        return new SensNameAddin();
    }
};

CClassConstSP const SensNameAddin::TYPE = CClass::registerClassLoadMethod(
    "SensNameAddin", typeid(SensNameAddin), load);

static Control* makeMegaControl(bool               writeToFile,
                                const string&      filename,
                                bool               defaultSens,
                                bool               defaultReq,
                                SensitivityArraySP extraSens,
                                const string&      assetPriceSource) {
    static const string method = "makeMegaControl";
    try {
        CStringArraySP       requests = OutputRequest::allRequests();
        OutputRequestArraySP req(new OutputRequestArray(0));
        for (int i = 0; i < requests->size() && defaultReq; i++) {
            req->push_back(OutputRequestSP(new OutputRequest((*requests)[i])));
        }

        SensitivityArraySP sens(new SensitivityArray(0));
        if (defaultSens) {
            FactoryHashtable::const_iterator iter = factories.begin();
            while (iter != factories.end()) {
                Lookup     lookup = iter->second;
                SensitivityFactory::ICreation* factory = lookup.factory.get();
                SensitivityFactory::IDefault*  defaultFactory =
                    dynamic_cast<SensitivityFactory::IDefault*>(factory);
                if (defaultFactory) {
                    sens->push_back(SensitivitySP(defaultFactory->createDefault()));
                }
                ++iter;
            }
        }

        CControlSP ctrl(new Control(sens,
                                    req,
                                    writeToFile,
                                    filename,
                                    assetPriceSource));

        if (!(!extraSens)) {
            for (int j = 0; j < extraSens->size(); j++) {
                ctrl->addSensitivity((*extraSens)[j]);
            }
        }
        return ctrl.release();
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

Control* SensitivityFactory::megaControl() {
    return makeMegaControl(false, "", true, true, SensitivityArraySP(   ), Control::ASSET_THEO);
}

Control* SensitivityFactory::megaControl(bool requests){
    return makeMegaControl(false, "", true, requests, SensitivityArraySP(   ), Control::ASSET_THEO);
}

// build a control with all output requests and as many sensitivities
// as can be defaulted
class MegaControlAddin: public CObject {
    static CClassConstSP const TYPE;

    bool               writeToFile;
    string             filename;
    bool               defaultSens;
    SensitivityArraySP extraSens;
    string             assetPriceSource;
    bool               shuffleSens;
    int                shuffleSeed;

    static IObjectSP megaCtrl(MegaControlAddin* params){
        static const string method = "MegaControlAddin::megaCtrl";
        try {

            int seed = params->shuffleSeed;
            if (seed == -999) { 
                seed = MegaShuffleInterface::createShuffleSeed();
            }

            Control *ctrl = makeMegaControl(params->writeToFile,
                params->filename,
                params->defaultSens,
                true,
                params->extraSens,
                params->assetPriceSource);

            if (params->shuffleSens) { 
                ctrl->shuffleSens(seed);
            }
            return IObjectSP(ctrl);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** for reflection */
    MegaControlAddin():  CObject(TYPE), assetPriceSource(Control::ASSET_THEO), shuffleSens(false), shuffleSeed(-999) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MegaControlAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultMegaControlAddin);
        FIELD(writeToFile, "writeToFile");
        FIELD(filename, "filename");
        FIELD(defaultSens, "defaultSens");
        FIELD(extraSens, "extraSens");
        FIELD_MAKE_OPTIONAL(extraSens);
        FIELD(assetPriceSource, "assetPriceSource");
        FIELD_MAKE_OPTIONAL(assetPriceSource);
        FIELD(shuffleSens, "shuffleSens");
        FIELD_MAKE_OPTIONAL(shuffleSens);
        FIELD(shuffleSeed, "shuffleSeed");
        FIELD_MAKE_OPTIONAL(shuffleSeed);

        Addin::registerClassObjectMethod("MEGA_CONTROL",
                                         Addin::RISK,
                                         "build a control with all output "
                                         "requests & sensitivities defaulted",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)megaCtrl);
    }

    static IObject* defaultMegaControlAddin(){
        return new MegaControlAddin();
    }
};

CClassConstSP const MegaControlAddin::TYPE = CClass::registerClassLoadMethod(
    "MegaControlAddin", typeid(MegaControlAddin), load);

DRLIB_END_NAMESPACE

