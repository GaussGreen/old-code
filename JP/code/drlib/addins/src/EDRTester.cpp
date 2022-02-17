//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EDRTester.hpp
//
//   Description : Models Library Interface Tester
//
//   Author      : Andrew J Swain
//
//   Date        : 15 June 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/EDRTester.hpp"
#include "edginc/Addin.hpp"
#include "edginc/EDRInterface.h"
#include "edginc/Library.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE



class EDRTesterHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(EDRTester, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultEDRTester);
    }

    static IObject* defaultEDRTester(){
        return new EDRTester();
    }
};

EDRTester::EDRTester(): CObject(TYPE) {}

CClassConstSP const EDRTester::TYPE = 
CClass::registerClassLoadMethod("EDRTester", 
                                typeid(EDRTester), 
                                EDRTesterHelper::load);
   
EDRTester::Atomic::Atomic(): CObject(TYPE) {}

IObjectSP EDRTester::Atomic::run(EDRTester::Atomic* params)  {
    static const string method = "EDRTester::Atomic::run";
    try {
        bool ok = false;

        EDRObject oi     = 0;
        EDRObject od     = 0;
        EDRObject os     = 0;
        EDRObject ob     = 0;
        EDRObject odm    = 0;
        EDRObject odt    = 0;
        EDRObject omatp  = 0;
        EDRObject omattp = 0; 
        EDRObject obm    = 0; 
        EDRObject oenum  = 0; 
        char*     exp    = 0;
        char*     str    = 0;
        int       cols;
        int       rows;

        EDRTester::Atomic* out = 0;

        Library::shutdown();
        if (!EdrStartup()) {
            goto done;
        }
    
        // build them up
        oi = EdrIntNew(params->i);
        od = EdrDoubleNew(params->d);
        os = EdrStringNew(params->s.c_str());
        ob = EdrBoolNew(params->b);
        if (params->enumType.get() && params->enumValue.get() &&
            !(oenum = EdrEnumNew(params->enumType->stringValue().c_str(),
                                 params->enumValue->stringValue().c_str()))){
            goto done;
        }
        if (!oi || !od || !os || !ob) {
            goto done;
        } 

        // do matrix 
        cols = params->dm.numCols();
        rows = params->dm.numRows();
        int i;
        int j;

        odm = EdrDoubleMatrixEmptyNew(cols, rows);
        if (!odm) {
            goto done;
        } 
        
        for (i = 0; i < cols; i++) {
            for (j = 0; j < rows; j++) {
                EdrDoubleMatrixSet(odm,
                                   i,
                                   j,
                                   params->dm[i][j]);

            }
        }

        odt    = EdrDateTimeNew(params->dt.c_str(), params->tm.c_str());
        omatp  = EdrMaturityPeriodNew(params->matp.c_str());
        omattp = EdrMaturityTimePeriodNew(params->matpdt.c_str(), params->matptm.c_str());
        obm    = EdrBenchmarkDateNew(params->bmdt.c_str(), params->bmtm.c_str());
        
        if (!odt || !omatp || !omattp || !obm) {
            goto done;
        } 

        // and pull them to bits
        out = new(EDRTester::Atomic);

        EdrIntGet(oi, &out->i);
        EdrDoubleGet(od, &out->d);
        str = EdrStringGet(os);
        out->s = string(str);
        free(str);
        EdrBoolGet(ob, &out->b);
        if (oenum){
            out->enumType.reset(CString::create(
                                    string(EdrObjectGetType(oenum))));
            out->enumValue.reset(CString::create(string(EdrEnumGet(oenum))));
        }
        EdrDoubleMatrixGetSize(odm, &cols, &rows);
        out->dm = DoubleMatrix(cols, rows);
        for (i = 0; i < cols; i++) {
            for (j = 0; j < rows; j++) {
                double value;
                EdrDoubleMatrixGet(odm, i, j, &value);
                out->dm[i][j] = value;
            }
        }

        exp = EdrDateTimeGet(odt);
        out->dt = string(exp);
        free(exp);

        exp = EdrExpiryGet(omatp);
        out->matp = string(exp);
        free(exp);

        exp = EdrExpiryGet(omattp);
        out->matpdt = string(exp);
        free(exp);

        exp = EdrExpiryGet(obm);
        out->bmdt = string(exp);
        free(exp);
           
        ok = true;
    done:
        EdrObjectFree(oi);
        EdrObjectFree(od);
        EdrObjectFree(os);
        EdrObjectFree(ob);
        EdrObjectFree(oenum);
        EdrObjectFree(odm);
        EdrObjectFree(odt);
        EdrObjectFree(omatp);
        EdrObjectFree(omattp);
        EdrObjectFree(obm);

        if (!EdrErrorCleared()) {
            char* error = EdrErrorStackTrace();
            ModelException e(method, error);
            free(error);
            throw e;
        }

        EdrShutdown();
        Library::startup();
        return IObjectSP(out);
    }
    catch (exception& e) {
        EdrShutdown();
        Library::startup();
        throw ModelException(e, method);
    }
}

class  EDRAtomicHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(EDRTester::Atomic, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultAtomic);

        FIELD(i, "int");
        FIELD(d, "double");
        FIELD(s, "string");
        FIELD(b, "bool");
        
        FIELD(dm, "double matrix");

        FIELD(dt, "date");
        FIELD(tm, "time");

        FIELD(matp, "maturity period");
        FIELD(matpdt, "mat time period (date)");
        FIELD(matptm, "mat time period (time)");
        FIELD(bmdt, "benchmark (date)");
        FIELD(bmtm, "benchmark (time)");
        FIELD(enumType, "Optional type of enum to create");
        FIELD_MAKE_OPTIONAL(enumType);
        FIELD(enumValue, "Optional value of enum to create");
        FIELD_MAKE_OPTIONAL(enumValue);
        Addin::registerClassObjectMethod("IFACE_ATOMIC_TESTER",
                                         Addin::UTILITIES,
                                         "tests EDR atomic methods",
                                         EDRTester::Atomic::TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)EDRTester::Atomic::run);


    }

    static IObject* defaultAtomic(){
        EDRTester::Atomic* p = new EDRTester::Atomic();
        return p;
        //return new EDRTester::Atomic();
    }
};

CClassConstSP const EDRTester::Atomic::TYPE = 
CClass::registerClassLoadMethod(
    "EDRTester::Atomic", typeid(EDRTester::Atomic), 
    EDRAtomicHelper::load);


EDRTester::Array::Array(): CObject(TYPE) {}

typedef vector<EDRObject> EDRObjectArray;

IObjectSP EDRTester::Array::run(EDRTester::Array* params)  {
    static const string method = "EDRTester::Array::run";
    try {
        bool ok = false;

        int length;
        int i;

        EDRObject oarray = 0;
        EDRObject od     = 0;

        DoubleArraySP  out;
        EDRObjectArray items;

        Library::shutdown();
        if (!EdrStartup()) {
            goto done;
        }
    
        length = params->dbls->size();

        // build up an array
        if (!(oarray = EdrArrayNew("Double", length))) {
            goto done;
        }

        for (i = 0; i < length; i++) {
            EDRObject elem = EdrDoubleNew((*params->dbls)[i]) ;
            if (!elem) {
                goto done;
            }
            items.push_back(elem);
            if (!EdrArraySet(oarray, i, elem)) {
                goto done;
            }
        }

        // stuff something on the end
        if (!(od = EdrDoubleNew(params->extra))) {
            goto done;
        }

        if (!EdrArrayAppend(oarray, od)) {
            goto done;
        }

        // now pull it to bits
        if (!EdrArrayLength(oarray, &length)) {
            goto done;
        }

        out = DoubleArraySP(new DoubleArray(length));

        for (i = 0; i < length; i++) {
            EDRObject elem = EdrArrayItem(oarray, i);
            if (!elem) {
                goto done;
            }
            double d;
            if (!(EdrDoubleGet(elem, &d))) {
                goto done;
            }

            (*out)[i] = d;

            EdrObjectFree(elem);
        }
           
        ok = true;
    done:
        EdrObjectFree(oarray);
        EdrObjectFree(od);
        for (i = 0; i < (int)items.size(); i++) {
            EdrObjectFree(items[i]);
        }

        if (!EdrErrorCleared()) {
            char* error = EdrErrorStackTrace();
            ModelException e(method, error);
            free(error);
            throw e;
        }

        EdrShutdown();
        Library::startup();
        return IObjectSP(out);
    }
    catch (exception& e) {
        EdrShutdown();
        Library::startup();
        throw ModelException(e, method);
    }
}

class  EDRArrayHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(EDRTester::Array, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultArray);

        FIELD(dbls, "list of doubles");
        FIELD(extra, "added to end of list");

        Addin::registerClassObjectMethod("IFACE_ARRAY_TESTER",
                                         Addin::UTILITIES,
                                         "tests EDR Array methods",
                                         EDRTester::Array::TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)EDRTester::Array::run);


    }

    static IObject* defaultArray(){
        return new EDRTester::Array();
    }
};

CClassConstSP const EDRTester::Array::TYPE = 
CClass::registerClassLoadMethod(
    "EDRTester::Array", typeid(EDRTester::Array), 
    EDRArrayHelper::load);

    
EDRTester::TypeInfo::TypeInfo(): CObject(TYPE) {}

typedef vector<EDRObject> EDRObjectTypeInfo;

IObjectSP EDRTester::TypeInfo::run(EDRTester::TypeInfo* params)  {
    static const string method = "EDRTester::TypeInfo::run";
    try {
        bool ok = false;

        int i;
        int numRows;

        const char *  typekey = params->typeName.c_str();
        ObjectArraySP output;

        Library::shutdown();
        if (!EdrStartup()) {
            goto done;
        }
        // is it an enum?
        EDRBool isEnum;
        if (!EdrTypeIsEnum(typekey, &isEnum)){
            goto done;
        }

        if (!isEnum){
            output.reset(new ObjectArray(4)); // 4 columns

            CStringArraySP paramName;
            CStringArraySP paramType;
            CStringArraySP paramOptional;
            CStringArraySP paramInfo;
            
            if (!EdrTypeNumItems(typekey, &numRows)) {
                goto done;
            }
            
            // construct our arrays
            paramName = CStringArraySP(new StringArray(numRows+1));
            paramType = CStringArraySP(new StringArray(numRows+1));
            paramOptional = CStringArraySP(new StringArray(numRows+1));
            paramInfo = CStringArraySP(new StringArray(numRows+1));
            
            // save in the output
            (*output)[0] = paramName;
            (*output)[1] = paramType;
            (*output)[2] = paramOptional;
            (*output)[3] = paramInfo;
            
            // bit of a cheat, recycle first element to store modifiers
            int modifiers;
            if (!EdrTypeGetModifiers(typekey, &modifiers)) {
                goto done;
            }
            
            if (modifiers & EDR_PUBLIC){
                (*paramName)[0] = "public";
            }
            if (modifiers & EDR_PROTECTED){
                (*paramType)[0] = "protected";
            }
            if (modifiers & EDR_ABSTRACT){
                (*paramOptional)[0] = "abstract";
            }
            if (modifiers & EDR_INTERFACE){
                (*paramInfo)[0] = "interface";
            }
            
            for (i = 1; i < numRows+1; i++){
                (*paramName)[i] = EdrTypeFieldName(typekey, i-1);
                (*paramType)[i] = EdrTypeFieldType(typekey, i-1);
                
                int optional;
                EdrTypeFieldIsOptional(typekey, i-1, &optional);
                (*paramOptional)[i] = optional ? "Optional": "Mandatory";
                (*paramInfo)[i] = EdrTypeFieldDescription(typekey, i-1);
            }
        } else {
            output.reset(new ObjectArray(2)); // 2 columns
            CStringArraySP enumValues = CStringArray::SP();
            CStringArraySP enumDescriptions = CStringArray::SP();
            for (int i = 0; i < 2; i++){
                const char** edrEnumValues = 
                    EdrEnumValues(typekey, i == 0? EDR_FALSE: EDR_TRUE);
                if (!edrEnumValues){
                    goto done;
                }
                for (const char** pos = edrEnumValues; *pos; pos++){
                    (i == 0? enumValues: enumDescriptions)->
                        push_back(string(*pos));
                }
                free(edrEnumValues);
            }
            // save in the output
            (*output)[0] = enumValues;
            (*output)[1] = enumDescriptions;
        }            
        ok = true;
    done:
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

class  EDRTypeInfoHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(EDRTester::TypeInfo, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultTypeInfo);
        FIELD(typeName, "type name");

        Addin::registerClassObjectMethod("IFACE_TYPEINFO_TESTER",
                                         Addin::UTILITIES,
                                         "tests EDR type info methods",
                                         EDRTester::TypeInfo::TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)EDRTester::TypeInfo::run);


    }

    static IObject* defaultTypeInfo(){
        return new EDRTester::TypeInfo();
    }
};

CClassConstSP const EDRTester::TypeInfo::TYPE = 
CClass::registerClassLoadMethod(
    "EDRTester::TypeInfo", typeid(EDRTester::TypeInfo), 
    EDRTypeInfoHelper::load);
  

class EDRTester::ClassDescription: public CObject {       
    string        typeName;
    StringArray   properties;
    StringArray   constructors;
    StringArray   fieldNames;
    StringArray   fieldTypes;
    StringArray   fieldFlags;
    StringArray   fieldDescriptions;
public:
    static CClassConstSP const TYPE;
    ClassDescription(const string& typeName): 
        CObject(TYPE), typeName(typeName) {
        static const string method = "ClassDescription::ClassDescription";
        try {
            // fill up the description
            const char* typekey = typeName.c_str();
            properties = StringArray(0);
            
            int modifiers;
            if (!EdrTypeGetModifiers(typekey, &modifiers)) {
                throw ModelException(method);
            }
            
            if (modifiers & EDR_PUBLIC){
                properties.push_back("public");
            }
            if (modifiers & EDR_PROTECTED){
                properties.push_back("protected");
            }
            if (modifiers & EDR_ABSTRACT){
                properties.push_back("abstract");
            }
            if (modifiers & EDR_INTERFACE){
                properties.push_back("interface");
            }
            
            // list constructors for protected types
            constructors = StringArray(0);
            const char** builders = 0;
            
            if (modifiers & EDR_PROTECTED){
                builders = EdrTypeListConstructorTypes(typekey);
                if (builders) {
                    for (int i = 0; builders[i]; i++){
                        constructors.push_back(builders[i]);
                    }            
                }
            }
            // list fields for public types
            int numRows;
            if (modifiers & EDR_PUBLIC){
                if (!EdrTypeNumItems(typekey, &numRows)) {
                    throw ModelException(method);
                }
            } else {
                numRows = 0;
            }

            fieldNames        = StringArray(numRows);
            fieldTypes        = StringArray(numRows);
            fieldFlags        = StringArray(numRows);
            fieldDescriptions = StringArray(numRows);
            
            for (int i = 0; i < numRows; i++){
                fieldNames[i] = EdrTypeFieldName(typekey, i);
                fieldTypes[i] = EdrTypeFieldType(typekey, i);
                
                int optional;
                EdrTypeFieldIsOptional(typekey, i, &optional);
                fieldFlags[i] = optional ? "Optional": "Mandatory";
                fieldDescriptions[i] = EdrTypeFieldDescription(typekey, i);
            }
            free(builders);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }


    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const{

        int i;
        string myprefix = prefix == "" ? "" : prefix + "_";
        string catenatedProps;
        for (i = 0; i < properties.size(); i++) {
            if (i > 0){
                catenatedProps += ".";
            }
            catenatedProps += properties[i];
        }
        stream << linePrefix << myprefix << typeName << "_" << 
                "properties: " << catenatedProps << endl;

        for (i = 0; i < constructors.size(); i++) {
            stream << linePrefix << myprefix << typeName << "_" << 
                "constructors[" << i << "]: " << constructors[i] << endl;
        }
        for (i = 0; i < fieldNames.size(); i++) {
            stream << linePrefix << myprefix << typeName << "_" << 
                "field_type[" << fieldNames[i] << "]: " << 
                fieldTypes[i] << endl;
            stream << linePrefix << myprefix << typeName << "_" << 
                "field_optional[" << fieldNames[i] << "]: " 
                   << fieldFlags[i] << endl;
        }
    }
    
private:
    ClassDescription(): CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ClassDescription, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultClassDescription);
        FIELD(typeName, "type name");
        FIELD(properties, "properties");
        FIELD(constructors, "constructors");
        FIELD(fieldNames, "fieldNames");
        FIELD(fieldTypes, "fieldTypes");
        FIELD(fieldFlags, "fieldFlags");
        FIELD(fieldDescriptions, "fieldDescriptions");
    }

    static IObject* defaultClassDescription(){
        return new ClassDescription();
    }
};


CClassConstSP const EDRTester::ClassDescription::TYPE = 
CClass::registerClassLoadMethod(
    "EDRTester::ClassDescription", typeid(EDRTester::ClassDescription), load);

// noddy wrapper to write out an array of ClassDescriptions
class ClassDescriptionArray: public CObject {       
public:
    static CClassConstSP const TYPE;

    ClassDescriptionArray(): CObject(TYPE), theArray(new ObjectArray(0)) {}
    
    void push_back(EDRTester::ClassDescription* cd) {
        theArray->push_back(IObjectSP(cd));
    }

    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const {
        for (int i = 0; i < theArray->size(); i++) {
            (*theArray)[i]->outputWrite(linePrefix, prefix, stream);
        }
    }
    
private:
    ObjectArraySP theArray;
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ClassDescriptionArray, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultClassDescriptionArray);
        FIELD(theArray, "list of descriptions");
    }

    static IObject* defaultClassDescriptionArray(){
        return new ClassDescriptionArray();
    }
};

CClassConstSP const ClassDescriptionArray::TYPE = 
CClass::registerClassLoadMethod(
    "ClassDescriptionArray", typeid(ClassDescriptionArray), load);



class EDRTester::PublicInterface: public CObject {
public:
    static CClassConstSP const TYPE;       
    static IObjectSP run(PublicInterface* params){
        static const string method = "EDRTester::PublicInterface::run";
        try {
            ClassDescriptionArray* output = new ClassDescriptionArray();
            
            Library::shutdown();
            if (!EdrStartup()) {
                throw ModelException(method);
            }
            
            const char** allTypes = EdrTypeList();
            
            if (!allTypes) {
                throw ModelException(method);
            }
            
            StringArray sorted(0);
            for (int i =0; allTypes[i]; i++){
                sorted.push_back(allTypes[i]);
            }
            
            // now sort them
            sort(sorted.begin(), sorted.end());
            
            for (int j = 0; j < sorted.size(); j++) {
                output->push_back(new ClassDescription(sorted[j]));
            }
            
            free(allTypes);
            EdrShutdown();
            Library::startup();
            return IObjectSP(output);
        }
        catch (exception& e) {
            EdrShutdown();
            Library::startup();
            throw ModelException(e, method);
        } 
    }

private:
    PublicInterface(): CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PublicInterface, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPublicInterface);

        Addin::registerClassObjectMethod("IFACE_PUBLIC_INTERFACE_TESTER",
                                         Addin::UTILITIES,
                                         "tests EDR PublicInterface methods",
                                         EDRTester::PublicInterface::TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)EDRTester::PublicInterface::run);
    }

    static IObject* defaultPublicInterface(){
        return new PublicInterface();
    }
};

CClassConstSP const EDRTester::PublicInterface::TYPE = 
CClass::registerClassLoadMethod(
    "EDRTester::PublicInterface", typeid(EDRTester::PublicInterface), load);
 

DRLIB_END_NAMESPACE
