
#include "edginc/config.hpp"
#define QLIB_OUTPUTNAME_CPP
#include "edginc/OutputName.hpp"
#include "edginc/Format.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/AtomicArray.hpp"
#include ext_hash_set

DRLIB_BEGIN_NAMESPACE

string const OutputName::DELIMITER = "_";

/** Destructor */
OutputName::~OutputName(){
    delete name2;
}

/** Creates a new output name using the given string */
OutputName::OutputName(const string& name1):
    CObject(TYPE), name1(name1), name2(0){
}

/** Creates a new output name using the given pair of strings */
OutputName::OutputName(const string& name1, const string& name2):
    CObject(TYPE), name1(name1), name2(new string(name2)){
}

/** Creates a new output name using the two output names. */
OutputName::OutputName(const OutputName* output1, const OutputName* output2):
    CObject(TYPE), name2(0){
    static const string routine("OutputName::OutputName");
    if (!output1 || !output2){
        throw ModelException(routine, "NULL inputs");
    }
    /* this reveals that the design is not yet complete. Are we going to
       restrict identities of objects (exluding results) to a single string? 
       If so, then we're okay as long as the results don't refer to more than
       two objects */
    if (output1->name2 || output2->name2){
        // hmm problems
        throw ModelException(routine, "Too many identifiers for output");
    }
    name1 = output1->name1;
    name2 = new string(output2->name1);
}
   
/** Creates a new output name using the two output names. */
OutputNameSP OutputName::SP(OutputNameConstSP name1, OutputNameConstSP name2) {
    return OutputNameSP(new OutputName(name1.get(), name2.get()));
}

/** Creates a new output name using the given string. */
OutputNameSP OutputName::SP(const string& name) {
    return OutputNameSP(new OutputName(name));
}

/** returns a hashcode for the object based upon the strings making up
    the name */
int OutputName::hashCode() const{
    if (name2){
        return (hash_string(name1) ^ hash_string(*name2));
    } else {
        return hash_string(name1);
    }
}
    
bool OutputName::equals(const OutputName* outputName) const{
    if (outputName == 0 || name1 != outputName->name1){
        return false;
    }
    if (name2 == 0 && outputName->name2 == 0){
        return true;
    }
    if (name2 == 0 || outputName->name2 == 0){
        return false;
    }
    
    return *name2 == *outputName->name2;
}

bool OutputName::equals(const string& val) const{
    if (name2 != 0 || !(val == name1)){
        return false;
    }
    return true;
}

/** Returns true if the this name identity is given by two strings
    which matche the one provided */
bool OutputName::equals(const string& val1, const string val2) const{
    if (name2 == 0 || !(val1 == name1) || !(val2 == *name2)){
        return false;
    }
    return true;
}

/** formatted version of outputname  */
string OutputName::toString() const{
    if (name2 != 0){
        return name1 + DELIMITER + *name2;
    } else {
        return name1;
    }
}

/** Returns true if output name is just "" */
bool OutputName::isEmpty() const{
    return (name1.empty() && (!name2 || name2->empty()));
}

/** Returns the number of strings which identifies this output name */
int OutputName::idCount() const{
    return (name2? 2: 1);
}

/** Returns the component string, identified by index, which
    identifies this output name */
const string& OutputName::idGet(int index) const{
    if (index == 0){
        return name1;
    }
    if (index == 1 && name2){
        return *name2;
    }
    throw ModelException("OutputName::idGet", 
                         "Index "+Format::toString(index)+" out of range");
}

/** boil an array of name pairs down into a list of single names */
OutputNameArraySP OutputName::singleNameArray(OutputNameArrayConstSP names) {
    static const string method = "OutputName::singleNameArray";
    try {
        int i;

        OutputNameArraySP singleNames(new OutputNameArray(0));

        for (i = 0; i < names->size(); i++) {
            singleNames->push_back(OutputNameSP(new OutputName((*names)[i]->name1)));
            if ((*names)[i]->idCount() > 1) {
                singleNames->push_back(OutputNameSP(new OutputName(*(*names)[i]->name2)));
            }               
        }

        return singleNames;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** is an array of names the same name again & again ? */
bool OutputName::uniqueNameArray(OutputNameArrayConstSP names) {
    bool unique = true;
    for (int i=1; i < names->size() && unique; ++i) {
        if (!(*names)[i]->equals((*names)[i-1].get())) {
            unique = false;
        }
    }
    return unique;
}

/** remove duplicates and empty names from an array. Order is otherwise
    unchanged */
OutputNameArraySP OutputName::trim(OutputNameArrayConstSP names) {
    int size = names->size();
    // performance improvement for common cases
    if (size == 1){
        const OutputNameSP& name = names->front();
        if (name->isEmpty()){
            return OutputNameArraySP(new OutputNameArray());
        }
        return OutputNameArraySP(new OutputNameArray(1, name));
    }
    OutputNameArraySP trimmed(new OutputNameArray());
    if (size == 0){
        return trimmed;
    }
    hash_set<OutputNameConstSP, OutputName::HashUtil, OutputName::HashUtil> s;
    for (int i = 0; i < size; i++){
        const OutputNameSP& name = (*names)[i];
        if (!name->isEmpty() && s.find(name) == s.end()){
            trimmed->push_back(name);
            s.insert(name);
        }
    }
    return trimmed;
}

/** Returns the intersection of the two arrays ie the names that appear
    in both lists */
OutputNameArraySP OutputName::intersection(OutputNameArrayConstSP set1,
                                           OutputNameArrayConstSP set2){
    hash_set<OutputNameConstSP, OutputName::HashUtil, 
        OutputName::HashUtil> s(set1->begin(), set1->end());
    OutputNameArraySP intersection(new OutputNameArray());
    // there may be an stl function to do this but by the time I spent 5 mins
    // looking for one it seemed stupid - it's not hard.
    for (int i = 0; i < set2->size(); i++){
        const OutputNameSP& name = (*set2)[i];
        if (s.find(name) != s.end()){
            intersection->push_back(name);
        }
    }      
    return intersection;
}
        
/** Returns the 'difference' of the two arrays ie the names that appear
    in the first list but not in the second */
OutputNameArraySP OutputName::difference(OutputNameArrayConstSP namesToInclude,
                                         OutputNameArrayConstSP namesToExclude){
    hash_set<OutputNameConstSP, OutputName::HashUtil, 
        OutputName::HashUtil> s(namesToExclude->begin(), namesToExclude->end());
    OutputNameArraySP difference(new OutputNameArray());
    for (int i = 0; i < namesToInclude->size(); i++){
        const OutputNameSP& name = (*namesToInclude)[i];
        if (s.find(name) == s.end()){
            // not excluded so add to the list
            difference->push_back(name);
        }
    }      
    return difference;
}
        
/** Appends to an array */
void OutputName::append(OutputNameArray& names,
                        OutputNameArrayConstSP namesToAppend) {
    for (int n = 0; n < namesToAppend->size(); ++n)
        names.push_back((*namesToAppend)[n]);
}

/** for reflection purposes */
OutputName::OutputName(): CObject(TYPE), name2(0) {
    // empty
}

bool OutputName::HashUtil::operator()(const OutputNameConstSP& name1, 
                                      const OutputNameConstSP& name2) const{
    return name1->equals(name2.get());
}
int OutputName::HashUtil::operator()(const OutputNameConstSP& name) const{
    return name->hashCode();
}


// for DR Interface purposes
class OutputNameProxy: public CObject {
public:
    static CClassConstSP const TYPE;

private:
    StringArray names;

    typedef smartPtr<OutputNameProxy> OutputNameProxySP;

    static IObjectSP toOutputName(const IObjectConstSP& obj) {
        const OutputNameProxy* proxy = dynamic_cast<const OutputNameProxy*>(obj.get());
        if (!proxy) {
            throw ModelException("OutputNameProxy::toOutputName",
                                 "object is not a OutputNameProxy");
        }
    
        if (proxy->names.size() > 2) {
            throw ModelException("OutputNameProxy::toOutputName", 
                                 "too many name supplied (max is 2)");
        }

        OutputNameSP outputName;
        if (proxy->names.size() == 2) {
            outputName = OutputNameSP(new OutputName(proxy->names[0],
                                                     proxy->names[1]));
        }
        else {
            string name = proxy->names.empty() ? "" : proxy->names[0];
            outputName = OutputNameSP(new OutputName(name));
        }
        return outputName;
    }

    static IObjectSP fromOutputName(const IObjectConstSP& obj) {
        const OutputName* outputName = dynamic_cast<const OutputName*>(obj.get());
        if (!outputName) {
            throw ModelException("OutputNameProxy::fromOutputName",
                                 "object is not a OutputName");
        }
    
        OutputNameProxySP proxy(new OutputNameProxy());
        for (int i = 0; i < outputName->idCount(); i++) {
            proxy->names.push_back(outputName->idGet(i));
        }
        return proxy;
    }

    OutputNameProxy() : CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); 
        REGISTER(OutputNameProxy, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultOutputNameProxy);
        registerObjectProxy(OutputName::TYPE,
                            OutputNameProxy::TYPE,
                            fromOutputName,
                            toOutputName);
        FIELD(names, "names");
    }

    static IObject* defaultOutputNameProxy(){
        return new OutputNameProxy();
    }  
};

CClassConstSP const OutputNameProxy::TYPE = CClass::registerClassLoadMethod(
    "OutputNameProxy", typeid(OutputNameProxy), load);

class OutputNameHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDRIProxyType(OutputNameProxy::TYPE); // use proxy for dri
        REGISTER(OutputName, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultOutputName);
        FIELD(name1, "main identifier");
        FIELD(name2, "secondary identifier");
        FIELD_MAKE_OPTIONAL(name2);
    }
    static IObject* defaultOutputName(){
        return new OutputName();
    }
};

CClassConstSP const OutputName::TYPE = CClass::registerClassLoadMethod(
    "OutputName", typeid(OutputName), OutputNameHelper::load);

DEFINE_TEMPLATE_TYPE(OutputNameArray);

DRLIB_END_NAMESPACE
