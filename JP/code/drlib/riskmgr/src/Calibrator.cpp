//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Calibrator.cpp
//
//   Description : Calibrator
//
//   Author      : regis Guichard
//
//   Date        : 21 May 2002
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/Format.hpp"
#include "edginc/Results.hpp"
#include "edginc/ObjectIteration.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Maths.hpp"
#include "edginc/NotApplicable.hpp"
#include "edginc/NonPricingModel.hpp"
#include <set>
#include ext_hash_map
#include "edginc/ClientRunnable.hpp"
#include "edginc/SpreadSheetMode.hpp"
#include "edginc/Timer.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE
static void loadAdjust(CClassSP& clazz){
    REGISTER_INTERFACE(IObject, clazz);
    EXTENDS(IObject);
}

// IADUSTABLE INTERFACE
CClassConstSP const Calibrator::IAdjustable::TYPE = 
CClass::registerInterfaceLoadMethod(
    "Calibrator::IAdjustable", typeid(IAdjustable), loadAdjust);

Calibrator::IAdjustable::~IAdjustable(){}
class Calibrator::IAdjustable::FieldInfo{
    /** Utility class for use in hashmap templates */
    struct Hash{
        size_t operator()(const CFieldConstSP p) const {return (size_t)p;}
    };
    friend class Calibrator_ObjectiveFunc;
    friend class Calibrator::IAdjustable;
public:
    typedef hash_map<CFieldConstSP, RangeSP, Hash> RangeForField;
    typedef hash_map<CFieldConstSP, TGetExpiriesMethod, Hash> GetExpiriesMethodForField;
private:
    static RangeForField adjFields;
    static GetExpiriesMethodForField getExpiriesMethods;
public:
};
Calibrator::IAdjustable::FieldInfo::RangeForField 
Calibrator::IAdjustable::FieldInfo::adjFields;
Calibrator::IAdjustable::FieldInfo::GetExpiriesMethodForField 
Calibrator::IAdjustable::FieldInfo::getExpiriesMethods;

/** Stores range in map of fields */
void Calibrator::IAdjustable::registerField(const CClassConstSP& clazz, 
                                            const string&        fieldName,
                                            Range*               range){
    CFieldConstSP field = clazz->getDeclaredField(fieldName);
    FieldInfo::adjFields[field] = RangeSP(range);
}

/** Stores range in map of fields */
void Calibrator::IAdjustable::registerBootstrappableField(const CClassConstSP& clazz, 
                                                          const string&        fieldName,
                                                          Range*               range,
                                                          TGetExpiriesMethod   getExpiriesMethod){
    CFieldConstSP field = clazz->getDeclaredField(fieldName);
    FieldInfo::adjFields[field] = RangeSP(range);
    FieldInfo::getExpiriesMethods[field] = getExpiriesMethod;
}

/** Returns the recorded range for the specified field */
const Range& Calibrator::IAdjustable::getRange(const CFieldConstSP& field){
    FieldInfo::RangeForField::const_iterator iter = 
        FieldInfo::adjFields.find(field);
    if (iter != FieldInfo::adjFields.end()){
        return (*iter->second.get());
    }
    throw ModelException("Calibrator::IAdjustable::getRange",
                         "Field '"+field->getName()+"' in type "+
                         field->getDeclaringClass()->getName()+
                         " not registered with Calibrator");
}

Calibrator::IAdjustable::TGetExpiriesMethod
Calibrator::IAdjustable::getGetExpiriesMethod(const CFieldConstSP& field){
    FieldInfo::GetExpiriesMethodForField::const_iterator iter = 
        FieldInfo::getExpiriesMethods.find(field);
    if (iter != FieldInfo::getExpiriesMethods.end() && iter->second != 0){
        return iter->second;
    }
    throw ModelException("Calibrator::IAdjustable::getGetExpiriresMethod",
                         "Field '"+field->getName()+"' in type "+
                         field->getDeclaringClass()->getName()+
                         " has no associated 'getExpirires' method");
}

bool Calibrator::IAdjustable::hasGetExpiriesMethod(const CFieldConstSP& field){
    FieldInfo::GetExpiriesMethodForField::const_iterator iter = 
        FieldInfo::getExpiriesMethods.find(field);
    return iter != FieldInfo::getExpiriesMethods.end() && iter->second != 0;
}

/** Returns all fields that are registered with the calibrator */
CFieldArray  Calibrator::IAdjustable::getFields(const CClassConstSP& clazz){
    // get all (inc parent) fields in the class
    CFieldArray allFields(Addin::getDataClassFields(clazz));
    CFieldArray fields;
    for (unsigned int i = 0; i < allFields.size(); i++){
        FieldInfo::RangeForField::const_iterator iter = 
            FieldInfo::adjFields.find(allFields[i]);
        if (iter != FieldInfo::adjFields.end()){
            fields.push_back(allFields[i]);
        }
    }
    return fields;
}

/** Checks fields (which can be calibrated) in object are within range */
void Calibrator::IAdjustable::checkRange(const IAdjustable* adjustable){
    static const string method("Calibrator::IAdjustable::checkRange");
    try{
        // get the relevant fields for this object
        CFieldArray fields(getFields(adjustable->getClass()));
        // turn pointer into SP
        IObjectConstSP obj(IObjectConstSP::attachToRef(adjustable));
        // then loop over them
        for (unsigned int i = 0; i < fields.size(); i++){
            const Range& range = getRange(fields[i]);
            const string& fieldName = fields[i]->getName();
            // only support doubles or double arrays at the moment
            if (fields[i]->getType() == CDouble::TYPE){
                double val = fields[i]->getDouble(obj);
                Range::checkVariableIsInRange(range, val, fieldName);
            } else if (fields[i]->getType() == DoubleArray::TYPE){
                IObjectConstSP dbArrayObj = fields[i]->constGet(obj);
                if (dbArrayObj.get()){
                    const DoubleArray& dbArray =
                        dynamic_cast<const DoubleArray&>(*dbArrayObj);
                    for (int i = 0; i < dbArray.size(); i++){
                        Range::checkVariableIsInRange(range, dbArray[i], 
                                                      fieldName);
                    }
                }   
            } else {
                throw ModelException(method, "Only doubles or DoubleArrays "
                                     "supported");
            }
        }
        
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

// OBJFUNC
CClassConstSP const Calibrator::ObjFunc::TYPE = 
CClass::registerClassLoadMethod("Calibrator::ObjFunc", 
                                typeid(Calibrator::ObjFunc), load);

void Calibrator::ObjFunc::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(Calibrator::ObjFunc, clazz);
    SUPERCLASS(CObject);
}

Calibrator::ObjFunc::ObjFunc(const CClassConstSP& clazz):
    CObject(clazz){}

void Calibrator::ObjFunc::Helper::getMarket(CInstrument* inst,
                                            IModel*      model,
                                            CControl*    control,
                                            MarketData*  market){
    try{
        // call model->getInstrumentAndModelMarket to initiate market data selection
        model->getInstrumentAndModelMarket(market, inst);

        // call getMarket to initiate market data selection
        control->getMarket(IModelConstSP::attachToRef(model),
                           MarketDataConstSP::attachToRef(market),
                           IInstrumentCollection::singleton(inst));
    }
    catch(exception& e){
        throw ModelException(e, "Calibrator::ObjFunc::Helper::getMarket");
    }
}

void Calibrator::ObjFunc::validate(){}
typedef Calibrator::ObjFuncArray CalibratorObjFuncArray; // MSVC7 bug
DEFINE_TEMPLATE_TYPE_WITH_NAME("Calibrator::ObjFuncArray", CalibratorObjFuncArray);

// OBJFUNC LEASTSQUARE
CClassConstSP const Calibrator::ObjFuncLeastSquare::TYPE = 
CClass::registerClassLoadMethod("Calibrator::ObjFuncLeastSquare", 
                                typeid(Calibrator::ObjFuncLeastSquare), load);

void Calibrator::ObjFuncLeastSquare::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(Calibrator::ObjFuncLeastSquare, clazz);
    SUPERCLASS(Calibrator::ObjFunc);
    FIELD(funcs, "");
    FIELD_MAKE_TRANSIENT(funcs);
}

Calibrator::ObjFuncLeastSquare::ObjFuncLeastSquare(const CClassConstSP& clazz):
    Calibrator::ObjFunc(clazz), funcs(0), objFuncNormalized(""){}


double Calibrator::ObjFuncLeastSquare::calcValue() const{
    try{
        int nbFuncs = getNbFuncs();

        if (!funcs){
            funcs = CDoubleArraySP(new CDoubleArray(nbFuncs));
        } else {
            if (nbFuncs != funcs->size()) {
                // nbFuncs may vary for the same objective function
                funcs = CDoubleArraySP(new CDoubleArray(nbFuncs));
            }
        }
        DoubleArray& funcvals = *funcs;
        calcValue(funcvals);
        double res = 0.0;
        int iFunc = 0;
        for (; iFunc < nbFuncs; ++iFunc){
            res += Maths::square(funcvals[iFunc]);
        }
        if (CString::equalsIgnoreCase(objFuncNormalized,Calibrator::ObjFuncLeastSquare::DO_NOT_NORMALIZE)){
            return res;
        }
        else{
            return (res / nbFuncs);
        }
    }
    catch(exception& e){
        throw ModelException(e, "Calibrator::ObjFuncLeastSquare::calcValue");
    }    
}

void Calibrator::ObjFuncLeastSquare::isObjFuncNormalized(string doNormalize){
    objFuncNormalized = doNormalize;
}

string Calibrator::ObjFuncLeastSquare::getNormalizedWeights(){
    return objFuncNormalized;
}

const string Calibrator::ObjFuncLeastSquare::DO_NOT_NORMALIZE  = "do_not_normalize"; 

// work around for msvc 7 bug
typedef Calibrator::ObjFuncLeastSquareArray CalibratorObjFuncLeastSquareArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("Calibrator::ObjFuncLeastSquareArray", CalibratorObjFuncLeastSquareArray);

// INSTANCEID
Calibrator::InstanceID::InstanceID(const CClassConstSP&  clazz,
                                   const string&  typeToCalibrate,
                                   const string&  name, 
                                   const string&  fieldName): 
    CObject(clazz),
    typeToCalibrate(typeToCalibrate), 
    name(name), 
    fieldName(fieldName),
    classToCalibrate(0), 
    field(0){}

// for default constructor
Calibrator::InstanceID::InstanceID(const CClassConstSP&  clazz): 
    CObject(clazz),
    classToCalibrate(0), 
    field(0){}

void Calibrator::InstanceID::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(InstanceID, clazz);
    SUPERCLASS(CObject);
    FIELD(typeToCalibrate,"Identifies type of object to calibrate");
    FIELD(name, "Identifies instance of object to calibrate");
    FIELD(fieldName, "Identifies field in object to calibrate");
}

void Calibrator::InstanceID::validatePop2Object(){
    static const string method("Calibrator::InstanceID:"
                               "validatePop2Object");
    classToCalibrate = CClass::forName(typeToCalibrate);
    if (!Calibrator::IAdjustable::TYPE->
        isAssignableFrom(classToCalibrate)){
        throw ModelException(method,
                             "Object of type "+
                             typeToCalibrate+" does not "
                             "support calibration");
    }
    CClassConstSP c = classToCalibrate;
    field = 0;
    do {
        field = c->hasDeclaredField(fieldName);
    } while (!field && (c = c->getSuperClass()) != 0);
    if (!field){
        throw ModelException(method,
                             fieldName+ " not found in "+ typeToCalibrate);
    }
}
    
IObject* Calibrator::InstanceID::clone() const{
    IObject* myCopy = CObject::clone();
    Calibrator::InstanceID& id = 
        dynamic_cast<Calibrator::InstanceID&>(*myCopy);
    id.classToCalibrate = classToCalibrate;
    id.field = field;
    return myCopy;
}

/** Applies adjustments to array of InstanceID's */
void Calibrator::InstanceID::applyAdjustment(const InstanceIDArray& ids,
                                             const IObjectSP&       adjGroup,
                                             const CDoubleArray&    x){
    int index = 0;
    for (int i = 0; i < ids.size(); i++){
        ids[i]->applyAdjustment(adjGroup, x, index);
        index += ids[i]->numVariables();
    }
}

/** Gets names from array of InstanceID's */
void Calibrator::InstanceID::getNames(const InstanceIDArray& ids,
                                      StringArray&           names){
    int index = 0;
    for (int j = 0; j < ids.size(); j++){
        ids[j]->getNames(names, index);
        index += ids[j]->numVariables();
    }
}

/** Gets ranges from array of InstanceID's */
void Calibrator::InstanceID::getRanges(const InstanceIDArray& ids,
                                       RangeArray&            ranges){
    int index = 0;
    for (int j = 0; j < ids.size(); j++){
        ids[j]->getRanges(ranges, index);
        index += ids[j]->numVariables();
    }
}

/** Gets initial values from array of InstanceID's */
void Calibrator::InstanceID::getInitialValues(const InstanceIDArray& ids,
                                              const IObjectSP&       adjGroup,
                                              DoubleArray&           x){
    int index = 0;
    for (int i = 0; i < ids.size(); i++){
        ids[i]->getInitialValues(x, index, adjGroup);
        index += ids[i]->numVariables();
    }
}

/** Gets current values from array of InstanceID's */
void Calibrator::InstanceID::getValues(const InstanceIDArray& ids,
                                       const IObjectSP&       adjGroup,
                                       DoubleArray&           x){
    int index = 0;
    for (int i = 0; i < ids.size(); i++){
        ids[i]->getValues(x, index, adjGroup);
        index += ids[i]->numVariables();
    }
}

/** Put the calibrated results into some sort of context -
    not clear what is the best way of doing this at the moment */
IObjectSP Calibrator::InstanceID::writeResults(const InstanceIDArray& ids,
                                               const DoubleArray&     calibratedVals){
    int numVars = ids.size();
    StringArraySP types(new StringArray(0));
    types->reserve(numVars);
    StringArraySP names(new StringArray(0));
    names->reserve(numVars);
    StringArraySP fieldNames(new StringArray(0));
    fieldNames->reserve(numVars);
    ObjectArraySP values(new ObjectArray(0));
    values->reserve(numVars);
    StringArraySP dbTypes(new StringArray(0));
    dbTypes->reserve(numVars);
    ObjectArraySP arrayIndices(new ObjectArray(0));
    arrayIndices->reserve(numVars);
    int index = 0;
    for (int i = 0; i < numVars; i++){
        types->push_back(ids[i]->typeToCalibrate);
        names->push_back(ids[i]->name);
        fieldNames->push_back(ids[i]->fieldName);
        values->push_back(ids[i]->writeResults(calibratedVals, index));
        dbTypes->push_back(ids[i]->isScalar() ? "SCALAR" : "VECTOR");
        arrayIndices->push_back(ids[i]->getArrayIndices());
        index += ids[i]->numVariables();
    }
    ObjectArraySP results(new ObjectArray(6));
    (*results)[OBJ_TYPE] = types;
    (*results)[NAME] = names;
    (*results)[FIELD_NAME] = fieldNames;
    (*results)[VALUE] = values;
    (*results)[DB_TYPE] = dbTypes;
    (*results)[ARRAY_INDEX] = arrayIndices;
    
    return results;
}
        
CClassConstSP const Calibrator::InstanceID::TYPE = 
CClass::registerClassLoadMethod(
    "Calibrator::InstanceID", typeid(Calibrator::InstanceID), load);

typedef Calibrator::InstanceIDArray CalibratorInstanceIDArray; // msvc7 bug
DEFINE_TEMPLATE_TYPE_WITH_NAME("Calibrator::InstanceIDArray", CalibratorInstanceIDArray);

// InstanceID for doubles
Calibrator::InstanceIDDb::InstanceIDDb(): InstanceID(TYPE), useOverride(false), 
                                          override(0.0), arrayIndex(0), rangeOverride(NULL),  skipTransient(true){}

IObject* Calibrator::InstanceIDDb::defaultCtor(){
    return new InstanceIDDb();
}

void Calibrator::InstanceIDDb::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(InstanceIDDb, clazz);
    SUPERCLASS(InstanceID);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(useOverride,"Whether to override initial value");
    FIELD(override, "Optional override for initial value");
    FIELD_MAKE_OPTIONAL(override);
    FIELD(arrayIndex, "For arrays, index of element to calibrate");
    FIELD_MAKE_OPTIONAL(arrayIndex);
	FIELD(rangeOverride,"Optional override for range of field to calibrate. Need open range for calibrator.");
	FIELD_MAKE_OPTIONAL(rangeOverride);
    FIELD(skipTransient, "whether skip transient fields");
    FIELD_MAKE_OPTIONAL(skipTransient);
}

void Calibrator::InstanceIDDb::validatePop2Object(){
    static const string method("Calibrator::InstanceIDDb");
    InstanceID::validatePop2Object();
    CClassConstSP type = field->getType();
    if (type != CDouble::TYPE && type != CDoubleArray::TYPE ){
        throw ModelException(method, "Field "+fieldName+
                             " is not of type double or DoubleArray");
    }
    if (arrayIndex < 0 || (arrayIndex > 0 && type != CDoubleArray::TYPE)){
        throw ModelException(method, "Array index must be >= 0 and is "
                             "only applicable for arrays");
    }
}

Calibrator::InstanceIDDb::InstanceIDDb(const string&  typeToCalibrate,
                                       const string&  name, 
                                       const string&  field,
                                       bool           useOverride,
                                       double         override,
                                       int            arrayIndex):
    InstanceID(TYPE, typeToCalibrate, name, field), 
    useOverride(useOverride), override(override), arrayIndex(arrayIndex), rangeOverride(NULL){
    validatePop2Object();
}
/** override constructor to allow rangeOverride */
Calibrator::InstanceIDDb::InstanceIDDb(const string&  typeToCalibrate,
	const string&  name, 
	const string&  field,
	bool           useOverride,
	double         override,
	int            arrayIndex,
	RangeSP		   rangeOverride):
InstanceID(TYPE, typeToCalibrate, name, field), 
	useOverride(useOverride), override(override), arrayIndex(arrayIndex), rangeOverride(rangeOverride)
{
		validatePop2Object();
	}

/** override constructor to allow skipTransient */
Calibrator::InstanceIDDb::InstanceIDDb(const string&  typeToCalibrate,
	const string&  name, 
	const string&  field,
	bool           useOverride,
	double         override,
	int            arrayIndex,
	bool		   skipTransient):
InstanceID(TYPE, typeToCalibrate, name, field), 
	useOverride(useOverride), override(override), arrayIndex(arrayIndex), 
    rangeOverride(NULL), skipTransient(skipTransient)
{
		validatePop2Object();
}

/** Must be called before any other methods - allows object to see the
    market data (eg length of arrays) */
void Calibrator::InstanceIDDb::initialise(const IObjectSP&  adjGroup){
}

/** Returns the number of variables */
int Calibrator::InstanceIDDb::numVariables() const{
    return 1;
}

/** Populates the names array from 'index' to 
    'index' + numVariables() */
void Calibrator::InstanceIDDb::getNames(StringArray& names, int index) const{
    names[index] = typeToCalibrate 
        + "_" + name
        + "_" + fieldName;
    if (field->getType() == CDoubleArray::TYPE){
        names[index] += "[" + Format::toString(arrayIndex) + "]";
    }
}

/** Populates the ranges array from 'index' to 
    'index' + numVariables() */
void Calibrator::InstanceIDDb::getRanges(RangeArray& ranges, int index) const
{
	// range is either rangeOverride or value defined in the field
	const Range& range = rangeOverride.get() == 0 ? 
		Calibrator::IAdjustable::getRange(field) : 
		*rangeOverride ;

    ranges[index] = RangeSP(new Range(range));
}

/** Returns true if scalar type; false if vector type */
bool Calibrator::InstanceIDDb::isScalar() const{
    CClassConstSP type = field->getType();
    return (type == CDouble::TYPE);
}

/** If vector type, returns indices of elts that got calibrated.
    Returns 0, otherwise. */
IObjectSP Calibrator::InstanceIDDb::getArrayIndices() const{
    return CIntSP(CInt::create(arrayIndex));
}
/** tiny class used for call back mechanism for initial
    values or for adjusting value */
class Calibrator::InstanceIDDb::ReadWriteVal:
    virtual public ObjectIteration::IAction
{
    bool                read;  // true: read values in object
    double              value; // source or destination
    const InstanceIDDb* instance;
    bool                found;
public:
    ReadWriteVal(bool                read,  // true: read value
                 double              value, // read from or written to
                 const InstanceIDDb* instance);
    
    //// returns value
    double getValue() const;
    
    //// was instance found
    bool instanceFound() const;
    
    /* method invoked by recurse routine each time we hit an object
       of the right type */
    bool invoke(ObjectIteration::State& state, IObjectConstSP obj);
};

Calibrator::InstanceIDDb::ReadWriteVal::ReadWriteVal(
    bool                read,  // true: read value
    double              value, // read from or written to
    const InstanceIDDb* instance): 
    read(read), value(value), instance(instance), 
    found(false) {}

//// returns value
double Calibrator::InstanceIDDb::ReadWriteVal::getValue() const{
    return value;
}

//// was instance found
bool Calibrator::InstanceIDDb::ReadWriteVal::instanceFound() const{
    return found;
}

/* method invoked by recurse routine each time we hit an object
   of the right type */
bool Calibrator::InstanceIDDb::ReadWriteVal::invoke(
    ObjectIteration::State& state, IObjectConstSP obj){
    const Calibrator::IAdjustable& adjustable = 
        dynamic_cast<const Calibrator::IAdjustable&>(*obj);
    const string& name = adjustable.getName();
    if (name == instance->name){
        // we respect the constness below, but need to duplicate code/
        // use a template to avoid the cast.
        IObjectSP nonConstObj = read? 
            IObjectSP::constCast(obj): state.getObject();
        CClassConstSP type = instance->field->getType();
        DoubleArray* dbArray;
        if (type == CDoubleArray::TYPE){
            IObjectSP dbArrayObj = instance->field->get(nonConstObj);
            dbArray = &dynamic_cast<DoubleArray&>(*dbArrayObj);
            if (dbArray->size() <= instance->arrayIndex){
                throw ModelException("InstanceIDDb::ReadWriteVal",
                                     "Array index out of bounds");
            }
        }
        if (read){
            if (type == CDouble::TYPE){
                value = instance->field->getDouble(nonConstObj);
            } else {
                value = (*dbArray)[instance->arrayIndex];
            }
        } else {
            if (type == CDouble::TYPE){
                instance->field->setDouble(nonConstObj, value);
            } else {
                (*dbArray)[instance->arrayIndex] = value;
            }
            // tell object its fields have changed
            // this may be called too many times and hit performance
            nonConstObj->fieldsUpdated(CFieldArray(1, instance->field));
        }
        found = true;
    }
    return true;
}

void Calibrator::InstanceIDDb::getInitialValues(
    DoubleArray&         values,
    int                  index,
    const IObjectSP&  adjGroup) const{
    if (useOverride){
        values[index] = override;
    } else {
        ObjectIteration iteration(classToCalibrate);
        ReadWriteVal readWriteVal(true, 0.0, this);
        // but skip transient fields
 //       iteration.setSkipTransient(true);
        iteration.setSkipTransient(skipTransient);
        iteration.recurse(readWriteVal, adjGroup); // go
        values[index] = readWriteVal.getValue();
    }
}

/** Populates the values array from 'index' to
    'index' + numVariables() */
void Calibrator::InstanceIDDb::getValues(
    DoubleArray&         values,
    int                  index,
    const IObjectSP&  adjGroup) const{
    ObjectIteration iteration(classToCalibrate);
    ReadWriteVal readWriteVal(true, 0.0, this);
    // but skip transient fields
//    iteration.setSkipTransient(true);
    iteration.setSkipTransient(skipTransient);
    iteration.recurse(readWriteVal, adjGroup); // go
    values[index] = readWriteVal.getValue();
}

/** Uses the values in x from 'index' to 
    'index' + numVariables() to adjust the field being calibrated */
void Calibrator::InstanceIDDb::applyAdjustment(const IObjectSP&  adjGroup,
                                               const CDoubleArray&  x,
                                               int                  index) const{
    ObjectIteration iteration(classToCalibrate);
    ReadWriteVal readWriteVal(false, x[index], this);
    // but skip transient fields
    //   iteration.setSkipTransient(true);
    iteration.setSkipTransient(skipTransient);
    iteration.recurse(readWriteVal, adjGroup); // go
    if (!readWriteVal.instanceFound()){
        throw ModelException("InstanceIDDb::applyAdjustment",
                             "Object with name "+name+" of type "+
                             typeToCalibrate+" not found");
    }
}

/** Return the calibratedVals[index] to 
    calibratedVals[index+numVariables] as an object */
IObjectSP Calibrator::InstanceIDDb::writeResults(const DoubleArray& calibratedVals,
                                                 int                index) const{
    return IObjectSP(CDouble::create(calibratedVals[index]));
}

/** Read off the 'values' object and populate calibratedVars array */
void Calibrator::InstanceIDDb::readResults(const IObject& values,
                                           DoubleArray&   calibratedVals,
                                           int            index) const{
    if (!CDouble::TYPE->isInstance(values)){
        throw ModelException("InstanceDb::readResults", "CDouble type expected");
    }
    const CDouble& val = dynamic_cast<const CDouble&>(values);
    calibratedVals[index] = val.doubleValue();
}

CClassConstSP const Calibrator::InstanceIDDb::TYPE = 
CClass::registerClassLoadMethod(
    "Calibrator::InstanceIDDb", typeid(Calibrator::InstanceIDDb), load);
 
// InstanceID for double arrays
class Calibrator::InstanceIDDbArray: public InstanceID,
                  public InstanceID::IBootstrappable{
private:
    int                size; // length of array - transient
    DoubleArraySP      override; // overrides initial value in object
	RangeSP			   rangeOverride; // overrides range in object

    class ReadWriteVal;
    friend class ReadWriteVal;
    class GetExpiriesAction;
    friend class GetExpiriesAction;

    InstanceIDDbArray(): InstanceID(TYPE), size(-1), rangeOverride(NULL){}
    
	static IObject* defaultCtor(){
        return new InstanceIDDbArray();
    }

    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(InstanceIDDbArray, clazz);
        SUPERCLASS(InstanceID);
        IMPLEMENTS(InstanceID::IBootstrappable);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD_NO_DESC(size);
        FIELD_MAKE_TRANSIENT(size);
        FIELD(override, "Optional override for initial value(s)");
        FIELD_MAKE_OPTIONAL(override);

		FIELD(rangeOverride,"Optional override for range of field to calibrate (do not use when bootstrapping)");
		FIELD_MAKE_OPTIONAL(rangeOverride);

    }
    /** validates we've been initialised */
    void check() const{
        if (size < 0){
            throw ModelException("Calibrator::InstanceIDDbArray", 
                                 "Not initialised");
        }
    }
public:
    static CClassConstSP const TYPE;

    virtual void validatePop2Object(){
        InstanceID::validatePop2Object();
        if (field->getType() != DoubleArray::TYPE){
            throw ModelException("Calibrator::InstanceIDDbArray", 
                                 "Field "+fieldName+
                                 " is not of type DoubleArray");
        }
    }

    InstanceIDDbArray(const string&         typeToCalibrate,
                      const string&         name, 
                      const string&         field,
                      const DoubleArraySP&  override): // can be null
        InstanceID(TYPE, typeToCalibrate, name, field), size(-1),
        override(override),
		rangeOverride(NULL)
		{
        validatePop2Object();
    }

    /** Returns the number of variables */
    virtual int numVariables() const{
        check();
        return size;
    }
    
    /** Populates the names array from 'index' to 
        'index' + numVariables() */
    virtual void getNames(StringArray& names, int index) const{
        check();
        for (int i = 0; i < size; i++){
            names[index+i] = typeToCalibrate 
                + "_" + name
                + "_" + fieldName
                + "[" + Format::toString(i) + "]";
        }
    }

    /** Populates the ranges array from 'index' to 
        'index' + numVariables() */
    virtual void getRanges(RangeArray& ranges, int index) const{
        check();

		// range is either rangeOverride or value defined in the field
		const Range& range = !rangeOverride ? 
			Calibrator::IAdjustable::getRange(field) : 
			*rangeOverride;
		
        for (int i = 0; i < size; i++){
            ranges[index+i] = RangeSP(new Range(range));
        }
    }

    /** Returns true if scalar type; false if vector type */
    virtual bool isScalar() const{
        return false;
    }

    /** If vector type, returns indices of elts that got calibrated.
        Returns 0, otherwise. */
    virtual IObjectSP getArrayIndices() const{
        IntArraySP indices(new IntArray(size));
        for (int i = 0; i < size; ++i){
            (*indices)[i] = i;
        }
        return indices;
    }

private:
    /** tiny class used for call back mechanism for initial
        values or for adjusting value */
    class ReadWriteVal: virtual public ObjectIteration::IAction{
        bool                     findLength; // true: find length of arrays
        bool                     read;   // true: read values in object
        vector<double>::iterator values; // source or destination
        const InstanceIDDbArray* instance;
        bool                     found;
    public:
        //// Constructor to find length of array
        ReadWriteVal(InstanceIDDbArray*   instance): 
            findLength(true), read(true), instance(instance), found(false){}

        ReadWriteVal(bool                     read,  // true: read value
                     vector<double>::iterator values, //read from or written to
                     const InstanceIDDbArray* instance): 
            findLength(false), read(read), values(values), 
            instance(instance), found(false) {}

        //// was instance found
        bool instanceFound() const{
            return found;
        }
        /* method invoked by recurse routine each time we hit an object
           of the right type */
        bool invoke(ObjectIteration::State& state, IObjectConstSP obj){
            const Calibrator::IAdjustable& adjustable = 
                dynamic_cast<const Calibrator::IAdjustable&>(*obj);
            const string& name = adjustable.getName();
            if (name == instance->name && instance->size != 0){
                // we respect the constness below, but need to duplicate code/
                // use a template to avoid the cast.
                IObjectSP nonConstObj = read? 
                    IObjectSP::constCast(obj): state.getObject();
                IObjectSP dbArrayObj = instance->field->get(nonConstObj);
                DoubleArray& dbArray = dynamic_cast<DoubleArray&>(*dbArrayObj);
                if (findLength){
                    // this saves a lot of code
                    const_cast<int&>(instance->size) = dbArray.size();
                } else  if (dbArray.size() != instance->size){
                    throw ModelException("InstanceIDDbArray::ReadWriteVal",
                                         "Internal error - array lengths"
                                         " differ");
                } else {
                    vector<double>::iterator pos = values;
                    for (int i = 0; i < instance->size; i++, pos++){
                        if (read){
                            *pos = dbArray[i];
                        } else {
                            dbArray[i] = *pos;
                        }
                    }
                    if (!read){
                        // tell object its fields have changed
                        //this may be called too many times and hit performance
                        nonConstObj->fieldsUpdated(
                            CFieldArray(1, instance->field));
                    }
                }
                found = true;
            }
            return true;
        }
    };
    
    /** tiny class used for call back mechanism for getting
        the expiries off the boostrappable object */
    class GetExpiriesAction: virtual public ObjectIteration::IActionConst{
        friend class InstanceIDDbArray;
        const InstanceIDDbArray* instance;
        Calibrator::IAdjustable::TGetExpiriesMethod getExpiriesMethod;
        bool                     found;
        ExpiryArraySP            expiries;

    public:
        GetExpiriesAction(const InstanceIDDbArray* instance): 
            instance(instance), found(false){
            static const string method("Calibrator::InstanceIDDbArray"
                                       "::GetExpiriesAction::GetExpiriesAction");
            try{
                getExpiriesMethod 
                    = Calibrator::IAdjustable::getGetExpiriesMethod(instance->field);
                if (!getExpiriesMethod){
                    throw ModelException(method, 
                                         "The returned getExpiriesMethod is null");
                }
            }
            catch(exception& e){
                throw ModelException(e, method);
            }
        }

        /* method invoked by recurse routine each time we hit an object
           of the right type */
        bool invoke(const ObjectIteration::State& state, IObjectConstSP obj){
            static const string method("Calibrator::InstanceIDDbArray"
                                       "::GetExpiriesAction::invoke");
            const Calibrator::IAdjustable& adjustable = 
                dynamic_cast<const Calibrator::IAdjustable&>(*obj);
            const string& name = adjustable.getName();
            if (name == instance->name){
                // get the expiries
                ExpiryArraySP expiries(getExpiriesMethod(obj.get()));
                if (!this->expiries){
                    this->expiries = expiries;
                }
                // if already there, check that thee expiries have the same
                // values in all the copies of the object that are in the scope
                else{
                    if (!Expiry::equals(expiries.get(), this->expiries.get())){
                        throw ModelException(
                            method, 
                            "Object of type "
                            + obj->getClass()->getName()
                            + " and market name "
                            + name
                            + " was found twice, but non-equal arrays "
                            + "of expiries were returned upon invokation "
                            + "of the getExpiryMethod");
                    }
                }
                found = true;
            }
            return true;
        }
    };

public:
    /** Must be called before any other methods - allows object to see the
        market data (eg length of arrays) */
    virtual void initialise(const IObjectSP&  adjGroup){
        // find length of array
        ObjectIteration iteration(classToCalibrate);
        ReadWriteVal readWriteVal(this);
        // but skip transient fields
        iteration.setSkipTransient(true);
        iteration.recurse(readWriteVal, adjGroup); // go
        if (override.get() && size != override->size()){
            throw ModelException("Calibrator::InstanceIDDbArray::initialise",
                                 "Length of override values differs from that "
                                 "in object ("+field->getName()+" "+name+")");
        }
    }

    /** Populates the values array from 'index' to
        'index' + numVariables() */
    virtual void getInitialValues(
        DoubleArray&         values,
        int                  index,
        const IObjectSP&  adjGroup) const{
        check();
        if (override.get()){
            for (int i = 0; i < size; i++){
                values[index+i] = (*override)[i];
            }
        } else {
            ObjectIteration iteration(classToCalibrate);
            ReadWriteVal readWriteVal(true, values.begin()+index, this);
            // but skip transient fields
            iteration.setSkipTransient(true);
            iteration.recurse(readWriteVal, adjGroup); // go
        }
    }
    
    /** Populates the values array from 'index' to
        'index' + numVariables() */
    virtual void getValues(DoubleArray&         values,
                           int                  index,
                           const IObjectSP&  adjGroup) const{
        check();
        ObjectIteration iteration(classToCalibrate);
        ReadWriteVal readWriteVal(true, values.begin()+index, this);
        // but skip transient fields
        iteration.setSkipTransient(true);
        iteration.recurse(readWriteVal, adjGroup); // go
    }

    /** Uses the values in x from 'index' to 
        'index' + numVariables() to adjust the field being calibrated */
    void applyAdjustment(const IObjectSP&  adjGroup,
                         const CDoubleArray&  x,
                         int                  index) const{
        check();
        ObjectIteration iteration(classToCalibrate);
        ReadWriteVal readWriteVal(false,
                                  // cast saves a lot of code
                                  const_cast<CDoubleArray&>(x).begin()+index,
                                  this);
        // but skip transient fields
        iteration.setSkipTransient(true);
        iteration.recurse(readWriteVal, adjGroup); // go
        if (!readWriteVal.instanceFound()){
            throw ModelException("InstanceIDDb::applyAdjustment",
                                 "Object with name "+name+" of type "+
                                 typeToCalibrate+" not found");
        }
    }
    
    /** Return the calibratedVals[index] to 
        calibratedVals[index+numVariables] as an object */
    virtual IObjectSP writeResults(const DoubleArray& calibratedVals,
                                   int                index) const{
        check();
        DoubleArraySP val(new DoubleArray(size));
        for (int i = 0; i < size; i++){
            (*val)[i] = calibratedVals[index+i];
        }
        return val;
    }

    /** Read off the 'values' object and populate calibratedVars array */
    virtual void readResults(const IObject& values,
                             DoubleArray&   calibratedVals,
                             int            index) const{
        check();
        if (!DoubleArray::TYPE->isInstance(values)){
            throw ModelException("InstanceDb::readResults", "CDouble type expected");
        }
        const CDoubleArray& vals = dynamic_cast<const CDoubleArray&>(values);
        for (int i = 0; i < size; i++){
            calibratedVals[index+i] = vals[i];
        }
    }

    /** Extracts the InstanceIDDb corresponding to the index idxMat */
    virtual InstanceIDSP getInstanceID(int idxMat) const{
        static const string method = "getInstanceID";
        // the index can't exceed the number of variables in the array
        if (idxMat < 0
            || idxMat >= numVariables()){
            throw ModelException("Calibrator::InstanceIDDbArray",
                                 "The maturity index ("
                                 + Format::toString(idxMat + 1)
                                 + ") should be in [1, "
                                 + Format::toString(numVariables())
                                 + "]");
        }       
        // return a InstanceIDDb object
        bool useOverride = !!override;
        return InstanceIDDbSP(
            new InstanceIDDb(typeToCalibrate, 
                             name, 
                             fieldName, 
                             useOverride, 
                             useOverride ? (*override)[idxMat]: 0.0,
                             idxMat));
    }

    virtual ExpiryArraySP getExpiries(const IObjectConstSP& adjGroup) const{
        static const string method("Calibrator::InstanceIDDbArray::getExpiries");
        try{
            check();
            ObjectIteration iteration(classToCalibrate);
            GetExpiriesAction action(this);
            // but skip transient fields
            iteration.setSkipTransient(true);
            iteration.recurse(action, adjGroup); // go
            if (!action.found){
                throw ModelException(method,
                                     "Object with name "+name+" of type "+
                                     typeToCalibrate+" not found");
            }
            if (!action.expiries){
                throw ModelException(method,
                                     "No expiries were found for field "
                                     + field->getType()->getName()
                                     + " in object with name "
                                     + name
                                     + " of type "
                                     + typeToCalibrate);
            }
            if (action.expiries->size() != numVariables()){
                throw ModelException(method,
                                     "The size ("
                                     + Format::toString(action.expiries->size())
                                     + ") of the array of expiries that was found is different from that ("
                                     + Format::toString(numVariables())
                                     + ") of the vector field of type "
                                     + field->getType()->getName()
                                     + " in object with name "
                                     + name
                                     + " of type "
                                     + typeToCalibrate);
            }
            return action.expiries;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }
};

CClassConstSP const Calibrator::InstanceIDDbArray::TYPE = 
CClass::registerClassLoadMethod(
    "Calibrator::InstanceIDDbArray", typeid(Calibrator::InstanceIDDbArray),
    load);



/** Put the calibrated results into some sort of context -
    not clear what is the best way of doing this at the moment */
pair<Calibrator::InstanceIDArray, DoubleArray> 
Calibrator::InstanceID::readResults(const IObjectConstSP& res){
    static const string method("Calibrator::InstanceID::readResults");
    try{
        // will assume res is in the right format
        // if it isn't, it shouldn't core but the error message might be
        // a bit cryptic
        if (!res){
            throw ModelException(method, "null res");
        }
        IObject& res1 = const_cast<IObject&>(*res);
        ObjectArray& res2 = dynamic_cast<ObjectArray&>(res1);
        if (res2.size() != Calibrator::InstanceID::NB_RES){
            throw ModelException(method, "res has wrong size");
        }
        // obj types
        if (!res2[Calibrator::InstanceID::OBJ_TYPE]){
            throw ModelException(method, "obj type is null");
        }
        StringArray& types = dynamic_cast<StringArray&>(*res2[Calibrator::InstanceID::OBJ_TYPE]);
        int numVars = types.size();
        // names
        if (!res2[Calibrator::InstanceID::NAME]){
            throw ModelException(method, "name is null");
        }
        StringArray& names = dynamic_cast<StringArray&>(*res2[Calibrator::InstanceID::NAME]);
        if (numVars != names.size()){
            throw ModelException(method, "name array is wrong size");
        }
        // field names
        if (!res2[Calibrator::InstanceID::FIELD_NAME]){
            throw ModelException(method, "field name is null");
        }
        StringArray& fieldNames = dynamic_cast<StringArray&>(*res2[Calibrator::InstanceID::FIELD_NAME]);
        if (numVars != fieldNames.size()){
            throw ModelException(method, "field name array is wrong size");
        }
        // values
        if (!res2[Calibrator::InstanceID::VALUE]){
            throw ModelException(method, "value is null");
        }
        ObjectArray& values = dynamic_cast<ObjectArray&>(*res2[Calibrator::InstanceID::VALUE]);
        if (numVars != values.size()){
            throw ModelException(method, "value array is wrong size");
        }
        // db types
        if (!res2[Calibrator::InstanceID::DB_TYPE]){
            throw ModelException(method, "db type is null");
        }
        StringArray& dbTypes = dynamic_cast<StringArray&>(*res2[Calibrator::InstanceID::FIELD_NAME]);
        if (numVars != dbTypes.size()){
            throw ModelException(method, "db type array is wrong size");
        }
        // array indices
        if (!res2[Calibrator::InstanceID::ARRAY_INDEX]){
            throw ModelException(method, "array index is null");
        }
        ObjectArray& arrayIndices = dynamic_cast<ObjectArray&>(*res2[Calibrator::InstanceID::ARRAY_INDEX]);
        if (numVars != arrayIndices.size()){
            throw ModelException(method, "array index array is wrong size");
        }
        // loop through each var
        InstanceIDArray ids(numVars);
        DoubleArray calibratedVals(numVars);
        int index = 0;
        for (int i = 0; i < numVars; i++){
            string typeToCalibrate = types[i];
            string name = names[i];
            string fieldName = fieldNames[i];
            string dbType = dbTypes[i];
            if (!arrayIndices[i]){
                throw ModelException(method, "null array index");
            }
            // InstanceDb type
            if (dbType == "SCALAR"
                || (dbType == "VECTOR" && CInt::TYPE->isInstance(arrayIndices[i]))){
                CInt& arrayIndex = dynamic_cast<CInt&>(*arrayIndices[i]);
                ids[i] = InstanceIDSP(
                    new InstanceIDDb(typeToCalibrate,
                                     name, 
                                     fieldName,
                                     false, // use override ?
                                     0.0,   // override
                                     arrayIndex.intValue()));
            }
            else if (dbType == "VECTOR"){
                ids[i] = InstanceIDSP(
                    new InstanceIDDbArray(typeToCalibrate,
                                          name, 
                                          fieldName,
                                          DoubleArraySP(   ))); // override
            }
            else{
                throw ModelException(method, "unexpected db type");
            }
            if (!values[i]){
                throw ModelException(method, "null value");
            }
            ids[i]->readResults(*values[i], calibratedVals, index);
            index += ids[i]->numVariables();
        }
        return pair<InstanceIDArray, DoubleArray>(ids, calibratedVals);
    }   
    catch(exception& e){
        throw ModelException(e, method);
    }
}

// OBJFUNC::IBOOTSTRAPPABLE INTERFACE

typedef Calibrator::ObjFunc::IBootstrappable Calibrator_ObjFunc_IBootstrappable;

void Calibrator_ObjFunc_IBootstrappable::load(CClassSP& clazz){
    REGISTER_INTERFACE(IObject, clazz);
    EXTENDS(IObject);
}

CClassConstSP const Calibrator_ObjFunc_IBootstrappable::TYPE =
CClass::registerInterfaceLoadMethod(
    "Calibrator::ObjFunc::IBootstrappable", 
    typeid(Calibrator_ObjFunc_IBootstrappable), load);

// ObjFunc::IGenericBootstrappable interface
void Calibrator::ObjFunc::IGenericBootstrappable::load(CClassSP& clazz){
	REGISTER_INTERFACE(IObject, clazz);
	EXTENDS(IObject);
}

CClassConstSP const Calibrator::ObjFunc::IGenericBootstrappable::TYPE =
	CClass::registerInterfaceLoadMethod(
    	"Calibrator::ObjFunc::IGenericBootstrappable", 
    	typeid(Calibrator::ObjFunc::IGenericBootstrappable), load);

// INSTANCEID::IBOOTSTRAPPABLE INTERFACE
void Calibrator::InstanceID::IBootstrappable::load(CClassSP& clazz){
    REGISTER_INTERFACE(IObject, clazz);
    EXTENDS(IObject);
}

CClassConstSP const Calibrator::InstanceID::IBootstrappable::TYPE = 
CClass::registerInterfaceLoadMethod(
    "Calibrator::InstanceID::IBootstrappable", 
    typeid(Calibrator::InstanceID::IBootstrappable), load);

// INSTANCEID::IBOOTSTRAPPABLEARRAY INTERFACE
DEFINE_TEMPLATE_TYPE(Calibrator::InstanceID::IBootstrappableArray);

// CALIBRATOR
class Calibrator_ObjectiveFunc: public MFunctionND{
    int                nbVars;
    int                nbFuncs;
    bool               isLeastSq;
    Calibrator::ObjFunc*            objFunc;
    Calibrator::ObjFuncLeastSquare* leastSqFunc;
    IObjectSP          adjGroup;    // need a SP as 'applyAdjustement' expects one
    Calibrator::InstanceIDArray instances;

public:
    virtual void operator()(const CDoubleArray&  x,
                            CDoubleArray&        f) const {
        static const string method = "Calibrator_ObjectiveFunc::operator()";

        try{
            /* Consistency check */
            if (x.size() != nbVars){
                throw ModelException(
                    method,
                    Format::toString("x's size should be %ld; got %ld",
                                     nbVars, x.size()));
            }
            if (f.size() != nbFuncs){
                throw ModelException(
                    method, Format::toString("f's size should be %ld; got %ld",
                                             nbFuncs, f.size()));
            }
            if (SpreadSheetMode::abort()){
                throw ModelException(method,
                                     "calibration was interupted by user");
            }
            Calibrator::InstanceID::applyAdjustment(instances, adjGroup, x);

            if (isLeastSq){
                leastSqFunc->calcValue(f);
            }
            else{
                f[0] = objFunc->calcValue();
            }
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    Calibrator_ObjectiveFunc(const RangeArray&  defRanges,
                             int                nbVars,
                             int                nbFuncs,
                             bool               isLeastSq,
                             Calibrator::ObjFunc* objFunc,
                             Calibrator::ObjFuncLeastSquare* leastSqFunc,
                             IObject*           adjGroup,
                             const Calibrator::InstanceIDArray&  instances):
        MFunctionND(nbVars, nbFuncs, defRanges),
        nbVars(nbVars), nbFuncs(nbFuncs),
        isLeastSq(isLeastSq),
        objFunc(objFunc),
        leastSqFunc(leastSqFunc),
        // should not in any way make a deep copy here
        // otherwise the wrong object will get 'adjusted'
        adjGroup(IObjectSP::attachToRef(adjGroup)),
        instances(instances){}
};


/** Run calibrator */
CResultsSP Calibrator::run(
    ObjFunc&               objFunc,
    const InstanceIDArray& instanceIDs) const {
    static const string method = "Calibrator::run";
    try {
        // start clock
        Timer timer;

        CResultsSP res(new Results());
        /** Get hold of market data names to calibrate */
        IObjectSP adjGroup(objFunc.getAdjustableGroup());
        
        // find number of variables to calibrate
        int nbVars = 0;
        for (int i = 0; i < instanceIDs.size(); i++){
            // initialise instanceIDs before we invoke any methods on them
            instanceIDs[i]->initialise(adjGroup);
            nbVars += instanceIDs[i]->numVariables();
        }

        res->storeScalarGreek(nbVars, 
                              "Calibrator",
                              OutputNameSP(new OutputName("nbVars")));

        if (nbVars > 0){
            /* Check given objective function supports given optimizer 
               and deduce nbFuncs */
            int nbFuncs;
            ObjFuncLeastSquare* leastSqFunc = 0;
            bool isLeastSq = 
                LeastSquareOptimizer::TYPE->isInstance(optimizer.get());
            if (isLeastSq) {
                if (!(leastSqFunc = dynamic_cast<ObjFuncLeastSquare*>(&objFunc))){
                    throw ModelException(method,
                                         objFunc.getClass()->getName() + " does"
                                         " not support Least-Square optimizers");
                }
                nbFuncs = leastSqFunc->getNbFuncs();
            }
            else{
                nbFuncs = 1;
            }

            // Get the range of definition of each shift. First set up empty array
            RangeArray defRanges(nbVars, RangeSP(new InfiniteRange()) /* bit bogus */);
            InstanceID::getRanges(instanceIDs, defRanges);

            /* Get initial guesses */
            DoubleArray xguess(nbVars);
            InstanceID::getInitialValues(instanceIDs, adjGroup, xguess);

            // get id names
            StringArray names(nbVars);
            InstanceID::getNames(instanceIDs, names);

            // Optimizers require open ranges, so
            // 1) convert ranges into open ranges, and 
            // 2) make sure intial guesses are in open ranges 
            int iVar = 0;
            for (; iVar < nbVars; ++iVar){
                (*(defRanges[iVar])).makeOpen();
                if (!Range::variableIsInRange(*(defRanges[iVar]), xguess[iVar])){
                    throw ModelException(method,
                                         "The initial guess for variable " 
                                         + names[iVar] + "\n"
                                         + "must belong to the open range "
                                         + (*(defRanges[iVar])).toString()
                                         + "; got " + Format::toString(xguess[iVar]));
                }
            }

            Calibrator_ObjectiveFunc targetFunc(defRanges,
                                                nbVars,
                                                nbFuncs,
                                                isLeastSq,
                                                &objFunc,
                                                leastSqFunc,
                                                adjGroup.get(),
                                                instanceIDs);

            /* Call optimizer */
            DoubleArray x(nbVars);
            optimizer->minimize(targetFunc,
                                xguess,
                                names,
                                x);

            // export calibrated values
            res->storeGreek(InstanceID::writeResults(instanceIDs, x), 
                            "Calibrator", 
                            OutputNameSP(new OutputName("calibratedVars")));

            // export additional stats info if least square
            LeastSquareOptimizer* leastSqOpt = dynamic_cast<LeastSquareOptimizer*>(optimizer.get());
            bool leastSquareStatus;
            if (leastSquareStatus = (leastSqOpt && leastSqOpt->statsIsAvailable())){
                res->storeScalarGreek(leastSqOpt->getVariance(), 
                                      "Calibrator",
                                      OutputNameSP(new OutputName("leastSquareStats", "variance")));
                res->storeGreek(IObjectSP(leastSqOpt->getCovarianceMatrix().clone()),
                                "Calibrator",
                                OutputNameSP(new OutputName("leastSquareStats", "covarianceMatrix")));
                res->storeGreek(IObjectSP(leastSqOpt->getUpperLimits().clone()),
                                "Calibrator",
                                OutputNameSP(new OutputName("leastSquareStats", "upperLimits")));
                res->storeGreek(IObjectSP(leastSqOpt->getCorrelationMatrix().clone()),
                                "Calibrator",
                                OutputNameSP(new OutputName("leastSquareStats", "correlationMatrix")));
                res->storeGreek(IObjectSP(leastSqOpt->getLowerLimits().clone()),
                                "Calibrator",
                                OutputNameSP(new OutputName("leastSquareStats", "lowerLimits")));
                res->storeGreek(IObjectSP(leastSqOpt->getStandardErrors().clone()),
                                "Calibrator",
                                OutputNameSP(new OutputName("leastSquareStats", "standardErrors")));
            }
            res->storeGreek(IObjectSP(CBool::create(leastSquareStatus)), 
                            "Calibrator",
                            OutputNameSP(new OutputName("leastSquareStats", "status")));
        }

        // calculate obj func and export obj func value
        double f = objFunc.calcValue();
        res->storeScalarGreek(f,
                              "Calibrator",
                              OutputNameSP(new OutputName("calibratedObjFuncValue")));

        // export calc time
        double calcTime = timer.calcTime();
        res->storeScalarGreek(calcTime, 
                              "Calibrator",
                              OutputNameSP(new OutputName("CALC_TIME")));   // CALC_TIME will be ignored by testcmp.pl

        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

//overloaded run 
CResultsSP Calibrator::run(ObjFunc&  objFunc,
const InstanceIDArray& instanceIDs,
DoubleArray xguess,
double objFuncVal) const {
    static const string method = "Calibrator::run";
    try {
        // start clock
        Timer timer;

        CResultsSP res(new Results());
        /** Get hold of market data names to calibrate */
        IObjectSP adjGroup(objFunc.getAdjustableGroup());

        // find number of variables to calibrate
        int nbVars = 0;
        for (int i = 0; i < instanceIDs.size(); i++){
            // initialise instanceIDs before we invoke any methods on them
            instanceIDs[i]->initialise(adjGroup);
            nbVars += instanceIDs[i]->numVariables();
        }

        res->storeScalarGreek(nbVars, 
            "Calibrator",
            OutputNameSP(new OutputName("nbVars")));

        if (nbVars > 0){
            /* Check given objective function supports given optimizer 
            and deduce nbFuncs */
            int nbFuncs;
            ObjFuncLeastSquare* leastSqFunc = 0;
            bool isLeastSq = 
                LeastSquareOptimizer::TYPE->isInstance(optimizer.get());
            if (isLeastSq) {
                if (!(leastSqFunc = dynamic_cast<ObjFuncLeastSquare*>(&objFunc))){
                    throw ModelException(method,
                        objFunc.getClass()->getName() + " does"
                        " not support Least-Square optimizers");
                }
                nbFuncs = leastSqFunc->getNbFuncs();
            }
            else{
                nbFuncs = 1;
            }

            // get id names
            StringArray names(nbVars);
            InstanceID::getNames(instanceIDs, names);

            // export calibrated values
            res->storeGreek(InstanceID::writeResults(instanceIDs, xguess), 
                "Calibrator", 
                OutputNameSP(new OutputName("calibratedVars")));

            // export additional stats info if least square
            LeastSquareOptimizer* leastSqOpt = dynamic_cast<LeastSquareOptimizer*>(optimizer.get());
            bool leastSquareStatus;
            if (leastSquareStatus = (leastSqOpt && leastSqOpt->statsIsAvailable())){
                res->storeScalarGreek(leastSqOpt->getVariance(), 
                    "Calibrator",
                    OutputNameSP(new OutputName("leastSquareStats", "variance")));
                res->storeGreek(IObjectSP(leastSqOpt->getCovarianceMatrix().clone()),
                    "Calibrator",
                    OutputNameSP(new OutputName("leastSquareStats", "covarianceMatrix")));
                res->storeGreek(IObjectSP(leastSqOpt->getUpperLimits(). clone()),
                    "Calibrator",
                    OutputNameSP(new OutputName("leastSquareStats", "upperLimits")));
                res->storeGreek(IObjectSP(leastSqOpt->getCorrelationMatrix().clone()),
                    "Calibrator",
                    OutputNameSP(new OutputName("leastSquareStats", "correlationMatrix")));
                res->storeGreek(IObjectSP(leastSqOpt->getLowerLimits().clone()),
                    "Calibrator",
                    OutputNameSP(new OutputName("leastSquareStats", "lowerLimits")));
                res->storeGreek(IObjectSP(leastSqOpt->getStandardErrors().clone()),
                    "Calibrator",
                    OutputNameSP(new OutputName("leastSquareStats", "standardErrors")));
            }
            res->storeGreek(IObjectSP(CBool::create(leastSquareStatus)), 
                "Calibrator",
                OutputNameSP(new OutputName("leastSquareStats", "status")));
        }

        // calculate obj func and export obj func value
        double f = objFuncVal;
        res->storeScalarGreek(f,
            "Calibrator",
            OutputNameSP(new OutputName("calibratedObjFuncValue")));

        // export calc time
        double calcTime = timer.calcTime();
        res->storeScalarGreek(calcTime, 
            "Calibrator",
            OutputNameSP(new OutputName("CALC_TIME")));   // CALC_TIME will be ignored by testcmp.pl

        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/** tiny class used for call back mechanism */
class Calibrator::NameCollector: virtual public ObjectIteration::IActionConst{
public:
    set<string> names; // use set to avoid duplicates
    /* method invoked by recurse routine each time we hit an object
       of the right type */
    bool invoke(const ObjectIteration::State& state, IObjectConstSP obj){
        const Calibrator::IAdjustable& adjustable = 
            dynamic_cast<const Calibrator::IAdjustable&>(*obj);
        names.insert(adjustable.getName());
        return true;
    }
};
//// Utility: Returns all distinct names of objects of specified type 
StringArray Calibrator::getNames(CClassConstSP   clazz,
                                 IObjectConstSP  adjGroup){
    NameCollector collector;
    ObjectIteration iteration(clazz);
    // but skip transient fields
    iteration.setSkipTransient(true);
    iteration.recurse(collector, adjGroup); // go
    StringArray namesToCalibrate;
    namesToCalibrate.reserve(collector.names.size());
    const set<string>& names = collector.names;
    for (set<string>::const_iterator iter = names.begin();
         iter != names.end(); ++iter){
        namesToCalibrate.push_back(*iter);
    }
    return namesToCalibrate;
}

Calibrator::Calibrator(const OptimizerNDSP& optimizer):
    CObject(TYPE), optimizer(optimizer){}

/** for reflection */
Calibrator::Calibrator(): CObject(TYPE){}

class CalibratorHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Calibrator, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(optimizer, "Optimizer");
    }

    static IObject* defaultCtor(){
        return new Calibrator();
    }
};

CClassConstSP const Calibrator::TYPE = CClass::registerClassLoadMethod(
    "Calibrator", typeid(Calibrator), CalibratorHelper::load);

// ADDINS
/* Addin class */
class CalibratorAddin : public CObject {
    OptimizerNDSP              optimizer;
    Calibrator::ObjFuncSP      objFunc;
    // general case is an array of these next 4 fields
    string                     typeToCalibrate;
    StringArray                fields; // optional (default: all)
    StringArray                namesToCalibrate; // optional (default: all)
    DoubleMatrix               guesses; // optional (default: current values)
    CMarketDataSP              market;
public:
    static CClassConstSP const TYPE;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CalibratorAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(optimizer, "N-dimensional optimizer");
        FIELD(objFunc, "Objective function");
        FIELD(typeToCalibrate, "Type of object to calibrate");
        FIELD(fields, "Fields to calibrate");
        FIELD_MAKE_OPTIONAL(fields);
        FIELD(namesToCalibrate, "Market name of objects to calibrate");
        FIELD_MAKE_OPTIONAL(namesToCalibrate);
        FIELD(guesses, "Initial values for calibrator");
        FIELD_MAKE_OPTIONAL(guesses);
        FIELD(market, "Market");
        FIELD_MAKE_OPTIONAL(market);

        // registration for addin function
        Addin::registerClassObjectMethod(
            "CALIBRATOR",
            Addin::RISK,
            "Calibrates an array of parameter"
            "shifts against an objective function",
            CalibratorAddin::TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)go);        
    }

//private: // shuts the compiler up
    CalibratorAddin(const CalibratorAddin& rhs);
    CalibratorAddin& operator=(const CalibratorAddin& rhs);

    static IObjectSP go(CalibratorAddin* params){
        return params->run();
    }
    IObjectSP run(){
        static const string method("CalibratorAddin::run");
        Calibrator calibrator(optimizer);
        // if there's a market, give objFunc a chance to get it
        // otherwise, trust (!) that objFunc already has all the market
        // data it needs
        if (market.get()){
            objFunc->getMarket(market.get());
        }
        // further validation irrespective of getMarket having been
        // called or not
        objFunc->validate();
        // now massage input data into form expected by calibrator
        CClassConstSP clazz = CClass::forName(typeToCalibrate);
        // do validation 
        bool useOverride = !guesses.empty();
        if (useOverride){
            if (namesToCalibrate.empty()){
                throw ModelException(method, "If 'guesses' supplied then "
                                     "names of the corresponding object(s) "
                                     "must be supplied too");
            }
            if (fields.empty()){
                throw ModelException(method, "If 'guesses' supplied then "
                                     "names of the corresponding fields(s) "
                                     "must be supplied too");
            }
            int numCols = guesses.numCols();
            if (numCols != namesToCalibrate.size()){
                throw ModelException(method, "Number of columns of 'guesses' "
                                     "supplied must match number of names");
            }
            int numRows = guesses.numRows();
            if (numRows != fields.size()){
                throw ModelException(method, "Number of rows of 'guesses' "
                                     "supplied must match number of fields");
            }
        }
        // create the InstanceIDs. Start by getting names of the instances
        if (namesToCalibrate.empty()){
            namesToCalibrate = Calibrator::getNames(clazz,
                                                    objFunc->getAdjustableGroup());
        }
        // validate that any supplied fields are valid
        CFieldArray realFields;
        if (!fields.empty()){
            // to do: move block somewhere central ?
            realFields.reserve(fields.size());
            for (int i = 0; i < fields.size(); i++){
                const string& field = fields[i];
                CClassConstSP c = clazz;
                CFieldConstSP realField;
                do {
                    realField = c->hasDeclaredField(field);
                } while (!realField && (c = c->getSuperClass()) != 0);
                if (!realField){
                    throw ModelException(method, field + " not found in "+
                                         typeToCalibrate);
                }
                realFields.push_back(realField);
            }
        } else {
            realFields = Calibrator::IAdjustable::getFields(clazz);
        }
        int numIDs = realFields.size() * namesToCalibrate.size();
        Calibrator::InstanceIDArray ids;
        ids.reserve(numIDs);
        for (int i = 0; i < namesToCalibrate.size(); i++){
            for (unsigned int j = 0; j < realFields.size(); j++){
                // currently all doubles at the moment
                ids.push_back(Calibrator::InstanceIDSP(
                                  new Calibrator::InstanceIDDb(typeToCalibrate, 
                                                               namesToCalibrate[i],
                                                               realFields[j]->getName(),
                                                               useOverride,
                                                               useOverride? 
                                                               guesses[i][j]: 0.0,
                                                               0)));
            }
        }
        // run the calibrator
        return calibrator.run(*objFunc, ids);
    }

    CalibratorAddin(): CObject(TYPE){}

    static IObject* defaultCtor(){
        return new CalibratorAddin();
    }

};

CClassConstSP const CalibratorAddin::TYPE = CClass::registerClassLoadMethod(
    "CalibratorAddin", typeid(CalibratorAddin), CalibratorAddin::load);



class AdjustableFieldsAddin : public CObject {
    string  objType;

public:
    static CClassConstSP const TYPE;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(AdjustableFieldsAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(objType, "Type of object");

        // registration for addin function
        Addin::registerClassObjectMethod(
            "ADJUSTABLE_FIELDS_GET",
            Addin::RISK,
            "Returns all adjustable fields for a given object type",
            AdjustableFieldsAddin::TYPE,
            false,
            Addin::expandSimple,
            (Addin::ObjMethod*)go);        
    }

//private: // shuts the compiler up
    AdjustableFieldsAddin(const AdjustableFieldsAddin& rhs);
    AdjustableFieldsAddin& operator=(const AdjustableFieldsAddin& rhs);

    static IObjectSP go(AdjustableFieldsAddin* params){
        CClassConstSP clazz(CClass::forName(params->objType));
        StringArraySP fieldNames;
        StringArraySP scalarTypes;
        if (Calibrator::IAdjustable::TYPE->
            isAssignableFrom(clazz)){
            CFieldArray fields(Calibrator::IAdjustable::getFields(clazz));
            int nbFields = fields.size();
            fieldNames = StringArraySP(new StringArray(nbFields));
            scalarTypes = StringArraySP(new StringArray(nbFields));
            int iField = 0;
            for (; iField < nbFields; ++iField){
                (*fieldNames)[iField] = fields[iField]->getName();
                (*scalarTypes)[iField] = (fields[iField]->typeIsArray() ? "VECTOR" : "SCALAR");
            }
        }
        ObjectArraySP res(new ObjectArray(2));
        (*res)[0] = fieldNames;
        (*res)[1] = scalarTypes;
        return res;
    }

    AdjustableFieldsAddin(): CObject(TYPE){}

    static IObject* defaultCtor(){
        return new AdjustableFieldsAddin();
    }
};

CClassConstSP const AdjustableFieldsAddin::TYPE = CClass::registerClassLoadMethod(
    "AdjustableFieldsAddin", typeid(AdjustableFieldsAddin), load);

/* Addin class */
class AdjustableFieldRangeAddin : public CObject {
    string  objType;
    string  fieldName;

public:
    static CClassConstSP const TYPE;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(AdjustableFieldRangeAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(objType, "Type of object");
        FIELD(fieldName, "Name of field");

        // registration for addin function
        Addin::registerClassStringMethod(
            "ADJUSTABLE_FIELD_RANGE_GET",
            Addin::RISK,
            "Returns the range of definition of a given field",
            AdjustableFieldRangeAddin::TYPE,
            (Addin::StringMethod*)go);        
    }

    //private: // shuts the compiler up
    AdjustableFieldRangeAddin(const AdjustableFieldRangeAddin& rhs);
    AdjustableFieldRangeAddin& operator=(const AdjustableFieldRangeAddin& rhs);

    static string go(AdjustableFieldRangeAddin* params){
        CClassConstSP clazz(CClass::forName(params->objType));
        if (!Calibrator::IAdjustable::TYPE->isAssignableFrom(clazz)){
            return "";
        }
        CClassConstSP c = clazz;
        CFieldConstSP realField;
        do {
            realField = c->hasDeclaredField(params->fieldName);
        } while (!realField && (c = c->getSuperClass()) != 0);
        if (!realField){
            return "";
        }
        return Calibrator::IAdjustable::getRange(realField).toString();
    }

    AdjustableFieldRangeAddin(): CObject(TYPE){}

    static IObject* defaultCtor(){
        return new AdjustableFieldRangeAddin();
    }
};

CClassConstSP const AdjustableFieldRangeAddin::TYPE = 
CClass::registerClassLoadMethod(
    "AdjustableFieldRangeAddin", typeid(AdjustableFieldRangeAddin), load);

// Helper func
template <class ObjectType>
void ObjFuncHelper_clearOffZeros(DoubleArray&                             first,
                                 array<smartPtr<ObjectType>, ObjectType>& second){
    typedef smartPtr<ObjectType> ObjectTypeSP;
    // remove those 'seconds' for which 'first' == 0.0
    vector<double>::iterator first_it(first.begin());
    typename vector<ObjectTypeSP>::iterator second_it(second.begin());
    while (first_it != first.end()){
        // locate zero 'firsts'
        typename vector<ObjectTypeSP>::iterator second_it_start(second_it);
        vector<double>::iterator first_it_start(first_it);
        bool remove = false;
        while (first_it != first.end() && Maths::isZero(*first_it)){
            remove = true;
            ++first_it;
            ++second_it;
        }
        // remove zero 'firsts', if necessary
        if (remove){
            first_it = first.erase(first_it_start, first_it);
            second_it = second.erase(second_it_start, second_it);
        }
        // otherwise, next
        else{
            ++first_it;
            ++second_it;
        }
    }
}


// PENALTY FUNC
typedef MarketWrapper<MarketObject> GenericMarketObjectWrapper;
// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(GenericMarketObjectWrapper);

void Calibrator::PenaltyFunc::validatePop2Object(){
    static const string method("Calibrator::PenaltyFunc::validatePop2Object");
    try{
        // at construction time, market wrappers get converted to strings
        // therefore if object is a string, interpret it as a market wrapper
        // of the same name
        CClassConstSP clazz(object->getClass());
        if (CString::TYPE->isAssignableFrom(clazz)){
            CString& name = dynamic_cast<CString&>(*object);
            object = IObjectSP(new GenericMarketObjectWrapper(name.stringValue()));
        }
        // if the object is an empty market object, make sure
        // the objectType is non-empty (as, then, there will be some market 
        // data to get)
        clazz = object->getClass();
        if (MarketObjectWrapper::TYPE->isAssignableFrom(clazz)){
            MarketObjectWrapper& mow = dynamic_cast<MarketObjectWrapper&>(*object);
            if (!mow.getMO() && objectType.empty()){
                throw ModelException(method,
                                     "object is an empty market wrapper ('" 
                                     + clazz->getName() + "').\n"
                                     + "objectType must be provided");
            }
        }

        if (ids.size() != weights.size()){
            throw ModelException(method,
                                 "Nber of ids and nber of weights should be the same;\ngot "
                                 + Format::toString(ids.size()) + " and "
                                 + Format::toString(weights.size()) + ", respectively");
        }

        Maths::checkNonNegative(weights, "weights");

        // remove those vars for which weight == 0.0
        ObjFuncHelper_clearOffZeros(weights, ids);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void Calibrator::PenaltyFunc::getMarket(MarketData* market){
    static const string method = "Calibrator::PenaltyFunc::getMarket";
    try{
        CClassConstSP clazz(object->getClass());
        // if object is a market object wrapper (eg VolBaseWrapper)
        // need to get market
        if (MarketObjectWrapper::TYPE->isAssignableFrom(clazz)){
            // populate the wrapper with object of type objectType (eg VolSVJ)
            CClassConstSP objClazz(CClass::forName(objectType));
            NonPricingModel dummyModel;
            MarketObjectWrapper& mow = dynamic_cast<MarketObjectWrapper&>(*object);
            mow.getData(&dummyModel, market, objClazz);
        }
        // otherwise, no data to get
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void Calibrator::PenaltyFunc::validate(){
    static const string method = "Calibrator::PenaltyFunc::validate";
    try{
        // find number of variables to calibrate
        nbVars = 0;
        int i = 0;
        for (; i < ids.size(); i++){
            // initialise instanceIDs before we invoke any methods on them
            ids[i]->initialise(object);
            nbVars += ids[i]->numVariables();
        }
        if (nbVars > 0){
            // get initial values
            initvals.resize(nbVars);
            currvals.resize(nbVars);
            InstanceID::getInitialValues(ids, object, initvals);
            // create 'actual size' weight array
            usedweights.resize(nbVars);
            int j = 0;
            for (i = 0; i < ids.size(); i++){
                double weight = sqrt(weights[i]);
                int kmax = ids[i]->numVariables();
                for (int k = 0; k < kmax; ++k, ++j){
                    usedweights[j] = weight;
                }
            }
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

IObjectSP Calibrator::PenaltyFunc::getAdjustableGroup(){
    return object;
}

int Calibrator::PenaltyFunc::getNbFuncs() const{
    return nbVars;
}

void Calibrator::PenaltyFunc::calcValue(CDoubleArray& funcvals) const{
    try{        
        InstanceID::getValues(ids, object, currvals);
        for(int iVar = 0; iVar < nbVars; ++iVar){
            funcvals[iVar] = usedweights[iVar] * (currvals[iVar] - initvals[iVar]);
        }
    }
    catch(exception& e){
        throw ModelException(e, "Calibrator::PenaltyFunc::calcValue");
    }
}

// for reflection
Calibrator::PenaltyFunc::PenaltyFunc(): 
    Calibrator::ObjFuncLeastSquare(TYPE),
    nbVars(0){}

class Calibrator_PenaltyFuncHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Calibrator::PenaltyFunc, clazz);
        SUPERCLASS(Calibrator::ObjFuncLeastSquare);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(ids, "Array of target variables (and their target values)");
        FIELD(weights, "Array of weights");
        FIELD(object, "Object containing all the target variables");
        FIELD(objectType, "Type of object");
        FIELD_MAKE_OPTIONAL(objectType);
        FIELD(nbVars, "");
        FIELD_MAKE_TRANSIENT(nbVars);
        FIELD(initvals, "");
        FIELD_MAKE_TRANSIENT(initvals);
        FIELD(currvals, "");
        FIELD_MAKE_TRANSIENT(currvals);
        FIELD(usedweights, "");
        FIELD_MAKE_TRANSIENT(usedweights);
    }

    static IObject* defaultCtor(){
        return new Calibrator::PenaltyFunc();
    }
};

CClassConstSP const Calibrator::PenaltyFunc::TYPE = CClass::registerClassLoadMethod(
    "Calibrator::PenaltyFunc", typeid(Calibrator::PenaltyFunc), Calibrator_PenaltyFuncHelper::load);

// OBJFUNC COMBO
void Calibrator::ObjFuncCombo::validatePop2Object(){
    static const string method("Calibrator::ObjFuncCombo::validatePop2Object");
    try{
        if (weights.size() == 0){
            int nbWeights = objFuncs.size(); 
            weights.resize(nbWeights);
            for (int iWeight = 0; iWeight < nbWeights; ++iWeight){
                weights[iWeight] = 1.0;
            }
        }
        if (objFuncs.size() != weights.size()){
            throw ModelException(method,
                                 "Nber of objFuncs and nber of weights should be the same;\ngot "
                                 + Format::toString(objFuncs.size()) + " and "
                                 + Format::toString(weights.size()) + ", respectively");
        }

        Maths::checkNonNegative(weights, "weights");

        // remove those funcs for which weight == 0.0
        ObjFuncHelper_clearOffZeros(weights, objFuncs);

        nbFuncs = objFuncs.size();
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void Calibrator::ObjFuncCombo::getMarket(MarketData* market){
    static const string method = "Calibrator::ObjFuncCombo::getMarket";
    try{
        for (int iFunc = 0; iFunc < nbFuncs; ++iFunc){
            objFuncs[iFunc]->getMarket(market);
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void Calibrator::ObjFuncCombo::validate(){
    static const string method = "Calibrator::ObjFuncCombo::validate";
    try{
        for (int iFunc = 0; iFunc < nbFuncs; ++iFunc){
            objFuncs[iFunc]->validate();
            adjGroupArray[iFunc] = objFuncs[iFunc]->getAdjustableGroup();
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

IObjectSP Calibrator::ObjFuncCombo::getAdjustableGroup(){
    return IObjectSP::attachToRef(&adjGroupArray);
}

double Calibrator::ObjFuncCombo::calcValue() const{
    try{        
        double val = 0.0;
        for (int iFunc = 0; iFunc < nbFuncs; ++iFunc){
            val += weights[iFunc] * objFuncs[iFunc]->calcValue();
        }
        return val;
    }
    catch(exception& e){
        throw ModelException(e, "Calibrator::ObjFuncCombo::calcValue");
    }
}

/** Makes additional validations for calibration with bootstrapping */
void Calibrator::ObjFuncCombo::validate(const InstanceID::IBootstrappableArray& ids) const{
    static const string method = "Calibrator::ObjFuncCombo::validate";
    try{        
        for (int iFunc = 0; iFunc < nbFuncs; ++iFunc){
            if (!Calibrator::ObjFunc::IBootstrappable::TYPE->isInstance(objFuncs[iFunc])){
                throw ModelException(method,
                                     "The "
                                     + Format::toString(iFunc + 1)
                                     + "-th function is not bootstrappable");
            }
            Calibrator::ObjFunc::IBootstrappableSP bootstrapFunc(
                Calibrator::ObjFunc::IBootstrappableSP::dynamicCast(objFuncs[iFunc]));
            bootstrapFunc->validate(ids);
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/** Updates the instrument for calibration with bootstrapping
    before each calibrator run */
void Calibrator::ObjFuncCombo::update(int idxMat){
    static const string method = "Calibrator::ObjFuncCombo::validate";
    try{        
        for (int iFunc = 0; iFunc < nbFuncs; ++iFunc){
            if (!Calibrator::ObjFunc::IBootstrappable::TYPE->isInstance(objFuncs[iFunc])){
                throw ModelException(method,
                                     "The "
                                     + Format::toString(iFunc + 1)
                                     + "-th function is not bootstrappable");
            }
            Calibrator::ObjFunc::IBootstrappableSP bootstrapFunc(
                Calibrator::ObjFunc::IBootstrappableSP::dynamicCast(objFuncs[iFunc]));
            bootstrapFunc->update(idxMat);
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/** Reset the instrument for calibration with bootstrapping
    after each calibrator run */
void Calibrator::ObjFuncCombo::reset(){
    static const string method = "Calibrator::ObjFuncCombo::validate";
    try{        
        for (int iFunc = 0; iFunc < nbFuncs; ++iFunc){
            if (!Calibrator::ObjFunc::IBootstrappable::TYPE->isInstance(objFuncs[iFunc])){
                throw ModelException(method,
                                     "The "
                                     + Format::toString(iFunc + 1)
                                     + "-th function is not bootstrappable");
            }
            Calibrator::ObjFunc::IBootstrappableSP bootstrapFunc(
                Calibrator::ObjFunc::IBootstrappableSP::dynamicCast(objFuncs[iFunc]));
            bootstrapFunc->reset();
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

// for reflection
Calibrator::ObjFuncCombo::ObjFuncCombo(): 
    Calibrator::ObjFunc(TYPE),
    nbFuncs(0){}

/** Invoked when Class is 'loaded' */
void Calibrator::ObjFuncCombo::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(Calibrator::ObjFuncCombo, clazz);
    SUPERCLASS(Calibrator::ObjFunc); 
    IMPLEMENTS(Calibrator::ObjFunc::IBootstrappable);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(objFuncs, "Array of objective functions");
    FIELD(weights, "Array of weights");
    FIELD_MAKE_OPTIONAL(weights);
    FIELD(nbFuncs, "");
    FIELD_MAKE_TRANSIENT(nbFuncs);
    FIELD(adjGroupArray, "");
    FIELD_MAKE_TRANSIENT(adjGroupArray);
}

IObject* Calibrator::ObjFuncCombo::defaultCtor(){
    return new Calibrator::ObjFuncCombo();
}

CClassConstSP const Calibrator::ObjFuncCombo::TYPE = CClass::registerClassLoadMethod(
    "Calibrator::ObjFuncCombo", typeid(Calibrator::ObjFuncCombo), load);

// OBJFUNC LEAST SQUARE COMBO
void Calibrator::ObjFuncLeastSquareCombo::validatePop2Object(){
    static const string method("Calibrator::ObjFuncLeastSquareCombo::validatePop2Object");
    try{
        if (weights.size() == 0){
            int nbWeights = objFuncs.size(); 
            weights.resize(nbWeights);
            for (int iWeight = 0; iWeight < nbWeights; ++iWeight){
                weights[iWeight] = 1.0;
            }
        }
        if (objFuncs.size() != weights.size()){
            throw ModelException(method,
                                 "Nber of objFuncs and nber of weights should be the same;\ngot "
                                 + Format::toString(objFuncs.size()) + " and "
                                 + Format::toString(weights.size()) + ", respectively");
        }

        Maths::checkNonNegative(weights, "weights");

        // remove those funcs for which weight == 0.0
        ObjFuncHelper_clearOffZeros(weights, objFuncs);

        // nb funcs
        nbFuncs = objFuncs.size();

        // needs to sqrt the weights, as least square optimizers take
        // the square of each func (and therefore the weight gets squared too)
        for (int iFunc = 0; iFunc < nbFuncs; ++iFunc){
            weights[iFunc] = sqrt(weights[iFunc]);
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void Calibrator::ObjFuncLeastSquareCombo::getMarket(MarketData* market){
    static const string method = "Calibrator::ObjFuncLeastSquareCombo::getMarket";
    try{
        for (int iFunc = 0; iFunc < nbFuncs; ++iFunc){
            objFuncs[iFunc]->getMarket(market);
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void Calibrator::ObjFuncLeastSquareCombo::validate(){
    static const string method = "Calibrator::ObjFuncLeastSquareCombo::validate";
    try{
        vals.resize(nbFuncs);
        adjGroupArray.resize(nbFuncs);
        for (int iFunc = 0; iFunc < nbFuncs; ++iFunc){
            // first call the 'validate' method,
            // as 'validate' also initialize transient fields
            objFuncs[iFunc]->validate();
            // then only you can invoke other methods
            int currNbFuncs = objFuncs[iFunc]->getNbFuncs();
            vals[iFunc].resize(currNbFuncs);
            totalNbFuncs += currNbFuncs;
            adjGroupArray[iFunc] = objFuncs[iFunc]->getAdjustableGroup();
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

IObjectSP Calibrator::ObjFuncLeastSquareCombo::getAdjustableGroup(){
    return IObjectSP::attachToRef(&adjGroupArray);
}

int Calibrator::ObjFuncLeastSquareCombo::getNbFuncs() const{
    return totalNbFuncs;
}

Calibrator::ObjFuncLeastSquareArray Calibrator::ObjFuncLeastSquareCombo::getObjFuncArray(){
    return objFuncs;
}

void Calibrator::ObjFuncLeastSquareCombo::calcValue(CDoubleArray& funcvals) const{
    try{        
        int index = 0;
        for (int iFunc = 0; iFunc < nbFuncs; ++iFunc){
            objFuncs[iFunc]->calcValue(vals[iFunc]);
            int currNbFuncs = objFuncs[iFunc]->getNbFuncs();
            for (int iCurrFunc = 0; iCurrFunc < currNbFuncs; ++iCurrFunc){
                funcvals[index + iCurrFunc] = weights[iFunc] * vals[iFunc][iCurrFunc];
            }
            index += currNbFuncs;
        }
    }
    catch(exception& e){
        throw ModelException(e, "Calibrator::ObjFuncLeastSquareCombo::calcValue");
    }
}

/** Makes additional validations for calibration with bootstrapping */
void Calibrator::ObjFuncLeastSquareCombo::validate(const InstanceID::IBootstrappableArray& ids) const{
    static const string method = "Calibrator::ObjFuncLeastSquareCombo::validate";
    try{        
        for (int iFunc = 0; iFunc < nbFuncs; ++iFunc){
            if (!Calibrator::ObjFunc::IBootstrappable::TYPE->isInstance(objFuncs[iFunc])){
                throw ModelException(method,
                                     "The "
                                     + Format::toString(iFunc + 1)
                                     + "-th function is not bootstrappable");
            }
            Calibrator::ObjFunc::IBootstrappableSP bootstrapFunc(
                Calibrator::ObjFunc::IBootstrappableSP::dynamicCast(objFuncs[iFunc]));
            bootstrapFunc->validate(ids);
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/** Updates the instrument for calibration with bootstrapping
    before each calibrator run */
void Calibrator::ObjFuncLeastSquareCombo::update(int idxMat){
    static const string method = "Calibrator::ObjFuncLeastSquareCombo::validate";
    try{        
        for (int iFunc = 0; iFunc < nbFuncs; ++iFunc){
            if (!Calibrator::ObjFunc::IBootstrappable::TYPE->isInstance(objFuncs[iFunc])){
                throw ModelException(method,
                                     "The "
                                     + Format::toString(iFunc + 1)
                                     + "-th function is not bootstrappable");
            }
            Calibrator::ObjFunc::IBootstrappableSP bootstrapFunc(
                Calibrator::ObjFunc::IBootstrappableSP::dynamicCast(objFuncs[iFunc]));
            bootstrapFunc->update(idxMat);
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/** Reset the instrument for calibration with bootstrapping
    before each calibrator run */
void Calibrator::ObjFuncLeastSquareCombo::reset(){
    static const string method = "Calibrator::ObjFuncLeastSquareCombo::validate";
    try{        
        for (int iFunc = 0; iFunc < nbFuncs; ++iFunc){
            if (!Calibrator::ObjFunc::IBootstrappable::TYPE->isInstance(objFuncs[iFunc])){
                throw ModelException(method,
                                     "The "
                                     + Format::toString(iFunc + 1)
                                     + "-th function is not bootstrappable");
            }
            Calibrator::ObjFunc::IBootstrappableSP bootstrapFunc(
                Calibrator::ObjFunc::IBootstrappableSP::dynamicCast(objFuncs[iFunc]));
            bootstrapFunc->reset();
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

// for reflection
Calibrator::ObjFuncLeastSquareCombo::ObjFuncLeastSquareCombo(): 
    Calibrator::ObjFuncLeastSquare(TYPE),
    nbFuncs(0),
    totalNbFuncs(0){}

/** Invoked when Class is 'loaded' */
void Calibrator::ObjFuncLeastSquareCombo::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(Calibrator::ObjFuncLeastSquareCombo, clazz);
    SUPERCLASS(Calibrator::ObjFuncLeastSquare);
    IMPLEMENTS(Calibrator::ObjFunc::IBootstrappable);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(objFuncs, "Array of objective functions");
    FIELD(weights, "Array of weights");
    FIELD_MAKE_OPTIONAL(weights);
    FIELD(nbFuncs, "");
    FIELD_MAKE_TRANSIENT(nbFuncs);
    FIELD(totalNbFuncs, "");
    FIELD_MAKE_TRANSIENT(totalNbFuncs);
    FIELD(vals, "");
    FIELD_MAKE_TRANSIENT(vals);
    FIELD(adjGroupArray, "");
    FIELD_MAKE_TRANSIENT(adjGroupArray);
}

IObject* Calibrator::ObjFuncLeastSquareCombo::defaultCtor(){
    return new Calibrator::ObjFuncLeastSquareCombo();
}

CClassConstSP const Calibrator::ObjFuncLeastSquareCombo::TYPE = CClass::registerClassLoadMethod(
    "Calibrator::ObjFuncLeastSquareCombo", typeid(Calibrator::ObjFuncLeastSquareCombo), load);

/* external symbol to allow class to be forced to be linked in */
bool CalibratorLinkIn(){
    return (Calibrator::TYPE != 0 && 
            CalibratorAddin::TYPE != 0);
}

/**
 * Helper class responsible for retrieving the 
 * object corresponding to a given InstanceID
 * This can be used to retrive the object to be bootstrapped
 * */
class FindBootstrappableAction:
    public virtual ObjectIteration::IAction
{
public:
    
    /** Constructor */
    FindBootstrappableAction(string instanceName):
        instanceName(instanceName), bootstrappable(0) {}

    /*
     * Method invoked by recurse routine each time we hit an object
     * of the right type
     * */
    virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj) {
        const Calibrator::IAdjustable& adjustable = 
            dynamic_cast<const Calibrator::IAdjustable&>(*obj);
        const string& name = adjustable.getName();
        if (name == instanceName) 
        {
             bootstrappable = obj;
            // stop recursion
            return false;
        } else 
        {
            // continue recursion
            return true;
        }
    }
    
    /** Returns the bootstrappable */
    IObjectConstSP getBootstrappable() const {
        if (bootstrappable.get() != 0) {
            return bootstrappable;
        } else {
            throw ModelException(
                "FindBootstrapperAction::getBootstrappable",
                "Internal error: " + instanceName + " not found");
        }
    }

private:
    /** Name of the InstanceID to */
    string instanceName;
    
    /** The object we are looking for */
    IObjectConstSP bootstrappable;
};

/** 
* Returns the object corresponding to typeToCalibrate for this InstanceID
*/
IObjectConstSP Calibrator::InstanceID::getObject(const IObjectSP& adjGroup) const
{
    ObjectIteration iteration(classToCalibrate);
    FindBootstrappableAction findBootstrappable(name);
    // skip transient fields
    iteration.setSkipTransient(true);
    iteration.recurse(findBootstrappable, adjGroup); // go
    return findBootstrappable.getBootstrappable();

}

/** access to field */
string Calibrator::InstanceID::getFieldName() const
{
    return fieldName;
}

// create InstanceID for the initial guess
Calibrator::InstanceIDSP Calibrator::InstanceID::dbArray(
	    const string&         typeToCalibrate,
        const string&         name, 
		const string&         field) {
	DoubleArraySP noOverrideThankYou;
	return InstanceIDSP(new Calibrator::InstanceIDDbArray(
		typeToCalibrate, name, field, noOverrideThankYou));
}

// get expiries for initial guess
ExpiryArraySP Calibrator::InstanceID::getExpiries(InstanceIDSP ids, ObjFuncSP objFunc){

	Calibrator::InstanceIDDbArraySP iddbArray = Calibrator::InstanceIDDbArraySP::dynamicCast(ids);
	IObjectSP adjGroup(objFunc->getAdjustableGroup());
	iddbArray->initialise(adjGroup);
	ExpiryArraySP expiries = iddbArray->getExpiries(adjGroup);
	return expiries;
}

// initialize ids for initial guess
void Calibrator::InstanceID::initialise(InstanceIDSP ids, ObjFuncSP objFunc){

	Calibrator::InstanceIDDbArraySP iddbArray = Calibrator::InstanceIDDbArraySP::dynamicCast(ids);
	IObjectSP adjGroup(objFunc->getAdjustableGroup());
	iddbArray->initialise(adjGroup);
}

/** method to extract values from results used in bootstrapping */
// copied form CalibratorAddin
void Calibrator::getResultValues(CResultsSP currentResults,  // (I)
                            DoubleArray & aggregCalibValues)
{
    static const string method = "Calibrator::getResultValues";

    // extract calibrated values (XXX this is a complete shit show !!!)
    IObjectConstSP calibVars0 =
        currentResults->retrieveGreek(
        "Calibrator",
        OutputNameSP(new OutputName("calibratedVars")));
    IObjectSP calibVars1 = IObjectSP::constCast(calibVars0);
    if (!calibVars1 || !ObjectArray::TYPE->isInstance(calibVars1)){
        throw ModelException(method, "format error for calibrated vars");
    }
    ObjectArraySP calibVars2 = ObjectArraySP::dynamicCast(calibVars1);
    if (calibVars2->size() != Calibrator::InstanceID::NB_RES 
        || !ObjectArray::TYPE->isInstance((*calibVars2)[Calibrator::InstanceID::VALUE])){
            throw ModelException(method, "format error for calibrated values");
        }
        ObjectArraySP vals0 = ObjectArraySP::dynamicCast((*calibVars2)[Calibrator::InstanceID::VALUE]);
        int nbVals = vals0->size();
        //DoubleArray vals1;
        int iVal;
        for (iVal = 0; iVal < nbVals; ++iVal){
            if (!(*vals0)[iVal]){
                throw ModelException(method, "format error for calibrated vars");
            }
            else if (CDouble::TYPE->isInstance((*vals0)[iVal])){
                CDoubleSP vals2 = CDoubleSP::dynamicCast((*vals0)[iVal]);
                //vals1.push_back(vals2->doubleValue());
                // aggregate the ids and the calibrated values
                aggregCalibValues.push_back(vals2->doubleValue());
            }
            else{
                // XXX should never get there (for now!)
                throw ModelException(method, "format error for calibrated vars");
            }
        }

}
DRLIB_END_NAMESPACE
