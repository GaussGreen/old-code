//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KPyComponent.cpp
//
//   Description : General python component where a serialized python 
//                 class contains the calc;
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KPyComponent.hpp"

DRLIB_BEGIN_NAMESPACE

class KPyComponentTree : public FDProduct {
public:
    KPyComponentConstSP inst;
    FDProductArray underlyings;
    

    /******************** methods ********************/
    KPyComponentTree(const KPyComponentConstSP &inst, FDModel* model);
    virtual DateTime getStartDate(void) const { return model->getDate(0);}
    virtual string getOutputName(void) const { return inst->outputName; }
    virtual void addModelResetDates(
        const DateTimeArray &modelResetDates, 
        const DateTimeArray &lateResetDates);

    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
        if (inst->recordOutputName)
            recordSliceToOutputName(ctrl, results, model, 
                inst->isTopComponent(), inst->outputName, "", getValue(0, model->getDate(0)));
        inst->recordExtraOutput(ctrl, results);
    }

    virtual void init(Control*) const {}
    virtual void initProd(void);
    virtual void update(int& step, UpdateType type);
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const;
    virtual const TreeSlice & getCashFlow(int step) const;

private:
    TreeSliceSP resultslice;

    // pointer to the de-serialized python functor/component constructed
    // and sent through the QLib interface 
    // (in instrument field serializedPyObject)
    PyObject* pythonObject;
};

/****************************** KPyComponentTree ********************************/

KPyComponentTree::KPyComponentTree(const KPyComponentConstSP &inst, FDModel* model) :
    FDProduct(model), inst(inst) {

    try {
        for (i = 0; i < underlyings.size(); i++)
            underlyings[i] = model->createProduct(inst->x0);
    }
    catch (exception& e){
        string errMsg = "Error creating underlying products of KPyComponent , "
                        "object outputName = " + inst->outputName;
        throw ModelException(e, __FUNCTION__, errMsg);
    }
}

// propagate any reset date events to underlying components
void KPyComponentTree::addModelResetDates(
    const DateTimeArray &modelResetDates, 
    const DateTimeArray &lateResetDates)
{
    try {
        
    }
    catch (exception& e){
        throw ModelException(e, __FUNCTION__, inst->outputName);
    }
}

/** initializing and setting product variables */
void KPyComponentTree::initProd(void){
    try {

        // de-serialize the pickle object into a PyObject
        pythonObject = load(serializedPyObject.c_str());

        // internal slice to store results of expression
        resultSlice = model->createSlice();
        resultSlice->name = inst->outputName+"_KPyComponent_slice";

        // call python object constructor at this point
        // ie. __init__(self, schedule, underlyingList, resultSlice):
        // to register references to the products internal resultSlice,
        // the supplied schedule and list of IProdCreators
        // - have to construct the appropriate types for these objects
        // as PyObjects, get the refcounts correct etc.
    }
    catch (exception& e){
        throw ModelException(e, __FUNCTION__, inst->outputName);
    }
}


const TreeSlice & KPyComponentTree::getValue(int step, DateTime eventDate) const {
    try {
        PyObject *pArgs, *pValue,  *pFunc;

        // get function reference for python object's getValue function
        pFunc = PyObject_GetAttrString(pythonObject, "getValue");

        // construct update function argument list
        // expected prototype = getValue(Date stepDate, Date eventDate)
        DateTime currentDate = model->getDate(step); // the current date
        
        int nbArgs = 2;
        pArgs = PyTuple_New(nbArgs);
        pValue = PyDate_FromString(currentDate.p());
        PyTuple_SetItem(pArgs, 0, pValue);  // add argument to argument list
        pValue = PyDate_FromString(eventDate.p());
        PyTuple_SetItem(pArgs, 0, pValue);  // add argument to argument list

        // call the actual python function with argument list
        pValue = PyObject_CallObject(pFunc, pArgs);

    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__, "timeStep = " + Format::toString(step) +
                             ", eventDate = " + eventDate.toString() + ", outputName = " + 
                             getOutputName());
    }
}


void KPyComponentTree::update(int& step, UpdateType type) {
    try {
        PyObject *pArgs, *pValue,  *pFunc;

        // get function reference for python object's update function
        pFunc = PyObject_GetAttrString(pythonObject, "update");

        // construct update function argument list
        // expected prototype = update(Date stepDate)
        DateTime currentDate = model->getDate(step); // the current date
        pValue = PyDate_FromString(currentDate.p());
        int nbArgs = 1;
        pArgs = PyTuple_New(nbArgs);
        PyTuple_SetItem(pArgs, 0, pValue);  // add argument to argument list

        // call the actual python function with argument list
        pValue = PyObject_CallObject(pFunc, pArgs);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

/*********************************** KPyComponent ***********************************/

void KPyComponent::addResetDates(const DateTimeArray &resetDatesP) {
    resetDates.insert(resetDates.end(), resetDatesP.begin(), resetDatesP.end());
}

void KPyComponent::validatePop2Object(void) {
    try {
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__, outputName);
    }
}

void KPyComponent::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        validatePop2Object();

        // propagate any reset dates from parents to underlying instruments
        for (i = 0; i < underlyings.size(); i++)
            underlyings[i]->addResetDates(resetDates);
            underlyings[i]->setup(model, market);
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__, outputName);
    }
}

// known values
double KPyComponent::getValue(DateTime date, CashflowInfo &cfi) const {
    try {
        throw ModelException("Not yet implemented - expression acting on "
            "past values");
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__, outputName);
    }
}


FDProductSP KPyComponent::createProduct(FDModel * model) const {
    if (!setupCalled) 
        throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
    return FDProductSP(new KPyComponentTree(KPyComponentConstSP(this), model));
}

void KPyComponent::load(CClassSP& clazz) {

    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Wrapper style Simple IR swap interface component");
    REGISTER(KPyComponent, clazz);
    SUPERCLASS(KComponent);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);
    // exported fields
    //FIELD(schedule, "Schedule of some type")
    //FIELD_MAKE_OPTIONAL(schedule);
    FIELD(serializedPyObject, "serialized python pricing object using pickle\n"
        "Pricing object must implement the following functions - based on FDProduct:\n"
        "update(update(int step, UpdateType type)\n
        "getValue(int step, DateTime eventDate)\n");
    FIELD(underlyings, "generic list of underlying IProdCreators")
    
    // transient
    FIELD(resetDates,""); 
    FIELD_MAKE_TRANSIENT(resetDates);
}

CClassConstSP const KPyComponent::TYPE = CClass::registerClassLoadMethod(
    "KPyComponent", typeid(KPyComponent), KPyComponent::load);

/******************************/
// for type linking
bool KPyComponentLoad(void){
    return (KPyComponent::TYPE != 0);
}

DRLIB_END_NAMESPACE

