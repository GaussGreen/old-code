//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KAccumulatedCF.cpp
//
//   Description : Tarn component
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KAccumulatedCF.hpp"
#include "edginc/RatesUtils.hpp"

DRLIB_BEGIN_NAMESPACE

class KAccumulatedCFTree : public FDProduct,
                           public TreeSliceLayer::StateSupport
{
    FDModel * model;
public:
    /************************ state variables support *********************/
	virtual bool gridChanges(int step) const;
	virtual void setGrid(int step, bool init = false);
    virtual void transition(
        const vector< const TreeSlice * > & srcSlices,
        const vector< double > & gridIn,
        vector< double > & gridOut ) const;
    /************************ end of state variables support *********************/

    /******************** methods ********************/
    KAccumulatedCFTree(const KAccumulatedCF* inst, FDModel* model);

    virtual DateTime getStartDate(void) const { return model->getDate(0);}
    virtual string getOutputName(void) const { return inst->outputName; }
    virtual void init(Control*) const;
    virtual void initProd(void);
    virtual void update(int& step, UpdateType type) {}
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const { return getGridSlice(); }
    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results);

    const KAccumulatedCF*  inst; // this overrides base
    FDProductSP underlying;
    TreeSliceSP value;
};

/****************************** KAccumulatedCFTree ********************************/

KAccumulatedCFTree::KAccumulatedCFTree(const KAccumulatedCF* inst, FDModel* model) :
    FDProduct(model),
    TreeSliceLayer::StateSupport( "KAccumulated", inst->todayValue, INTERP_LINEAR ),
    model(model),
    inst(inst)
{
    try {
        underlying = model->createProduct(inst->underlier);
    }
    catch (exception& e){ 
        throw ModelException(e,"KAccumulatedCFTree::KAccumulatedCFTree"); 
    }
}

void KAccumulatedCFTree::recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
    if (inst->recordOutputName)
        recordSliceToOutputName(ctrl, results, model, 
            inst->isTopComponent(), inst->outputName, "", getValue(0, model->getDate(0)));
    inst->recordExtraOutput(ctrl, results);
}

void KAccumulatedCFTree::init(Control*) const {
    model->addCritDates(inst->cfDates);
}

/** initializing and setting product variables */
void KAccumulatedCFTree::initProd(void){
    const string curve = inst->discount->getName();

    model->registerStateSupport(this);

    // value slice just for getValue(). no DEV
    value = model->createSlice(curve);
    value->name = inst->outputName+"_value";
    *value = 0.0;
    startDEV(value);
}

/************************ state variables operation *********************/
void KAccumulatedCFTree::setGrid(int step, bool init) {
    currGrid.resize(inst->numState);
    DateTime stepDate = model->getDate(step);
    double loVal = inst->todayValue;
    double hiVal = inst->boundaryValues->interpolate(stepDate);
    double slope = (hiVal-loVal)/(inst->numState-1.0);
    for (int i=0; i<inst->numState; i++){ 
        currGrid[i] = slope * i + loVal;
    }

    if( init )
        prevGrid = currGrid;
}

// compute grid levels for state variables
void KAccumulatedCFTree::transition(
    const vector< const TreeSlice * > & srcSlices,
    const vector< double > & gridIn,
    vector< double > & gridOut ) const
{
    int nbGrid = gridOut.size();
    ASSERT( nbGrid == (int)gridIn.size() );

    double cashFlow = srcSlices[ 0 ]->calc();

    for( int i = 0; i < nbGrid; ++i )
        gridOut[ i ] = gridIn[ i ] + cashFlow;
}

bool KAccumulatedCFTree::gridChanges(int step) const {
    return RatesUtils::happensNow(inst->cfDates, model->getDate(step));
}

/********************************** KAccumulatedCF **********************************/
/** implement FDModel product interface */

void KAccumulatedCF::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        underlier->setup(model, market);
        
        KComponentSP comp = KComponentSP::dynamicCast(underlier);
        comp->getCfDates(cfDates);
        DateTime::doSortUniq(cfDates);
    }
    catch (exception& e) {
        throw makeException(e,__FUNCTION__);
    }
}

FDProductSP KAccumulatedCF::createProduct(FDModel * model) const {
    if (!setupCalled) 
        throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
        
    return FDProductSP(new KAccumulatedCFTree(this, model));
}

void KAccumulatedCF::load(CClassSP& clazz)
{
    try {
        clazz->setPublic();
        REGISTER(KAccumulatedCF, clazz)
        SUPERCLASS(KComponent);
        IMPLEMENTS(FDModel::IIntoProduct);
        EMPTY_SHELL_METHOD(KAccumulatedCF::defaultConstructor);
        //IMPLEMENTS(LastSensDate);
        FIELD(numState, "number of states i.e. size of state grid");
        FIELD(todayValue, "current values for the state variable grid. Defaults to 0.0.");
        FIELD_MAKE_OPTIONAL(todayValue);
        FIELD(boundaryValues, "high values for the state variable grid");
        FIELD(underlier,"");
        Addin::registerConstructor(Addin::UTILITIES, KAccumulatedCF::TYPE);
    }
    catch (exception& e) {
        throw ModelException(e,__FUNCTION__);
    }
}

CClassConstSP const KAccumulatedCF::TYPE = CClass::registerClassLoadMethod(
    "KAccumulatedCF", typeid(KAccumulatedCF), KAccumulatedCF::load);

/**********************************/
bool KAccumulatedCFLoad()
{
    return KAccumulatedCF::TYPE !=0;
}

DRLIB_END_NAMESPACE
