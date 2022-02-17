//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KKnockIO.hpp
//
//   Description : one touch KO swap component
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KKnockIO.hpp"
#include "edginc/RatesUtils.hpp"
#include "edginc/TreeSliceOperExt.hpp"
#include "edginc/BarrierLevel.hpp"
#include "edginc/OutputRequestUtil.hpp"

DRLIB_BEGIN_NAMESPACE

/****************************** KKnockIOTree ********************************/

class KKnockIOTree : public FDProduct {

public:
    /******************** methods ********************/
    KKnockIOTree(const KKnockIO* inst, FDModel* model);

    virtual DateTime getStartDate(void) const { return model->getDate(0);}
    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results);
    virtual void init(Control*) const;
    virtual void initProd(void);
    virtual void update(int& step, UpdateType type);
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const;

protected:
    void createSliceForDev(TreeSliceSP &s, string const& name);
    void createSliceArray(TreeSliceArray &a, string const& name, int n);

    KKnockIO const *  inst;
    FDProductArray   underlyings;
    FDProductArray   indicators;
    FDProductArray   rebates;
    int nbDates;
    TreeSliceSP mainValueSl;
    TreeSliceSP rebateValueSl;
    TreeSliceArray undDelayedSlAr;
    TreeSliceArray rebateDelayedSlAr;
    TreeSliceSP cumulUndSl;
    TreeSliceSP getValueSl;
    string discountCurveName;
    bool multiMode;
};

static void createProducts(FDModel* model, IProdCreatorArray const& instArray, FDProductArray &prodArray, DateTimeArray &resetDates) {
    try {
        int nbDates = resetDates.size();
        prodArray.resize(nbDates);
        if (instArray.size()==1) {
            FDProductSP prod = model->createProduct(instArray[0]);
            for (int i=0; i<nbDates; ++i) {
                prodArray[i] = prod;
            }
        }
        else {
            for (int i=0; i<nbDates; ++i) {
                prodArray[i] = model->createProduct(instArray[i]);
            }
        }
        for (int i=0; i<nbDates; ++i) {
            prodArray[i]->addModelResetDate(resetDates[i]);
        }
    }
    catch (exception& e){
        throw ModelException(e, __FUNCTION__);
    }
}

KKnockIOTree::KKnockIOTree(const KKnockIO* inst, FDModel* model)
: FDProduct(model), inst(inst), nbDates(inst->sched->nbDates) 
{
    try {
        multiMode = (inst->underlyings.size()>1);
        discountCurveName = inst->getActualDiscount()->getName();
        createProducts(model, inst->indicators, indicators, inst->sched->notifDate);
        createProducts(model, inst->underlyings, underlyings, inst->sched->exerciseDate);
        createProducts(model, inst->rebates, rebates, inst->sched->exerciseDate);
    }
    catch (exception& e){
        throw makeException(e, __FUNCTION__);
    }
}

void KKnockIOTree::init(Control*) const{
    try {
        model->addCritDates(inst->sched->notifDate);
        model->addCritDates(inst->sched->exerciseDate);
    }
    catch (exception& e){
        throw makeException(e, __FUNCTION__);
    }
}

void KKnockIOTree::createSliceForDev(TreeSliceSP &s, string const& name) {
    s = model->createSlice(discountCurveName);
    s->name = inst->outputName + "_"+name;
    *s = 0.;
    startDEV(s);
}

void KKnockIOTree::createSliceArray(TreeSliceArray &a, string const& name, int n) {
    a.resize(n);
    for (int i=0; i<n; ++i) {
        a[i] = model->createSlice(discountCurveName);
        a[i]->name = inst->outputName + "_"+name+"["+Format::toString(i)+"]";
    }
}

void KKnockIOTree::initProd(void){
    try {
        if (!multiMode && inst->knockIOType == KKnockIO::KnockIOType::K_OUT) {
            createSliceForDev(rebateValueSl, "rebateValueSl");
            createSliceArray(rebateDelayedSlAr, "rebateDelayedSlAr", nbDates);
        }
        if (multiMode && inst->knockIOType == KKnockIO::KnockIOType::K_IN) {
            createSliceForDev(cumulUndSl, "cumulUndSl");
        }
        createSliceArray(undDelayedSlAr, "undDelayedSlAr", nbDates);
        createSliceForDev(mainValueSl, "mainValueSl");
    }
    catch (exception& e){
        throw makeException(e, __FUNCTION__);
    }
}

void KKnockIOTree::update(int& step, UpdateType type){
    try {
        DateTime curDate = model->getDate(step);

        // on exercise date, store value of the underlying
        for (int i=0; i<nbDates; ++i) {
            if (inst->sched->exerciseDate[i] == curDate) {
                const TreeSlice &undSl = underlyings[i]->getValue(step, curDate);
                *undDelayedSlAr[i] = undSl;
                startDEV(undDelayedSlAr[i]);
                if (inst->knockIOType == KKnockIO::KnockIOType::K_OUT) {
                    const TreeSlice &rebateSl = rebates[i]->getValue(step, curDate);
                    *rebateDelayedSlAr[i] = rebateSl;
                    startDEV(rebateDelayedSlAr[i]);
                }
            }
        }

        for (int i=0; i < nbDates; ++i) {
            if (inst->sched->notifDate[i]!=curDate) {
                continue;
            }
            const TreeSlice &undDelayedSl = *undDelayedSlAr[i];
            const TreeSlice &indicSl = indicators[i]->getValue(step, curDate);
            if (multiMode) {
                if (inst->knockIOType == KKnockIO::KnockIOType::K_IN) {
                    // mainValue contains knock-in
                    *cumulUndSl += undDelayedSl;
                    *mainValueSl = indicSl * *cumulUndSl + (1 - indicSl) * *mainValueSl;
                }
                else {
                    // mainValue contains knock-out + rebate
                    const TreeSlice &rebateDelayedSl = *rebateDelayedSlAr[i];
                    *mainValueSl = indicSl * (undDelayedSl + *mainValueSl) + (1. - indicSl) * rebateDelayedSl;
                }
            }
            else {
                // mainValueSl contains knock-in
                *mainValueSl = indicSl * undDelayedSl + (1 - indicSl) * *mainValueSl;
                // rebate is stored separately
                if (inst->knockIOType == KKnockIO::KnockIOType::K_OUT) {
                    // update rebateValueSl
                    const TreeSlice &rebateDelayedSl = *rebateDelayedSlAr[i];
                    *rebateValueSl = indicSl * *rebateValueSl + (1. - indicSl) * rebateDelayedSl;
                }
            }
            if (inst->knockIOType == KKnockIO::KnockIOType::K_OUT) {
                stopDEV(rebateDelayedSlAr[i]);
                rebateDelayedSlAr[i].reset(); // set pointer to NULL to free memory
            }
            stopDEV(undDelayedSlAr[i]);
            undDelayedSlAr[i].reset(); // set pointer to NULL to free memory
        }
    }
    catch (exception& e){
        throw makeException(e, __FUNCTION__);
    }
}

const TreeSlice & KKnockIOTree::getValue(int step, DateTime eventDate) const {
    try {
        if (historicalValueAvailable(*inst, *getValueSl, eventDate))
            return *getValueSl;

        if (eventDate != model->getDate(step)) {
            throw ModelException("Cannot be valued at a date != currentDate");
        }
        if (!multiMode && inst->knockIOType == KKnockIO::KnockIOType::K_OUT) {
            *getValueSl = underlyings[0]->getValue(step, eventDate) - *mainValueSl + *rebateValueSl;
        }
        else *getValueSl = *mainValueSl;

        return *getValueSl;
    }
    catch (exception& e){
        throw makeException(e, __FUNCTION__);
    }
}

void KKnockIOTree::recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
    try {
        if (inst->recordOutputName) {
            recordSliceToOutputName(ctrl, results, model,
                inst->isTopComponent(), inst->outputName, "", getValue(0, model->getDate(0)));
        }
        inst->recordExtraOutput(ctrl, results);
    }
    catch (exception& e){
        throw makeException(e, __FUNCTION__);
    }
}

/*********************************** KKnockIO ***********************************/


static void setupProducts(const IModel* model, const MarketData* market, IProdCreatorArray const& instArray, DateTimeArray &resetDates)
{
    try {
        if (instArray.size()==1) {
            instArray[0]->addResetDates(resetDates);
        }
        else {
            if (instArray.size() != resetDates.size()) {
                throw ModelException("instArray.size() != resetDates.size()");
            }
            for (int i=0; i<instArray.size(); ++i) {
                instArray[i]->addResetDate(resetDates[i]);
            }
        }
        for (int i=0; i < instArray.size(); ++i) {
            instArray[i]->setup(model, market);
        }
    }
    catch (exception& e){
        throw ModelException(e, __FUNCTION__);
    }
}


void KKnockIO::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        setupProducts(model, market, indicators, sched->notifDate);
        setupProducts(model, market, underlyings, sched->exerciseDate);
        setupProducts(model, market, rebates, sched->exerciseDate);
        sched->setup(model, market);
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

FDProductSP KKnockIO::createProduct(FDModel * model) const {
    if (!setupCalled)
        throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
    return FDProductSP(new KKnockIOTree(this, model));
}

void KKnockIO::load(CClassSP& clazz) {

    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Knock In/Out component based on SY component");
    REGISTER(KKnockIO, clazz);
    SUPERCLASS(KComponent);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(knockIOType,"");
    FIELD(indicators,"");
    FIELD(sched, "");
    FIELD(rebates, "");
    Addin::registerConstructor(Addin::UTILITIES, KKnockIO::TYPE);
}

CClassConstSP const KKnockIO::TYPE = CClass::registerClassLoadMethod(
    "KKnockIO", typeid(KKnockIO), KKnockIO::load);

START_PUBLIC_ENUM_DEFINITION(KKnockIO::KnockIOType::Enum, "");
ENUM_VALUE_AND_NAME(KKnockIO::KnockIOType::K_IN, "IN", "");
ENUM_VALUE_AND_NAME(KKnockIO::KnockIOType::K_OUT, "OUT", "");
END_ENUM_DEFINITION(KKnockIO::KnockIOType::Enum);


/******************************/
// for type linking
bool KKnockIOLoad(void){
    return (KKnockIO::TYPE != 0);
}

DRLIB_END_NAMESPACE

