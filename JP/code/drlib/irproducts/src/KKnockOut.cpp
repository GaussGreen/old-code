//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KKnockOut.cpp
//
//   Description : KnockOut component
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KKnockOut.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/TreeSliceOperExt.hpp"
#include "edginc/Addin.hpp"

#define  BARRIER_TOL   1E-4   /* Error tolerance for barriers */

DRLIB_BEGIN_NAMESPACE

/*************/

class KKnockOutTree : public FDProduct {
public:

    /************************ variables ************************/
    const KKnockOut*     inst;    

    FDProductSP prodFalse;
    FDProductSP prodTrue;

    FDProductArray knockoutProds;

    /************************ methods ************************/

    KKnockOutTree(const KKnockOut* inst, FDModel* model);

    DateTime getStartDate(void) const {return model->getValueDate();}
    virtual string getOutputName(void) const { return inst->outputName; }

    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results);
    void addModelResetDates(
        const DateTimeArray &modelResetDates, 
        const DateTimeArray &lateResetDates);
    virtual void init(Control*) const {}
    virtual void update(int& step, UpdateType type) {}
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const;
    virtual void initProd(void);

private:
    TreeSliceSP slice;
};


void koCombinedBinary(TreeSlice &result, 
                const TreeSlice &value1, double lo1, double hi1,
                const TreeSlice &value2, double lo2, double hi2,
                const TreeSlice &valForTrue,
                const TreeSlice &valForFalse,
                KKnockOut::Barrier::Enum style,
                vector<bool>& barrierBelongsToOutside,
                bool smooth);


/******************************** KKnockOutTree ********************************/

KKnockOutTree::KKnockOutTree(const KKnockOut* inst, FDModel* model) :
    FDProduct(model), inst(inst)
{   
    try {
        int i;
        knockoutProds.resize(inst->knockouts.size());

        for (i=0; i<inst->knockouts.size(); i++) {
            knockoutProds[i] = model->createProduct(inst->knockouts[i]->und);
        }

        prodFalse   = model->createProduct(inst->undFalse);
        prodTrue  = model->createProduct(inst->undTrue); 
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

void KKnockOutTree::addModelResetDates(
    const DateTimeArray &modelResetDates, 
    const DateTimeArray &lateResetDates)
{
    try {
        for (size_t i=0; i<knockoutProds.size(); i++) {
            knockoutProds[i]->addModelResetDates(modelResetDates, lateResetDates);
        }
        prodFalse->addModelResetDates(modelResetDates, lateResetDates);
        prodTrue->addModelResetDates(modelResetDates, lateResetDates);
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}


void KKnockOutTree::recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) 
{
    try {
        OutputRequest* request;
        request = ctrl->requestsOutput(OutputRequest::BARRIER_LEVEL);
        if (!request)
            return;

        DateTime valueDate = model->getValueDate();
        DateTime upperDate = BarrierLevel::barrierWindow(valueDate);
        upperDate = valueDate.rollDate(100000); // to ease debug
        for (int knk=0; knk<inst->knockouts.size(); ++knk) {
            Schedule* s;
            bool isUp;

            for (int b=0; b<2; ++b) // for up and low barriers
            {
                if (b==0) {
                    s = inst->knockouts[knk]->loBarrier.get();
                    isUp = false;
                } else {
                    s = inst->knockouts[knk]->hiBarrier.get();
                    isUp = true;
                }

                CashFlowArraySP subset(s->subset(valueDate, upperDate));
                if (!subset->empty()) {
                    BarrierLevelArray levels;
                    for (int i = 0; i < subset->size(); i++) {
                        BarrierLevel bl(isUp, (*subset)[i].date, (*subset)[i].amount);
                        levels.push_back(bl);
                    }
                    OutputRequestUtil::recordBarrierLevels(
                        ctrl, results, inst->outputName+"_"+Format::toString(knk), &levels);
                }
            }
        }
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

/** initialising and setting product variables */
void KKnockOutTree::initProd(void) {
    try {
        slice = model->createSlice();
        slice->name = "KKnockOut_"+inst->outputName;
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

static bool happensNow(const DateTimeArray& array, DateTime stepDate, int &pos) 
{
    DateTimeArray::const_iterator it;
    it = find(array.begin(), array.end(), stepDate);
    if (it==array.end()) return false;
    pos = it - array.begin();
    return true;
/* // more classic version
    for (int i=0; i < array.size(); ++i) {
        if (array[i] == stepDate) {
            pos = i;
            return true;
        }
    }
    return false;
    */
}

/** update products slice at each step */
const TreeSlice & KKnockOutTree::getValue(int step, DateTime eventDate) const { 
    try {
        DateTime stepDate = model->getDate(step);

        if (eventDate < model->getToday()) {
            CashflowInfo cfi;
            *slice = inst->getValue(eventDate, cfi);
            return *slice;
        }

        double loBar[2], hiBar[2]; // barriers today
        TreeSliceSP index[2];      // weighted observation index

        int barIdx;
        if (!happensNow(inst->knockouts[0]->loBarrier->getDateArray(), eventDate, barIdx))
            throw ModelException("Cannot be valued on " + eventDate.toString()+
                                 " because barrier level is not defined at this date.");

        vector<bool> barrierBelongsToOutside;

        for (size_t i=0; i<knockoutProds.size(); ++i) {
            loBar[i] = inst->knockouts[i]->loBarrier->getValueArray()[barIdx];
            hiBar[i] = inst->knockouts[i]->hiBarrier->getValueArray()[barIdx];

            if( loBar[i] > hiBar[i] ) {
                throw ModelException("At date " + eventDate.toString() +
                    + " the lower barrier " + Format::toString(loBar[i]) +
                    " is greater than the upper barrier " + Format::toString(hiBar[i]) +
                    " for KOEvent["+Format::toString(i)+"]");
            }

            const TreeSlice &obsIndex = knockoutProds[i]->getValue(step, eventDate);
            index[i] = obsIndex.clone(false);
            *index[i] = obsIndex * inst->knockouts[i]->idxWeight;

            { 
                bool isIn = KKnockOut::barrierIsIn(inst->barrierType,0);
                bool barrierBelongsToOutsideL = 
                    (inst->touchBarrierChangesState ? isIn : !isIn);
                barrierBelongsToOutside.push_back(barrierBelongsToOutsideL);
            }
        }

        const TreeSlice& sliceFalse = prodFalse->getValue(step, eventDate);
        const TreeSlice& sliceTrue  = prodTrue->getValue(step, eventDate);
        bool smooth = inst->smoothing == KKnockOut::Smoothing::SMOOTHING_STEP;

        if (knockoutProds.size() == 1) {
            double tol = (barrierBelongsToOutside[0] ? BARRIER_TOL : -BARRIER_TOL);
            if (smooth) {
                tol = 0.;
            }
            bool isIn = KKnockOut::barrierIsIn(inst->barrierType,0);
            const TreeSlice *sliceUp = (isIn ? &sliceFalse : &sliceTrue);
            const TreeSlice *sliceMid = (isIn ? &sliceTrue : &sliceFalse);
            const TreeSlice *sliceDown = (isIn ? &sliceFalse : &sliceTrue);

            koTernary(*slice,
                *index[0], loBar[0]*(1.+tol), hiBar[0]*(1.-tol),
                *sliceUp, *sliceMid, *sliceDown,
                barrierBelongsToOutside[0], // sic: see hyb3_koopt.c: Hyb3_KoOption_t
                smooth);
        }
        else {
            koCombinedBinary(*slice, 
                *index[0], loBar[0], hiBar[0],
                *index[1], loBar[1], hiBar[1],
                sliceTrue, sliceFalse, 
                inst->barrierType,
                barrierBelongsToOutside,
                smooth);
        }

        return *slice; 
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }

}

/********************************** KKnockOut::KOEvent *************************/

void KKnockOut::KOEvent::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("KOEvent component");
    REGISTER(KKnockOut::KOEvent, clazz)
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(und, "")
    FIELD(loBarrier, "")
    FIELD(hiBarrier, "") 
    FIELD(idxWeight, "weight assigned to the index, defaults to 1.0");
    FIELD_MAKE_OPTIONAL(idxWeight);
}

CClassConstSP const KKnockOut::KOEvent::TYPE = CClass::registerClassLoadMethod(
    "KKnockOut::KOEvent", typeid(KKnockOut::KOEvent), KKnockOut::KOEvent::load);

typedef KKnockOut::KOEventArray KKnockOutKOEventArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("KKnockOut::KOEventArray", KKnockOutKOEventArray);

/********************************** KKnockOut **********************************/

void KKnockOut::addResetDates(const DateTimeArray &resetDates) {
    try {
        for (int i=0; i<knockouts.size(); i++) {
            knockouts[i]->und->addResetDates(resetDates);
        }
        undFalse->addResetDates(resetDates);
        undTrue->addResetDates(resetDates);
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

void KKnockOut::validatePop2Object(void) {
    try {
        if ( knockouts.size()==0 )
            throw ModelException("Empty knockout list");
        if ( knockouts.size()>2 )
            throw ModelException("At most 2 knockout underliers are supported currently");
        if (!undFalse)
            undFalse = upUnd;
        if (!undTrue)
            undTrue = midUnd;
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

bool KKnockOut::barrierIsIn(Barrier::Enum barrier, int idx) {
    ASSERT(idx==0 || idx==1);
    switch (barrier){
        case Barrier::ONE_IN:  return true;
        case Barrier::ONE_OUT: return false;
        case Barrier::IN_IN:   return true;
        case Barrier::IN_OUT:  return idx==0;
        case Barrier::OUT_IN:  return idx==1;
        case Barrier::OUT_OUT: return false;
        case Barrier::IN:      return true;
        case Barrier::OUT:     return false;
        default: ASSERT(0);
    }
}

bool KKnockOut::barrierBoolValue(int koIdx, int barIdx, CashflowInfo &cfiLoc, DateTime stepDate) const {
                                 
    double lo = knockouts[koIdx]->loBarrier->getValueArray()[barIdx];
    double hi = knockouts[koIdx]->hiBarrier->getValueArray()[barIdx];

    if (lo > hi) {
        throw ModelException("In knockouts[" + Format::toString(koIdx)
            + "] the lower barrier is greater than the upper barrier");
    }
    double index;
    {
        double obsIndex = knockouts[koIdx]->und->getValue(stepDate, cfiLoc);
        index = obsIndex * knockouts[koIdx]->idxWeight;
    }

    bool barIn = barrierIsIn(barrierType, koIdx);
    bool barrierBelongsToOutside = (touchBarrierChangesState ? barIn : !barIn);
    if (index < lo || (index == lo && barrierBelongsToOutside) )
        return !barIn;
    if (hi < index || (hi == index && barrierBelongsToOutside) )
        return !barIn;
    return barIn;
}

// the IGetKnownValue interface for retrieving a single sample
double KKnockOut::getValue(DateTime stepDate, CashflowInfo &cfi) const
{
    try {
        CashflowInfo cfiLoc; // local CashflowInfo, in case we want to cancel the computation
        double valFalse = undFalse->getValue(stepDate, cfiLoc);
        double valTrue = undTrue->getValue(stepDate, cfiLoc);
        int barIdx;

        if (!happensNow(knockouts[0]->loBarrier->getDateArray(), stepDate, barIdx))
            throw ModelException("Cannot be valued on "+stepDate.toString());

        double result;

        if (knockouts.size() == 1) {
            if (barrierBoolValue(0, barIdx, cfiLoc, stepDate))
                result = valTrue;
            else result = valFalse;
        }
        else {
            bool isTrue[2];
            for (int i=0; i<2; ++i) {
                isTrue[i] = barrierBoolValue(i, barIdx, cfiLoc, stepDate);
            }
            if (barrierType==Barrier::ONE_IN
            ||  barrierType==Barrier::ONE_OUT) {
                result = ( isTrue[0] ||  isTrue[1] ? valTrue : valFalse);
            }
            else result = ( isTrue[0] &&  isTrue[1] ? valTrue : valFalse);
        }

        if (cfiLoc.amountType != CashflowInfo::AmountType::KNOWN) {
            // we need to know the exact value, an estimation is not good enough.
            cfi.updateAmountType(CashflowInfo::AmountType::UNKNOWN);
            return 0.;
        }
        cfi.merge(cfiLoc, outputName);

        return result;
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

void KKnockOut::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        if (knockouts.size() > 2) throw ModelException(
            "At most 2 knockout indexes are supported at this time.");

        if (knockouts.size() == 0) throw ModelException(
            "At least 1 knockout index is needed.");

        for (int i=0; i<knockouts.size(); ++i) {
            const DateTimeArray& d0 = knockouts[0]->loBarrier->getDateArray();
            const DateTimeArray& d1 = knockouts[i]->loBarrier->getDateArray();
            const DateTimeArray& d2 = knockouts[i]->hiBarrier->getDateArray();
            
            if (!(d1==d2)) throw ModelException("Barrier dates for hi&lo must the"
                " same for knockouts["+Format::toString(i)+"]");

            if (!(d0==d1)) throw ModelException("Barrier dates for knockouts[0]"
                " and knockouts["+Format::toString(i)+"] must be the same");

            knockouts[i]->und->setup(model, market);
        }
        if (knockouts.size()==1 && !downUnd) {
            throw ModelException("You should provide "
                "\"downUnd\" when there is only one \"knockouts\"");
        }
        undFalse->setup(model, market);
        undTrue->setup(model, market);

        // check consistency of barrierType with the number of "knockouts"

        bool isSimpleBarrier = (barrierType==KKnockOut::Barrier::IN
                             || barrierType==KKnockOut::Barrier::OUT);

        if (knockouts.size()==1) {
            if (!isSimpleBarrier) {
                /*
                throw ModelException("When nb of observable provided (knockouts) is 1 "
                "the barrier type can only be IN or OUT");
                */
                // for compatibility
                if (barrierIsIn(barrierType, 0)) {
                    barrierType=KKnockOut::Barrier::IN;
                }
                else barrierType=KKnockOut::Barrier::OUT;
            }
        }
        else {
            if (isSimpleBarrier) {
                throw ModelException("When nb of observable provided (knockouts) is 2 "
                "the barrier type can be neither IN nor OUT. "
                "Use a combined type (IN_IN,...).");
            }
        }
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

/** implement FDModel product interface */
FDProductSP KKnockOut::createProduct(FDModel * model) const {    
    if (!setupCalled) throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
    return FDProductSP(new KKnockOutTree(this, model));
}

void KKnockOut::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("KKnockOut composite component");
    REGISTER(KKnockOut, clazz)
    SUPERCLASS(KComponent);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(KKnockOut::defaultConstructor);
    FIELD(knockouts, "")
    FIELD(smoothing, "Type of smoothing");
    FIELD(barrierType, "");
    FIELD_MAKE_OPTIONAL(barrierType);
    FIELD(touchBarrierChangesState, "");
    FIELD_MAKE_OPTIONAL(touchBarrierChangesState);
    FIELD(upUnd, "index value returned when underlyings are considered to exceed the hiBarrier. Also value for FALSE for combined barriers (IN_IN,...)");
    FIELD(midUnd, "index value returned when underlyings are considered to be in the barriers. Also value for TRUE for combined barriers (IN_IN,...)");
    FIELD(downUnd, "index value returned when underlyings are considered to be below the loBarrier. Not used for combined barriers (IN_IN,...)");
    FIELD_MAKE_OPTIONAL(downUnd);
    FIELD(undFalse, "");  FIELD_MAKE_OPTIONAL(undFalse);
    FIELD(undTrue, "");  FIELD_MAKE_OPTIONAL(undTrue);
    Addin::registerConstructor(Addin::UTILITIES, KKnockOut::TYPE);
}
   
CClassConstSP const KKnockOut::TYPE = CClass::registerClassLoadMethod(
    "KKnockOut", typeid(KKnockOut), KKnockOut::load);


START_PUBLIC_ENUM_DEFINITION(KKnockOut::Barrier::Enum, "");
ENUM_VALUE_AND_NAME(KKnockOut::Barrier::ONE_IN, "ONE_IN", "");
ENUM_VALUE_AND_NAME(KKnockOut::Barrier::ONE_OUT, "ONE_OUT", "");
ENUM_VALUE_AND_NAME(KKnockOut::Barrier::IN_IN, "IN_IN", "");
ENUM_VALUE_AND_NAME(KKnockOut::Barrier::IN_OUT, "IN_OUT", "");
ENUM_VALUE_AND_NAME(KKnockOut::Barrier::OUT_IN, "OUT_IN", "");
ENUM_VALUE_AND_NAME(KKnockOut::Barrier::OUT_OUT, "OUT_OUT", "");
ENUM_VALUE_AND_NAME(KKnockOut::Barrier::IN, "IN", "");
ENUM_VALUE_AND_NAME(KKnockOut::Barrier::OUT, "OUT", "");
END_ENUM_DEFINITION(KKnockOut::Barrier::Enum);

START_PUBLIC_ENUM_DEFINITION(KKnockOut::Smoothing::Enum, "");
ENUM_VALUE_AND_NAME(KKnockOut::Smoothing::SMOOTHING_NONE, "SMOOTHING_NONE", "");
ENUM_VALUE_AND_NAME(KKnockOut::Smoothing::SMOOTHING_STEP, "SMOOTHING_STEP", "");
END_ENUM_DEFINITION(KKnockOut::Smoothing::Enum);

/**********************************************************************************************************/
/****************************************** Knock Out operations ******************************************/
/**********************************************************************************************************/

void koCombinedBinary(TreeSlice &result, 
                const TreeSlice &value0, double lo0, double hi0,
                const TreeSlice &value1, double lo1, double hi1,
                const TreeSlice &valForTrue,
                const TreeSlice &valForFalse,
                KKnockOut::Barrier::Enum style,
                vector<bool>& barrierBelongsToOutside,
                bool smooth)
{
    double tol[2];
    double valIn[2];
    double valOut[2];
    bool boundaryIn[2];
    TreeSliceSP tmp[2];

    for (int i=0; i<2; ++i) {
        if (KKnockOut::barrierIsIn(style,i)) {
            valIn[i] = 1;
            valOut[i] = 0;
        }
        else {
            valIn[i] = 0;
            valOut[i] = 1;
        }
        boundaryIn[i] = barrierBelongsToOutside[i]; // see hyb3_koopt.c: Hyb3_KoOption_t
        tol[i] = (barrierBelongsToOutside[i] ? BARRIER_TOL : -BARRIER_TOL);
        if (smooth) {
            tol[i] = 0.;
        }
        tmp[i] = result.clone(false);
    }

    koTernary(*tmp[0], value0, lo0*(1+tol[0]), hi0*(1-tol[0]), valOut[0], valIn[0], valOut[0], boundaryIn[0], smooth);
    koTernary(*tmp[1], value1, lo1*(1+tol[1]), hi1*(1-tol[1]), valOut[1], valIn[1], valOut[1], boundaryIn[1], smooth);

    switch (style) {
        case KKnockOut::Barrier::IN_IN  :
        case KKnockOut::Barrier::IN_OUT :
        case KKnockOut::Barrier::OUT_IN :
        case KKnockOut::Barrier::OUT_OUT:
            result = *tmp[0] * *tmp[1]  * (valForTrue - valForFalse) + valForFalse; break;
        case KKnockOut::Barrier::ONE_IN : 
        case KKnockOut::Barrier::ONE_OUT : 
            result = (*tmp[0] + (1-*tmp[0])* *tmp[1]) * (valForTrue - valForFalse) + valForFalse; break; // (*tmp[1] || *tmp[2])
        default: throw ModelException(__FUNCTION__, "Unsupported knock-out style");
    }
    return;
}

/**********************************/
bool KKnockOutLoad(void) {
    return (KKnockOut::TYPE !=0);
}

DRLIB_END_NAMESPACE
