//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KOption.cpp
//
//   Description : option component
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KOption.hpp"
#include "edginc/PhysicalSettlement.hpp"
#include "edginc/EventResults.hpp"
#include "edginc/TreeSliceOperExt.hpp"
#include "edginc/Addin.hpp"


DRLIB_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
//  KOptionTree
//  FDModel's pricing product/calc 
//  for the KOption component
//-----------------------------------------------------------------------------

KOptionTree::KOptionTree(const KOptionConstSP &inst, FDModel* model) :
    FDProduct(model), inst(inst), 
    optionType(inst->optionType)
{
    try {
        int nbExer = inst->sched->exerciseDate.size();
        undProd.resize(nbExer);
        if (inst->undList.size()>1) {
            // different underlyings for each exercise
            for (int i=0; i<nbExer; ++i) {
                undProd[i] = model->createProduct(inst->undList[i]);
            }
        }
        else {
            // same underlying for all exercises
            FDProductSP und = model->createProduct(inst->undList[0]);
            for (int i=0; i<nbExer; ++i) {
                undProd[i] = und;
            }
        }
        for (int i=0; i<nbExer; ++i) {
            DateTime exerDate = inst->sched->exerciseDate[i];
            undProd[i]->addModelResetDate(exerDate);
        }
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

void KOptionTree::recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
    try {
        string ccy = inst->discount->getCcy();
        const TreeSlice &price = getValue(0, model->getDate(0));

        recordSliceToOutputName(ctrl, results, model, inst->isTopComponent(), 
            inst->outputName, ccy, price);

        // record actual option price - != KOption price if callable/puttable/underlying
        recordSliceToOutputName(ctrl, results, model, false, 
            inst->outputName+"_OPTION_PRICE", ccy, *optionPrice);
        
        if (inst->isPhysical) {
            // If this is the first option from the top (isPhysical), 
            // report the value into the official output name as well
            recordSliceToOutputName(ctrl, results, model, false, 
                "OPTION_PRICE", ccy, *optionPrice);
        }

        if (keepValueIdx>=0) {
            recordSliceToOutputName(ctrl, results, model, false, 
                inst->outputName+"_SKIP_NEXT_EXERCISE", ccy, *optSkipExerPrice);
        }
        inst->recordExtraOutput(ctrl, results);
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

void KOptionTree::init(Control*) const{
    model->addCritDates(inst->sched->notifDate);
    model->addCritDates(inst->sched->exerciseDate);

    // register zeros to discount the strikes from exercise to notification
    DateTime today = model->getToday();
    for (int i=0; i<inst->sched->nbDates; ++i) {
        if (today <= inst->sched->notifDate[i]) {
            model->registerZero(inst->sched->notifDate[i], inst->sched->exerciseDate[i], inst->discount->getName());
        }
    }
}

void KOptionTree::initProd(void){
    try {
        const string curve = inst->discount->getName();

        optionPrice = model->createSlice(curve);
        optionPrice->name = inst->outputName + "_optionPrice";
        *optionPrice = 0.;
        startDEV(optionPrice);

        optSkipExerPrice = model->createSlice(curve);
        optSkipExerPrice->name = inst->outputName + "_optSkipExerPrice";

        undExer = model->createSlice(curve);
        undExer->name = inst->outputName + "_undExer";
        *undExer = 0.;
        startDEV(undExer);

        getValuePrice = model->createSlice();

        // calc keepValueIdx
        DateTime today = model->getToday();
        keepValueIdx=-1;
        DateTimeArray const &notifDates = inst->sched->notifDate.getDates();
        for (int i=0; i<notifDates.size(); ++i) {
            if (today <= notifDates[i]) {
                if (i+1==notifDates.size()) {
                    break;
                }
                keepValueIdx = i+1;
                break;
            }
        }
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

/** update products slice at each step */
void KOptionTree::update(int& step, UpdateType type) {
    try {
        DateTime currentDate = model->getDate(step); // the current date
        DateTime today = model->getToday();
        DateTimeArray const &exerDates = inst->sched->exerciseDate.getDates();
        DateTimeArray const &notifDates = inst->sched->notifDate.getDates();
        int nbExer = exerDates.size();

        int exerId;

        for (exerId=0; exerId<nbExer; ++exerId) {
            if (exerDates[exerId]==currentDate) {
                if (!undExer->isZero()) {
                    throw ModelException("Notification/Exercise periods overlap "
                        "(exercise #"+Format::toString(exerId)+")");
                }
                *undExer = undProd[exerId]->getValue(step, currentDate);
            }
        }
        // on notification date, update option price
        for (exerId=0; exerId<nbExer; ++exerId) {

            /* [begin] select which exercise we will be working on */

            if (notifDates[exerId]!=currentDate) 
                continue;

            if (currentDate < today)
                continue;

             /* [end] select which exercise we will be working on */

            double k1 = (*inst->strikes)[exerId];
            double k2 = (inst->strikesHi.get() ? (*inst->strikesHi)[exerId] : 0.);
            TreeSliceSP tmpValue = optionPrice->clone(false);
            TreeSliceSP zero;
            model->getZero(currentDate, exerDates[exerId], inst->discount->getName(), zero);
            TreeSlice &z = *zero;
            const TreeSlice &s = *undExer;

            switch (optionType) {
                case KOption::Type::CALL:
                case KOption::Type::CALLABLE:
                    *tmpValue = smax(s - k1*z, 0.);
                    *optionPrice = smax(*tmpValue, *optionPrice); 
                    break;
                case KOption::Type::PUT:
                case KOption::Type::PUTTABLE:
                    *tmpValue = smax(k1*z - s, 0.);
                    *optionPrice = smax(*tmpValue, *optionPrice); 
                    break;
                case KOption::Type::UNDERLYING:
                    break;
                case KOption::Type::CALL_SPREAD: 
                    *tmpValue = smin(s - k1*z, (k2 - k1)*z); 
                    *optionPrice = smax(*tmpValue, *optionPrice); 
                    break;
                case KOption::Type::PUT_SPREAD:  
                    *tmpValue = smin(k2*z - s, (k2 - k1)*z); 
                    *optionPrice = smax(*tmpValue, *optionPrice); 
                    break;
                case KOption::Type::STRANGLE:
                    *tmpValue = smax(k1*z - s, s - k2*z);
                    *optionPrice = smax(*tmpValue, *optionPrice); 
                    break;
                case KOption::Type::COLLAR: 
                    *tmpValue = smax(s - k2*z, 0.) 
                           - smax(k1*z - s, 0.);
                    *optionPrice = cond(*tmpValue > 0.,
                        smax(*tmpValue, *optionPrice),
                        smin(*tmpValue, *optionPrice)
                        );
                    break;
                case KOption::Type::BUTTERFLY: 
                    *tmpValue = smin(s - k1*z, k2*z - s);
                    *optionPrice = smax(*tmpValue, *optionPrice); 
                    break;
                default:
                    throw ModelException("optionType not supported");
            }
            *undExer = 0.; // undExer is now available to be reused for next exercise
            if (exerId==keepValueIdx) {
                *optSkipExerPrice = *optionPrice;
                startDEV(optSkipExerPrice);
            }
        }
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

const TreeSlice & KOptionTree::getValue(int step, DateTime eventDate) const {
    try {
        if (model->getDate(step)!=eventDate)
            throw ModelException("Does not support eventDate!=current date");
        // only support callable, puttable and underlying if timeStep = 0 at this stage
        // so it is a simple callable = underlying - option , puttable = underlying + option
        double los;
        if (inst->longOrShort == KOption::LongOrShort::OPTION_LONG)
            los = 1.0;
        else
            los = -1.0;
        
        switch (optionType) {
            case KOption::Type::CALLABLE: 
            case KOption::Type::PUTTABLE:
            case KOption::Type::UNDERLYING: {
                if (step != 0) {
                    throw ModelException( 
                        "CALLABLE/PUTTABLE/UNDERLYING only supported when component "
                        "calculates it as a deterministic value "
                        "ie. at timestep = 0.  Cannot be calculated at "
                        "tree timestep = " + Format::toString(step) + ", date = " + 
                        model->getDate(step).toString());
                }
                const TreeSlice &underlying = undProd[0]->getValue(step, eventDate);

                if (optionType==KOption::Type::CALLABLE) {
                    *getValuePrice = (underlying - *optionPrice)*los;
                }
                else if (optionType==KOption::Type::PUTTABLE) {
                    *getValuePrice = (underlying + *optionPrice)*los;
                }
                else *getValuePrice = underlying; // optionType is UNDERLYING

                return *getValuePrice;
                break;
            }
            default:
                *getValuePrice = *optionPrice * los;
                return *getValuePrice;
                break;
        }
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

//-----------------------------------------------------------------------------
//  KOption instrument
//-----------------------------------------------------------------------------
void KOption::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        /***** [begin] validity checks & calc transient */
        if (optionType == Type::CALL
        ||  optionType == Type::PUT) {
            // do not report KnownCashflows and PAYMENT_DATES for the underlyings
            undIsPhysical=false; 
        }
        // american
        if (exerType==ExerType::AMERICAN) {// exerType should be removed soon, american not supported
            throw ModelException("American style not supported");
        }

        if (isNotified)
            throw ModelException("isNotified flag is not supported for the moment");

        // discount
        if (discount.getName().empty()) 
            throw ModelException("\"discount\" not provided");

        // underlying
        if (und.get()) {
            if (undList.size())
                throw ModelException("Either \"und\" or \"undList\" must be provided, not both");
            undList.push_back(und);
            und.reset(0);
        }
        else if (undList.size()) {
            if (undList.size() != 1 
            && undList.size() != sched->exerciseDate.size()) {
                throw ModelException("\"undList\" must contain one entry or "
                "as many entries as there are coupons ("
                +Format::toString(sched->exerciseDate.size())+")");
            }
        }
        else throw ModelException("\"und\" or \"undList\" must be provided");

        if (undList.size()>1 
        && (optionType == Type::CALLABLE 
            || optionType == Type::PUTTABLE
            || optionType == Type::UNDERLYING)) {
            throw ModelException("Multiple underlying mode (undList.size()!=0)"
            " does not work with CALLABLE, PUTTABLE or UNDERLYING options");
        }

        // isNotified, dateExercised
        if (isNotified && dateExercised.empty()) {
            throw ModelException("isNotified but dateExercised is empty");
        }
        if (!dateExercised.empty()) {
            if (find(sched->exerciseDate.begin(), sched->exerciseDate.end(), dateExercised)
            == sched->exerciseDate.end()) {
                throw ModelException("dateExercised "+dateExercised.toString()
                +" is not part of the exercise schedule");
            }
        }

        // instSettle
        if (!instSettle.get()) {
            instSettle.reset(new PhysicalSettlement());
        }
        // exercise date schedule
        int nbExer = sched->exerciseDate.size();
        if (nbExer == 0) 
            throw ModelException("Exercise schedule is empty.");

        switch (optionType) {
            case Type::CALL:
            case Type::PUT:
            case Type::CALLABLE:
            case Type::PUTTABLE:
                if (strikesHi.get())
                    throw ModelException("strikeHi not needed for CALL, PUT, CALLABLE, PUTTABLE");
                break;
            case Type::CALL_SPREAD:
            case Type::PUT_SPREAD:
            case Type::STRANGLE:
            case Type::COLLAR:
            case Type::BUTTERFLY:
                if (!strikesHi)
                    throw ModelException("strikeHi needed for  CALL_SPREAD, PUT_SPREAD, STRANGLE, COLLAR");
                break;
            case Type::UNDERLYING:
                break; // do not care
            default: ASSERT(0);
        }
        if (strikes->size() != nbExer) {
            throw ModelException(
                "Number of strikes " + Format::toString(strikes->size()) +
                " must equal number of exercise/notification dates in schedule " +
                Format::toString(nbExer));
        }
        if (strikesHi.get()) {
            if (strikes->size() != strikesHi->size()) {
                throw ModelException("strikes->size() != strikesHi->size()");
            }
            for (int i=0; i<strikes->size(); ++i) {
                if ((*strikes)[i]>(*strikesHi)[i]) {
                    throw ModelException("strikes[i]>strikesHi[i], i="
                    +Format::toString(i));
                }
            }
        }

        if (sched->notifDate.size() != nbExer) {
            throw ModelException(
                "Number of notification dates " 
                + Format::toString(sched->notifDate.size()) +
                " must equal number of exercise dates in schedule " +
                Format::toString(nbExer));
        }

        // check that [notif;exercise] are ordered and non-overlapping periods
        // and that "dateExercised" is a valid exercise date
        bool dateExercisedFound=false;
        for (int i=0; i<nbExer; ++i) {
            if (sched->notifDate[i] > sched->exerciseDate[i]) {
                throw ModelException(
                    Format::toString(i+1) + "-th notification date ("
                    + sched->notifDate[i].toString()+") is after "
                    "the corresponding exercise date ("
                    + sched->exerciseDate[i].toString()+")");
            }
            if (i) {
                if (sched->notifDate[i] <= sched->exerciseDate[i-1]) {
                    throw ModelException(
                        Format::toString(i+1) + "-th notification date ("
                        + sched->notifDate[i].toString()+") is not after "
                        "the previous exercise date ("
                        + sched->exerciseDate[i-1].toString()+")");
                }
            }
            dateExercisedFound = dateExercisedFound || (sched->exerciseDate[i]==dateExercised);
        }
        if (isNotified && !dateExercisedFound) {
            throw ModelException("\"dateExercised\" "+dateExercised.toString()
                +" not found in the exercise date schedule");
        }
        /***** [end] validity checks & calc transient */
        {
            DateTimeArray const &exerciseDate = sched->exerciseDate.getDates();
            for (int i=0; i<undList.size(); ++i) {
                undList[i]->addResetDates(exerciseDate);
                undList[i]->setup(model, market);
            }
        }
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

bool KOption::isDead(DateTime valueDate, double *price) const {
    try {
        return false;
        if (isNotified) {
            if (valueDate > dateExercised) {
                if (price) 
                    *price = 0.;
                return true;
            }
            // have to price the underlying
            return false;
        }
        // here !isNotified
        if (getToday() >= sched->notifDate.back()) {
            // cannot be excercised anymore
            if (price) 
                *price = 0;
            return true;
        }
        // need a real pricing
        return false;
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

/** implement FDModel product interface */
FDProductSP KOption::createProduct(FDModel * model) const {
    try {
        if (!setupCalled) 
            throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
        return FDProductSP(new KOptionTree(KOptionConstSP(this), model));
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

static const string& typeEnumToCallabilityStr(KOption::Type::Enum x) {
    static string const composite = "!Composite"; // maybe add it to Callability
    switch (x) {
        case KOption::Type::CALL:     return Callability::CALL;
        case KOption::Type::PUT:      return Callability::PUT;
        case KOption::Type::CALLABLE: return Callability::CALLABLE;
        case KOption::Type::PUTTABLE: return Callability::PUTTABLE;
        default: return composite;
    }
}

// Callability::IEventHandler interface
void KOption::getEvents(const Callability*, 
                        IModel*            model, 
                        const DateTime&    eventDate,
                        EventResults*      events) const
{
    try {
        if (!isPhysical)
            return;
        // returns ALL callability dates
        // no comment is made on the desirability of calling
        for (int i = 0; i < sched->notifDate.size(); i++) {
            events->addEvent( new Callability(sched->notifDate[i], 
                sched->exerciseDate[i],
                typeEnumToCallabilityStr(optionType),
                (*strikes)[i]));
        }
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

KComponent::CashflowNodeSP KOption::reportCashflowsTree(bool amountsNeeded) const {
    try {
        CashflowNodeSP cfn(new CashflowNode);
        cfn->comp = KComponentConstSP(this);
        reportCashFlows(*cfn->cashflows, amountsNeeded);
        if (undIsPhysical) {
            KComponent *und = dynamic_cast<KComponent*>(undList[0].get());
            if (und) {
                CashflowNodeSP leaf = und->reportCashflowsTree(amountsNeeded);
                if (leaf.get()) {
                    cfn->underlyings.push_back(leaf);
                }
            }
        }
        return cfn;
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

void KOption::reportCashFlows(CashflowInfoArray &cashflowInfos, bool amountsNeeded ) const {
    try {
        cashflowInfos.clear();

        if (cfType == CashflowInfo::CfType::UNSET)
            return;

        int nbExer = sched->nbDates;
        // optional: gives a hint on how many cells will be needed in the array
        cashflowInfos.reserve(nbExer);

        // populate "cashflowInfos" array
        for (int i=0; i<nbExer; ++i) {
            CashflowInfoSP cfi(new CashflowInfo);

            cfi->componentName = outputName;
            cfi->date = sched->exerciseDate[i];
            cfi->amount = 0;
            cfi->cfType = cfType;
            cfi->amountType = CashflowInfo::AmountType::UNKNOWN;

            cashflowInfos.push_back(cfi);
        }
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

void KOption::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Option component");
    REGISTER(KOption, clazz)
    SUPERCLASS(KComponent);
    EMPTY_SHELL_METHOD(defaultConstructor);
    IMPLEMENTS(FDModel::IIntoProduct);
    IMPLEMENTS(Callability::IEventHandler);
    FIELD(und, "underlying index for option payoff (old way)")
    FIELD_MAKE_OPTIONAL(und);
    FIELD(undList, "underlying index for option payoff (new way)")
    FIELD_MAKE_OPTIONAL(undList);
    FIELD(optionType, "Type of option");
    FIELD(exerType, "Can exercise any time from first exercise date to last");
    FIELD(longOrShort, "Default value = OPTION_LONG");
    FIELD_MAKE_OPTIONAL(longOrShort);
    FIELD(smoothing, "Payoff smoothing flag");
    FIELD_MAKE_OPTIONAL(smoothing);
    FIELD(isNotified, "if exercise notification served");
    FIELD_MAKE_OPTIONAL(isNotified);
    FIELD(dateExercised, "date exercise takes effect");
    FIELD_MAKE_OPTIONAL(dateExercised);
    FIELD(levelExercised, "underlying level at which option has been exercised");
    FIELD_MAKE_OPTIONAL(levelExercised);
    FIELD(instSettle, "instrument settlement");
    FIELD_MAKE_OPTIONAL(instSettle);
    FIELD(sched, "Object containing option notification and exercise dates");
    FIELD(strikes, "Array of strikes values (one per date entry in the option schedule)");
    FIELD(strikesHi, "Array of high strikes values (one per date entry in the option schedule) if needed by optionType");
    FIELD_MAKE_OPTIONAL(strikesHi);
    FIELD(cfType, "Type of cashflow to report, if any wanted");
    FIELD_MAKE_OPTIONAL(cfType);
    Addin::registerConstructor(Addin::UTILITIES, KOption::TYPE);
}

CClassConstSP const KOption::TYPE = CClass::registerClassLoadMethod(
    "KOption", typeid(KOption), KOption::load);

START_PUBLIC_ENUM_DEFINITION(KOption::Type::Enum, "");
ENUM_VALUE_AND_NAME(KOption::Type::CALL, "CALL", "");
ENUM_VALUE_AND_NAME(KOption::Type::PUT, "PUT", "");
ENUM_VALUE_AND_NAME(KOption::Type::CALLABLE, "CALLABLE", "");
ENUM_VALUE_AND_NAME(KOption::Type::PUTTABLE, "PUTTABLE", "");
ENUM_VALUE_AND_NAME(KOption::Type::UNDERLYING, "UNDERLYING", "");
ENUM_VALUE_AND_NAME(KOption::Type::CALL_SPREAD, "CALL_SPREAD", "");
ENUM_VALUE_AND_NAME(KOption::Type::PUT_SPREAD, "PUT_SPREAD", "");
ENUM_VALUE_AND_NAME(KOption::Type::STRANGLE, "STRANGLE", "");
ENUM_VALUE_AND_NAME(KOption::Type::COLLAR, "COLLAR", "");
ENUM_VALUE_AND_NAME(KOption::Type::BUTTERFLY, "BUTTERFLY", "");
END_ENUM_DEFINITION(KOption::Type::Enum);


START_PUBLIC_ENUM_DEFINITION(KOption::ExerType::Enum, "");
ENUM_VALUE_AND_NAME(KOption::ExerType::EUROPEAN, "EUROPEAN", "");
ENUM_VALUE_AND_NAME(KOption::ExerType::AMERICAN, "AMERICAN", "");
END_ENUM_DEFINITION(KOption::ExerType::Enum);

START_PUBLIC_ENUM_DEFINITION(KOption::Smoothing::Enum, "");
ENUM_VALUE_AND_NAME(KOption::Smoothing::SMOOTHING_NONE, "SMOOTHING_NONE", "");
ENUM_VALUE_AND_NAME(KOption::Smoothing::SMOOTHING_RATES, "SMOOTHING_RATES", "");
END_ENUM_DEFINITION(KOption::Smoothing::Enum);

START_PUBLIC_ENUM_DEFINITION(KOption::LongOrShort::Enum, "");
ENUM_VALUE_AND_NAME(KOption::LongOrShort::OPTION_LONG, "OPTION_LONG", "");
ENUM_VALUE_AND_NAME(KOption::LongOrShort::OPTION_SHORT, "OPTION_SHORT", "");
END_ENUM_DEFINITION(KOption::LongOrShort::Enum);

bool KOptionLoad() {return KOption::TYPE !=0;}

DRLIB_END_NAMESPACE
