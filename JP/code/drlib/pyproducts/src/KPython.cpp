//----------------------------------------------------------------------------
//
//   Group       : Global QR
//
//   Filename    : KPython.cpp
//   Date        : 1 Nov 2006
//   Description : python flexible payoff component for tree/fd
//----------------------------------------------------------------------------

#include "edginc/PyInterface.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
//  KPythonTree
//  FDModel's pricing product/calc 
//  for the KPython component
//-----------------------------------------------------------------------------
KPythonTree::KPythonTree(const KPythonConstSP &inst, FDModel* model) :
    FDProduct(model), inst(inst)
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

        //if (inst->srcCode.size()==0)
        //    throw ModelException("no python src code supplied");

        useSrcCode = inst->srcCode.size()>0;
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

void KPythonTree::recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
    try {
        string ccy = inst->discount->getCcy();

        //!!! const TreeSlice &price = getValue(0);

        //!!! recordSliceToOutputName(ctrl, results, model, inst->isTopComponent(), 
        //!!! inst->outputName, ccy, price);

        // record actual option price - != KPython price if callable/puttable/underlying
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

void KPythonTree::init(Control*) const{
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

void KPythonTree::initProd(void){
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
void KPythonTree::update(int& step, UpdateType type) {
    try {
        if (useSrcCode){ // add ! when we actually store code text, only for testing performance
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

                // [begin] select which exercise we will be working on

                if (notifDates[exerId]!=currentDate) 
                    continue;

                if (currentDate < today)
                    continue;

                // [end] select which exercise we will be working on

                double k1 = (*inst->strikes)[exerId];
                TreeSliceSP zero;
                model->getZero(currentDate, exerDates[exerId], inst->discount->getName(), zero);
                TreeSlice &z = *zero;
                const TreeSlice &s = *undExer;

                *optionPrice = smax(s - k1*z, *optionPrice); 

                *undExer = 0.; // undExer is now available to be reused for next exercise
                if (exerId==keepValueIdx) {
                    *optSkipExerPrice = *optionPrice;
                    startDEV(optSkipExerPrice);
                }
            }
        }
        else{

            try {
                pyInterfaceInit();
                boost::python::object module( boost::python::handle<>( PyImport_ImportModule( "KPython" ) ) );
                boost::python::object updateFunction = module.attr( "update" );
                KPythonTreePy treePy( *this );
                updateFunction( treePy, step, type );
            } catch ( boost::python::error_already_set ) {
                PyErr_Print();
            }
        }
    } catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

const TreeSlice & KPythonTree::getValue(int step) const {
      return *getValuePrice;
}

//-----------------------------------------------------------------------------
//  KPython instrument
//-----------------------------------------------------------------------------
void KPython::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        /***** [begin] validity checks & calc transient */
        // discount
        if (discount.getName().empty()) 
            throw ModelException("\"discount\" not provided");

        // underlying
        if (undList.size()==0) 
            throw ModelException("\"und\" or \"undList\" must be provided");

        // instSettle
        if (!instSettle.get()) {
            instSettle.reset(new PhysicalSettlement());
        }
        // exercise date schedule
        int nbExer = sched->exerciseDate.size();
        if (nbExer == 0) 
            throw ModelException("Exercise schedule is empty.");

        if (strikes->size() != nbExer) {
            throw ModelException(
                "Number of strikes " + Format::toString(strikes->size()) +
                " must equal number of exercise/notification dates in schedule " +
                Format::toString(nbExer));
        }
        if (sched->notifDate.size() != nbExer) {
            throw ModelException(
                "Number of notification dates " 
                + Format::toString(sched->notifDate.size()) +
                " must equal number of exercise dates in schedule " +
                Format::toString(nbExer));
        }

        // check that [notif;exercise] are ordered and non-overlapping periods
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

/** implement FDModel product interface */
FDProductSP KPython::createProduct(FDModel * model) const {
    try {
        if (!setupCalled) 
            throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
        return FDProductSP(new KPythonTree(KPythonConstSP(this), model));
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

// Callability::IEventHandler interface
void KPython::getEvents(const Callability*, 
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
                "CALL", // dummy !!!
                (*strikes)[i]));
        }
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

void KPython::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Option component");
    REGISTER(KPython, clazz)
    SUPERCLASS(KComponent);
    EMPTY_SHELL_METHOD(defaultConstructor);
    IMPLEMENTS(FDModel::IIntoProduct);
    IMPLEMENTS(Callability::IEventHandler);
    FIELD(undList, "underlying index for option payoff (new way)")
    FIELD_MAKE_OPTIONAL(undList);
    FIELD(instSettle, "instrument settlement");
    FIELD_MAKE_OPTIONAL(instSettle);
    FIELD(sched, "Object containing option notification and exercise dates");
    FIELD(strikes, "Array of strikes values (one per date entry in the option schedule)");
    FIELD(srcCode,     "source code text");
    FIELD_MAKE_OPTIONAL(srcCode);
    Addin::registerConstructor(Addin::UTILITIES, KPython::TYPE);
}

CClassConstSP const KPython::TYPE = CClass::registerClassLoadMethod(
    "KPython", typeid(KPython), KPython::load);

/******************************/
// for type linking
bool PYPRODUCTS_DLL KPythonLoad(){
    return (KPython::TYPE != 0);
}

DRLIB_END_NAMESPACE
