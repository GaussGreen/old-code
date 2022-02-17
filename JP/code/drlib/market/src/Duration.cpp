//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : Duration.cpp
//
//   Description : Defines an interface for MarketObjects that can return
//                 a duration calculation and a ClientRunnable means of
//                 performing the calculation
//
//   Author      : Gordon Stephens
//
//   Date        : 12 April 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/Duration.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/VectorShift.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/Addin.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/RiskQuantityEvaluator.hpp"
#include "edginc/Control.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/Results.hpp"
#include "edginc/NamedRiskQuantity.hpp"

DRLIB_BEGIN_NAMESPACE

CClassConstSP const Duration::TYPE = 
    CClass::registerClassLoadMethod("Duration", typeid(Duration), load);

CClassConstSP const Duration::IParHandler::TYPE =
    CClass::registerInterfaceLoadMethod("Duration::IParHandler", 
                                        typeid(Duration::IParHandler), 
                                        load);

template<> CClassConstSP const MarketWrapper<Duration::IParHandler>::TYPE = 
    CClass::registerClassLoadMethod("MarketWrapper<Duration::IParHandler>", 
                                    typeid(MarketWrapper<Duration::IParHandler>),
                                    load);

CClassConstSP const Duration::IParHandlerWithClosedForm::TYPE =
    CClass::registerInterfaceLoadMethod(
              "Duration::IParHandlerWithClosedForm",
              typeid(Duration::IParHandlerWithClosedForm), 
              load);

template<> CClassConstSP const MarketWrapper<Duration::IParHandlerWithClosedForm>::TYPE = 
    CClass::registerClassLoadMethod(
          "MarketWrapper<Duration::IParHandlerWithClosedForm>", 
          typeid(MarketWrapper<Duration::IParHandlerWithClosedForm>),
          load);

CClassConstSP const Duration::IParHandlerWithoutClosedForm::TYPE =
    CClass::registerInterfaceLoadMethod(
              "Duration::IParHandlerWithoutClosedForm",
              typeid(Duration::IParHandlerWithoutClosedForm), 
              load);

CClassConstSP const Duration::IParHandlerWithoutClosedForm_VectorShift::TYPE =
    CClass::registerInterfaceLoadMethod(
              "Duration::IParHandlerWithoutClosedForm_VectorShift",
              typeid(Duration::IParHandlerWithoutClosedForm_VectorShift), 
              load);

template<> CClassConstSP const MarketWrapper<Duration::IParHandlerWithoutClosedForm>::TYPE = 
    CClass::registerClassLoadMethod(
          "MarketWrapper<Duration::IParHandlerWithoutClosedForm>", 
          typeid(MarketWrapper<Duration::IParHandlerWithoutClosedForm>),
          load);


// Returns the discount curve (if any)
YieldCurveWrapper Duration::getDiscount() const {
    return discount;
}

// Returns the benchmark expiries (can be an empty SP)
const ExpiryArrayConstSP Duration::getBenchmarks() const {
    return benchmarks;
}

// Returns the model
IModelSP Duration::getModel() const {
    return model;
}

// Returns the market
MarketDataSP Duration::getMarket() const {
    return market;
}


IObjectSP Duration::run() {
    static const string method = "Duration::run";

    try {
        //get the required class of market data, if specified
        CClassConstSP clazz = type.empty() ? IParHandler::TYPE : CClass::forName(type);
        //get data into the wrapper
        name.getData(model.get(), market, clazz);
        //and extract
        MarketObjectSP mosp(name.getMO());
        
	    //build storage for the duration results
        ExpiryResultArraySP durCalcs(new ExpiryResultArray(0));

        //get data into discount if necessary
        if (!discount.isEmpty()) {
            discount.getData(model.get(), market);
        }

        if (usedClosedForm) {
            // Check if we have an IParHandlerWithClosedForm, and if so get
            // the Duration 
            if (IParHandlerWithClosedForm::TYPE->isInstance(name.get())) {
                IParHandlerWithClosedForm* moPar = 
                    dynamic_cast<IParHandlerWithClosedForm*>(mosp.get());
                
                return moPar->getDuration(this);
            } 
            else {
                throw ModelException(method, 
                                     "Closed form duration requested but not supported");
            }
        } 
        else {
            // Check if we have an IParHandlerWithoutClosedForm, and if so get
            // the Duration 
            if (IParHandlerWithoutClosedForm::TYPE->isInstance(name.get())) {
                IParHandlerWithoutClosedForm* moPar = 
                    dynamic_cast<IParHandlerWithoutClosedForm*>(mosp.get());

                //get the points to tweak if not already specified (optional parameter)
                if (!benchmarks) {
                    benchmarks = moPar->getParBenchmarks();
                }

                //get the pointwise tweaker
                smartConstPtr<PerNameRiskPropertySensitivity<ExpiryWindow> > pTwk(
                    moPar->getPointwiseTweaker());

                //build up control
                SensitivityArraySP sens(new SensitivityArray(0));
                OutputRequestArraySP reqs(new OutputRequestArray(0));

                //for each tweakPoint
                for (int i=0; i<benchmarks->size(); i++) {
                    //get par instrument for this benchmark
                    ExpirySP bMark((*benchmarks)[i]);
                    InstrumentSP parImnt = moPar->getParInstrument(bMark);
                    //get the market data into the instrument
                    model->getInstrumentAndModelMarket(market.get(), parImnt.get());

                    RiskQuantityEvaluator::ValueArraySP values =
                        RiskQuantityEvaluator().values(
                            IRiskQuantityFactoryArray::SP(
                                1, pTwk->withQualifier(ExpiryWindow::SP(ExpirySP(),
                                                                        bMark,
                                                                        ExpirySP()))),
                            MultiTweakGroup::SP(
                                IInstrumentCollection::singleton(parImnt),
                                model)
                        );

                    if (values->size() != 1) {
                        throw ModelException(method, "unique tweak name not returned");
                    }

                    double duration = - (*values)[0]->value() *
                                        (*values)[0]->riskQuantity->unit *
                                        10000;  // due to tweak size of 1bp

                    //create result element
                    ExpiryResult bmarkDuration((*benchmarks)[i], duration);
                    //and store
                    durCalcs->push_back(bmarkDuration);
                }
            }

            /**
             * This is a stopgap to provide compatibility with RhoPointwise until we
             * get round to porting it to the IRiskQuantityFactory sensitivities
             * framework
             */

            //@{

            else if (IParHandlerWithoutClosedForm_VectorShift::TYPE->isInstance(name.get())) {
                IParHandlerWithoutClosedForm_VectorShift* moPar = 
                    dynamic_cast<IParHandlerWithoutClosedForm_VectorShift*>(mosp.get());

                //get the points to tweak if not already specified (optional parameter)
                if (!benchmarks) {
                    benchmarks = moPar->getParBenchmarks();
                }

                //get the pointwise tweaker
                VectorShiftSP pTwk(moPar->getPointwiseTweaker());

                //build up control
                SensitivityArraySP sens(new SensitivityArray(0));
                OutputRequestArraySP reqs(new OutputRequestArray(0));
                Control ctrl(sens, reqs, false, "");

                //Add the sensitivity
                ctrl.addSensitivity(pTwk);

                //for each tweakPoint
                for (int i=0; i<benchmarks->size(); i++) {
                    //get par instrument for this benchmark
                    ExpirySP bMark((*benchmarks)[i]);
                    InstrumentSP parImnt = moPar->getParInstrument(bMark);
                    //get the market data into the instrument
                    model->getInstrumentAndModelMarket(market.get(), parImnt.get());

                    pTwk->setExpiryToTweak(bMark);

                    //calculate price & sensitivities
                    CResults rmResults;
                    ctrl.calculate(model.get(), parImnt.get(), &rmResults);

                    //determine the name used for this sensitivity result
                    OutputNameArrayConstSP names(pTwk->names(mosp.get()));
                    //ensure that names returned 1 and only 1 value
                    if (names->size() != 1) {
                        throw ModelException(method, "unique tweak name not returned");
                    }
                    //extract pTwk results & convert to useable form
                    IObjectConstSP durResults(rmResults.retrieveGreek(
                        pTwk->getSensOutputName(), 
                        (*names)[0]));

                    const ExpiryResultArray* durVector = 
                        dynamic_cast<const ExpiryResultArray*>(durResults.get());

                    if (!durVector || durVector->size() > 1) {
                        throw ModelException(method, 
                                             "failed to retrieve tweak results correctly");
                    }

                    double duration = ((*durVector)[0]).getResult();

                    //scaled by 10000 due to tweak size of 1bp
                    duration = -duration * 10000;

                    //create result element
                    ExpiryResult bmarkDuration((*benchmarks)[i],duration);
                    //and store
                    durCalcs->push_back(bmarkDuration);
                }
            }

            //@}

            else {
                throw ModelException(method, 
                                     "Non-closed form duration requested but not supported");
            }
        }

        return durCalcs;
    }
    catch(exception& e) {
        throw ModelException(e, method);
    }
}

void Duration::load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(Duration, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultDuration);
        FIELD       (name, "the market object to get durations from");
        FIELD              (market, "the market cache containing mktObject");
        FIELD              (model, "Model used to price par instrument / "
                            "retrieve market data (default value is ClosedForm)");
        FIELD_MAKE_OPTIONAL(model);
        FIELD       (discount, "discount curve");
        FIELD_MAKE_OPTIONAL(discount);
        FIELD              (benchmarks, 
                            "the benchmarks for which durations are required");
        FIELD_MAKE_OPTIONAL(benchmarks);
        FIELD       (type, "the explicit type of mktWrapper");
        FIELD_MAKE_OPTIONAL(type);
        FIELD       (usedClosedForm, 
                            "return a closed form result if the data supports it");
        FIELD_MAKE_OPTIONAL(usedClosedForm);
        Addin::registerObjectMethod("DURATION",
                                    Addin::RISK,
                                   "Returns term structure of durations",
                                    false,
                                    Addin::expandMulti,
                                    &Duration::run);
}


void Duration::IParHandler::load(CClassSP& clazz) {
    REGISTER_INTERFACE(Duration::IParHandler, clazz);
    EXTENDS(IObject);
}

void Duration::IParHandlerWithClosedForm::load(CClassSP& clazz) {
    REGISTER_INTERFACE(Duration::IParHandlerWithClosedForm, clazz);
    EXTENDS(IParHandler);
}

void Duration::IParHandlerWithoutClosedForm::load(CClassSP& clazz) {
    REGISTER_INTERFACE(Duration::IParHandlerWithoutClosedForm, clazz);
    EXTENDS(IParHandler);
}

void Duration::IParHandlerWithoutClosedForm_VectorShift::load(CClassSP& clazz) {
    REGISTER_INTERFACE(Duration::IParHandlerWithoutClosedForm_VectorShift, clazz);
    EXTENDS(IParHandler);
}

IObject* Duration::defaultDuration() {
    return new Duration();
}

Duration::Duration(): CObject(TYPE), 
                      model(IModelSP(new ClosedForm())), 
                      usedClosedForm(false) 
{}

DRLIB_END_NAMESPACE
