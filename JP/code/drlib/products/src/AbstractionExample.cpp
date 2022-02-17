//----------------------------------------------------------------------------
//
//   Group       : EDG Quant Research
//
//   Filename    : AbstractionExample.cpp
//
//   Description : Tests Pyramid client support for abstraction
//
//   Author      : Andrew J Swain
//
//   Date        : 5 June 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GenericNFactor.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/FRIfaces.hpp"
#include "edginc/ObservationBuilder.hpp"
#include "edginc/Results.hpp"

DRLIB_BEGIN_NAMESPACE

class AbstractionExample: public GenericNFactor, 
                          virtual public CClosedFormLN::IIntoProduct {
public:
    static CClassConstSP const TYPE;

    class Subtempl: public CObject {
    public:
        static CClassConstSP const TYPE;  

        // data
        IObservationBuilderSP  builder;   // has an abstraction inside
        InstrumentSettlementSP settles;   // is an abstraction

        Subtempl():CObject(TYPE) {}; 
   
    private:
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz){
            clazz->setPublic(); // make visible to EAS/spreadsheet
            REGISTER(AbstractionExample::Subtempl, clazz);
            SUPERCLASS(CObject);
            EMPTY_SHELL_METHOD(defaultSubtemplate);
            FIELD(builder, "date builder");
            FIELD(settles, "settles");
        }
        
        static IObject* defaultSubtemplate(){
            return new AbstractionExample::Subtempl();
        }
    };

private:
    /// fields ////////

    // 2 day counts (empty classes typically) at top level
    DayCountConventionSP     fixDCC;
    DayCountConventionSP     floatDCC;

    // a "subtemplate" with abstractions
    Subtempl                 subtemplate;

    // an array of abstractions - one per equity so can use "Equity" template
    FRIfaces::ILValueExpressionArraySP variables;

public:
    friend class AbsEgClosedForm;

    // validation
    void validatePop2Object(){
        if (assets->NbAssets() != variables->size()) {
            throw ModelException("AbstractionExample::validatePop2Object",
                                 "need one variable entry per asset");
        }
    }
   
    /** Implementation of ClosedFormLN::IntoProduct interface - the
        implementation of this is below */
    virtual CClosedFormLN::IProduct* createProduct(
        CClosedFormLN* model) const;

    void price(CResults* results) const {
        static const string method = "AbstractionExample::price";
        try {
            double value = 42.0;
            value *= notional;
            results->storePrice(value, discount->getCcy());
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }


private:
    AbstractionExample(): GenericNFactor(TYPE) {} // for reflection
    AbstractionExample(const AbstractionExample& rhs); // not implemented
    AbstractionExample& operator=(const AbstractionExample& rhs); // not implemented

    static IObject* defaultAbstractionExample(){
        return new AbstractionExample();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AbstractionExample, clazz);
        SUPERCLASS(GenericNFactor);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        EMPTY_SHELL_METHOD(defaultAbstractionExample);
        FIELD(fixDCC,             "fixDCC");
        FIELD(floatDCC,           "floatDCC");
        FIELD(subtemplate, "subtemplate");
        FIELD(variables,          "variables");
    }
};

/** private class */
class AbsEgClosedForm: public CClosedFormLN::IProduct{
private:
    const AbstractionExample* ae; // a reference

public:
    AbsEgClosedForm(const AbstractionExample* ae): ae(ae){}

    void price(CClosedFormLN* model,
               Control*       control, 
               CResults*      results) const{
        ae->price(results);
    }
};


CClosedFormLN::IProduct* AbstractionExample::createProduct(
    CClosedFormLN* model) const {
    return new AbsEgClosedForm(this);
}

CClassConstSP const AbstractionExample::TYPE = CClass::registerClassLoadMethod(
    "AbstractionExample", typeid(AbstractionExample), AbstractionExample::load);

CClassConstSP const AbstractionExample::Subtempl::TYPE = CClass::registerClassLoadMethod(
    "AbstractionExample::Subtempl", 
    typeid(AbstractionExample::Subtempl), 
    AbstractionExample::Subtempl::load);


// * for class loading (avoid having header file) */
bool AbstractionExampleLoad() {
    return (AbstractionExample::TYPE != 0);
}

DRLIB_END_NAMESPACE

