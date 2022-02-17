//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KOption.hpp
//
//   Description : option komponent
//
//----------------------------------------------------------------------------
#ifndef _KOPTION_HPP
#define _KOPTION_HPP

#include "edginc/KComponent.hpp"
#include "edginc/OptionSched.hpp"
#include "edginc/InstrumentSettlement.hpp"

DRLIB_BEGIN_NAMESPACE

class KOption : public KComponent,
                virtual public FDModel::IIntoProduct,
                virtual public Callability::IEventHandler,
                virtual public DeadPricer
{
public:
    static CClassConstSP const TYPE;
    friend class KOptionTree;

    struct Type {
        enum Enum {CALL, PUT, CALLABLE, PUTTABLE, UNDERLYING, CALL_SPREAD, PUT_SPREAD, STRANGLE, COLLAR, BUTTERFLY };
    };
    typedef BoxedEnum<Type::Enum> TypeBoxedEnum;

    struct ExerType {
        enum Enum {EUROPEAN, AMERICAN};
    };
    typedef BoxedEnum<ExerType::Enum> ExerTypeBoxedEnum;

    struct Smoothing {
        enum Enum {SMOOTHING_NONE, SMOOTHING_RATES};
    };
    typedef BoxedEnum<Smoothing::Enum> SmoothingBoxedEnum;

    struct LongOrShort {
        enum Enum {OPTION_LONG, OPTION_SHORT};
    };
    typedef BoxedEnum<LongOrShort::Enum> LongOrShortBoxedEnum;

    /****************** exported fields ************/
public: // can be modified by shell instruments
    IProdCreatorSP  und;            // underlying index
    IProdCreatorArray undList;      // list of underlying indexes
protected:
    Type::Enum      optionType;     // call, put
    ExerType::Enum  exerType;       // American or European
    Smoothing::Enum   smoothing;      // currently only used by rates products.
                                    // defaults to SMOOTHING_NONE.
    LongOrShort::Enum longOrShort;    // multiplies price by 1 (long) or short (-1)
    bool            isNotified;     // if exercise notification served
    DateTime        dateExercised;  // settlement date of exercise can be derived from this
    double          levelExercised; // value of the underlying at exercise

    OptionSchedDatesSP sched;
    DoubleArraySP      strikes;
    DoubleArraySP      strikesHi;
    InstrumentSettlementSP  instSettle;       // instrument settlement

    CashflowInfo::CfType::Enum cfType;

    /****************** methods ************/
public:

    /* Callability::IEventHandler */
    void getEvents( const Callability*, IModel* model, 
        const DateTime&  eventDate, EventResults* events) const;

protected:
    /* KComponent:: */
    virtual void setup(const IModel* model, const MarketData* market);
    virtual void reportCashFlows(CashflowInfoArray &cashflowInfos, bool amountsNeeded ) const;
    virtual CashflowNodeSP reportCashflowsTree(bool amountsNeeded) const;

    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;
    /* DeadPricer:: */
    virtual bool isDead(DateTime valueDate, double *price) const;

    KOption(CClassConstSP const &type) 
    :   KComponent(type), smoothing(Smoothing::SMOOTHING_NONE), longOrShort(LongOrShort::OPTION_LONG),
        isNotified(false), dateExercised(DateTime(0,0)),  levelExercised(0.0),
        cfType(CashflowInfo::CfType::UNSET) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KOption(TYPE); }
};

typedef smartPtr<KOption> KOptionSP;
typedef smartConstPtr<KOption> KOptionConstSP;

class KOptionTree : public FDProduct {
    /************************ variables ************************/
private:
    KOptionConstSP    inst;
    KOption::Type::Enum optionType; // copied from inst->optionType for convenience

protected:
    FDProductArray    undProd;
    bool              isNotified; // copied from inst->isNotified for convenience
    int               keepValueIdx; // index of the exercice at which the value of 
                                    // option must be saved in optSkipExerPrice

public:
    /************************ methods ************************/
    KOptionTree(const KOptionConstSP &inst, FDModel* model);

    virtual void init(Control*) const;
    virtual void initProd();
    virtual void update(int& step, UpdateType type);
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const;
    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results);
    virtual string getOutputName(void) { return inst->outputName; }

protected:
    TreeSliceSP optionPrice;
    TreeSliceSP optSkipExerPrice;
    TreeSliceSP getValuePrice;    // used to return results from getValue function
    TreeSliceSP undExer;          // value of the underlying at exercise date
};

DRLIB_END_NAMESPACE

#endif


