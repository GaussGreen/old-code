#ifndef QLIB_KPYTHON_HPP
#define QLIB_KPYTHON_HPP

#include "edginc/config.hpp"
#include "edginc/KComponent.hpp"
#include "edginc/OptionSched.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/PhysicalSettlement.hpp"
#include "edginc/EventResults.hpp"

DRLIB_BEGIN_NAMESPACE

class PYPRODUCTS_DLL KPython : public KComponent,
                               virtual public FDModel::IIntoProduct,
                               virtual public Callability::IEventHandler
{
public:
    friend class KPythonTreePy;
    static CClassConstSP const TYPE;
    friend class KPythonTree;

    /* Callability::IEventHandler */
    void getEvents( const Callability*, IModel* model, 
        const DateTime&  eventDate, EventResults* events) const;

    /****************** exported fields ************/
    // can be modified by shell instruments
    IProdCreatorArray undList;      // list of underlying indexes
protected:
    OptionSchedDatesSP sched;
    DoubleArraySP      strikes;
    InstrumentSettlementSP  instSettle;       // instrument settlement

    CStringArray      srcCode; // python source code

    /* KComponent:: */
    virtual void setup(const IModel* model, const MarketData* market);

    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;

    KPython(CClassConstSP const &type) : KComponent(type){}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KPython(TYPE); }
};

typedef smartPtr<KPython> KPythonSP;
typedef smartConstPtr<KPython> KPythonConstSP;

class KPythonTree : public FDProduct {
    /************************ variables ************************/
private:
    friend class KPythonTreePy;
    KPythonConstSP    inst;

protected:
    FDProductArray    undProd;
    bool              isNotified; // copied from inst->isNotified for convenience
    int               keepValueIdx; // index of the exercice at which the value of 
                                    // option must be saved in optSkipExerPrice
    TreeSliceSP optionPrice;
    TreeSliceSP optSkipExerPrice;
    TreeSliceSP getValuePrice;    // used to return results from getValue function
    TreeSliceSP undExer;          // value of the underlying at exercise date

    bool        useSrcCode; // for testing performance, to be removed

public:
    /************************ methods ************************/
    KPythonTree(const KPythonConstSP &inst, FDModel* model);

    virtual void init(Control*) const;
    virtual void initProd();
    virtual void update(int& step, UpdateType type);
    virtual const TreeSlice & getValue(int step) const;
    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results);
    virtual string getOutputName(void) { return inst->outputName; }
};

typedef smartPtr<KPythonTree> KPythonTreeSP;
typedef smartConstPtr<KPythonTree> KPythonTreeConstSP;

DRLIB_END_NAMESPACE

#endif