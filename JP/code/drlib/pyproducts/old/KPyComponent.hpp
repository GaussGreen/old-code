//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KPyComponent.hpp
//
//   Description : General python component where a serialized python 
//                 class contains the functor/calc
//
//----------------------------------------------------------------------------
#ifndef _KPYCOMPONENT_HPP
#define _KPYCOMPONENT_HPP

#include "edginc/Schedule.hpp"
#include "edginc/KComponent.hpp"

DRLIB_BEGIN_NAMESPACE

class KPyComponent : public KComponent,
                     virtual public FDModel::IIntoProduct
{
public:
    static CClassConstSP const TYPE;
    friend class KPyComponentTree;

    /****************** exported fields ************/
protected:
    //ScheduleSP schedule;  // optional schedule
    string serializedPyObject;  // "pickled" serialized python object or maybe fileName of script source
    CModel::IProdCreatorArray underlyings;
    

    /****************** transient fields ************/
    DateTimeArray resetDates;  // propagate any dates from parent

    /****************** methods ************/
public:

    KPyComponent() : KComponent(TYPE) {}

    virtual void validatePop2Object(void);
    virtual void addResetDates(const DateTimeArray &resetDates);

protected:
    /* KComponent:: */
     virtual void setup(const IModel* model, const MarketData* market);
    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;

    virtual double getValue(DateTime date, CashflowInfo &cfi ) const;

    KPyComponent(CClassConstSP const &type) : KComponent(type) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KPyComponent(TYPE); }
};

typedef smartPtr<KPyComponent> KPyComponentSP;
typedef smartConstPtr<KPyComponent> KPyComponentConstSP;

DRLIB_END_NAMESPACE

#endif


