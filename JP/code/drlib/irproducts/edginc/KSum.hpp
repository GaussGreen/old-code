//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KSum.hpp
//
//   Description : sum compoenent
//
//----------------------------------------------------------------------------
#ifndef _KSUM_HPP
#define _KSUM_HPP

#include "edginc/Schedule.hpp"
#include "edginc/KComponent.hpp"
#include "edginc/IndexSpecFX.hpp"

DRLIB_BEGIN_NAMESPACE

class KSum : public KComponent,
             virtual public FDModel::IIntoProduct
{
public:
    static CClassConstSP const TYPE;
    friend class KSumTree;

    /****************** exported fields ************/
protected:
    CModel::IProdCreatorArray listK;
    DoubleArray weights;
    double constant;
    ScheduleSP constSchedule;

    /****************** transient fields ************/
	IndexSpecFXArray fx;
    DateTimeArray resetDates;

    /****************** methods ************/
public:
    void add(IProdCreatorSP el, double weight);
    KSum() : KComponent(TYPE), constant(0.0) {}

    virtual void validatePop2Object(void);
	virtual void addResetDates(const DateTimeArray &resetDates);
    CModel::IProdCreatorArray getListK() const {return listK;}
    KComponent::CashflowNodeSP reportCashflowsTree(bool amountsNeeded) const;

protected:
    /* KComponent:: */
     virtual void setup(const IModel* model, const MarketData* market);
    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;

    /* IProdCreator */
    virtual double getValue(DateTime date, CashflowInfo &cfi ) const;

    KSum(CClassConstSP const &type) : KComponent(type), constant(0.0) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KSum(TYPE); }
};

typedef smartPtr<KSum> KSumSP;
typedef smartConstPtr<KSum> KSumConstSP;

DRLIB_END_NAMESPACE

#endif


