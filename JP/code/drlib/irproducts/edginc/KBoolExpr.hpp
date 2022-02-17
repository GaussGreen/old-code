//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Description : Simple boolean expression (!)A (&& (!)B)
//
//----------------------------------------------------------------------------
#ifndef _KBOOLEXPR_HPP
#define _KBOOLEXPR_HPP

#include "edginc/KComponent.hpp"
#include "edginc/IIndicCreator.hpp"
#include "edginc/Schedule.hpp"

DRLIB_BEGIN_NAMESPACE

class KBoolExpr : public KComponent,
                  virtual public IIndicCreator,
                  virtual public FDModel::IIntoProduct
{
public:
    static CClassConstSP const TYPE;
    friend class KBoolExprTree;

    /****************** exported fields ************/
protected:
    IIndicCreatorSP arg1, arg2;
    bool not1, not2;

    /****************** methods ************/
public:
    /* KComponent:: */
     virtual void setup(const IModel* model, const MarketData* market);
    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;

    /* IProdCreator */
    virtual double getValue(DateTime date, CashflowInfo &cfi ) const;
	virtual void addResetDates(const DateTimeArray &resetDates);

    KBoolExpr(CClassConstSP const &type) : KComponent(type), not1(true), not2(true) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KBoolExpr(TYPE); }
};
DECLARE(KBoolExpr);

DRLIB_END_NAMESPACE

#endif


