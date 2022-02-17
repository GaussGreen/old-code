//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KKOSwap.hpp
//
//   Description : one touch KO swap component
//
//----------------------------------------------------------------------------
#ifndef _KKOSwap_HPP
#define _KKOSwap_HPP

#include "edginc/Schedule.hpp"
#include "edginc/KComponent.hpp"
#include "edginc/FlexDates.hpp"

DRLIB_BEGIN_NAMESPACE

class KKOSwap : public KComponent,
                virtual public FDModel::IIntoProduct
{
public:
    static CClassConstSP const TYPE;
    friend class KKOSwapTree;

    //virtual void validatePop2Object(void);

    bool isObsDate(const DateTime& date) const;

    // get ko level
    double getKOLevel(const DateTime& date) const;

protected:
    /****************** exported fields ************/
    IProdCreatorSP  underlier;
    IProdCreatorSP  obsUnd;
    FlexDates       obsDates;
    bool            isUpOut;
    ScheduleSP      barrier;

    /****************** methods ************/
    /* KComponent:: */
    virtual void setup(const IModel* model, const MarketData* market);
    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;

    /* IProdCreator */
    //virtual double getValue(DateTime date, CashflowInfo &cfi ) const;

    KKOSwap(CClassConstSP const &type = TYPE) : KComponent(type), isUpOut(true) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KKOSwap(TYPE); }
};

typedef smartPtr<KKOSwap> KKOSwapSP;
typedef smartConstPtr<KKOSwap> KKOSwapConstSP;

DRLIB_END_NAMESPACE

#endif


