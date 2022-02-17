//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KKnockIO.hpp
//
//   Description : one touch KO swap component
//
//----------------------------------------------------------------------------
#ifndef _KKnockIO_HPP
#define _KKnockIO_HPP

#include "edginc/KComponent.hpp"
#include "edginc/FDModel.hpp"
#include "edginc/OptionSched.hpp"

DRLIB_BEGIN_NAMESPACE

/*
This KKnockIO components can work in two different modes:
1) with one underlying
Knock-In, you get the underlying as soon as the indicator reaches 1
Knock-Out is: underlying minus knock-in.
2) with more than one underlying (one each exercise date)
Knock-In gets all the underlyings from the date the indicator reaches 1
Knock-Out gets all the underlyings until the date the indicator reaches 1
*/

class KKnockIO : public KComponent,
                 virtual public FDModel::IIntoProduct
{
public:
    static CClassConstSP const TYPE;
    friend class KKnockIOTree;

    /****************** exported fields ************/
public:
    struct KnockIOType {
        enum Enum {K_IN, K_OUT};
    };
    typedef BoxedEnum<KnockIOType::Enum> KnockIOTypeBoxedEnum;

protected:
    KnockIOType::Enum knockIOType;
    IProdCreatorArray indicators;
    IProdCreatorArray underlyings;
    IProdCreatorArray rebates;
    OptionSchedDatesSP sched;

    /****************** methods ************/
    /* KComponent:: */
    virtual void setup(const IModel* model, const MarketData* market);
    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;

    /* IProdCreator */
    //virtual double getValue(DateTime date, CashflowInfo &cfi ) const;

    KKnockIO(CClassConstSP const &type = TYPE) : KComponent(type) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KKnockIO(TYPE); }
};

typedef smartPtr<KKnockIO> KKnockIOSP;
typedef smartConstPtr<KKnockIO> KKnockIOConstSP;

DRLIB_END_NAMESPACE

#endif


