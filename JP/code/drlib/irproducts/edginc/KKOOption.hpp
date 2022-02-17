//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KKOOption.hpp
//
//   Description : KO option component - KKOOption, a subclass of KOption
//
//----------------------------------------------------------------------------
#ifndef _KKOOPTION_HPP
#define _KKOOPTION_HPP

#include "edginc/Schedule.hpp"
#include "edginc/KOption.hpp"

DRLIB_BEGIN_NAMESPACE


class KKOOption : public KOption {
public:
    static CClassConstSP const TYPE;
    friend class KKOOptionTree;

    bool isObsDate(const DateTime& date) const;

    // get ko level
    double getKOLevel(const DateTime& date) const;

    /****************** exported fields ************/
protected:
    IProdCreatorSP  obsUnd;
    FlexDates       obsDates;
    bool            isUpOut;
    ScheduleSP      barrier;
    bool            isKnockIn;

    /****************** methods ************/
    
protected:
    /* KComponent:: */
    virtual void setup(const IModel* model, const MarketData* market);
    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;
    /* DeadPricer:: comment out for now*/
    //virtual bool isDead(DateTime valueDate, double *price) const;

    KKOOption(CClassConstSP const &type) 
        : KOption(type), isUpOut(true), isKnockIn(false) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KKOOption(TYPE); }
};

typedef smartPtr<KKOOption> KKOOptionSP;
typedef smartConstPtr<KKOOption> KKOOptionConstSP;

DRLIB_END_NAMESPACE

#endif


