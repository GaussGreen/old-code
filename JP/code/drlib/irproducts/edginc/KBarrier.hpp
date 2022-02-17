//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Description : Barrier component
//
//----------------------------------------------------------------------------
#ifndef _KBarrier_HPP
#define _KBarrier_HPP

#include "edginc/Schedule.hpp"
#include "edginc/KComponent.hpp"
#include "edginc/IIndicCreator.hpp"

DRLIB_BEGIN_NAMESPACE

class KBarrier : public KComponent,
                  virtual public IIndicCreator,
                  virtual public FDModel::IIntoProduct
{
public:
    static CClassConstSP const TYPE;

    struct Smoothing {
        enum Enum {NONE, STEP};
    };
    typedef BoxedEnum<Smoothing::Enum> SmoothingBoxedEnum;
	
    /* IProdCreator:: */
	virtual void addResetDates(const DateTimeArray &resetDates);
    virtual void setup(const IModel* model, const MarketData* market);
    virtual double getValue(DateTime date, CashflowInfo &cfi ) const;

    virtual FDProductSP createProduct(FDModel * model) const;    

protected:
    bool barrierBoolValue(int koIdx, int barIdx, CashflowInfo &cfiLoc, DateTime currentDate) const;

public:
    friend class KBarrierTree;

    /****************** exported fields ************/

    IProdCreatorSP  und;             // underlying index
    ScheduleSP      levels;          // low barrier
    ScheduleSP      levelsHi;        // optional high barrier
    bool            onBarrierIsUp;
    bool            upIsTrue;
    Smoothing::Enum smoothing;

    /****************** methods ************/

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KBarrier(); }

    KBarrier(void) : KComponent(TYPE), onBarrierIsUp(true), upIsTrue(true), smoothing(Smoothing::NONE) {}
};
DECLARE(KBarrier);

DRLIB_END_NAMESPACE

#endif


