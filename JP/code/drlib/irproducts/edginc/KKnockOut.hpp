//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KKnockOut.hpp
//
//   Description : KnockOut component
//
//----------------------------------------------------------------------------
#ifndef _KKnockOut_HPP
#define _KKnockOut_HPP

#include "edginc/Schedule.hpp"
#include "edginc/KComponent.hpp"

DRLIB_BEGIN_NAMESPACE

const string KKnockOutBarrierValues();
    
class KKnockOut : public KComponent,
                  virtual public FDModel::IIntoProduct
{
public:
    // IN/OUT represent the normal state before the event occurs.
    // When used to monitor a knock-out, use ONE_IN (IN is the "normal" state until k-o)
    struct Barrier {
        enum Enum { ONE_IN, ONE_OUT, IN_IN, IN_OUT, OUT_IN, OUT_OUT, IN, OUT };
    };
    typedef BoxedEnum<Barrier::Enum> BarrierBoxedEnum;

    struct Smoothing {
        enum Enum {SMOOTHING_NONE, SMOOTHING_STEP};
    };
    typedef BoxedEnum<Smoothing::Enum> SmoothingBoxedEnum;
	
    /****************** KOEvent *************/
    class KOEvent : public CObject {
    public:
        static CClassConstSP const TYPE;

        friend class KKnockOut;
        friend class KKnockOutTree;

        IProdCreatorSP  und;             // underlying index
        double          idxWeight;
        ScheduleSP      loBarrier;
        ScheduleSP      hiBarrier;

        /**** methods ****/

    protected:
        KOEvent(CClassConstSP const &type) : CObject(type), idxWeight(1.0) {}

    private:
        static IObject* defaultConstructor(void) { return new KOEvent(TYPE); }
        static void load(CClassSP& clazz);
    };
    typedef smartPtr<KOEvent> KOEventSP;
    typedef array<KOEventSP, KOEvent> KOEventArray;
    /**********/

    /* KKnockOut public variables and methods */

    static CClassConstSP const TYPE;

    void validatePop2Object(void);
 
    /* IProdCreator:: */
	virtual void addResetDates(const DateTimeArray &resetDates);
    virtual void setup(const IModel* model, const MarketData* market);
    virtual double getValue(DateTime date, CashflowInfo &cfi ) const;

    static bool barrierIsIn(Barrier::Enum barrier, int idx);

protected:
    bool barrierBoolValue(int koIdx, int barIdx, CashflowInfo &cfiLoc, DateTime stepDate) const;


public:
    friend class KKnockOutTree;


    /****************** exported fields ************/
    KOEventArray    knockouts;
    Smoothing::Enum smoothing;
    Barrier::Enum   barrierType;
    bool            touchBarrierChangesState;
    // old fields   
    IProdCreatorSP  upUnd;  // value for FALSE for combined barriers (IN_IN,...)
    IProdCreatorSP  midUnd; // value for TRUE
    IProdCreatorSP  downUnd;
    // new fields   
    IProdCreatorSP  undFalse;  // value for FALSE
    IProdCreatorSP  undTrue;  // value for TRUE

    /****************** methods ************/

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KKnockOut(); }

    KKnockOut(void) : KComponent(TYPE), barrierType(Barrier::ONE_IN), touchBarrierChangesState(true) {}

    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;    

};

typedef smartPtr<KKnockOut> KKnockOutSP;

DRLIB_END_NAMESPACE

#endif


