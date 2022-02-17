//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KTurb.hpp
//
//   Description : turbo leg component
// 
//   Generic Model Turbo interface/data structure
// 
//   This instrument interface class is used directly in 
//   pricing by the models
//    - all other IR cashflow legs should construct this
//      structure for pricing using different/simplified 
//      constructors if required
// 
//----------------------------------------------------------------------------
#ifndef KTURBO_HPP
#define KTURBO_HPP

#include "edginc/config.hpp"
#include "edginc/KComponent.hpp"
#include "edginc/CouponSched.hpp"
#include "edginc/IndexSpecFX.hpp"

DRLIB_BEGIN_NAMESPACE

class KTurbo : public KComponent,
               virtual public FDModel::IIntoProduct
{
public:
	class Leg : public CObject {
    public:
		static CClassConstSP const TYPE;
        friend class KTurbo;
        friend class KTurboTree;

        /************ exported fields ************/
    public:
        YieldCurveWrapper  discount;      /**< Currency of the leg */
    protected:
        DoubleArraySP  notionals; /**< Notional/coupon (in currency) */
        DoubleArraySP  spreads;   /**< Fixed spread/coupon (as %) */
        DoubleArraySP  weights;   /**< Weight/coupon for floating index */
        
        bool payInitial;         /**< Pay initial notional */
        bool payPrincipal;       /**< Pay final notional */
        bool payExercise;        /**< Pay principals on exercise event */

        RateType::Enum rateType;   /**< Simple/Continuous rate */
        DayCountConventionSP dcc;/**< day count convention for fwd turbo stub calculation */
        CModel::IProdCreatorArray indexes;

        /******** transient ********/
        DoubleArray  dcfs;        /**< Day count fraction/coupon */
        IndexSpecFXSP fx;

        /******************** methods *******************/
    public:
        void validatePop2Object(void);

    protected:
        Leg(CClassConstSP const &type) 
            : CObject(type), payInitial(false), payPrincipal(false), 
            payExercise(false) {}

    private:
        static IObject* defaultConstructor(void) { return new Leg(TYPE); }
        static void load(CClassSP& clazz);
    };
    typedef smartPtr<Leg> LegSP;
    typedef array<LegSP, Leg> LegArray;

    /****************************** variables ********************************/
    static CClassConstSP const TYPE;
    friend class KTurboTree;

public: // can be modified by shell instruments
    LegArray legs;

protected:
    /******** swap schedule ********/
    CouponSchedDatesSP sched;

    /* cap/floor/other fields.
       must be one column/coupon.
       Floor = strikeRate[0][...], Cap = strikeRate[1][...], others vary */
    DoubleArraySP  floor, cap, notional;
    DayCountConventionSP dcc;
    
    /******** transient fields ********/
    DoubleArray dcfs;
    
    /******** methods ********/
public:
    virtual DateTime getLastDate(void) const;

    /* KComponent:: */
    void reportEvents(const KnownCashflows*, IModel* model,
    	    const DateTime& eDate, EventResults* events ) const;
    /* KComponent:: */
    virtual void getAccEndDates(DateTimeArray &accEndDates) const;

protected:
    /* KComponent:: */
    virtual void setup(const IModel* model, const MarketData* market);
    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;
    
    KTurbo(CClassConstSP const &type) : KComponent(type) {};
private:
    static void load(CClassSP& clazz); 
    static IObject* defaultConstructor(void) { return new KTurbo(TYPE); }
};
typedef smartPtr<KTurbo> KTurboSP;
typedef smartConstPtr<KTurbo> KTurboConstSP;
typedef array<KTurboSP, KTurbo> KTurboArray;
typedef smartPtr<KTurboArray> KTurboArraySP;

DRLIB_END_NAMESPACE
#endif
