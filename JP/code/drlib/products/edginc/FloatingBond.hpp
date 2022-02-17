//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FloatingBond.hpp
//
//   Description : Floating Bond implementation 
//
//   Author      : André Segger
//
//   Date        : 30 April 2002
//
//
//----------------------------------------------------------------------------

#ifndef FLOATING_BONDCASHFLOWS_HPP
#define FLOATING_BONDCASHFLOWS_HPP

#include "edginc/Bond.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/Theta.hpp"
#include "edginc/Asset.hpp"
#include "edginc/Settlement.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/RateParallel.hpp"
#include "edginc/RatePointwise.hpp"

DRLIB_BEGIN_NAMESPACE
/** A Bond implementation where you input the cashflows */

class PRODUCTS_DLL FloatingBond:     public Bond,
                        public Theta::IShift
{
public:
    static CClassConstSP const TYPE;
    friend class FloatingBondHelper;

    virtual void fieldsUpdated(const CFieldArray& fields);

    /** change the settlement of the bond if applicable */
    virtual void setSettlement(SettlementSP newSettle);

    virtual ~FloatingBond();
    
    /** get the accrued interest at a date */
    virtual double getAccruedAtDate(const DateTime &aiDate) const;

    /** returns only the cashFlows that occur after the start date */
    virtual CashFlowArraySP getCashFlows(const DateTime &startDate) const;
    
    /** returns all the cashFlows */
    virtual CashFlowArraySP getCashFlows() const;

    virtual DateTime getMaturityDate() const;

    virtual DateTime getUnadjMaturityDate() const;

    virtual DateTimeArraySP getExCouponDates() const;
    
    virtual double getFaceValue() const;

    virtual double getRedemption() const;

    virtual double getMaturityCoupon() const;

    virtual void validatePop2Object();

    virtual double priceFromYield(double yield, bool clean, DateTime yDate) const;

    virtual Bond* getBondToPut(const DateTime &putDate, double putLevel) const;

    void    getMarket(const IModel* model, const MarketData* market);

    void exportLiborRates(const string& header);

    void calculateLiborRates(bool isTheta);

    // get all cashflows known today
    virtual CashFlowArraySP getKnownCashFlows() const;

    // Theta
    bool sensShift(Theta* shift);

    virtual DateTime getAccrualStartDate() const;

    // when does bond settle?
    virtual DateTime settles(const DateTime& tradeDate) const;

    //////// Payment Class  //////////////
    class PRODUCTS_DLL FloatingBondPayment: public CObject{
    public:
        static CClassConstSP const TYPE;
        friend class FloatingBond;
        friend class FloatingBondPaymentHelper;
    
        FloatingBondPayment();

        FloatingBondPayment(const DateTime& refixDate,
                            const DateTime& paymentDate,
                            const double&   fixing,
                            const bool      isFixed);

        /** overrides default */
        virtual void validatePop2Object();
 
        /** access method for refix date */ 
        const DateTime& getRefixDate()    { return refixDate;   }

        /** access method for payment date */
        const DateTime& getPaymentDate()  { return paymentDate; }

        /** access method for fixing level */
        double          getFixing()       { return fixing;      }

        /** access method for fixed flag */
        bool            getIsFixed()      { return isFixed;     }        

    private:
        FloatingBondPayment(const FloatingBondPayment &rhs);
        FloatingBondPayment& operator=(const FloatingBondPayment& rhs);

        DateTime refixDate;
        DateTime paymentDate;
        double   fixing;
        bool     isFixed;
    };
        
    //////// End of Payment Class  //////////////

    typedef smartPtr<FloatingBondPayment>                     FloatingBondPaymentSP;
    typedef array<FloatingBondPaymentSP, FloatingBondPayment> FloatingBondPaymentArray;
    typedef smartPtr<FloatingBondPaymentArray>                FloatingBondPaymentArraySP;
    typedef smartConstPtr<FloatingBondPayment>                FloatingBondPaymentConstSP;

    FloatingBond(const double                      faceValue,
                 const DateTime                    datedDate,
                 const double                      redemptionPct,
                 const MaturityPeriodSP            liborPeriod,
                 const string                      liborDCCString,
                 const double                      liborSpread,
                 const YieldCurveWrapper&          liborReferenceCurve,
                 const DateTime                    valueDate,
                 const FloatingBondPaymentArray&   fixings);

    FloatingBond(const double                      faceValue,
                 const DateTime                    datedDate,
                 const double                      redemptionPct,
                 const MaturityPeriodSP            liborPeriod,
                 const string                      liborDCCString,
                 const double                      liborSpread,
                 const YieldCurveWrapper&          liborReferenceCurve,
                 const DateTime                    valueDate,
                 const FloatingBondPaymentArray&   fixings,
                 const bool                        deferredCoupon,
                 const CAssetWrapper               asset,
                 const DateTime                    couponDeferralEndDate,
                 const double                      couponDeferralThreshold);

private:
    FloatingBond();
    virtual void throwIfNotValid() const;
    virtual void validateInputs() const;
    

    double                      faceValue;
    DateTime                    datedDate;
    double                      redemptionPct;            // Optional
    MaturityPeriodSP            liborPeriod;
    string                      liborDCCString;
    double                      liborSpread;
    YieldCurveWrapper           liborReferenceCurve;
    DateTime                    valueDate;
    FloatingBondPaymentArray    fixings;

    // parameters for coupon deferral
    bool                        deferredCoupon;
    CAssetWrapper               asset;
    DateTime                    couponDeferralEndDate;
    double                      couponDeferralThreshold;

    // cap/floor parameters
    bool                        capFloatingCoupons;
    double                      floatingCouponCap;
    bool                        floorFloatingCoupons;
    double                      floatingCouponFloor;

    SettlementSP                settle;           // Optional

    // derived
    bool                 isValid;
    string               errorMessage; // $unregistered
    CashFlowArraySP      cashFlowsProcessed;
    DateTimeArraySP      exCouponDatesProcessed;
    DayCountConventionSP accruedDCC; // $unregistered
    DayCountConventionSP liborDCC;
};

typedef smartConstPtr<FloatingBond> FloatingBondConstSP;
typedef smartPtr<FloatingBond> FloatingBondSP;

DRLIB_END_NAMESPACE
#endif
