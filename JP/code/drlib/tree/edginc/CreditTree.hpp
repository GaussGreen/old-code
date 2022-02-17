//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : CreditTree.hpp
//
//   Description : top level model class for trees modeling credit
//
//----------------------------------------------------------------------------

#ifndef QR_CREDITTREE_HPP
#define QR_CREDITTREE_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/RateTree.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(CreditTree)
FORWARD_DECLARE(RiskyZeroBond)
FORWARD_DECLARE(ProtectionLeg)

class CreditTree : public RateTree {
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    /* METHODS TO HANDLE PRIMITIVE INSTRUMENTS */
    /**Add a risky zero-coupon bond for pricing.*/
    virtual void insertRiskyZero(const RiskyZeroBond& riskyZeroBond) = 0;
    /**Get a reference to a specific ZC bond's price at a given time-step.*/
    virtual const TreeSlice& getRiskyZeroValue(const RiskyZeroBond& riskyZeroBond, int step) = 0;
    /**Add a protection leg for pricing.*/
    virtual void insertProtectionLeg(const ProtectionLeg& protLeg) = 0;
    /**Get a reference to a specific protection leg's price at a given time step.*/
    virtual const TreeSlice& getProtectionLegValue(const ProtectionLeg& protLeg, int step) = 0;

    /** override FDModel::makeProduct to retrieve IndexSpec product */
    virtual FDProductSP makeProduct(const IProdCreatorSP & creator);

    virtual ~CreditTree();

protected:

    CreditTree(CClassConstSP type = TYPE);
};

/**Unit risky zero-coupon bond: instrument that pays 1 times outstanding principal at payment date.
   If (kLo>0 || kHi<1) tranche, pays max(L,kHi) - max(L,kLo).
   isMarketRecovery and digitalRecovery determine how losses (L) are related to number of defaults.
   if (!isMarketRecovery), L = (1-digitalRecovery)*defaultedNotional.
   if (zeroPaymentIfAllDefaulted) payout is zero if all notional has defaulted, regardless
   of the actual losses.
   */
class RiskyZeroBond : public CObject {
public:
    static CClassConstSP const TYPE;
    DateTime startDate;
    DateTime payDate;
    double kLo;
    double kHi;
    bool isMarketRecovery; // losses are determined by market recovery rate
    double digitalRecovery; // if (!isMarketRecovery), this is the rate
    bool zeroPaymentIfAllDefaulted; // if true, pays nothing if all names defaulted (no recovery).
    string discountCurveName;
    string creditCurveName;

    RiskyZeroBond(const string& discountCurve, const string& creditCurve,
        const DateTime& startDt, const DateTime& payDt, double attachmentStrike=0.0,
        double detachmentStrike=1.0) 
        : CObject(TYPE), startDate(startDt), payDate(payDt), 
          kLo(attachmentStrike), kHi(detachmentStrike), isMarketRecovery(true),
          digitalRecovery(0.0), zeroPaymentIfAllDefaulted(true) {}
    virtual ~RiskyZeroBond() {}
};

/**Unit protection leg. Pays each loss that occurs between startDate and endDate. Note that payment
   is only for losses that occur after startDate - no payments are made for losses before,
   so settlement is conditional. 
   */
class ProtectionLeg : public CObject {
public:
    DateTime startDate;
    DateTime endDate;
    double kLo;
    double kHi;
    bool isMarketRecovery; // losses are determined by market recovery rate
    double digitalRecovery; // if (!isMarketRecovery), this is the rate
    bool payAtEndDate; // true if all losses paid at end date, else pay as occur.
    string discountCurveName;
    string creditCurveName;

    ProtectionLeg(const string& discountCurve, const string& creditCurve,
        const DateTime& startDt, const DateTime& endDt,
        bool payAtEndDt = false,
        double attachmentStrike=0.0,
        double detachmentStrike=1.0) 
        : CObject(TYPE), startDate(startDt), endDate(endDt), kLo(attachmentStrike),
        kHi(detachmentStrike), isMarketRecovery(true),
          digitalRecovery(0.0), payAtEndDate(payAtEndDt) {}
    virtual ~ProtectionLeg() {}
};

/* ?? SHOULD WE INCLUDE ACCRUED RECOVERY AS A PRIMITIVE, TOO? OR JUST USE 
   APPROXIMATION BASED ON PROTECTION? EITHER AI = PROTECTION/2 OR
   AI = PROTECTION *(-1/lnZ) - Z/(1-Z)). PROBABLY APPROX GOOD ENOUGH, AS INCLUDING
   AS SEPARATE PRIMITIVE WILL INCREASE NUMBER OF SLICES A LOT. */

DRLIB_END_NAMESPACE

#endif
