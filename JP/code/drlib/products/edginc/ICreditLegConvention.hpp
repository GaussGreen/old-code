//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : ICreditLegConvention.hpp
//
//   Description : Interface for "credit legs" conventions (see below)
//
//   Author      : Antoine Gregoire
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#ifndef I_CREDIT_LEG_CONVENTION_HPP
#define I_CREDIT_LEG_CONVENTION_HPP

#include "edginc/YieldCurve.hpp"
#include "edginc/ICreditFeeLeg.hpp"
#include "edginc/ICreditContingentLeg.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Interface for "credit legs" conventions:
 * ICreditLegConvention objects must be able to generate
 * a fee leg and a contingent leg.
 *
 * Concrete classes are:
 * - CDSLegConvention: conventions for single name CDS, CDS on indices, CDO
 * - CDSLegConventionIMM: same as CDSLegConvention, but using standard IMM dates
 *
 * */
class PRODUCTS_DLL ICreditLegConvention : public virtual IGetMarket {
public:
    /** Virtual destructor */
    virtual ~ICreditLegConvention();

    /** Returns the default start date given the trade date */
    virtual DateTime startDate(const DateTime& tradeDate) const = 0;

    /** Generates the fee leg given the end date */
    virtual ICreditFeeLegSP generateFeeLeg(
        const DateTime& startDate,
        const DateTime& endDate,
        double coupon,
        double upfrontPayment,
        double notional = 1.0) const = 0;

    /** Generates the contingent leg given the end date */
    virtual ICreditContingentLegSP generateContingentLeg(
        const DateTime& startDate,
        const DateTime& endDate,
        double coupon,
        double upfrontPayment,
        double notional = 1.0) const = 0;

    /** Access to discount curve name */
    virtual string getDiscountName() const = 0;

    /** Access to discount curve */
    virtual YieldCurveConstSP getDiscount() const = 0;

    /** access to recoverNotional flag */
    virtual bool getRecoverNotional() const = 0;

    /** TYPE */
    static CClassConstSP const TYPE;
};

// support for ref count pointers 
typedef smartPtr<ICreditLegConvention> ICreditLegConventionSP;
//support for wrapper
typedef MarketWrapper<ICreditLegConvention> ICreditLegConventionWrapper;

DRLIB_END_NAMESPACE

#endif /*I_CREDIT_LEG_CONVENTION_HPP*/

