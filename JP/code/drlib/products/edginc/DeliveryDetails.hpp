//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Filename    : DeliveryDetails.hpp
//
//   Description : Class used to contain the information regarding a delivery
//                 for physically settled defaults
//
//   Author      : Jose Hilera
//
//   Date        : May 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_DELIVERYDETAILS_HPP
#define QLIB_DELIVERYDETAILS_HPP

#include "edginc/Atomic.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


FORWARD_DECLARE(IBadDayAdjuster);

class DeliveryDetails : public CObject {
public:
    friend class PhysicalSettlementOverrideName;
    static CClassConstSP const TYPE;

    ~DeliveryDetails();
    
    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    /** Returns the delivery date of this delivery */
    DateTime getDeliveryDate() const;

    /** Returns the recovery rate for this delivery - the parameter
     * is the default recovery rate, in case it is not overriden here */
    double getDeliveryRecoveryRate(const double recoveryRate) const;

    /** Returns the (potentially estimate) calculation date for this
     * recovery, as used for CDO tranches */
    DateTime deliveryCalcDate(const DateTime& creditEventDate,
                              const DateTime& determinationDate,
                              const DateTime& valueDate,
                              IBadDayAdjusterConstSP bda) const;

    /** Returns the % notional expected to arrive on this delivery */
    double getNotionalFractionInDelivery() const;
    
private:
    // For reflection
    static void load (CClassSP& clazz);
    DeliveryDetails();
    static IObject* defaultDeliveryDetails();
    
    // FIELDS
    DateTime deliveryDate; // Date of this delivery
    CDoubleArraySP deliveryAmountFractions; // Array of notional fractions being delivered
    CDoubleSP recovery;   // Overall RR on this delivery

    /** An estimate for the delay between the credit event and the calculation 
     * date corresponding to this delivery */
    CIntSP defaultToCalculationDelay;

    /** Date when the recovery rate is determined for the obligations in this 
     * delivery (if the deliveries DO settle separately) */
    DateTime calculationDate;
};

DECLARE(DeliveryDetails);

DRLIB_END_NAMESPACE

#endif
