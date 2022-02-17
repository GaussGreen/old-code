//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StepDownBond.hpp
//
//   Description   Step down bond instrument. 
//
//
//----------------------------------------------------------------------------

#ifndef EDG_STEPDOWN_BOND_HPP
#define EDG_STEPDOWN_BOND_HPP

#include "edginc/Generic1Factor.hpp"
#include "edginc/PayStream.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Tree1f.hpp"
#include "edginc/FD1F.hpp"
#include "edginc/SwapLegIntFace.hpp"
#include "edginc/PastObservations.hpp"

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL CStepDownBond : public Generic1Factor, 
					public FDModel::IIntoProduct,
                    //public Theta::IShift,
                    //public ISensitiveStrikes,
                    public LastSensDate
{
public:
    static CClassConstSP const TYPE;

    virtual void Validate();
 
    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;

    // below are copied from Vanilla
    virtual DateTime getValueDate() const;

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
        returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

    virtual DateTime endDate(const Sensitivity* sensControl) const;

    virtual void recordOutputRequests(Control* control, Results* results, double fairValue) const;

    CVolRequestConstSP makeVolRequest() const;

    /** get historical max or min value arrays. If it's not sampled yet, return 0 */
    CashFlowArray getHistMinOrMaxValues() const; 

    void addOutputRequests(Control* control,
                           Results* results,
                           const double& fairValue,
                           const double& indVol) const;
    
    bool sensShift(Theta* shift);

private:
    friend class CStepDownBondHelper;
	friend class CStepDownBondFDProd;

    static void load(CClassSP& clazz);
    CStepDownBond();
    CStepDownBond(CClassConstSP clazz);

    // ladder dates, level and coupon rates
    CDoubleArray    LadderLevelArray; // this should be matrix but 1D array for generic type. size= number of coupon periods * rungNumber
    CDoubleArray    LadderCouponArray; // this should be matrix but 1D array for generic type. size= number of coupon periods * (rungNumber+1)
    int             RungNumber; // number of rungs
    bool            IntraDayMonitor; // false means once in a day
    double          CurrentMinOrMax;  // current min (or max if step up!) level reached in this monioring period, MUST update this !!! 
    bool            isMonitorOnMax; // default=false, step down, true=step up
    bool            LinearCoupon; // default=false, step coupon, true=coupon rates linearly interpolated between ladder levels
    int             LinearGridSize; // size of grid for linear coupon use, default=61
    double          VolStrike;  // vol strike for LN case only, 
    bool            SwitchPlainSwap;  //True for price with PlainSwap(default), False retuns only LadderCoupon without Plain Swap.

    PastObservationsSP pastValues;          //Historic extreame value for each Ladder Observation time.    

    // libor leg
    FixedLegSP      fixedLeg;       // fixed leg
    
    LiborLegSP      liborLeg;       // libor leg

    int                     DEBUG_num_segments;
    
};

DRLIB_END_NAMESPACE
#endif
