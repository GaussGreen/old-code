//----------------------------------------------------------------------------
//
// Group       : Credit Hybrids QR
//
// Description : View into instrument as required for FtDs
//
// Date        : October 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_FTDPRODUCT_HPP
#define QLIB_FTDPRODUCT_HPP

//#include "edginc/Results.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/GeneralisedCDO.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE(Control);
FORWARD_DECLARE(Results);
FORWARD_DECLARE(MarketFactor);
FORWARD_DECLARE(AbstractCashFlow);
FORWARD_DECLARE(EffectiveCurve);
FORWARD_DECLARE(YieldCurve);
FORWARD_DECLARE(NToDefaultLossConfig);

class PRODUCTS_DLL FtDProduct : 
    public virtual IGeneralisedConvolutionProduct
{
public:
    /** Constructor for an FtDProduct that prices using a 
        one-factor Gaussian copula */
    FtDProduct(ConvolutionEngineConstSP convolutionEngine,
               const GeneralisedCDO* inst);

    //++++++++++++  IGeneralisedConvolutionProduct methods
    //
    /** Invoked by model to compute the price.
        Prices the FTD product using a semi-analytical closed form model
        and stores the result in the "results" object */
    virtual void price(const IConvolutionModel* convolutionModel,
                       Control*                 control, 
                       Results*                 results,
                       IForwardRatePricerSP     model) const;

    /** Returns the object that defines the losses for this product */
    virtual ICreditLossConfigConstSP getLossConfig() const;

    /** Returns the last observation date for credit related data */
    virtual DateTime lastObservationDate() const;

    /** Returns the max of maxDate and the last date from when a yield
        curve is used (eg for discounting or rate estimation) by
        payoff method.  This method does not have to worry about eg what
        dates are used when the clean spread curves are built. Used to
        control when to stop tweaking */
    virtual DateTime lastYCSensDate(const DateTime& maxDate) const;

    /** Returns the [riskless] curve for discounting (could remove this if
        we had better methods on ICDSParSpreads - only used for calculating
        durations) */
    virtual YieldCurveConstSP getDiscount() const;
    //
    //------------  IGeneralisedConvolutionProduct methods

    static NToDefaultLossConfigConstSP getFtDLossCfg(const GeneralisedCDO* inst);

private:
    class FtDOutputs;

    /** Generates market factors that are one dimension and normal(0,1) */
    void generateMarketFactors(MarketFactorArray& integrands) const;

    /** Computes accrual at a certain date; Should be a member of CashFlow object */
    double getAccrualAmount(const AbstractCashFlow& cashflow, 
                            const DateTime& atDate) const;

    void storeOutputRequests(FtDOutputs myOutput, 
                             Control*   control, 
                             Results*   results) const;

    // FIELDS
    ConvolutionEngineConstSP    convolutionEngine;
    const GeneralisedCDO*       inst;
};

DRLIB_END_NAMESPACE

#endif
