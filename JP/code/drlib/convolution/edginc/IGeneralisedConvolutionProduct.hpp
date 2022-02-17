//----------------------------------------------------------------------------
//
// Group       : CH Quantitative Research
//
// Description : Interface that all products created by a GeneralisedCDO when
//               using a ConvolutionEngine must implement

//
// Date        : October 2006
//
//----------------------------------------------------------------------------

#ifndef QR_IGENERALISEDCONVOLUTIONPRODUCT_HPP
#define QR_IGENERALISEDCONVOLUTIONPRODUCT_HPP

#include "edginc/Control_forward.hpp"
#include "edginc/Results_forward.hpp"
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

class IConvolutionModel;
FORWARD_DECLARE(IForwardRatePricer);
FORWARD_DECLARE(ICreditLossConfig);
FORWARD_DECLARE(YieldCurve);


class CONVOLUTION_DLL IGeneralisedConvolutionProduct : 
    public virtual VirtualDestructorBase 
{
public:
    virtual ~IGeneralisedConvolutionProduct();

    /** Invoked by model to compute the price */
    virtual void price(const IConvolutionModel* convolutionModel,
                       Control*                 control, 
                       Results*                 results,
                       IForwardRatePricerSP     model) const = 0;

    /** Returns the object that defines the losses for this product */
    virtual ICreditLossConfigConstSP getLossConfig() const = 0;

    /** Returns the last observation date for credit related data */
    virtual DateTime lastObservationDate() const = 0;

    /** Returns the max of maxDate and the last date from when a yield
        curve is used (eg for discounting or rate estimation) by
        payoff method.  This method does not have to worry about eg what
        dates are used when the clean spread curves are built. Used to
        control when to stop tweaking */
    virtual DateTime lastYCSensDate(const DateTime& maxDate) const = 0;

    /** Returns the [riskless] curve for discounting (could remove this if
        we had better methods on ICDSParSpreads - only used for calculating
        durations) */
    virtual YieldCurveConstSP getDiscount() const = 0;


protected:
    IGeneralisedConvolutionProduct();

private:
    // Not defined, so do not use
    IGeneralisedConvolutionProduct(const IGeneralisedConvolutionProduct& rhs);
    IGeneralisedConvolutionProduct& operator=(const IGeneralisedConvolutionProduct& rhs);
};

DRLIB_END_NAMESPACE

#endif
