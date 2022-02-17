//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : BootstrappableCDSParSpreads.hpp
//
//   Description : Functionality common to all CDSParSpreads implementing the
//                 ICDSBootstrappable interface.
//                 Note this class is still abstract since it does not provide
//                 an implementations for all methods defined in the 
//                 ICDSBootstrappable interface.
//   
//   Date        : Oct 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_BOOTSTRAPPABLE_CDSPARSPREADS_HPP
#define QLIB_BOOTSTRAPPABLE_CDSPARSPREADS_HPP

#include "edginc/Class.hpp"
#include "edginc/Duration.hpp"
#include "edginc/ICDSBootstrappable.hpp"
#include "edginc/CDSPricer.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart) 
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(CreditIndexBasis)

class RiskyDurationCalculator;
typedef smartPtr<RiskyDurationCalculator> RiskyDurationCalculatorSP;

// For defining the gradient of par spreads with respect to clean spreads
struct CDSGradientResult;
typedef refCountPtr<CDSGradientResult> CDSGradientResultSP;


/** 
 * Functionality common to all CDSParSpreads implementing the ICDSBootstrappable
 * interface. Note this class is still abstract since it does not provide an 
 * implementations for all methods defined in ICDSBootstrappable.
 * */
class MARKET_DLL BootstrappableCDSParSpreads :
    public MarketObject,
    public virtual ICDSBootstrappable
{
public:
    static CClassConstSP const TYPE;

    virtual ~BootstrappableCDSParSpreads();
        
    /** Returns cash flow dates corresponding to expiry = index
     * (before bad days and holidays adjustment) */
    virtual const DateTimeArray getCashFlowDates(int index) const;

    /** Returns cash flows corresponding to expiry = index */
    virtual CashFlowArraySP getCashFlowArray(int index) const;

    /** Returns cash flows corresponding to expiry = index, adjusting
     * the spreads using the index basis */
    virtual CashFlowArraySP getCashFlowArray(CreditIndexBasisConstSP indexBasis,
                                             int index) const;

    /** Returns the RiskyDurationCalculator associated to this ICDSParSpreads */
    virtual RiskyDurationCalculatorSP getRiskyDurationCalculator() const;

    // required by CDSHelper : would have it private otherwise
    CashFlowArraySP getCashFlowArray(int index, 
                                 DateTime today,
                                 const BadDayConvention* bdc, 
                                 const Holiday* hols) const;

    // allow the spread to be specified
    CashFlowArraySP getCashFlowArray(int index, 
                                 DateTime today,
                                 const BadDayConvention* bdc, 
                                 const Holiday* hols,
                                 const double spread) const;

    /** Calculates the maximum tweak sizes */
    virtual DoubleArrayConstSP calculateTweakSizes() const;

    /** Returns the upfront fee, payable today, for a CDS
        with a particular maturity and particular fixed fee */
    /** Not expected to be overridden */
    double impliedUpfrontFee(const DateTime&   valueDate,       // (I)
                             YieldCurveConstSP discount,        // (I)
                             const DateTime&   maturity,        // (I)
                             const double      fee,             // (I)
                             const DateTime&   prevMaturity,    // (I) we've already bootstrapped to here
                             DefaultRatesSP&   defRates) const; // (M) the bootstrapping so far

protected:
    BootstrappableCDSParSpreads(CClassConstSP clazz);

private:
    class RDCalculator;

    /** Constructor */
    BootstrappableCDSParSpreads();
    
    /** Invoked when Class is 'loaded', used for reflection */
    static void load(CClassSP& clazz);
  
    /** Copy constructor : don't use */  
    BootstrappableCDSParSpreads(
        const BootstrappableCDSParSpreads& adjustedCDSParSpreads);

    /** Override '=' operator : don't use */
    BootstrappableCDSParSpreads& operator=(
        const BootstrappableCDSParSpreads& adjustedCDSParSpreads);

    void getCashFlowDates(int             index, 
                          const DateTime& today, 
                          DateTimeArray&  dates) const;

    // Allows for margin of error in computing max tweak sizes
    static const double CDS_ADJUSTMENT;        
    
    /** Calculates the gradient of par spreads with respect to clean spreads */
    CDSGradientResultSP calcMatrixGradients() const;

};


DECLARE(BootstrappableCDSParSpreads);
typedef MarketWrapper<BootstrappableCDSParSpreads> BootstrappableCDSParSpreadsWrapper;

DRLIB_END_NAMESPACE

#endif
