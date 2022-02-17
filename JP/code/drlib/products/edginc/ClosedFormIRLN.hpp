//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ClosedFormIRLN.hpp
//
//   Description : Closed Form Log Normal Algorithm (asks instrument to do it)
//                 Interest rate flavour
//
//   Author      : Mark A Robson/Andrew J Swain
//
//   Date        : 28 February 2005
//
//
//----------------------------------------------------------------------------

#ifndef CLOSED_FORM_IR_LN_HPP
#define CLOSED_FORM_IR_LN_HPP
#include "edginc/Model.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/IRVegaPointwise.hpp"

DRLIB_BEGIN_NAMESPACE

/** Implementation of Model algorithm where the algorithm is specific to
    the instrument and yet is based on a log normal methodology */
class PRODUCTS_DLL ClosedFormIRLN: public CModel,
                      public virtual IRVegaPointwise::ISensitivePoints{
public:
    static CClassConstSP const TYPE;
    friend class ClosedFormIRLNHelper;

    /** the class that the product must be able to create */
    class PRODUCTS_DLL IProduct{
    public:
        virtual void price(ClosedFormIRLN* model,
                           Control*        control, 
                           CResults*       results) const = 0;
        virtual ~IProduct() {};
    };

    /** interface that the instrument must implement */
    class PRODUCTS_DLL IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        friend class ClosedFormIRLNHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(ClosedFormIRLN* model) const = 0;
    };

    IObject* clone() const;

    /** Create a MarketDataFetcher which will be used for retrieving market data etc */
    virtual MarketDataFetcherSP createMDF() const;

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);

    /** Essentially relies on instrument implementing ISensitiveIRVolPoints.
        If not returns null. */
    virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
        OutputNameConstSP  outputName,
        const CInstrument* inst) const;

    // a version of black that uses 2Q smile unless smileStyle
    // is empty or set to NO_IR_Q_SMILE in which case it
    // just does the (log)normal thing
    double black(
        bool             isCall, 
        double           fwd, 
        double           strike, 
        double           pv, 
        double           variance,
        const IRVolBase* vol);

    ClosedFormIRLN();

    /** returns an ExposureHighlighter - a "model" that does everything except
        actually price, so you get to see what market data it uses 
        Default implementation supplied */
    virtual ExposureHighlighter* exposureHighlighter();

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant, because we price using a plain old vol
     * surface.  See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

private:
    static const string NO_IR_Q_SMILE;

    class MarketDataFetcherIR;
    ClosedFormIRLN(const ClosedFormIRLN &rhs);
    ClosedFormIRLN& operator=(const ClosedFormIRLN& rhs);

    string smileStyle;
};

DRLIB_END_NAMESPACE
#endif
