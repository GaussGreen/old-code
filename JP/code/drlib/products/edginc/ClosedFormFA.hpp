//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : ClosedFormFA.hpp
//
//   Description : Closed form model for taking firm asset info as input
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : September 3, 2001
//
//
//----------------------------------------------------------------------------

#ifndef CLOSEDFORMFA_HPP
#define CLOSEDFORMFA_HPP
#include "edginc/Model.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/IForwardRatePricer.hpp"

DRLIB_BEGIN_NAMESPACE

/** Credit default swaps can be priced closed form in a few different ways. 
    One way is using CDS par spreads as market input. Another is with credit 
    DR's asset-based closed-form model. The methods would take different 
    market data and so each is getting its own model object. This one is
    for credit DR's asset-based model.
*/

class PRODUCTS_DLL ClosedFormFA: public CModel,
                                 virtual public  IHasForwardRatePricer {
public:
    static CClassConstSP const TYPE;
    friend class ClosedFormFAHelper;

    /** the class that the product must be able to create */
    class PRODUCTS_DLL IProduct{
    public:
        virtual void price(ClosedFormFA* model,
                           Control*    control, 
                           CResults*   results) const = 0;
    };

    /** interface that the instrument must implement */
    class PRODUCTS_DLL IIntoProduct{
    public:
        friend class ClosedFormFAHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(ClosedFormFA* model) const = 0;
    };

    /** Simple constructor */
    ClosedFormFA();

    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);
    
    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant, because we price using a plain old vol
     * surface.  See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    /** Key method providing access to the pricer */
    virtual IForwardRatePricerSP getForwardRatePricer() const;
private:
    ClosedFormFA(const ClosedFormFA &rhs);
    ClosedFormFA& operator=(const ClosedFormFA& rhs);

    string  volType;
    mutable MarketDataFetcherLNSP mdf; // $unregistered
};

DRLIB_END_NAMESPACE
#endif
