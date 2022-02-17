//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : ClosedFormCDSPSandFA.hpp
//
//   Description : Closed form model for taking CDS par spreads as input plus additional
//                 calculation of E2C specific sensitivities
//
//   Author      : Andre Segger
//
//   Date        : September 13, 2002
//
//
//----------------------------------------------------------------------------

#ifndef CLOSEDFORMCDSPSFA_HPP
#define CLOSEDFORMCDSPSFA_HPP
#include "edginc/Model.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/E2CModel.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/IForwardRatePricer.hpp"

DRLIB_BEGIN_NAMESPACE

/** Credit default swaps can be priced closed form in a few different ways. 
    One way is using CDS par spreads as market input. Another is with credit 
    DRs' asset-based closed-form model. The methods would take different 
    market data and so each is getting its own model object. This one is
    for the CDS par spreads which calculates additional sensitivities based 
    on the E2C model.
*/

class MARKET_DLL ClosedFormCDSPSandFA: public CModel,
                            virtual public  IE2CModel,
                            virtual public  IHasForwardRatePricer {
public:
    static CClassConstSP const TYPE;
    friend class ClosedFormCDSPSandFAHelper;

    /** the class that the product must be able to create */
    class MARKET_DLL IProduct{
    public:
        virtual void price(ClosedFormCDSPSandFA* model,
                           Control*    control, 
                           CResults*   results) const = 0;
    };

    /** interface that the instrument must implement */
    class MARKET_DLL IIntoProduct{
    public:
        friend class ClosedFormCDSPSandFAHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(ClosedFormCDSPSandFA* model) const = 0;
    };

    /** Simple constructor */
    ClosedFormCDSPSandFA();

    virtual ~ClosedFormCDSPSandFA() {}

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);
 
    virtual void setParSpreadPricing(const bool setPar);

    bool isPricePar() const;

    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;
    
    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant, because we don't price using a parametric
     * (latent, non-observable) pdf that needs risk mapping.  See
     * IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    /** Key method providing access to the pricer */
    virtual IForwardRatePricerSP getForwardRatePricer() const;

private:
    string  volType;
    bool    pricePar;
    mutable MarketDataFetcherLNSP mdf; // $unregistered

    ClosedFormCDSPSandFA(const ClosedFormCDSPSandFA &rhs);
    ClosedFormCDSPSandFA& operator=(const ClosedFormCDSPSandFA& rhs);
};

DRLIB_END_NAMESPACE
#endif
