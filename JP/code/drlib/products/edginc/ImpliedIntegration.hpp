//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ImpliedIntegration.hpp
//
//   Description : Implied Integration Algorithm 
//
//   Author      : Andrew J Swain
//
//   Date        : 4 November 2002
//
//
//----------------------------------------------------------------------------

#ifndef IMPLIEDINTEGRATION_HPP
#define IMPLIEDINTEGRATION_HPP
#include "edginc/Model.hpp"
#include "edginc/ImpliedSample.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/Theta.hpp"
#include "edginc/Delta.hpp"
#include <map>

DRLIB_BEGIN_NAMESPACE 

/** Implied Integration Algorithm */
class PRODUCTS_DLL ImpliedIntegration: public CModel,
                          virtual public Theta::RestorableShift
{
public:
    static CClassConstSP const TYPE;

    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;

    /** the class that the product must be able to create */
    class PRODUCTS_DLL IProduct {
    public:
        /** invoke the pricing */
        virtual void price(ImpliedIntegration* model,
                           Control*            control, 
                           CResults*           results) = 0;

        /** following methods drive the Implied integration */
        virtual LinearImpliedSampler* sampler(
            const DateTime&  date,
            const PDFParams& params) const = 0;

        /** given a spot level, return the payoff */
        virtual Function1DDouble* payoff() const = 0;

        /** given a spot level, return integral of the payoff */
        virtual Function1DDouble* indefiniteIntegral() const = 0;

        /** when are we getting the distribution for */
        virtual DateTime time() const = 0;

        /** when is now ? */
        virtual DateTime today() const = 0;

        virtual ~IProduct() {};
    };

    /** interface that the instrument must implement */
    class PRODUCTS_DLL IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        friend class ImpliedIntegrationHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(ImpliedIntegration* model) const = 0;
    };

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       Control*      control, 
                       CResults*     results);

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant.  See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;
 
    /** integrate payoff */
    virtual double integrate(ImpliedIntegration::IProduct* product) const;

    /** what are the integration limits ? */
    StrikesPartition* limits(ImpliedIntegration::IProduct* product) const;

    /** for VEGA_MATRIX - chose some semi-arbitrary strikes to report against */
    /** mid/tailWidth -> determine hwo many strikes are used in mid/tail 
        regions -> get steps/width strikes
    */
    DoubleArraySP sensitiveStrikes(ImpliedIntegration::IProduct* product,
                                   int                           midWidth,
                                   int                           tailWidth);

    DoubleArraySP sensitiveStrikes(const DateTime& today,
                                   const DateTime& maturity,
                                   const CAsset*   asset,
                                   int             midWidth,
                                   int             tailWidth);

    /** Constructor takes type of vol to use */
    ImpliedIntegration(const string& volType);
   
    /** Does special things for Theta-type tweaks */
    virtual bool sensShift(Theta* shift);
    virtual void sensRestore(Theta* shift);
    
    virtual void flush(); 
    IObject* clone() const;

private:
    friend class ImpliedIntegrationHelper;
    ImpliedIntegration();
    ImpliedIntegration(const ImpliedIntegration &rhs);
    ImpliedIntegration& operator=(const ImpliedIntegration& rhs);

    DoubleArraySP sensitiveStrikes(
        StrikesPartition* bounds,
        int               midWidth,
        int               tailWidth);

    string volType;
    int    midSteps;
    int    tailSteps;
    double tailStdDev;

    // transients
    bool                           isPricing;
    bool                           isTimeShift;    
    bool                           isDeltaShift;

    //mutable LinearImpliedSamplerSP pricingSampler;
    
    typedef map<DateTime, LinearImpliedSamplerSP> SamplerMap;
    mutable SamplerMap samplerMap; // $unregistered

    mutable MarketDataFetcherSP mdf; // $unregistered
};

typedef smartPtr<ImpliedIntegration> ImpliedIntegrationSP;

DRLIB_END_NAMESPACE
#endif
