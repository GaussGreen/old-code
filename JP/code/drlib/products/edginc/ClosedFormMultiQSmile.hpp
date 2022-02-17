//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : ClosedFormMultiQSmile.hpp
//
//   Description : Closed-form European options using a Multi-Q 
//                 forward distribution
//
//   Author      : Charles Morcom
//
//   Date        : 16 December 2005
//
//
//----------------------------------------------------------------------------
#ifndef QR_CLOSEDFORMMULTIQSMILE_HPP
#define QR_CLOSEDFORMMULTIQSMILE_HPP

#include "edginc/Model.hpp"
#include "edginc/Results.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IForwardRatePricer.hpp"
DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ClosedFormMultiQSmile)
FORWARD_DECLARE(CInstrument)
FORWARD_DECLARE(Control)
FORWARD_DECLARE(MarketDataFetcher)

/** Pricing and tweaks are with respect to the BS implied vol smile at fixed strikes. */
  class PRODUCTS_DLL ClosedFormMultiQSmile: public CModel,
					    public virtual IHasForwardRatePricer {
public:
    static CClassConstSP const TYPE;
    friend class ClosedFormMultiQSmileHelper;

    /** the class that the product must be able to create */
    class PRODUCTS_DLL IProduct{
    public:
        virtual void price(ClosedFormMultiQSmile* model,
                           Control*        control, 
                           CResults*       results) const = 0;
        virtual ~IProduct() {};
    };

    /** interface that the instrument must implement */
    class PRODUCTS_DLL IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        friend class ClosedFormMultiQSmileHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(ClosedFormMultiQSmile* model) const = 0;
    };

    /** Create a MarketDataFetcher which will be used for retrieving market data etc */
    virtual MarketDataFetcherSP createMDF() const;

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);

    ClosedFormMultiQSmile();
    ClosedFormMultiQSmile(CClassConstSP clazz);
    virtual ~ClosedFormMultiQSmile();
  
  /** Key method providing access to the pricer */
    virtual IForwardRatePricerSP getForwardRatePricer() const;
  
private:
  
  IForwardRatePricerSP forwardRateModel;    /* used for determining fee leg cashflows */
    //class MarketDataFetcherBSSmile;
    //ClosedFormMultiQSmile(const ClosedFormMultiQSmile &rhs);
    //ClosedFormMultiQSmile& operator=(const ClosedFormMultiQSmile& rhs);
};

DRLIB_END_NAMESPACE
#endif
