//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : MultiQQuasiSmile.hpp
//
//   Description : Numerical European options using a Multi-Q 
//                 forward distribution
//
//   Author      : Mehdi Chaabouni
//
//   Date        : Feb 2006
//
//
//----------------------------------------------------------------------------
#ifndef QR_MultiQQuasiSmile_HPP
#define QR_MultiQQuasiSmile_HPP

#include "edginc/Model.hpp"
#include "edginc/Results.hpp"
#include "edginc/IForwardRatePricer.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(MultiQQuasiSmile);
FORWARD_DECLARE(CInstrument);
FORWARD_DECLARE(Control);
FORWARD_DECLARE(MarketDataFetcher);
    
/** Pricing and tweaks are with respect to the BS implied vol smile at fixed strikes. */
class MultiQQuasiSmile: public CModel,
                        public virtual IHasForwardRatePricer
{
public:
    static CClassConstSP const TYPE;
    friend class MultiQQuasiSmileHelper;

    /** the class that the product must be able to create */
    class IProduct{
    public:
        virtual void price(MultiQQuasiSmile* model,
                           Control*        control, 
                           CResults*       results) const = 0;
        virtual ~IProduct() {};
    };

    /** interface that the instrument must implement */
    class IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        friend class MultiQQuasiSmileHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(MultiQQuasiSmile* model) const = 0;
    };

    /** Create a MarketDataFetcher which will be used for retrieving market data etc */
    virtual MarketDataFetcherSP createMDF() const;

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);

    MultiQQuasiSmile();
    MultiQQuasiSmile(CClassConstSP clazz);
    virtual ~MultiQQuasiSmile();

    //---------------------------------
    // IHasForwardRatePricer methods
    //---------------------------------
    virtual IForwardRatePricerSP getForwardRatePricer() const;

private:
    IForwardRatePricerSP forwardRateModel;

    //class MarketDataFetcherBSSmile;
    //MultiQQuasiSmile(const MultiQQuasiSmile &rhs);
    //MultiQQuasiSmile& operator=(const MultiQQuasiSmile& rhs);
};

DRLIB_END_NAMESPACE
#endif
