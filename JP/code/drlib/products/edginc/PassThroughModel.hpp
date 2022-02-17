//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PassThroughModel.hpp
//
//   Description : Model that wraps another model - 
//                 similar to ClosedForm in that asks instrument to drive the
//                 process, but routes pricing etc through wrapped model
//
//   Author      : Andrew J Swain
//
//   Date        : 19 August 2003
//
//
//----------------------------------------------------------------------------

#ifndef PASSTHROUGHMODEL_HPP
#define PASSTHROUGHMODEL_HPP
#include "edginc/Model.hpp"

DRLIB_BEGIN_NAMESPACE

/** Model that wraps another model - 
    similar to ClosedForm in that asks instrument to drive the
    process, but routes pricing etc through wrapped model
*/
class PRODUCTS_DLL PassThroughModel: public CModel {
public:
    static CClassConstSP const TYPE;

    // price using wrapped model - gateway for instrument to get access
    void calculate(CInstrument* instrument, 
                   Control*     control, 
                   Results*     results);

    /** the class that the product must be able to create */
    class PRODUCTS_DLL IProduct {
    public:
        virtual void price(PassThroughModel* model,
                           Control*          control, 
                           Results*          results) const = 0;
        virtual ~IProduct() {};
    };

    /** interface that the instrument must implement */
    class PRODUCTS_DLL IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        friend class PassThroughModelHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(PassThroughModel* model) const = 0;
    };

    /** calculate single price and store result in Results */
    virtual void Price(CInstrument* instrument, 
                       Control*     control, 
                       Results*     results);
   
    /** Returns a [deep] copy of the market data with supplied name
        and type from the given market data cache. This gives the
        model a chance to choose a specific type of market data rather
        than just a general instance. For example, the method could
        request a Black-Scholes Vol rather than just any old vol. The
        default implementation provided by CModel just asks the market
        data for the object of the given type. Here if the type is
        derived from CVolBase we override the type with our specific
        type of volatility */
    virtual MarketObjectSP GetMarket(const MarketData*    market,
                                     const string&        name,
                                     const CClassConstSP& type) const;
 
    /** Invoked after instrument has got its market data. Allows model to
    get any extra data required. Default implementation does nothing */
    virtual void getMarket(const MarketData* market, IInstrumentCollectionSP instruments);

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Delegates to the wrapped model.  See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

protected:
    PassThroughModel(CClassConstSP clazz);

private:
    friend class PassThroughModelHelper;
    PassThroughModel();
    PassThroughModel(const PassThroughModel &rhs);
    PassThroughModel& operator=(const PassThroughModel& rhs);

    // fields
    IModelSP model;
};

DRLIB_END_NAMESPACE
#endif
