//----------------------------------------------------------------------------
//
//   File        : MQCMS.hpp
//
//   Description : A closed form model capable
//                 of generating and pricing deterministic forward rates
//
//----------------------------------------------------------------------------

#ifndef QLIB_MQCMSPRICER_HPP
#define QLIB_MQCMSPRICER_HPP

#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/Object.hpp"
#include "edginc/IForwardRatePricer.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_WRAPPER(MQQuasiIRVolCMS);


class MARKET_DLL MQCMSPricer : public CObject, 
                               public IForwardRatePricer
{
public:
    
    static CClassConstSP const TYPE;
    
    //allow public construction
    MQCMSPricer();
    
    //forward
    virtual void Forward(CObject *instrument, double *result);
    
    //getMarket
    virtual void getMarket(IModel *model, const MarketData *market);


    //get mqVol
    MQQuasiIRVolCMSSP getMQVol();
    double getBeta() const;
    bool getIsCashSettled() const;
    double getRichardsonOrder() const;

    class IProduct
    {
    public:
        //the method to return the forward rate
        virtual double fwd(MQCMSPricer* model) const = 0;
    };
    
    //Interface that target objects must implement
    class IIntoProduct
    {
    public:
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(MQCMSPricer* pricer) const = 0;
    };
    
private:
    MQQuasiIRVolCMSWrapper   mqVol;
    bool                     isCashSettled;     //(I), VNFM
    double                   beta;              //(I), VNFM, mean reversion
    int                      RichardsonOrder;   //(I), FAMQ

    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
};

DECLARE(MQCMSPricer);

DRLIB_END_NAMESPACE

#endif
