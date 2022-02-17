//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : MQQuasiIRVolCMS.cpp
//
//   Description : class of IRVol for MultiQ dealing with CMS
//
//   Author      : Keith Jia
//
//   Date        : August 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/Addin.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/MultiQDistribution.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/Maths.hpp"
#include "edginc/MQQuasiIRVolCMS.hpp"
#include "edginc/VolRequestMQCMS.hpp"
#include "edginc/VolProcessedMQCMS.hpp"
#include "edginc/FAMQDistributionCMS.hpp"
#include "edginc/VNFMCMS.hpp"
#include "edginc/MDFUtil.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/MarketDataFetcher.hpp"

DRLIB_BEGIN_NAMESPACE

/* Make sure the class links */
bool  MQQuasiIRVolCMSLoad() {
    return (MQQuasiIRVolCMS::TYPE != 0);
}

//-----------------------------------------------------
// MQQuasiIRVolCMS methods
//-----------------------------------------------------
MQQuasiIRVolCMS::MQQuasiIRVolCMS(CClassConstSP clazz) 
    : MQQuasiIRVol(clazz)
{}

MQQuasiIRVolCMS::~MQQuasiIRVolCMS()
{}

void MQQuasiIRVolCMS::validatePop2Object()
{
    static const string method = "MQQuasiIRVolCMS::validatePop2Object";

    try
    {
        //no additional validatePop2Object here.  Call parent.
        MQQuasiIRVol::validatePop2Object();
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}


//-----------------------------------------------------
// IIRVol methods
//-----------------------------------------------------

/** Combines market and instrument data togeter to tive a processed vol
 */ 
CVolProcessed* MQQuasiIRVolCMS::getProcessedVol(const CVolRequest* volRequest,
                                                const MarketObject*  yc) const
{
    /** volRequest gives parSwapSprd, maturity date(resetDate), and all the necessary
     ** zero rates for vnfm estimation */

    static const char method[] = "MQQuasiIRVolCMS::getProcessedVol";

    if(volRequest->getClass() == VolRequestMQCMS::TYPE)
    {
        const VolRequestMQCMS *rq 
            = DYNAMIC_CAST(VolRequestMQCMS, volRequest);

        ExpiryConstSP tenor = rq->getTenor();
        double parSwap = rq->getParFwd();
        DateTime resetDate = rq->getResetDate();

        
        MultiQDistributionSP mq = MultiQDistributionSP(CreateMQ(tenor, parSwap, resetDate));

        VNFMSP vnfm = VNFMSP(new VNFMCMS(rq->getVolStart(),
                                         rq->getReset(),
                                         rq->getSwapStart(),
                                         rq->getSwapMat(),
                                         rq->getFreq(),
                                         rq->getBeta(),
                                         rq->getIsCashSettled(),
                                         parSwap,
                                         atmBSVol,
                                         rq->getFwdAnnuity()));
                                    
        double convexA, convexP, convexT;
        double convexFreq = rq->getFreq();

        vnfm->Q3VNFMZero2Swap(rq->getSwapMat(), rq->getSwapMatZR(),
                              convexA, convexP, convexT);


        double delayA, delayP, delayT; 
        double delayFreq = rq->getFreq();

        vnfm->Q3VNFMZero2Swap(rq->getPayDate(), rq->getPayZR(),
                              delayA,
                              delayP,
                              delayT);


        FAMQDistributionCMSSP famq = FAMQDistributionCMSSP
            (new FAMQDistributionCMS(mq,
                                     convexA, 
                                     convexP, 
                                     convexT, 
                                     convexFreq,
                                     delayA, 
                                     delayP, 
                                     delayT, 
                                     delayFreq,
                                     rq->getRichardsonOrder()));


        VolProcessedMQCMS *volPred = new VolProcessedMQCMS(name, metric, mq.get(), famq.get());
        
        return volPred;
    
    } else {
        throw ModelException(method, "VolRequest of type " +
                             volRequest->getClass()->getName() +
                             " not supported");
    }


}

//-----------------------------------------------------
// Invoked when class is loaded
//-----------------------------------------------------

void MQQuasiIRVolCMS::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(MQQuasiIRVolCMS, clazz);
    SUPERCLASS(MQQuasiIRVol);
    EMPTY_SHELL_METHOD(defaultIRVolCMS);
}

IObject* MQQuasiIRVolCMS::defaultIRVolCMS(){
    return new MQQuasiIRVolCMS();
}

//-----------------------------------------------------
// Static variables
//-----------------------------------------------------

CClassConstSP const MQQuasiIRVolCMS::TYPE 
= CClass::registerClassLoadMethod("MQQuasiIRVolCMS", typeid(MQQuasiIRVolCMS), load);

DEFINE_TEMPLATE_TYPE(MQQuasiIRVolCMSWrapper);


/************************************************
 ** MQQuasiIRVolAddin
 ************************************************
 */

class MQQuasiIRVolCMSAddin : public CObject
{
public:
    static CClassConstSP const TYPE;

private:
    MQQuasiIRVolCMSSP irVol;
    ExpirySP tenor;
    double parFwdRate;
    DateTime resetDate;
    MarketDataSP market;
    IModelSP model;

    static IObject* defaultAddin()
    {
        return new MQQuasiIRVolCMSAddin();
    }

    MQQuasiIRVolCMSAddin() : CObject(TYPE)
    {}

    IObjectSP createMQ()
    {
        static const string method = "MQQuasiIRVolCMSAddin::createMQ";

        MarketDataSP mkt(this->market);
        IModelSP mdl(this->model);

        if(!mdl)
        {
            mdl.reset(new NonPricingModel());
            MarketDataFetcherSP mdf = mdl->getMDF();
            MDFUtil::setUseCurrencyBasis(*mdf, true);
        }

        try
        {
            irVol->getMarket(mdl.get(), mkt.get());
            IObject* obj = irVol->CreateMQ(tenor, parFwdRate, resetDate);
            return IObjectSP(obj);
        }
        catch(exception &e)
        {
            throw ModelException(e, method);
        }
    }

    static void load(CClassSP& clazz)
    {
        clazz->setPublic();
        REGISTER(MQQuasiIRVolCMSAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultAddin);

        FIELD(irVol, "MQQuasiIRVolCMS object");
        FIELD(market, "market data");
        FIELD(model, "model");
        FIELD_MAKE_OPTIONAL(model);
        FIELD(tenor, "");
        FIELD(parFwdRate, "");
        FIELD(resetDate, "");
        
        Addin::registerObjectMethod
            ("MQ_IN_MQQUASIIRVOL",
             Addin::MARKET,
             "Retrive MQDistribution in MQQuasiIRVol",
             true,
             Addin::returnHandle,
             &MQQuasiIRVolCMSAddin::createMQ);
    }
};


CClassConstSP const MQQuasiIRVolCMSAddin::TYPE 
= CClass::registerClassLoadMethod("MQQuasiIRVolCMSAddin", typeid(MQQuasiIRVolCMSAddin), load);
         



DRLIB_END_NAMESPACE

