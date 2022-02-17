//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : MQQuasiIRVolCMS.hpp
//
//   Description : Class of IR Vol for MultiQ quasi pricing for CMS
//
//   Author      : Keith Jia
//
//   Date        : July 2006
//
//
//----------------------------------------------------------------------------

#ifndef _MQQUASIIRVOLCMS_HPP
#define _MQQUASIIRVOLCMS_HPP

#include "edginc/MQQuasiIRVol.hpp"

DRLIB_BEGIN_NAMESPACE



class MARKET_DLL MQQuasiIRVolCMS: public MQQuasiIRVol {
public:
    static CClassConstSP const TYPE;

    virtual ~MQQuasiIRVolCMS();

    void validatePop2Object();

    //--------------------------------
    // IIRVol methods
    //--------------------------------
    virtual CVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                           const MarketObject*  yc) const;

protected:
    MQQuasiIRVolCMS(CClassConstSP clazz = TYPE);

private:
    static void load(CClassSP& clazz);
    static IObject* defaultIRVolCMS();

    friend class MQQuasiIRVolCMSAddin;
};


DECLARE(MQQuasiIRVolCMS);
typedef MarketWrapper<MQQuasiIRVolCMS> MQQuasiIRVolCMSWrapper;

DRLIB_END_NAMESPACE

#endif
