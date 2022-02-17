//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CDSVol.hpp
//
//   Description : Interface defining what CDS Vols can do 
//
//   Author      : Mark A Robson
//
//   Date        : November 26, 2004
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_ICDSVOL_HPP
#define QLIB_ICDSVOL_HPP

#include "edginc/MarketObject.hpp"

DRLIB_BEGIN_NAMESPACE
class IVolProcessed;
class CVolRequest;
class ICDSParSpreads;

/** Interface defining what CDS Spot Vols can do (here 'spot' is in
    the C&R sense - the vol used in the diffusion equation
    and not the Black vols) */
class MARKET_DLL ICDSVol: public virtual IObject{
public:
    static CClassConstSP const TYPE; // in CDSParSpreads.cpp

    /** Combines market and instrument data together to give a
        Processed Vol */
    virtual IVolProcessed* getProcessedVol(const CVolRequest*    volRequest,
                                           const ICDSParSpreads* cds) const = 0;

    ~ICDSVol();
private:
    static void load(CClassSP& clazz);  // in CDSParSpreads.cpp
};

typedef MarketWrapper<ICDSVol> ICDSVolWrapper; // in CDSParSpreads.cpp
typedef smartConstPtr<ICDSVol> ICDSVolConstSP;
typedef smartPtr<ICDSVol>      ICDSVolSP;
#ifndef QLIB_ICDSVOL_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<ICDSVol>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<ICDSVol>);
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<ICDSVol>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<ICDSVol>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<ICDSVol>);
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<ICDSVol>);
#endif

DRLIB_END_NAMESPACE

#endif




