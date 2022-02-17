//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRVolBase.hpp
//
//   Description : Abstract base for IR vols
//
//   Author      : Andrew J Swain
//
//   Date        : 11 January 2002
//
//
//----------------------------------------------------------------------------

#ifndef _IRVOLBASE_HPP
#define _IRVOLBASE_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/IRGridPoint.hpp"
#include "edginc/OutputName.hpp"
#include "edginc/VolRequest.hpp"
#include "edginc/VolProcessed.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/IRGridPoint.hpp"

DRLIB_BEGIN_NAMESPACE
class CInstrument;

/** This class just allows classes to insist upon an IRVol if they should so
    desire.
 */
class MARKET_DLL IRVolBase: public MarketObject{
public:
    static CClassConstSP const TYPE;

    virtual ~IRVolBase();

    /** Returns name of vol */
    virtual string getName() const = 0;

    /** Combines market and instrument data together to give a
        Processed Vol */
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const MarketObject*  yc) const = 0;

protected:
    IRVolBase(const CClassConstSP& clazz);

private:
    static void load(CClassSP& clazz);
    IRVolBase(const IRVolBase& irv);
    IRVolBase& operator=(const IRVolBase& irv);
};

typedef smartConstPtr<IRVolBase> IRVolBaseConstSP;
typedef smartPtr<IRVolBase> IRVolBaseSP;
#ifndef QLIB_IRVOLBASE_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<IRVolBase>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IRVolBase>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<IRVolBase>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IRVolBase>);
#endif

// support for wrapper class
typedef MarketWrapper<IRVolBase> IRVolBaseWrapper;
#ifndef QLIB_IRVOLBASE_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<IRVolBase>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<IRVolBase>);
#endif


class MARKET_DLL IRVolCommon: public IRVolBase{
public:
    static CClassConstSP const TYPE;

    virtual ~IRVolCommon();

protected:
    IRVolCommon(const CClassConstSP& clazz);

private:
    static void load(CClassSP& clazz);
    IRVolCommon(const IRVolCommon& irv);
    IRVolCommon& operator=(const IRVolCommon& irv);
    };

typedef smartConstPtr<IRVolCommon> IRVolCommonConstSP;
typedef smartPtr<IRVolCommon> IRVolCommonSP;
#ifndef QLIB_IRVOLBASE_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<IRVolCommon>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IRVolCommon>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<IRVolCommon>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IRVolCommon>);
#endif

// support for wrapper class
typedef MarketWrapper<IRVolCommon> IRVolCommonWrapper;
#ifndef QLIB_IRVOLBASE_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<IRVolCommon>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<IRVolCommon>);
#endif

DRLIB_END_NAMESPACE

#endif
