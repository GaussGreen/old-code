//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolRequestRaw.hpp
//
//   Description : Vol request for accessing the 'raw' vol for cases where
//                 you know the type of vol that you are expecting and the
//                 vol implements IVolProcessed
//
//   Author      : Mark A Robson
//
//   Date        : 16 May 2005
//
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VOLREQUESTRAW_HPP
#define EDR_VOLREQUESTRAW_HPP

#include "edginc/VolProcessed.hpp"
#include "edginc/VolRequest.hpp"
#include "edginc/VolBase.hpp"

DRLIB_BEGIN_NAMESPACE
class IMultiFactors;
class CAsset;

/** Vol request for accessing the 'raw' vol for cases where you know
    the type of vol that you are expecting and the vol implements
    IVolProcessed. Note that currently all implementations just return "this"
    in the getProcessedVol() method - however this returns a non-const 
    IVolProcessed so the returned object must not be modified */
class MARKET_DLL VolRequestRaw: public CVolRequest {
public:
    static CClassConstSP const TYPE;

    virtual ~VolRequestRaw();
    VolRequestRaw();

    /** Little utility lets you grab a copy of a volbase that's in a
        MultiFactors */
    static CVolBaseSP copyVolBase(const IMultiFactors& mAsset,
                                  int                  iAsset);
    /** Little utility lets you grab a copy of the volbase that's in an asset */
    static CVolBaseSP copyVolBase(const CAsset& asset);
private:
    VolRequestRaw(const VolRequestRaw &rhs);
    VolRequestRaw& operator=(const VolRequestRaw& rhs);
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
};

// smart pointers for VolRequestRaw
typedef smartConstPtr<VolRequestRaw> VolRequestRawConstSP;
typedef smartPtr<VolRequestRaw> VolRequestRawSP;

DRLIB_END_NAMESPACE

#endif
