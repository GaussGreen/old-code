//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolProcessedDispatch.hpp
//
//   Description : Manages system for 'double dispatch' method invocation for
//                 getProcessedVol method ie allows methods to be registered
//                 for particular types of VolRequests
//
//   Author      : Mark A Robson
//
//   Date        : 25 May 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VOL_PROCESSED_DISPATCH_HPP
#define EDR_VOL_PROCESSED_DISPATCH_HPP
#include "edginc/VolFunctorMember.hpp"
#include "edginc/VolFunctorStatic.hpp"
#include "edginc/VolFunctorAssetMember.hpp"
#include "edginc/VolFunctorAssetStatic.hpp"

DRLIB_BEGIN_NAMESPACE
class CAsset;
class IVolProcessed;
class CVolBase;

/** Manages system for 'double dispatch' method invocation for
    getProcessedVol method ie allows methods to be registered
    for particular types of VolRequests. The idea is rather than do a series
    of switches in getProcessedVol (either on assets or vols) you just write
    a series of getProcessedVol functions for each type of VolRequest that
    you support. Then your getProcessedVol() method should then just call the
    appropriate dispatch method below. You, of course, also need to call
    the appropriate registerXXXMethod below for each of your specialised
    getProcessedVol methods. Note that each of these specialised 
    getProcessedVol methods can either be static or member functions. For
    examples see VolBaseParamSurface.cpp, VolSurface.cpp and XCB.cpp */
class MARKET_DLL VolProcessedDispatch {
public:
    /** Register a class member [method] which is a 'getProcessedVol' method
        for Vol V, request R, and returns IVolProcessed P */
    template <class P, class V, class R> static void registerVolMethod(
        P* (V::* method)(const R*, const CAsset*) const){ /* holds function
                                                             pointer */
        VolFunctor* func = new VolFunctorMember<P, V, R, CAsset>(method);
        addToHashtable(V::TYPE, R::TYPE, func);
    }

    /** Register a static class method which is a 'getProcessedVol' method
        for Vol V, request R, and returns IVolProcessed P */
    template <class P, class V, class R> static void registerVolMethod(
        P* (method)(const V*, const R*, const CAsset*)){ /* holds C style
                                                            function pointer */
        VolFunctor* func = new VolFunctorStatic<P, V, R, CAsset>(method);
        addToHashtable(V::TYPE, R::TYPE, func);
    }

    /** Register a class member [method] which is a 'getProcessedVol' method
        for asset A, request R, and returns IVolProcessed P */
    template <class P, class A, class R> static void registerAssetMethod(
        P* (A::* method)(const R*) const){ /* holds function pointer */
        VolFunctor* func = new VolFunctorAssetMember<P, A, R>(method);
        addToHashtable(A::TYPE, R::TYPE, func);
    }

    /** Register a static class method which is a 'getProcessedVol' method
        for asset A, request R, and returns IVolProcessed P */
    template <class P, class A, class R> static void registerAssetMethod(
        P* (method)(const A*, const R*)){ /* holds C style function pointer */
        VolFunctor* func = new VolFunctorAssetStatic<P, A, R>(method);
        addToHashtable(A::TYPE, R::TYPE, func);
    }

    /** Invoke the method registered for this type of vol and this type of
        request */
    static IVolProcessed* dispatch(const CVolBase*      vol,
                                   const CVolRequest*   request,
                                   const CAsset*        asset);

    /** Invoke the method registered for this type of asset and this type of
        request */
    static IVolProcessed* dispatch(const CAsset*        asset,
                                   const CVolRequest*   request);

private:
    static void addToHashtable(CClassConstSP          volType,
                               CClassConstSP          requestType,
                               VolFunctor*            func);

    VolProcessedDispatch();
};

DRLIB_END_NAMESPACE

#endif
