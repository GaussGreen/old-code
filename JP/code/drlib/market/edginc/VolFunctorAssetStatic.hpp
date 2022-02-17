//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolFunctorAssetStatic.hpp
//
//   Description : Templated class which wraps C++ pointers to static member
//                 functions targeting the getProcessedVol() method on assets.
//                 The requirement being that these static functions take 
//                 two parameters (the asset and the Request). 
//                 NB the volRequest type is templated.
//
//   Author      : Mark A Robson
//
//   Date        : 25 May 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VOL_FUNCTOR_ASSET_STATIC_HPP
#define EDR_VOL_FUNCTOR_ASSET_STATIC_HPP
#include "edginc/VolFunctor.hpp"

DRLIB_BEGIN_NAMESPACE
class CAsset;
class IVolProcessed;

/** Templated class which wraps C++ pointers to static member
    functions targeting the getProcessedVol() method on assets.
    The requirement being that these static functions take 
    two parameters (the asset and the Request). 
    NB the volRequest type is templated. */
template<class P, class A,
         class R> class VolFunctorAssetStatic: public VolFunctor{
    // C pointer to static function
    typedef P* (Func)(const A*, const R*);
    Func* func;      // Pointer to static function
public:
    // constructor
    VolFunctorAssetStatic(Func function): func(function) {}
    ~VolFunctorAssetStatic(){} // and destructor
    //// the getProcessedVol for vol of type V with volRequest of type R
    //// the getProcessedVol for vol of type V with volRequest of type R.
    //// the extraParam should be ignored - it just allows us to derive from
    //// VolFunctor
    virtual IVolProcessed* getProcessedVol(const IObject*     asset,
                                           const CVolRequest* volRequest, 
                                           const void*        extraParam) const{
        const A* myAsset = (const A*) A::TYPE->staticCast(asset);
        const R* myVolRequest = static_cast<const R*>(volRequest);
        return func(myAsset, myVolRequest);
    }
};

DRLIB_END_NAMESPACE

#endif
