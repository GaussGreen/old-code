//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolFunctorAssetMember.hpp
//
//   Description : Templated class which wraps C++ pointers to member functions
//                 targeting the getProcessedVol() method on 'assets'
//                 The requirement being that these member functions take 
//                 one parameter namely the VolRequest. NB the 
//                 volRequest type is templated.
//
//   Author      : Mark A Robson
//
//   Date        : 25 May 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VOL_FUNCTOR_ASSET_MEMBER_HPP
#define EDR_VOL_FUNCTOR_ASSET_MEMBER_HPP
#include "edginc/VolFunctor.hpp"

DRLIB_BEGIN_NAMESPACE
class CAsset;
class IVolProcessed;

/** Templated class which wraps C++ pointers to member functions
    targeting the getProcessedVol() method.
    The requirement being that these member functions take 
    two parameters (the VolRequest and the extraParm (eg CAsset)). NB the 
    volRequest type is templated.
    Here P is a class derived from IVolProcessed, here A is the class that 
    contains the member function whilst R is the VolRequest */
template<class P, class A, class R> class VolFunctorAssetMember: public VolFunctor{
public:
    //C++ pointer to member function
    typedef P* (A::* Func)(const R*) const;
    // constructor
    VolFunctorAssetMember(Func func): func(func) {}
    ~VolFunctorAssetMember(){} // and destructor
    //// the getProcessedVol for vol of type V with volRequest of type R.
    //// the extraParam should be ignored - it just allows us to derive from
    //// VolFunctor
    virtual IVolProcessed* getProcessedVol(const IObject*     asset,
                                           const CVolRequest* volRequest,
                                           const void*        extraParam) const{
        const A* myAsset = (const A*) A::TYPE->staticCast(asset);
        const R* myVolRequest = static_cast<const R*>(volRequest);
        return (myAsset->*func)(myVolRequest);
    }
private:
    Func func;      // Pointer to member function
};

DRLIB_END_NAMESPACE

#endif
