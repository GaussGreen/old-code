//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolFunctorStatic.hpp
//
//   Description : Templated class which wraps C++ pointers to static member
//                 functions targeting the getProcessedVol() method.
//                 The requirement being that these static functions take 
//                 three parameters (the vol, the Request and the asset). 
//                 NB the volRequest type is templated.
//
//   Author      : Mark A Robson
//
//   Date        : 25 May 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VOL_FUNCTOR_STATIC_HPP
#define EDR_VOL_FUNCTOR_STATIC_HPP
#include "edginc/VolFunctor.hpp"

DRLIB_BEGIN_NAMESPACE
class CAsset;
class IVolProcessed;

/** Templated class which wraps pointers to C style static functions
    targeting the getProcessedVol() method.
    The requirement being that these functions take 
    three parameters (the Vol, the VolRequest and the extraParam (eg CAsset)).
    NB the volRequest type is templated.
    Here P is a class derived from IVolProcessed, V is the vol, 
    R is the VolRequest and A is the type of the extra parameter */
template<class P, class V, 
         class R, class A> class VolFunctorStatic: public VolFunctor{
    // C pointer to static function
    typedef P* (Func)(const V*, const R*, const A*);
    Func* func;      // Pointer to static function
public:
    // constructor
    VolFunctorStatic(Func function): func(function) {}
    ~VolFunctorStatic(){} // and destructor
    //// the getProcessedVol for vol of type V with volRequest of type R
    virtual IVolProcessed* getProcessedVol(const IObject*     vol,
                                           const CVolRequest* volRequest, 
                                           const void*        extraParam) const{
        const V* myVol = (const V*) V::TYPE->staticCast(vol);
        const R* myVolRequest = static_cast<const R*>(volRequest);
        const A* myExtraParam = reinterpret_cast<const A*>(extraParam);
        return func(myVol, myVolRequest, myExtraParam);
    }
};

DRLIB_END_NAMESPACE

#endif
