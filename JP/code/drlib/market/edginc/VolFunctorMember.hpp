//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolFunctorMember.hpp
//
//   Description : Templated class which wraps C++ pointers to member functions
//                 targeting the getProcessedVol() method.
//                 The requirement being that these member functions take 
//                 two parameters (the VolRequest and the asset). NB the 
//                 volRequest type is templated.
//
//   Author      : Mark A Robson
//
//   Date        : 25 May 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VOL_FUNCTOR_MEMBER_HPP
#define EDR_VOL_FUNCTOR_MEMBER_HPP
#include "edginc/VolFunctor.hpp"

DRLIB_BEGIN_NAMESPACE
class CAsset;
class IVolProcessed;

/** Templated class which wraps C++ pointers to member functions
    targeting the getProcessedVol() method.
    The requirement being that these member functions take 
    two parameters (the VolRequest and the extraParm (eg CAsset)). NB the 
    volRequest type is templated.
    Here P is a class derived from IVolProcessed, here T is the class that 
    contains the member function whilst R is the VolRequest and A is the
    type of the extra parameter */
template<class P, class V, 
         class R, class A> class VolFunctorMember: public VolFunctor{
public:
    //C++ pointer to member function
    typedef P* (V::* Func)(const R*, const A*) const;
    // constructor
    VolFunctorMember(Func func): func(func) {}
    ~VolFunctorMember(){} // and destructor
    //// the getProcessedVol for vol of type V with volRequest of type R
    virtual IVolProcessed* getProcessedVol(const IObject*     vol,
                                           const CVolRequest* volRequest, 
                                           const void*        extraParam) const{
        const V* myVol = (const V*) V::TYPE->staticCast(vol);
        const R* myVolRequest = static_cast<const R*>(volRequest);
        const A* myExtraParam = reinterpret_cast<const A*>(extraParam);
        return (myVol->*func)(myVolRequest, myExtraParam);
    }
private:
    Func func;      // Pointer to member function
};

DRLIB_END_NAMESPACE

#endif
