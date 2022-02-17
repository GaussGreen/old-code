//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolFunctor.hpp
//
//   Description : Base class for templated class which wraps C++ pointers
//                 to member functions
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

#ifndef EDR_VOL_FUNCTOR_HPP
#define EDR_VOL_FUNCTOR_HPP
#include "edginc/VirtualDestructorBase.hpp"

DRLIB_BEGIN_NAMESPACE
class IVolProcessed;
class CVolRequest;

/** Base class for templated class which wraps C++ pointers to member functions
    targeting the getProcessedVol() method on both 'asset's and vols.
    The requirement being that these member functions take 
    one or two parameters (the VolRequest and, for getProcessedVol on the vol,
    the 'asset'). NB the volRequest type is templated.
    Here T is the class that contains the member function whilst R
    is the VolRequest. Note that when the member function takes only one 
    parameter (eg CAsset::getProcessedVol) then extraParam is null */
class MARKET_DLL VolFunctor: public virtual VirtualDestructorBase{
public:
    ~VolFunctor(); // destructor - in VolProcessedDispatch.cpp
    //// the getProcessedVol for vol of type V with volRequest of type R
    virtual IVolProcessed* getProcessedVol(
        const IObject*     mainObject,
        const CVolRequest* volRequest, 
        const void*        extraParam) const = 0;
protected:
    VolFunctor();  // constructor - in VolProcessedDispatch.cpp
};

typedef smartPtr<VolFunctor> VolFunctorSP;
DRLIB_END_NAMESPACE

#endif
