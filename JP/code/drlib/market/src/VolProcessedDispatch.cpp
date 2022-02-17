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

#include "edginc/config.hpp"
#if defined(_MSC_VER)
// disable warning truncated decorated names information
#pragma warning(disable : 4503)
#endif
#include "edginc/VolProcessedDispatch.hpp"
#include "edginc/VolRequest.hpp"
#include "edginc/Asset.hpp"
#include ext_hash_map

DRLIB_BEGIN_NAMESPACE
VolFunctor::~VolFunctor(){}
VolFunctor::VolFunctor(){}

// hash table with <vol request class, VolFunctorSP> entries
typedef hash_map<CClassConstSP, VolFunctorSP, CClass::Hash> VolFunctorTable;

// hash table with <vol class, hash of VolFunctor by vol request class> entries
typedef hash_map<CClassConstSP, VolFunctorTable, CClass::Hash > DispatchTable;

static DispatchTable dispatchTable;

void VolProcessedDispatch::addToHashtable(CClassConstSP volType,
                                          CClassConstSP requestType,
                                          VolFunctor*   func){
    VolFunctorSP& functor = dispatchTable[volType][requestType];
    if (functor.get()){
        throw ModelException("VolProcessedDispatch::addToHashtable",
                             "Method already registered for vol of "
                             "type "+volType->getName()+" and request of type "+
                             requestType->getName());
    }
    functor.reset(func);
}

/** Look up method registered for this type of vol or asset and this type of
    request */
static VolFunctorSP lookUp(const IObject*       mainObject,
                           const CVolRequest*   request){
    static const string method("VolProcessedDispatch::dispatch");
    CClassConstSP mainType = mainObject->getClass();
    /*
      Algorithm for looking for method:
      For current type, look for method for type of given VolRequest.
      Repeat for parent types of VolRequest until you hit VolRequest itself.
      Then move to parent asset/vol and repeat process.
      One you identify method, call it
    */
    DispatchTable::const_iterator mainIter;
    do {
        DispatchTable::const_iterator mainIter = dispatchTable.find(mainType);
        if (mainIter != dispatchTable.end()){
            const VolFunctorTable& functorTable = mainIter->second;
            VolFunctorTable::const_iterator functorIter;
            CClassConstSP requestType = request->getClass();
            while (requestType != 0 &&
                   (functorIter = functorTable.find(requestType)) == 
                   functorTable.end()){
                // not found for this type of VolRequest. Try parent (it doesn't
                // really matter if we go past VolRequest we're going to fail)
                requestType = requestType->getSuperClass();
            }
            if (requestType){
                // found!
                try{
                    return functorIter->second;
                } catch (exception& e){
                    throw ModelException(e, method, "Object of type "+
                                         mainObject->getClass()->getName()+
                                         " failed to process volatility request"
                                         " of type "+
                                         request->getClass()->getName());
                }
            }
        }
        // not found for this type of mainObject - try parent
    } while ((mainType = mainType->getSuperClass()) != 0);
    // no method
    throw ModelException(method, "Object of type "+
                         mainObject->getClass()->getName()+" does not support "
                         "request of type "+request->getClass()->getName());
}

/** Invoke the method registered for this type of vol and this type of
    request */
IVolProcessed* VolProcessedDispatch::dispatch(const CVolBase*      vol,
                                              const CVolRequest*   request,
                                              const CAsset*        asset){
    VolFunctorSP functor(lookUp(vol, request));
    return functor->getProcessedVol(vol, request, asset);
}

/** Invoke the method registered for this type of asset and this type of
    request */
IVolProcessed* VolProcessedDispatch::dispatch(const CAsset*        asset,
                                              const CVolRequest*   request){
    VolFunctorSP functor(lookUp(asset, request));
    return functor->getProcessedVol(asset, request, 0 /* not used */);
}

DRLIB_END_NAMESPACE

