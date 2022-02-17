//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StateVariableCollector.hpp
//
//   Description : Class for 'collecting' 'true' state variables (actually 
//                 the generators) - sort of a glorified array
//
//   Date        : March 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_STATEVARIABLECOLLECTOR_HPP
#define EDR_STATEVARIABLECOLLECTOR_HPP

#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/ElemStateVariableGen.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for classes that can 'collect' 'true' state variables (actually
    the generators) - sort of a glorified array. We have made this an 
    interface in case any path generators want to have their own
    handling. */
class RISKMGR_DLL IStateVariableCollector: public virtual VirtualDestructorBase {
public:
    IStateVariableCollector();
    virtual ~IStateVariableCollector();

    /** Appends 'underlying' state variable to the list of 'true' (ie
        non-derived) variables. Note that only references are taken so
        the supplied object must not go out of scope. Also, when retrieving
        IStateVar from path generators the same pointer must be used as
        supplied to this function */
    virtual void append(const IStateVariableGen* svGen) = 0;

    /** Appends an elementary state variable to the list of 'true' (ie
        non-derived) variables. This method should be invoked by state
        variable generators that are both IMCStateVarClient (users) and
        IElemStateVariableGen when they want to append their IElemStateVariableGen
        part. The simple append method should be called when they append their
        IMCStateVarClient part. */
    virtual void appendElementary(const IElemStateVariableGen* elemGen) = 0;

    
    /** Accessor function to array of elementary state variable generators */
    virtual const IElemStateVariableGenArray& getElemStateVarGens() const = 0;
};

typedef smartPtr<IStateVariableCollector> IStateVariableCollectorSP;
#ifndef QLIB_STATEVARIABLECOLLECTOR_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<IStateVariableCollector>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<IStateVariableCollector>);
#endif



/** Default implementation of IStateVariableCollector. */
class RISKMGR_DLL StateVariableCollector: virtual public IStateVariableCollector {
public:
    /** Appends 'underlying' state variable to the list of 'true' (ie
        non-derived) variables. Note that only references are taken so
        the supplied object must not go out of scope. Also, when retrieving
        IStateVar from path generators the same pointer must be used as
        supplied to this function */
    virtual void append(const IStateVariableGen* svGen);
    
    /** Appends an elementary state variable to the list of 'true' (ie
        non-derived) variables. This method should be invoked by state
        variable generators that are both IMCStateVarClient (users) and
        IElemStateVariableGen when they want to append their IElemStateVariableGen
        part. The simple append method should be called when they append their
        IMCStateVarClient part. */
    virtual void appendElementary(const IElemStateVariableGen* elemGen);

    /** Accessor function to array of elementary state variable generators */
    virtual const IElemStateVariableGenArray& getElemStateVarGens() const;

private:
    IElemStateVariableGenArray elemStateVarGens;  //!< Elementary state variable generators
};

typedef smartPtr<StateVariableCollector> StateVariableCollectorSP;


/** Erases from the input vector of base pointers all entries
    that can be casted to derived pointers and returns them as
    a vector casted derived pointers. 
    Assumes there is no NULL pointers in the input or that it's safe to remove them. 
*/

template<class _TBasePtr, class _TDerivedPtr>
vector<_TDerivedPtr> filterPointers(vector<_TBasePtr>& basePtrs) {
    vector<_TDerivedPtr> castedPtrs;
    
    typedef typename vector<_TBasePtr>::iterator Iterator;
    
    for(Iterator iter = basePtrs.begin(); iter != basePtrs.end(); ++iter) {
        _TBasePtr basePtr = *iter;
        _TDerivedPtr castedPtr = dynamic_cast<_TDerivedPtr>(basePtr);
        if (castedPtr != NULL) {
            castedPtrs.push_back(castedPtr); 
            *iter = NULL; // mark elements for deletion; assumes that there is no NULLs in the input or it's ok to erase
        }
    }
    // remove all elements that were marked for deletion (successfully casted) -- keeps the order
    basePtrs.erase(remove(basePtrs.begin(), basePtrs.end(), (_TBasePtr) NULL), basePtrs.end()); 

    return castedPtrs;
}

/** Wrapper around state variables filtering */
template<class _TDerivedSV>
vector<const _TDerivedSV*> filterStateVars(IElemStateVariableGenArray& stateVarGens) {
    return filterPointers<const IElemStateVariableGen*, const _TDerivedSV*>(stateVarGens);
}

DRLIB_END_NAMESPACE

#endif
