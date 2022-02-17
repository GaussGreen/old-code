//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StateVariableGen.hpp
//
//   Description : Interfaces of generators and users of SVs
//
//   Date        : March 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_STATEVARIABLEGEN_HPP
#define EDR_STATEVARIABLEGEN_HPP

#include "edginc/smartPtr.hpp"
#include "edginc/Object.hpp"
#include "edginc/StateVariable.hpp"
#include <map>

DRLIB_BEGIN_NAMESPACE

/** Interface implemented by all classes that can generate state variables
    (whether they be 'elementary' or 'derived' state variables). Here
    derived state variables are those state variables which, as far as the
    path generator is concerned, are actually a function of some other
    state variable(s). Elementary state variables are, on the other hand,
    have no such dependency.
    Note we do not any support for reference counting here. This is so
    that objects derived from IObject can also derive from this interface
    without causing problems */
class RISKMGR_DLL IStateVariableGen {
public:
    /** Interface that models implement if they support state variables */
    class RISKMGR_DLL IStateGen {
    public:
        /** Virtual destructor */
        virtual ~IStateGen() {};

        /** Creates a state variable corresponding to a generator */
        virtual IStateVariableSP create(const IStateVariableGen* svGen) = 0;

        /** Returns true if this path generator is being used to 'simulate'
        the past */
        virtual bool doingPast() const = 0;
    };
    
    //IStateVariableGen(): refCount(0){}
    virtual ~IStateVariableGen() {};

    /** Create the corresponding State Variable for this State Variable
        Generator (NB Implies one state variable per generator).
        The previous IStateVariableSP (may be null) should be passed in.
        The return object may or may not be the same as oldStateVar. 
        Note that whether or not oldStateVar is null does NOT imply
        the value of pathGen->doingPast(). In particular, if there is no
        past to do then the first time the create method will be called the
        return value of pathGen->doingPast() will be false and oldStateVar
        will be null. */
    virtual IStateVariableSP create(IStateVariableSP oldStateVar,
                                 IStateGen*     stateGen) const = 0;
};


/** Container for pairs of state variable generators and state variables*/
class RISKMGR_DLL StateVarDBase {
public:
    /** Appends a state variable to the database */
    void append(const IStateVariableGen* svGen, IStateVariableSP sv);

    /** Obtains the state variable corresponding to a generator */
    IStateVariableSP find(const IStateVariableGen* svGen) const;

    /** Retrieve all advanceable state variables */
    void filterAdvanceables(vector<IAdvanceableStateVariable*>& store);

    /** Destructor */
    ~StateVarDBase() {}

//private:
    /** The collection of StateVarGen, StateVar */
    typedef map<const IStateVariableGen*, IStateVariableSP> DBase;
    DBase svDBase;  //!< Collection of Generators + statevars
};

DRLIB_END_NAMESPACE

#endif
