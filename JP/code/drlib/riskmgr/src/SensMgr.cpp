
#include "edginc/config.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/TweakQualifierID.hpp"
#include "edginc/ObjectIteration.hpp"
#include "edginc/SensControlPerName.hpp"
#include "edginc/NoSubjectsFoundException.hpp"
#include ext_hash_set


DRLIB_BEGIN_NAMESPACE
ITweakID::ITweakID(){}
ITweakID::~ITweakID(){}
ITweakNameListID::ITweakNameListID(){}
ITweakNameListID::~ITweakNameListID(){}
ITweakNameResolver::~ITweakNameResolver(){}
ITweakNameResolver::ITweakNameResolver(){}
ITweakNameResolution::ITweakNameResolution(){}
ITweakNameResolution::~ITweakNameResolution(){}
ITweakOptID::ITweakOptID(){}
ITweakOptID::~ITweakOptID(){}
ITweakQualifierID::ITweakQualifierID(){}
ITweakQualifierID::~ITweakQualifierID(){}
ITweakTypeID::ITweakTypeID(){}
ITweakTypeID::~ITweakTypeID(){}

// hash function for IObjectSP, using only pointer value
struct ObjectPtrHash {
    size_t operator()(IObjectConstSP p) const{
        return (size_t)p.get();
    }
};
// set of objects using only pointer value for hash/equals
typedef hash_set<IObjectConstSP, ObjectPtrHash> ObjectHashSet;

/** Implements generic 'external' tweaking of instrument and/or market data. */

/** Support for collecting all the names of the different instances of
    an object. The invoke method below is called for each object of the
    right type */
class SensMgrConst::CollectName: virtual public ObjectIteration::IActionConst{
public:
    ITweakNameListID*   sens;
    OutputNameArraySP   namesList;
    ObjectHashSet       foundObjects; // avoid reporting anything twice

    CollectName(ITweakNameListID* sens): sens(sens){
        namesList.reset(new OutputNameArray());
    }
    
    /** method invoked by recurse routine to record name of object */
    bool invoke(const ObjectIteration::State& state, IObjectConstSP obj){
        if (foundObjects.find(obj) == foundObjects.end()){
            sens->appendName(*namesList, obj);
            foundObjects.insert(obj);
        }
        return true; // keep iterating
    }
};

/** Support for collecting the qualifier (eg benchmard dates) for a
    given sensitivity and a given output name. The invoke method below
    is called for each object of the right type */
class SensMgrConst::GetQualifier: virtual public ObjectIteration::IActionConst{
public:
    ITweakQualifierID*           sens;
    ITweakNameResolver*          nameResolver; // may be null
    IObjectConstSP               qualifier; 
    OutputNameConstSP            nameToMatch; // may be null

    /** Constructor */
    GetQualifier(ITweakQualifierID* sens): 
        sens(sens), nameResolver(sens->nameResolver()){
        if (nameResolver){
            nameToMatch = nameResolver->getMarketDataName();
        }
    }

    /** method invoked by recurse routine each time we hit an object
        of the right type */
    bool invoke(const ObjectIteration::State& state, IObjectConstSP obj){
        if (!nameToMatch || nameResolver->nameMatches(*nameToMatch, obj)){
            // get qualifier
            qualifier = sens->qualifier(obj);
            /* stop searching if non null - assume that first
               occurence will give us all information */
            if (qualifier.get()){
                state.quitRecursion(true /* exit immediately */);
                return false; // for clarity (return value here irrelevant)
            }
        }
        return true; // keep iterating
    }
};
/** Support for shifting an invidual object. The invoke method below
    is called for each object of the right type. This class does not make
    use of any Restorable shift methods that the objects might implement */
class SensMgr::Shift: virtual public ObjectIteration::IAction{
public:
    ////// fields /////
    ITweakID*           sens;
    ITweakNameResolver* nameResolver; // may be null
    OutputNameConstSP   nameToMatch;
    bool                shifted;
    ObjectHashSet       shiftedObjects; // avoid shifting anything twice
    /** Constructor */
    Shift(ITweakID*         sens,
          OutputNameConstSP overrideName): // may be null
        sens(sens), nameResolver(sens->nameResolver()), shifted(false){
        if (nameResolver){
            nameToMatch = !overrideName? 
                nameResolver->getMarketDataName():overrideName;
        } else if (overrideName.get()){
            throw ModelException("SensMgr::Shift", "Override name "
                                 "specified but but nameResolver is null");
        }
    }

    /** method invoked by recurse routine to shift object */
    bool invoke(ObjectIteration::State& state, IObjectConstSP obj){
        bool keepShifting = true;
        if (!nameToMatch || nameResolver->nameMatches(*nameToMatch, obj)){
            // get access to non const object
            IObjectSP theObject(state.getObject());
            // and then shift, if we haven't shifted it yet
            if (shiftedObjects.find(theObject) == shiftedObjects.end()){
                keepShifting = sens->shift(theObject);
                shiftedObjects.insert(theObject);
            } else {
                keepShifting = false; // object (and components) already tweaked
            }
            shifted = true;
        }
        return keepShifting;
    }
};

/** takes a reference to the given object to allow tweaking etc */
SensMgrConst::SensMgrConst(const IObject*  topObject): 
    topObject(IObjectSP::attachToRef(const_cast<IObject*>(topObject))){}

/** takes a reference to the given object to allow tweaking etc */
SensMgrConst::SensMgrConst(IObjectConstSP  topObject): 
    topObject(IObjectSP::constCast(topObject)){}

SensMgrConst::~SensMgrConst(){}

/** Returns true if object (or any of its components) derives from or implements
    the supplied target class. */
bool SensMgrConst::dependsUpon(IObjectConstSP   object,
                               CClassConstSP    targetClass){
    // need to fix ObjectIteration to avoid const cast
    return ObjectIteration::dependsUpon(IObjectSP::constCast(object),
                                        targetClass, true);
}

/** Returns array of output names which can be tweaked for given
    SensControl - in particular, if the SensControl has its own names
        they are ignored */
OutputNameArrayConstSP SensMgrConst::allNames(ITweakNameListID*   sens){
    try {
        /* find the interface objects need to support */
        CClassConstSP requiredIface = sens->shiftInterface();
        // create our Action class - handles each object of desired type
        CollectName action(sens);
        // create the ObjectIteration class
        ObjectIteration iteration(requiredIface);
        // Don't tweak transient fields
        iteration.setSkipTransient(true);
        // then recurse over all the relevant objects
        iteration.recurse(action, topObject);
        return action.namesList;
    } catch (exception& e){
        throw ModelException(e, "SensMgr::allNames");
    }
}

/** Returns ITweakQualifierID specific object needed for qualifing the
    greeks for a given ITweakQualifierID eg benchmark dates for pointwise
    vega. Will not return null object */
IObjectConstSP SensMgrConst::qualifier(ITweakQualifierID*   sens){
    static const string routine = "SensMgr::qualifier";
    try{
        /* find the interface objects need to support */
        CClassConstSP requiredIface = sens->shiftInterface();

        // create our Action class - handles each object of desired type
        GetQualifier action(sens);
        // create the ObjectIteration class
        ObjectIteration iteration(requiredIface);
        // Don't tweak transient fields
        iteration.setSkipTransient(true);
        
        // then recurse over all the relevant objects
        iteration.recurse(action, topObject);
        
        // check we actually got something
        if (!action.qualifier){
            throw ModelException(routine,
                                 "Could not find qualifier when looking for "
                                 "objects of type "+requiredIface->getName());
        }
        return action.qualifier;
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** Support for first() */
class SensMgrConst::CollectObject: virtual public ObjectIteration::IActionConst{
public:
    const ITweakNameResolver*    nameResolver; // may be null
    IObjectConstSP               it; 
    OutputNameConstSP            nameToMatch; // may be null

    /** Constructor */
    CollectObject(const ITweakNameResolver* nameResolver): 
        nameResolver(nameResolver){
        if (nameResolver){
            nameToMatch = nameResolver->getMarketDataName();
        }
    }

    /** method invoked by recurse routine each time we hit an object
        of the right type */
    bool invoke(const ObjectIteration::State& state, IObjectConstSP obj){
        if (!nameToMatch || const_cast<ITweakNameResolver*>(nameResolver)->
            nameMatches(*nameToMatch, obj)){
            // get qualifier
            it = obj;
            /* stop searching if non null - assume that first
               occurence will give us all information */
            if (it.get()){
                state.quitRecursion(true /* exit immediately */);
                return false; // for clarity (return value here irrelevant)
            }
        }
        return true; // keep iterating
    }
};

IObjectConstSP SensMgrConst::theFirst(CClassConstSP type,
                                      const ITweakNameResolver* id) const {
    static const string routine = "SensMgr::theFirst";
    try{
        /* find the interface objects need to support */
        // create our Action class - handles each object of desired type
        CollectObject action(id);
        
        // create the ObjectIteration class
        ObjectIteration iteration(type);
        // Don't tweak transient fields
        iteration.setSkipTransient(true);
        // then recurse over all the relevant objects
        iteration.recurse(action, topObject);
        
        // check we actually got something
        if (!action.it){
            throw NoSubjectsFoundException(string() +
                "Found no objects of type '" + type->getName() + "' " +
                (!action.nameToMatch ?
                     "" : "called " + action.nameToMatch->toString() + " ") +
                "in the supplied '" + topObject->getClass()->getName() + "'");
        }
        return action.it;
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** Support for all() */
class SensMgrConst_CollectObjects: virtual public ObjectIteration::IActionConst{
public:
    const ITweakNameResolver*    nameResolver; // may be null
    ObjectArraySP                them; 
    OutputNameConstSP            nameToMatch; // may be null

    /** Constructor */
    SensMgrConst_CollectObjects(const ITweakNameResolver* nameResolver): 
        nameResolver(nameResolver), them(new ObjectArray()) {
        if (nameResolver){
            nameToMatch = nameResolver->getMarketDataName();
        }
    }

    /** method invoked by recurse routine each time we hit an object
        of the right type */
    bool invoke(const ObjectIteration::State& state, IObjectConstSP obj){
        if (!nameToMatch || const_cast<ITweakNameResolver*>(nameResolver)->
            nameMatches(*nameToMatch, obj)){
            // we return an ObjectArrayConstSP (although that still not right)
            them->push_back(IObjectSP::constCast(obj));
        }
        return true; // keep iterating
    }
};

ObjectArrayConstSP SensMgrConst::all(CClassConstSP type,
                                     const ITweakNameResolver* id) const {
    static const string routine = "SensMgr::all";
    try{
        /* find the interface objects need to support */
        // create our Action class - handles each object of desired type
        SensMgrConst_CollectObjects action(id);
        
        // create the ObjectIteration class
        ObjectIteration iteration(type);
        // Don't tweak transient fields
        iteration.setSkipTransient(true);
        // then recurse over all the relevant objects
        iteration.recurse(action, topObject);
        
        return action.them;
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** takes a reference to the given object to allow tweaking etc */
SensMgr::SensMgr(IObject*  topObject): 
    SensMgrConst(topObject),  shifted(false){}

/** takes a reference to the given object to allow tweaking etc */
SensMgr::SensMgr(IObjectSP topObject): SensMgrConst(topObject), shifted(false){}

SensMgr::~SensMgr(){}

IObjectSP SensMgr::theFirst_mutable(CClassConstSP type,
                                    const ITweakNameResolver* id) const {
    return IObjectSP::constCast(theFirst(type, id));
}

/** Shifts data inside object as indicated by sens. There is no ability
    to restore the object once shifted */
IObjectSP SensMgr::shift(ITweakID*         sens){
    return shift(sens, OutputNameConstSP());
}

/** Shifts data inside object as indicated by sens. There is no ability
    to restore the object once shifted */
IObjectSP SensMgr::shift(ITweakID*         sens,
                         OutputNameConstSP overrideName){ // may be null
    try{
        // clear out the sens control
        sens->reset();
        /* find the interface objects need to support */
        CClassConstSP requiredIface = sens->shiftInterface();
        // create our Action class
        Shift action(sens, overrideName);
        // create the ObjectIteration class
        ObjectIteration iteration(requiredIface);
        // Dont tweak transient fields
        iteration.setSkipTransient(true);
        // then recurse over all the relevant objects
        iteration.recurse(action, topObject);
        shifted = action.shifted;
        return topObject;
    } catch (exception& e){
        throw ModelException(&e, "SensMgr::shift");
    }
}
                                
    
/** Indicates whether the object supplied was actually changed by the
    shift request */
bool SensMgr::getShiftStatus() const{
    return shifted;
}

/////////////// SensMgrOpt class //////////////////

/** simple subclass for handling exceptions we can deal with */
class ShiftException: public ModelException{
public:
    ShiftException(exception& e): ModelException(e){}
};

/** Support for shifting an individual object. The invoke method below
    is called for each object of the right type. The end result is the
    same as the Shift class above except this implementation is
    optimised to use any restore methods implemented by the objects so
    as to keep copying to a minimum */
class SensMgrOpt::OptimalShift: virtual public ObjectIteration::IAction{
public:
    ////// fields /////
    ITweakID*           sens;
    ITweakOptID*        optSens;
    ITweakNameResolver* nameResolver; // may be null
    OutputNameConstSP   nameToMatch;
    bool                shifted; // true if object shifted at all
    bool                inDeepCopy; // true if within a deep copied object
    int                 deepCopyDepth; // when we make the deep copy
    vector<bool>        shiftReturnValues;
    ObjectArray         clonedObjs; /* copies of objects which do not implement
                                       restore method */
    bool                inError;    // if a shift failed
    ModelException      theError;   // what the error was
    ObjectHashSet       shiftedObjects; // avoid shifting anything twice

    /** Constructor */
    OptimalShift(ITweakID*         sens, 
                 ITweakOptID*      optSens,       // may be null
                 OutputNameConstSP overrideName): // may be null
        sens(sens), optSens(optSens), nameResolver(sens->nameResolver()),
        shifted(false), inDeepCopy(false), deepCopyDepth(0), inError(false){
        if (nameResolver){
            nameToMatch = !overrideName? 
                nameResolver->getMarketDataName():overrideName;
        } else if (overrideName.get()){
            throw ModelException("SensMgrOpt::OptimalShift", "Override name "
                                 "specified but but nameResolver is null");
        }
    }

    /** method invoked by recurse routine to shift object */
    bool invoke(ObjectIteration::State& state, IObjectConstSP obj){
        static const string routine = "OptimalShift::invoke";
        bool keepShifting = true;
        try{
            if (!nameToMatch || nameResolver->nameMatches(*nameToMatch,obj)){
                if (state.recursingDownward()){
                    // get non-const access to object to change/replace etc
                    IObjectSP object(state.getObject());
                    // extend count of shifts applied (in case of failure)
                    shiftReturnValues.push_back(false);                
                    // avoid shifting anything twice
                    if (shiftedObjects.find(object) == shiftedObjects.end()){
                        /* need to shift object but first must decide
                           whether need to take a copy or not */
                        if (!inDeepCopy && 
                            (!optSens || !optSens->restorableShift(object))){
                            IObjectSP clone;
                            try{  // take deep copy
                                clone = IObjectSP(object.clone());
                            } catch (exception& e){
                                throw ShiftException(e);
                            }
                            deepCopyDepth = state.getDepth();
                            // prefer to modify copy rather than original -
                            // this preserves references to internal object
                            clonedObjs.push_back(object); // save original
                            object = clone; // switch to clone
                            state.setObject(clone); // replace with copy
                            inDeepCopy = true;
                        }
                        try{  // then shift the object
                            keepShifting = sens->shift(object);
                            shiftedObjects.insert(object);
                        } catch (exception& e){
                          throw ShiftException(e);
                        }
                    } else {
                        keepShifting = false; /* object (and components) 
                                                 already tweaked */
                    }
                    shiftReturnValues.back() = keepShifting; // save value
                    shifted = true;
                } else {
                    /* as we come back up through the objects turn off the 
                       inDeepCopy flag when we get the to the original depth */
                    if (inDeepCopy && deepCopyDepth == state.getDepth()){
                        inDeepCopy = false;
                    }
                }
            }
        } catch (ShiftException& e){
            theError = ModelException(e, routine);
            inError = true; // flag to restore method that we had a failure
            state.quitRecursion(true); // exit recursion immediately
        } catch (exception& e){
            // can't handle this - give up - need to throw stronger error
            throw ModelException(e, routine);
        }
        return keepShifting;
    }
};

/** Support for restoring an object which has been shifted using the
    above OptimalShift class */
class SensMgrOpt::OptimalRestoreShift: public OptimalShift{
public:
    ////// fields /////
    unsigned int       posShift; // pos in shiftReturnValues array
    int                posClone; // pos in clonedObjs array

    /** Constructor */
    OptimalRestoreShift(OptimalShift* optimalShift): 
        OptimalShift(*optimalShift), posShift(0), posClone(0) {
        inDeepCopy = false; // reset inherited fields 
        shiftedObjects.clear(); // empty contents of hash_set
    }

    /** method invoked by recurse routine to restore an object after an
        'optimal shift' */
    bool invoke(ObjectIteration::State& state, IObjectConstSP obj){
        static const string routine = "OptimalRestoreShift::invoke";
        bool keepShifting = true;
        if (!nameToMatch || nameResolver->nameMatches(*nameToMatch, obj)){
            try{
                if (state.recursingDownward()){
                    // first sort out deep copy flag
                    if (!inDeepCopy && 
                        (!optSens || !optSens->restorableShift(obj))){
                        // we took deep copy in OptimalShift method
                        deepCopyDepth = state.getDepth();
                        inDeepCopy = true;
                    }
                    /* return the same value as when we were shifting */
                    keepShifting = shiftReturnValues[posShift];
                    posShift++;
                    if (inError && posShift == shiftReturnValues.size()){
                        // where we got to in shifting - can now
                        // exit by working our way back up only
                        state.quitRecursion(false /* 'weak' quit */);
                    }
                } else {  // restore the objects as we come back up
                    bool hitErrorObj = 
                        inError && posShift == shiftReturnValues.size();
                    if (!inDeepCopy){
                        if (!optSens || !optSens->restorableShift(obj)){
                            throw ModelException(routine, "Internal error 1");
                        }
                        if (!hitErrorObj){
                            // easy, just get access non const object 
                            // and then get it to restore itself
                            IObjectSP theObject(state.getObject());
                            // but only if we haven't restored it already
                            if (shiftedObjects.find(theObject) == 
                                shiftedObjects.end()){
                                optSens->restore(theObject);
                                shiftedObjects.insert(theObject);
                            }
                        }
                    } else if (deepCopyDepth == state.getDepth()){
                        inDeepCopy = false;
                        // replace object with deep copy (only if clone worked)
                        if (posClone < clonedObjs.size()){
                            state.setObject(clonedObjs[posClone]);
                            posClone++;
                        }
                    }
                    if (hitErrorObj){
                        inError = false; // clear for remainder
                    }
                }
            } catch (exception& e){
                throw ModelException(e, routine);
            }
        }
        return keepShifting;
    }
};
  
/** takes a reference to the given object to allow tweaking etc */
SensMgrOpt::SensMgrOpt(IObject*  topObject, bool optimalShift): 
    SensMgr(topObject), doOptimalShift(optimalShift), optShift(0){
    // empty
}

/** Same as above but takes smart pointer */
SensMgrOpt::SensMgrOpt(IObjectSP topObject, bool optimalShift):
    SensMgr(topObject), doOptimalShift(optimalShift), optShift(0){
    // empty
}

/** Equivalent to SensMgrOpt(topObject, true) */
SensMgrOpt::SensMgrOpt(IObject*  topObject): 
    SensMgr(topObject), doOptimalShift(true), optShift(0){
    // empty
}

/** Equivalent to SensMgrOpt(topObject, true) */
SensMgrOpt::SensMgrOpt(IObjectSP  topObject): 
    SensMgr(topObject), doOptimalShift(true), optShift(0){
    // empty
}

SensMgrOpt::~SensMgrOpt(){
    delete optShift;
}

/** Shifts data inside object as indicated by sens. Shifted object may
    be copy of original or may be modified version of original. Use
    restore to return object to its original state */
IObjectSP SensMgrOpt::shift(ITweakOptID*    sens){
    return shift(sens, OutputNameConstSP());
}

//// same as above but switch on run time instance of sens
IObjectSP SensMgrOpt::shift(ITweakID*    sens){
    return shift(sens, OutputNameConstSP());
}

/** same as above shift method for ITweakID but name for
    overrideName parameter is used rather than a call to
    ITweakNameResolver. If overrideName is null then that is
    treated the same nameResolver returning null */
IObjectSP SensMgrOpt::shift(ITweakID* sens, OutputNameConstSP overrideName){
    sens->reset();
    if (!doOptimalShift){
        return copyAndShift(sens, overrideName);
    }
    ITweakOptID* sensOpt = dynamic_cast<ITweakOptID*>(sens);
    return optimalShift(sens, sensOpt /* may be null */, overrideName);
}    

/** Same as above but avoids the internal cast */
IObjectSP SensMgrOpt::shift(ITweakOptID* sens, OutputNameConstSP overrideName){
    // clear out the sens control
    sens->reset();
    return doOptimalShift? 
        optimalShift(sens, sens, overrideName): copyAndShift(sens,overrideName);
}

/** Shifts data inside object as indicated by sens. Does minimal copying
    of data to allow object to be restored to its original state */
IObjectSP SensMgrOpt::optimalShift(
    ITweakID*         sens,
    ITweakOptID*      optSens,       // may be null
    OutputNameConstSP overrideName){ // may be null
    static const string routine = "SensMgrOpt::optimalShift";
    try{
        /* find the interface objects need to support */
        CClassConstSP requiredIface = sens->shiftInterface();
        // create our Action class
        optShift = new OptimalShift(sens, optSens, overrideName);
        // create the ObjectIteration class
        ObjectIteration iteration(requiredIface);
        // Don't tweak transient fields
        iteration.setSkipTransient(true);
        // ask for action method to be invoked in each direction of iteration
        iteration.invokeTwice();
        // then recurse over all the relevant objects
        try{
            topObject = iteration.recurse(*optShift, topObject);
        } catch (exception& e){
            // NB any exceptions caused by shift/clone methods are handled
            // explicitly and aren't propagated to here
            // throw derived exception? (terminal - abort pricing)
            throw ModelException(e, routine, "Failed to restore object");
        }
        if (optShift->inError){
            try{
                // need to complete restoration of object
                optimalShiftRestore();
            } catch (exception& e){
                optShift->inError = true; // flag that we had a failure
                // throw derived exception? (terminal - abort pricing)
                throw ModelException(e, routine, 
                                     "Failed to restore object");
            }
            optShift->inError = true; // flag that we had a failure
            throw optShift->theError;
        }
        shifted = optShift->shifted;
        return topObject;
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}
 
/** Shifts data inside object as indicated by sens. Shifted object is
    copy of original. The restore method in this case performs no actions */
IObjectSP SensMgrOpt::copyAndShift(
    ITweakID*         sens,
    OutputNameConstSP overrideName){ // may be null
    try{
        IObjectSP topCopy = IObjectSP(topObject->clone());
        /* find the interface objects need to support */
        CClassConstSP requiredIface = sens->shiftInterface();
        // create our Action class
        Shift action(sens, overrideName);
        // create the ObjectIteration class
        ObjectIteration iteration(requiredIface);
        // Don't tweak transient fields
        iteration.setSkipTransient(true);
        // then recurse over all the relevant objects
        iteration.recurse(action, topCopy);
        shifted = action.shifted;
        return topCopy;
    } catch (exception& e){
        throw ModelException(e, "SensMgrOpt::copyAndShift");
    }
}
/** restores an object after a tweak */
void SensMgrOpt::restore(){
    // if optShift->inError is true then have restored object already or 
    // there was an unrecoverable error
    if (doOptimalShift && optShift && !optShift->inError){
        optimalShiftRestore();
    }
}

/** restores an object after a tweak */
void SensMgrOpt::optimalShiftRestore(){
    try{
        /* find the interface objects need to support */
        CClassConstSP requiredIface = optShift->sens->shiftInterface();
        // create our Action class using the data gathered when shifting
        OptimalRestoreShift action(optShift);
        // create the ObjectIteration class
        ObjectIteration iteration(requiredIface);
        // Don't tweak transient fields
        iteration.setSkipTransient(true);
        // ask for action method to be invoked in
        // each direction of the iteration
        iteration.invokeTwice();
        // then recurse over all the relevant objects
        iteration.recurse(action, topObject);
    } catch (exception& e){
        throw ModelException(&e, "SensMgrOpt::restore");
    }
}

DRLIB_END_NAMESPACE
