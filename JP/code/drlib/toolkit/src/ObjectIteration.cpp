#include "edginc/config.hpp"
#if defined(_MSC_VER)
// disable warning truncated decorated names information
#pragma warning(disable : 4503)
#endif

#include "edginc/ObjectIteration.hpp"
#include "edginc/Array.hpp"
#include ext_hash_map
#include ext_hash_set

DRLIB_BEGIN_NAMESPACE

ObjectIteration::IAction::~IAction(){}
ObjectIteration::IActionConst::~IActionConst(){}
ObjectIteration::IOverride::~IOverride(){}

CClassConstSP const ObjectIteration::IOverride::TYPE =
CClass::registerInterfaceLoadMethod(
    "ObjectIteration::IOverride", typeid(ObjectIteration::IOverride), 0);

// hash function for CClassConstSP
struct CClassPtrHash {
    size_t operator()(const CClassConstSP p) const{
        return (size_t)p;
    }
};

////////////// type definitions ///////////
// hash table with <class, array of fields> entries
typedef hash_map<CClassConstSP, CFieldArray, 
    CClassPtrHash> CFieldArrayHashTable;

// hash table with <class, hash of array of fields by class> entries
typedef hash_map<CClassConstSP, CFieldArrayHashTable, 
    CClassPtrHash > TweaksHashTable;

// set of classes
typedef hash_set<CClassConstSP, CClassPtrHash> ClassHashSet;
    
// hash table with <class, set of classes> entries
typedef hash_map<CClassConstSP, ClassHashSet, 
    CClassPtrHash> SuperClassHashTable;

//// hash table containing <Class, int> pairs
typedef hash_map<CClassConstSP, int, CClassPtrHash> ClassIntHash;

//// hash table containing <Class, ClassIntHash> pairs
typedef hash_map<CClassConstSP, ClassIntHash, 
    CClassPtrHash> ComponentsHashTable;

////////////// end of type definitions ///////////

// keep implementation out of header file (and keep stl includes out of 
// header file
class ObjectIteration::Imp{
public:
    ////// static class values //////

    /* for each target type, superClasses contains a set of types which
       are superclasses of the target type ie given an instance, obj, of the
       target type, it contains all types that would satisfy 
       type->instanceOf(obj) is true. */
    static SuperClassHashTable  superClasses;

    /* for each target type, tweakable contains a hash. The hash holds
       pairs <type, int> and the int says whether or not an instance
       of the type might contain (directly or indirectly via its
       components) an instance of the target type */
    static ComponentsHashTable   components;

    /* for each target type, tweaks contains a hash table which holds
       the list of fields that need to be tweaked for each type with
       respect to target type. */
    static TweaksHashTable     tweaks;

    // the ints below are the various values in hashes inside
    // ComponentsHashTable
    /** class contains, directly or indirectly, a component which needs to be
        iterated over */
    static const int COMPONENT_YES; // = 0
    /** class does not contain, directly or indirectly, a component
        which needs to be iterated over */
    static const int COMPONENT_NO;  // = 1
    /** indicates state is yet to be finally calculated - to cope 
        for recursion */
    static const int COMPONENT_NO_SO_FAR;  // = 2
    /** equivalent to no entry */
    static const int COMPONENT_UNKNOWN;  // = 3

    /////////////// Fields /////////////////
    State*   state; /* reference back to containing class (needed as we pass
                       State object to clients from here */
    IAction* action; /* what to do when we hit an object of the target type
                        (this is for objects that can be modified) */
    IActionConst* actionConst; /* as above, but for const objects. Only one
                                  should be non-null */

    // what type to look for - use IObject to catch all types
    CClassConstSP         targetClass;

    ClassHashSet*         superClassesSet;  /* optional in conjunction with 
                                               fieldsCache. */
    CFieldArrayHashTable* fieldsCache; /* optional table used to
                                          override and/or cache fields */
    ClassIntHash*         componentsHash; /* optional in conjunction with 
                                             fieldsCache */
    int    depth;  // holds count of level of recursion

    bool   skipTransient; // default = false;

    // Holds information on what is currently being 'Actioned' on 
    CFieldConstSP currentField;
    IObjectSP currentObject; // holds the current object

    // Do we invoke action on top level object supplied
    bool actionTopObj; // default = true;

    bool quit;  /* default = false. Set to true to terminate iteration by
                   doing no more downward recursing - however invoke is
                   still called on the upward pass (if requested) */
    bool quitNow; /* default = false. Set to true to terminate entire iteration
                     immediately */
    /** Iteration recurses down over objects and then back up. This flag is
        true during the initial downwards recurse */
    bool recursingDown;

    /* true: call invoke() on downward and upward iteration. Default false */
    bool invokeTwiceFlag; 

    bool useFieldCache; /* default is true - use cache of which fields to 
                           iterate over */
    bool performSet; // true: client has invoked setObject()
    bool objectHasChanged; // true: client has altered object
    IObjectSP objToSet; // when performSet is true, holds object to set

    

    /** do we skip the given field. If override is non null then use recurse
        method to determine whether or not to skip */
    bool skipField(const IOverride*     override,
                   const CFieldConstSP& field) const{
        bool skipThisField;
        if (override){
            skipThisField = !override->recurse(field, targetClass);
        } else if (skipTransient && field->isTransientForIteration()){
            // skip transient fields if skipTransient = true
            skipThisField = true;
        } else {
            skipThisField = false;
        }
        return skipThisField;
    }

    /** This is a work around to a problem with purify 5.3 (prototype
        8) - without this if an exception occurs then the code goes
        into an infinite loop repeatedly calling action->invoke */
    bool invokeWrapperPurifyFix(const IObjectSP& cmpt) {
        currentObject = cmpt; // save, for use by getObject()

        // Dummy try/catch to ""solve"" vc6.opt crashing

		try {
            return action? 
                action->invoke(*state, cmpt):
                actionConst->invoke(*state, cmpt);
		}
		catch (...) {
            throw;
		}
    }

    /** simple wrapper around invoke to give detailed error information */
    bool invokeWrapper(const IObjectSP&    cmpt,
                       const CFieldConstSP field,  /* may be null */
                       bool                recursingDown){

        // cerr << string(depth * 2, '·') << " " << (!field ? string("no field") : field->getDeclaringClass()->getName() + "." + field->getName()) << " on " << cmpt->getClass()->getName() << " @ " << cmpt.get() << "\n";

        bool returnVal;
        this->currentField = field;
        this->recursingDown = recursingDown;
        try{
            returnVal = invokeWrapperPurifyFix(cmpt);
        } catch (exception& e){
          static const string method("ObjectIteration::invokeWrapper");
            if (field){
                throw ModelException(e, method,
                                     "Action method failed for "
                                     "field '"+field->getName()+"' of type "+
                                     cmpt->getClass()->getName()+
                                     " in object of type "+
                                     field->getDeclaringClass()->getName());
            }
            // top level object
            throw ModelException(e, method,
                                 "Action method failed for "
                                 "top level object of type "+
                                 cmpt->getClass()->getName());
        }

        return returnVal;
    }

    /** invoke method on action object and recurse. Returns true if the
        supplied object has changed - callers should check original and new
        values of cmpt.get() to see if the original object needs to be replaced
        with the new value */
    bool invokeMethod(IObjectSP& cmpt, CFieldConstSP field /* may be null */){
        if (cmpt.get()){
            bool isInstance   = targetClass->isInstance(cmpt);
#if 0
            // only here for debugging
            cout << "Checking " << cmpt->getClass()->getName() << 
                " for target " << targetClass->getName() << endl;
#endif 
            bool recurse      = true;
            bool objectChanged = false;
            if (isInstance){
                // call action->invoke()
                // Dummy try/catch to ""solve"" vc6.opt crashing
				try {
					recurse = invokeWrapper(cmpt, field, true);
                }
                catch (...) {
                    throw;
                }
                objectChanged = clearObjectChangedFlag(cmpt);
            }
            if (!quit && recurse && !isPrimitive(cmpt)) {
                if (internalRecurse(cmpt)){  // recurse over components
                    objectChanged = true; /* return value indicates whether
                                             object has changed */
                }
            }
            // now work our way back up
            if (!quitNow && invokeTwiceFlag && isInstance){
                /* call action->invoke(), ignore return val */
                invokeWrapper(cmpt, field, false);
            }
            return (clearObjectChangedFlag(cmpt) || objectChanged);
        }
        return false;
    }

    /** Clears the 'performSet' flag (which indicates whether the
        client wants to replace the current object with an alternative
        object) as well as the objectHasChanged has changed flag. True
        is returned if either flag is true. Callers should watch to
        see if obj.get() changes to see whether the original object
        needs to be replaced with the new value */
    bool clearObjectChangedFlag(IObjectSP& obj){
        if (objectHasChanged){
            objectHasChanged = false; // clear
            if (performSet){
                performSet = false;
                obj = objToSet;
            }
            return true;
        }
        return false;
    }

    /** recurse over each element in array invoking action method for
        each non null component. Skips over 'atomic' fields or arrays
        of atomic fields. The {@link Action action} class drives what
        is invoked for each component. The return flag indicates
        whether any components of this object were altered */
    bool internalRecurseArray(const IArraySP& arrayObj){
        static const string routine("ObjectIteration::internalRecurseArray");
        int  arrayLen = arrayObj->getLength();
        CFieldConstSP field = currentField; // save
        bool hasChanged = false;
        try{
            for (int i = 0; i < arrayLen && !quit; i++){
                // cerr << string(depth * 2, '·') << " " << "[" << i << "]\n";
                IObjectSP cmpt(arrayObj->get(i));
                IObject*  origPtr = cmpt.get();
                if (invokeMethod(cmpt, field)){
                    hasChanged = true;
                    // cmpt has changed - do we need to do a set?
                    if (origPtr != cmpt.get()){
                        // replace object as it's a different object
                        try{
                            arrayObj->set(i, cmpt);
                        } catch (exception& e){
                            throw ModelException(e, routine, "Failed to set "
                                                 "value in array");
                        }
                    }
                }
            }
        } catch (exception& e) {
            throw ModelException(e, routine, "Failed to iterate over array "
                                 "of type " +arrayObj->getClass()->getName());
        }
        return hasChanged;
    }

    /** recurse over each component invoking action method for each non
        null component. Skips over 'atomic' fields or arrays of 
        atomic fields. The {@link Action action} class drives what is
        invoked for each component. The return flag indicates
        whether any components of this object were altered */
    bool internalRecurseObject(const IObjectSP& obj){
        static const string routine("ObjectIteration::internalRecurse");
        try{
            const IOverride* override = 0;
            bool             doneCast = false;
            CFieldArray      updatedFields; // fields that 'have changed'
            // iterate over the type description, 
            CClassConstSP c = obj->getClass();
            // pull data for each item out
            do {
                /* get the list of fields we need to recurse
                   over. Note that when using the cache, parent fields
                   are returned as well */
                const CFieldArray& fields = fieldsCache? 
                    getFieldsCache(c): c->getDeclaredFields();
#if 0
                // only here for debugging
                cout << "For objects of type " << c->getName() << 
                    " when looking for " 
                     << targetClass->getName() << ": " << endl;
                for (unsigned int i2 = 0; i2 < fields.size(); i2++){
                    cout << fields[i2]->getName() << ": " << 
                        fields[i2]->getType()->getName() << endl;
                }
#endif
                // cerr << string(depth * 2, '·') << " - " << c->getName() << " @ " << obj.get() << "\n";
                for (unsigned int i = 0; i < fields.size() && !quit; i++) {
                    if (!doneCast){
                        // need to see if object implements IOverride interface
                        if (IOverride::TYPE->isInstance(obj.get())){
                            override = &dynamic_cast<const IOverride&>(*obj);
                        }
                        doneCast = true;
                    }
                    CFieldConstSP field = fields[i];
                    // cerr << string(depth * 2, '·') << " " << field->getDeclaringClass()->getName() << "." << field->getName() << "\n";
                    if (!skipField(override, field)) {
                        IObjectSP cmpt(field->get(obj));
                        IObject*  origPtr = cmpt.get();
                        if (invokeMethod(cmpt, field)){
                            updatedFields.push_back(field); // save
                            // cmpt has changed - do we need to do a set?
                            if (origPtr != cmpt.get()){
                                // replace object as it's a different object
                                try{
                                    field->set(obj, cmpt);
                                } catch (exception& e){
                                    throw ModelException(e, routine,
                                                         "Failed to set"
                                                         " value in object");
                                }
                            }
                        }
                    }
                }
            } while (!fieldsCache && (c = c->getSuperClass()) && !quit);
            if (!updatedFields.empty()){
                // inform object
                obj->fieldsUpdated(updatedFields);
                return true; // object has changed
            }
            return false; // object unchanged
        } catch (exception& e){
            throw ModelException(e, routine, 
                                 "Failed to iterate over object of type "
                                 +obj->getClass()->getName());
        }
    }

    
    /** recurse over each component invoking action method for each
        non null component. Skips over 'atomic' fields or arrays of atomic
        fields. The {@link Action action} class drives what is invoked for
        each component. The return flag indicates whether any components of
        this object were altered */
    bool internalRecurse(const IObjectSP& obj){
        static const string routine("ObjectIteration::internalRecurse");
        CClassConstSP c = obj->getClass();
        depth++; // record the level of recursion
        bool hasChanged = false;
        if (c->isArray()){     // arrays need special handling
            // look at the type of the components in the array
            // can skip whole array if the components never need to be actioned
            if (!fieldsCache || classContainsComponent(c)){
                IArraySP array = IArraySP::dynamicCast(obj);
                hasChanged = internalRecurseArray(array);
            }
        } else {
            hasChanged = internalRecurseObject(obj);
        }
        depth--;
        return hasChanged;
    }
    
    
    /** utility method - is the object atomic. The definition of which
        means do you stop recursing through the components of an object
        once you hit an atomic one. This includes array of atomic objects */
    bool isPrimitive(IObjectSP& o){
        return o->getClass()->isPrimitive();
    }

    /** Creates an ObjectIteration object which controls the iteration over
        an object's components. The Action class controls what is done to
        each componenet found */
        Imp(State*         state, // NB not fully defined at this stage
             CClassConstSP  targetClass):
        state(state), action(0), actionConst(0), targetClass(targetClass), 
        superClassesSet(0), fieldsCache(0), depth(0),
        skipTransient(false), currentField(0), actionTopObj(true), quit(false),
        quitNow(false),
        recursingDown(true), invokeTwiceFlag(false), useFieldCache(true),
        performSet(false), objectHasChanged(false) {
    }

    /** recurse over each component invoking action method for each non
        null component. Skips over 'atomic' fields or arrays of atomic
        fields. Fails if obj or action are null.  Return value is same as
        supplied obj unless a 'set' operation is performed. In this case
        the return object will be different if the 'set' operation is
        applied to the top level object ie obj. This method must only be called
        on a new instance of Imp. */
    IObjectSP recurse(const IObjectSP& obj){
        static const string routine("ObjectIteration::recurse");
        try{
            if (useFieldCache && !fieldsCache){
                initialiseClassCache(); // initialises superClassesSet field
            }
            if (!obj.get()){
                throw ModelException(routine, "Object is null");
            }
            IObjectSP shiftedObj(obj); /* methods below take non const refs to
                                          shiftedObj */
            if (actionTopObj){
                invokeMethod(shiftedObj, 0);
            }else{
                internalRecurse(shiftedObj);
            }
            return shiftedObj;
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** Initialises the superClassesSet field of the class. This set
        contains all the superclasses of any types which implement the
        target class */
    void initialiseClassCache(){
        // retrieve existing set/create new set if not there - Note using
        // operator[] is v expensive if the element already exists
        SuperClassHashTable::iterator iter1 = superClasses.find(targetClass);
        if (iter1 == superClasses.end()){
            superClassesSet = &superClasses[targetClass];
        } else {
            superClassesSet = &iter1->second;
        }
        /* pull out entry from cache creating if it doesn't exist */
        TweaksHashTable::iterator iter2 = tweaks.find(targetClass);
        if (iter2 == tweaks.end()){
            fieldsCache = &tweaks[targetClass];
        } else {
            fieldsCache = &iter2->second;
        }
        // same for componentsHash
        ComponentsHashTable::iterator iter3 = components.find(targetClass);
        if (iter3 == components.end()){
            componentsHash = &components[targetClass];
        } else {
            componentsHash = &iter3->second;
        }

        // if the set is empty it must have just been created
        if (superClassesSet->empty()){
            /* get list of all classes. Search through for all that are
               instances of targetclass. Add each instance, together with all
               their superclasses to the set */

            // get a list of all our known types
            const CClassVec& classes = CClass::allClasses();

            // look through all types to see what implements 
            for (unsigned int i = 0; i < classes.size(); i++){
                // skip the class if it's in our set already
                ClassHashSet::const_iterator it =superClassesSet->find(classes[i]);

                if (it == superClassesSet->end() && 
                    targetClass->isAssignableFrom(classes[i])){
                    // not in set and implements/derived from targetClass
                    // then loop over this class & superclasses inserting them
                    CClassConstSP clazz = classes[i];
                    do{
                        superClassesSet->insert(clazz);
                    } while ((clazz = clazz->getSuperClass()));
                }
            }
        }
    }


/*
  have superClassesSet - contains all the superclasses of any types
  which implement the target class. Keep this. 

  Also want componentSet - contains all classes that have a field
  which is/implements the target class or contains a field which
  satisfies this condition. Key question: does this includes classes
  whose derived instance implement this eg does set include YieldCurve
  (for targetClass Rho)? Better performance if it does

  How to build this set? Given a type have to decide whether or not to
  put it into the set. Clearly if in superClassesSet then can put in
  directly. Otherwise, need to look at its fields and ask if they are
  in this set. Then if none of them are then (if answer to above
  question is yes) look over all derived instances. 

  One issue though is that we need to know if an element is not in the
  set because we can't fill all the elements in the set in one go
  without asking the same question repeatedly. So it needs to be a
  hash table rather than a set. The stored value reflects whether the
  field would belong in the set or not.

  This point leads on to the tricky question of recursion - we ask to
  see the value in this hash table when we are in the midst of trying
  to calculate the value. eg targetClass = delta. Start by asking
  whether asset is in the set or not, having no fields we move onto
  derived types. We start, say, by examining an XCB. It has fields,
  one of which is an array of assets - and hence the
  recursion. Clearly the answer at this stage is no but we need to
  know that this is only indicative at this stage - so do we need to
  be able to return a value meaning no_so_far?
  
  Have three routines. One which calculates the value and the other
  two which look up the value in the hash table. The latter of these
  forces a calculation if the value is no_so_far whilst the other one
  doesn't. Both routines call the one to calculate the value if it's
  not in the hash (before calling the routine they initialise the hash
  with a value of no_so_far).

*/

/** Returns true if the supplied clazz (or any derived classes) is or
    contains (directly or indirectly) an instance of target class */
    bool classContainsComponent(CClassConstSP clazz){
        // create space for storing unknown states - this is faster than 
        // iterating through the hash table looking for them
        vector<ClassIntHash::iterator> unknownStates;
        // call classContainsComponentState if NO_SO_FAR set
        // value in hash to NO
        int state = getClassCCState(clazz, true /* recalc any NO_SO_FAR */,
                                    unknownStates);
        if (state == COMPONENT_NO_SO_FAR){
            while (clazz->isArray()){
                clazz = clazz->getComponentType();
            }
            ClassIntHash::iterator iter = componentsHash->find(clazz);
            if (iter == componentsHash->end()){
                throw ModelException("ObjectIteration::classContainsComponent",
                                     "internal error");
            }
            // having finished the search for this type can now safely 
            // switch from NO to NO_SO_FAR
            iter->second = COMPONENT_NO;
            state = COMPONENT_NO;
        }

        // now we must switch any 'no so far' to 'unknown' as the 'no so far'
        // is only valid with respect to this particular clazz. Using a vector
        // is quicker than iterating through the hash
        for (unsigned int i = 0; i < unknownStates.size(); i++){
            if (unknownStates[i]->second == COMPONENT_NO_SO_FAR){
                // if COMPONENT_NO is the result for this clazz then it
                // must be true also for any components within
                unknownStates[i]->second = state == COMPONENT_NO? 
                    COMPONENT_NO: COMPONENT_UNKNOWN;
            }
        }

        return (state == COMPONENT_YES? true: false);
    }

    /** Returns the current state of the calculation whose value is
        returned by classContainsComponentState (CC = containsComponent).
        If forceCalc is false any value already in the hash is used
        otherwise a state of NO_SO_FAR forces a recalculation of the value */
    int getClassCCState(CClassConstSP                   clazz, 
                        bool                            forceCalc,
                        vector<ClassIntHash::iterator>& unknownStates){
        // if an array, recurse until no longer an array
        while (clazz->isArray()){
            clazz = clazz->getComponentType();
        }
        // if primitive the answer is always no
        if (clazz->isPrimitive()){
            return COMPONENT_NO;
        }
        // look up in hash table
        ClassIntHash::value_type newEntry(clazz, COMPONENT_NO_SO_FAR);
        pair<ClassIntHash::iterator, bool> entry = 
            componentsHash->insert(newEntry); // insert new entry if not there
        bool entryExists = !entry.second;
        // an entry of UNKNOWN is equivalent to no entry - so must switch
        // to NO_SO_FAR to avoid infinite recursion
        if (entry.first->second == COMPONENT_UNKNOWN){
            entry.first->second = COMPONENT_NO_SO_FAR;
            entryExists = false; // ie treat it as if it wasn't there
        }
        int state = entry.first->second; // for ease
        // if there, return the value (unless NO_SO_FAR and require a recalc) 
        if (entryExists &&
            (!forceCalc || state == COMPONENT_NO || state == COMPONENT_YES)){
            // entry already exists
            return state;
        }
        if (entryExists && forceCalc){
            // forceCalc flag is redundant - entry should not exist if forceCalc
            // set (unless it's 'unknown')
            throw ModelException("ObjectIteration::getClassCCState",
                                 "Internal error");
        }
        // if not there entry has been created and set to NO_SO_FAR (this 
        // addresses the problem of recursion)
        // call method to work out value and store it
        state = evaluateClassCCState(clazz, unknownStates);
        entry.first->second = state;
        if (state == COMPONENT_NO_SO_FAR){
            // record those for which we don't know the answer
            unknownStates.push_back(entry.first);
        }
        return state;
    }

    /** Determines if the supplied clazz (or any derived classes) is or
        contains (directly or indirectly) an instance of target class. It
        does this by examining the registered fields of that class */
    int evaluateClassCCState(CClassConstSP                   clazz,
                             vector<ClassIntHash::iterator>& unknownStates){
        // if clazz is in superClassesSet return YES
        if (superClassesSet->find(clazz) != superClassesSet->end()){
            return COMPONENT_YES;
        }
        // get list of derived instances of this class
        const CClassVec& derived = clazz->getAssignableClasses();
        int state = COMPONENT_NO; // initial value
        // loop over this list of classes
        for (unsigned int i = 0; i < derived.size(); i++){
            // for each class, first check to see if it implements the
            // target class
            CClassConstSP c = derived[i];
            if (superClassesSet->find(c) != superClassesSet->end()){
                return COMPONENT_YES;
            }
            // and then get fields for this class and loop over then
            do {
                const CFieldArray&  allFields = c->getDeclaredFields();
                for (unsigned int j = 0; j < allFields.size(); j++){
                    // evaluate this class (but avoid infinite recursion)
                    int fieldState = getClassCCState(allFields[j]->getType(), 
                                                     false,
                                                     unknownStates);
                    if (fieldState == COMPONENT_YES){
                        // can exit immediately
                        return COMPONENT_YES;
                    }
                    if (fieldState == COMPONENT_NO_SO_FAR){
                        // must record weaker condition
                        state = fieldState;
                    }
                }
                // for special case where we are examining clazz, need to
                // consider fields from parents
            } while (derived[i] == clazz && (c = c->getSuperClass()));
        }
        return state;
    }
 
    /** Retrieve (or create if it doesn't exists) the array of fields
     for supplied clazz, ie determine the set of fields in the
     supplied type (including parents) that need to be iterated over
     for the given target class. Here iterated over means that the
     iteration needs to recurse over the contents of the field eg if
     an asset has a field which is an equity then the list of fields
     of the asset that need to be iterated over for Delta is [at
     least] the field containing the equity */
    const CFieldArray& getFieldsCache(
        const CClassConstSP&   clazz){       // eg XCB
        
        if (clazz->isArray()){ // more of an assert really
            throw ModelException("ObjectIteration::getFieldsCache", 
                                 "Method not valid for array class");
        }
        CFieldArrayHashTable::value_type newEntry(clazz, CFieldArray(0));
        // insert if and only if there is nothing there already
        pair<CFieldArrayHashTable::iterator, bool> entry = 
            fieldsCache->insert(newEntry); // insert empty array
        CFieldArray&    fields = entry.first->second; /* our new array or
                                                         existing entry */
        // if the entry already exists or is primitive do nothing 
        if (entry.second && !clazz->isPrimitive()){
            /* then loop fields - see if they have a list > 0 in length */
            const CFieldArray&  allFields = clazz->getDeclaredFields();
            for (unsigned int i = 0; i < allFields.size(); i++){
                CClassConstSP fieldType = allFields[i]->getType();
                if (classContainsComponent(fieldType)){
                    fields.push_back(allFields[i]);
                }
            }
            
            /* finally have to include parents fields as well -
               but this function handles parents so we only have
               to call it once for the immediate parent */
            CClassConstSP c = clazz->getSuperClass();
            if (c != 0){
                // then add fields from parents
                const CFieldArray& allFields = getFieldsCache(c);
                fields.insert(fields.end(), 
                              allFields.begin(), allFields.end());
            }
        }
        return fields;
    }
    
};

/** Implements generic iteration over the fields of an object and
    then iteration their fields and so on. A strategy pattern is
    used to control what happens for each internal object iterated
    over */

// definition of static fields in class
SuperClassHashTable  ObjectIteration::Imp::superClasses;
ComponentsHashTable  ObjectIteration::Imp::components;
TweaksHashTable      ObjectIteration::Imp::tweaks;
const int ObjectIteration::Imp::COMPONENT_YES = 0;
const int ObjectIteration::Imp::COMPONENT_NO = 1;
const int ObjectIteration::Imp::COMPONENT_NO_SO_FAR = 2;
const int ObjectIteration::Imp::COMPONENT_UNKNOWN = 3; // equivalent to no entry

ObjectIteration::State::~State(){}
#if defined(_MSC_VER)
// disable warning about using 'this' in base member initialisation 
#pragma warning(disable : 4355)
#endif
ObjectIteration::State::State(CClassConstSP  targetClass):
    my(new Imp(this, targetClass)){}
#if defined(_MSC_VER)
#pragma warning(default : 4355) // undo above
#endif

/** Gives non-const access to the object passed in the 
    IAction::invoke method. */
IObjectSP ObjectIteration::State::getObject(){
    my->objectHasChanged = true;
    return my->currentObject;
}

/** Whilst recursing, a call to this function will replace the
    current object with the supplied one. The supplied object must
    of the correct type. No clone of the supplied object is made
    before the Field::set() method is used to replace the existing
    method. */
void ObjectIteration::State::setObject(IObjectSP object){
    my->performSet = true;
    my->objToSet = object;
    // set currentObject field? Scrap objToSet?
    my->objectHasChanged = true;
}

/** Signals that the entire recursion should be terminated (ie can
    be called after recurse() is invoked to signal that recurse()
    should exit). If quitImmediately then the recursion is
    terminated immediately otherwise the upward stage of the
    recursion will continue and the invoke method will be called
    if requested (but no more downward recursion will take place) */
void ObjectIteration::State::quitRecursion(bool quitImmediately) const{
    my->quit = true; // const property is lost in auto_ptr ...
    my->quitNow = quitImmediately;
}

/** Returns the {@link Field} where
    the current object being 'Actioned' resides. If called outside
    an action it returns the Field of the last object
    'Actioned'. If the top level object passed to {@link #recurse
    } is being actioned then null is returned */
CFieldConstSP ObjectIteration::State::getCurrentField() const{
    return my->currentField;
}

/** See invokeTwice(). Returns true if we iteration is on initial
    downwards pass. False if iteration is on second upwards pass */
bool ObjectIteration::State::recursingDownward() const{
    return my->recursingDown;
}

/** Returns the current level of recursion. In particular, returns
    0 when the top object supplied to recurse() is being
    actioned. Returns 1 when a field within the top object is
    being actioned. Returns 2 when a field within a field of the
    top object is being actioned and so on */
int ObjectIteration::State::getDepth() const{
    return my->depth;
}

/** Internal implementation of ObjectIteration class */
class ObjectIteration::InternalImp{
public:
    // what type to look for - use IObject to catch all types
    CClassConstSP         targetClass;
    bool useFieldCache;/* default is true - use cache of which fields to 
                           iterate over */
    bool skipTransient;// default = false;
    // Do we invoke action on top level object supplied
    bool actionTopObj;// default = true
    /* true: call invoke() on downward and upward iteration. Default false */
    bool invokeTwiceFlag;
    InternalImp(CClassConstSP targetClass): 
        targetClass(targetClass),
        useFieldCache(true), skipTransient(false), 
        actionTopObj(true), invokeTwiceFlag(false){}
};

/** Creates an ObjectIteration object which controls the iteration over
    an object's components. See the recurse methods */
ObjectIteration::ObjectIteration(CClassConstSP  targetClass):
    my(new InternalImp(targetClass)){}

/** Switches off the use of the cache of known 'sensitive'
    fields. When on this avoids descending down classes which do not
    contain any fields which are derived from the targetClass. Switching
    off would be useful when for example when the targetClass is IObject.*/
void ObjectIteration::noFieldCache(){
    my->useFieldCache = false;
}

/** Set this to true to avoid tweaking transient fields. 
    Default is false ie dont skip transient fields */
void ObjectIteration::setSkipTransient(bool skipTransient){
    my->skipTransient = skipTransient;
}

/** Set this to true in order to invoke the action on the top
    level object passed to {@link #recurse recurse} should that
    object be of the desired type. Default is false */
void ObjectIteration::actionTopObject(bool actionTopObject){
    my->actionTopObj = actionTopObject;
}

/** Iteration over the components of an object initially descends
    downwards recursing over the fields of each
    component. Eventually the iteration terminates as atomic
    objects (eg doubles, strings) are reached. Calling this method
    causes objects to be iterated over again but this time in the
    reverse order of the initial pass. The invoke method is again
    called for objects matching the target type. A client (within
    its invoke method) can determine whether the iteration is in
    the initial downwards iteration or the second upwards
    iteration by calling recursingDownward() */
void ObjectIteration::invokeTwice(){
    my->invokeTwiceFlag = true;
}

//// copy over setting from this to State object
void ObjectIteration::populateState(State& state) const{
    state.my->useFieldCache = my->useFieldCache;
    state.my->skipTransient = my->skipTransient;
    state.my->actionTopObj = my->actionTopObj;
    state.my->invokeTwiceFlag = my->invokeTwiceFlag;
}

/** recurse over each component invoking action method for each non
    null component. Skips over 'atomic' fields or arrays of atomic
    fields. Fails if obj is null.  Return value is same as
    supplied obj unless a 'set' operation is performed. In this case
    the return object will be different if the 'set' operation is
    applied to the top level object ie obj. Note that 
    (i) it is permissible to have recursive calls to recurse
    (ii) the various methods like invokeTwice(), setSkipTransient() etc
    only apply to future calls of recurse. */
IObjectSP ObjectIteration::recurse(IAction& action, IObjectSP obj){
    State state(my->targetClass);
    state.my->action = &action;
    populateState(state);
    return state.my->recurse(obj);
}

/** Same as above recurse method except that the supplied object (and
    implicitly any of its components) cannot be altered. */
void ObjectIteration::recurse(IActionConst& action, IObjectConstSP obj){
    State state(my->targetClass);
    state.my->actionConst = &action;
    populateState(state);
    state.my->recurse(IObjectSP::constCast(obj)); /* we respect the object's
                                                     constness internally */
}

/** Forces immediate caching of tweaking information for
    products. This is a one off event so nicer if not in main
    pricing loop */
void ObjectIteration::cacheClassInformation(
    CClassConstSP    baseClass,
    const CClassVec& targetClasses){
    // it doesn't particulary matter if this method does not get called
    // It just avoids confusing the quanttify analysis with one off events
    const CClassVec& allClasses = baseClass->getAssignableClasses();
    Imp iteration(0, baseClass); // create dummy one
    for (unsigned int j = 0; j < targetClasses.size(); j++){
        // loop over 'sensitivities'
        iteration.targetClass = targetClasses[j];
        iteration.initialiseClassCache();
        for (unsigned int i = 0; i < allClasses.size(); i++){
            // loop over 'instruments'
            iteration.getFieldsCache(allClasses[i]);
        }
    }
}

/** Returns true if object (or any of its components) derives from or implements
    the supplied target class. The skipTransientFields has the same meaning as
    the setSkipTransient method. */
bool ObjectIteration::dependsUpon(IObjectConstSP   object,
                                  CClassConstSP    targetClass,
                                  bool             skipTransientFields){
    //// a very simple implementation of the action class
    class Action: public virtual IActionConst{
    public: // class is private to this method anyway
        bool             dependsUpon;
    public:
        virtual ~Action(){}
        //// called when we find the object we are looking for
        virtual bool invoke(const State& state, IObjectConstSP obj){
            dependsUpon = true;
            state.quitRecursion(true /* quit immediately */);
            return false; // return parameter irrelevant having said quit now
        }
        Action(): dependsUpon(false) {}
    };
    // just create our Action
    Action action;
    // then create iteration
    ObjectIteration iteration(targetClass);
    iteration.setSkipTransient(skipTransientFields);
    iteration.recurse(action, object);
    return action.dependsUpon;
}

/** Print fields needed to iterate over when looking for
        targetClass in object of type clazz */
void ObjectIteration::printInfo(CClassConstSP    targetClass,
                                CClassConstSP    clazz){
    Imp iteration(0, targetClass); // create dummy one
    iteration.initialiseClassCache();
    const CFieldArray& fields = iteration.getFieldsCache(clazz);
    cout << "For objects of type " << clazz->getName() << 
        " when looking for " 
         << targetClass->getName() << ": " << endl;
    for (unsigned int i = 0; i < fields.size(); i++){
        cout << fields[i]->getName() << ": " << 
            fields[i]->getType()->getName() << endl;
    }
}

ObjectIteration::~ObjectIteration(){}

DRLIB_END_NAMESPACE
