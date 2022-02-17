// 	mrobson $

#ifndef EDG_OBJECTITERATION_H
#define EDG_OBJECTITERATION_H

#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

/** Implements generic recursive iteration over the fields of an object and
    then iteration over their fields and so on. A strategy pattern is
    used to control what happens for each internal object iterated
    over ie users supply a class implementing IAction. When an appropriate
    object is found the relevant method on IAction is invoked. One of the
    paraemeters to this method is 'State' which can be used to get additional
    information about the recursion. */
class TOOLKIT_DLL ObjectIteration{
    class Imp;
public:
    /** Captures the state of the ObjectIteration as presented to the client.
        An instance of this class is given to the client when an object of the
        required type is found */
    class TOOLKIT_DLL State{
    public:
        /** Gives non-const access to the object passed in the 
            IAction::invoke method. Note that this causes the containing
            object to be notified via the IObject::fieldsUpdated() method */
        IObjectSP getObject();

        /** Whilst recursing, a call to this function will replace the
            current object with the supplied one. The supplied object must
            of the correct type. No clone of the supplied object is made
            before the Field::set() method is used to replace the existing
            method. */
        void setObject(IObjectSP object);

        /** Signals that the entire recursion should be terminated (ie can
            be called after recurse() is invoked to signal that recurse()
            should exit). If quitImmediately then the recursion is
            terminated immediately otherwise the upward stage of the
            recursion will continue and the invoke method will be called
            if requested (but no more downward recursion will take place) */
        void quitRecursion(bool quitImmediately) const;

        /** Returns the {@link Field} where
            the current object being 'Actioned' resides. If called outside
            an action it returns the Field of the last object
            'Actioned'. If the top level object passed to {@link #recurse
            } is being actioned then null is returned */
        CFieldConstSP getCurrentField() const;

        /** See invokeTwice(). Returns true if we iteration is on initial
            downwards pass. False if iteration is on second upwards pass */
        bool recursingDownward() const;

        /** Returns the current level of recursion. In particular, returns
            0 when the top object supplied to recurse() is being
            actioned. Returns 1 when a field within the top object is
            being actioned. Returns 2 when a field within a field of the
            top object is being actioned and so on */
        int getDepth() const;
    private:
        State(CClassConstSP  targetClass);
        ~State();
        State(const State&); // make private
        State& operator=(const State& rhs); // make private
        State* operator&(); // make private
        friend class ObjectIteration;
        auto_ptr<Imp> my; // hide implementation
    };

    /** This class defines the interface that must be implelemented to use the
        ObjectIteration methods. They essentially form a set of call backs
        which define the action needed to be performed for each object */
    class TOOLKIT_DLL IAction{
    public:
        virtual ~IAction();
        /** invoke is called for each object which meets the criteria
            supplied to the ObjectIteration class. Returning false
            from this routine on the way down stops recursion on this
            object. Note to get non-const access you must use the getObject()
            method on State (this will inform dependent objects that the
            object has been altered). */
        virtual bool invoke(State& state, IObjectConstSP obj) = 0;
    };

    /** Same as IAction but does not allow any component object to be
        modified ie can be used with [top level] const objects. */
    class TOOLKIT_DLL IActionConst{
    public:
        virtual ~IActionConst();
        /** invoke is called for each object which meets the criteria
            supplied to the ObjectIteration class. Returning false
            from this routine on the way down stops recursion on this
            object. Note to get non-const access you must use an IAction
            object instead. Do not cast away constness */
        virtual bool invoke(const State& state, IObjectConstSP obj) = 0;
    };
  
     /** This class defines an interface that classes can optionally choose
         to implement. It allows greater (including runtime control) of what
         fields within are class are 'recursed' on. This overrides any
         default behaviour of transient fields */
    class TOOLKIT_DLL IOverride{
    public:
        static CClassConstSP const TYPE;

         /** When the ObjectIteration is recursing down through the
            different fields of a class, then it will invoke this
            method supplying the field which is about to be "recursed"
            on, together with the class of what is being looked
            for. The return value indicates whether or not the field
            should be "recursed" over. For example, a class with a field
            called "myAsset" which is of type CAsset. When doing the
            delta tweak, this method will be invoked with the field
            corresponding to the "myAsset" field and targetClass will
            be Delta::Shift or Delta::RestorableShift. Note that
            CFieldConstSP and CClassConstSP are in fact constant
            pointers so direct comparison can be used for
            performance. For an example see FXVol.cpp */
        virtual bool recurse(const CFieldConstSP& field,
                             const CClassConstSP& targetClass) const = 0;
        virtual ~IOverride();
    };

    /** Creates an ObjectIteration object which controls the iteration over
        an object's components. See the recurse methods */
    ObjectIteration(CClassConstSP  targetClass);

    /** Switches off the use of the cache of known 'sensitive'
        fields. When on this avoids descending down classes which do not
        contain any fields which are derived from the targetClass. Switching
        off would be useful when for example when the targetClass is IObject.
        Only has an affect for future calls to recurse. */
    void noFieldCache();

    /** Set this to true to avoid tweaking transient fields.
        Default is false ie tweak transient fields.
        Only has an affect for future calls to recurse. */
    void setSkipTransient(bool skipTransient);
    
    /** Set this to false in order to avoid invoking the action on the top
        level object passed to {@link #recurse} (should that
        object be of the desired type). Default is true. 
        Only has an affect for future calls to recurse. */
    void actionTopObject(bool actionTopObject);

    /** Iteration over the components of an object initially descends
        downwards recursing over the fields of each
        component. Eventually the iteration terminates as atomic
        objects (eg doubles, strings) are reached. Calling this method
        causes objects to be iterated over again but this time in the
        reverse order of the initial pass. The invoke method is again
        called for objects matching the target type. A client (within
        its invoke method) can determine whether the iteration is in
        the initial downwards iteration or the second upwards
        iteration by calling recursingDownward().
        Only has an affect for future calls to recurse. */
    void invokeTwice();

    /** recurse over each component invoking action method for each non
        null component. Skips over 'atomic' fields or arrays of atomic
        fields. Fails if obj is null.  Return value is same as
        supplied obj unless a 'set' operation is performed. In this case
        the return object will be different if the 'set' operation is
        applied to the top level object ie obj. Note that 
        (i) it is permissible to have recursive calls to recurse
        (ii) the various methods like invokeTwice(), setSkipTransient() etc
        only apply to future calls of recurse. */
    IObjectSP recurse(IAction& action, IObjectSP obj);

    /** Same as above recurse method except that the supplied object (and
        implicitly any of its components) cannot be altered. */
    void recurse(IActionConst& action, IObjectConstSP obj);

    /** Forces immediate caching of tweaking information for
        products. This is a one off event so nicer if not in main
        pricing loop */
    static void cacheClassInformation(CClassConstSP    baseClass,
                                      const CClassVec& targetClasses);

    /** Print fields needed to iterate over when looking for targetClass in
        object of type clazz */
    static void printInfo(CClassConstSP    targetClass,
                          CClassConstSP    clazz);

    /** Returns true if object (or any of its components) derives from
        or implements the supplied target class. The
        skipTransientFields has the same meaning as the
        setSkipTransient method. */
    static bool dependsUpon(IObjectConstSP   object,
                            CClassConstSP    targetClass,
                            bool             skipTransientFields);

    // needed as class Imp isn't defined here
    ~ObjectIteration();

private:
    // don't use
    ObjectIteration(const ObjectIteration& rhs);
    ObjectIteration& operator=(const ObjectIteration& rhs);

    void populateState(State& state) const;
    class InternalImp;
    auto_ptr<InternalImp> my; // hide implementation
};

DRLIB_END_NAMESPACE

#endif
