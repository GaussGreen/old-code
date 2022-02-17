
#ifndef EDG_SENS_MGR_H
#define EDG_SENS_MGR_H
#include "edginc/OutputName.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE
class ITweakNameListID;
class ITweakNameResolver;
class ITweakID;
class ITweakOptID;
class ITweakQualifierID;
class SensControlPerName;
class Results;

/** The methods of SensMgr which do not change the supplied object */
class RISKMGR_DLL SensMgrConst{
public:
    virtual ~SensMgrConst();

    /** takes a reference to the given object to allow tweaking etc */
    SensMgrConst(const IObject*  topObject);

    /** takes a reference to the given object to allow tweaking etc */
    SensMgrConst(IObjectConstSP topObject);

    /** Returns true if object (or any of its components) derives from
        or implements the supplied target class. */
    static bool dependsUpon(IObjectConstSP   object,
                            CClassConstSP    targetClass);
 
   /** Returns array of output names which can be tweaked for given
        SensControl - in particular, if the SensControl has its own names
        they are ignored */
    OutputNameArrayConstSP allNames(ITweakNameListID*   sens);

    /** Returns ITweakIDQualifier specific object needed for qualifing the
        greeks for a given SensControl eg benchmark dates for pointwise
        vega */
    IObjectConstSP qualifier(ITweakQualifierID* sens);

    /** Returns the first object of the given @a type found inside the
        topObject which has the specified @name; throws
        NoSubjectsFoundException if none is found */
    IObjectConstSP theFirst(CClassConstSP type,
                            const ITweakNameResolver* name) const;

    /** Returns all objects of the given @a type found inside the
        topObject which have the specified @name */
    ObjectArrayConstSP all(CClassConstSP type,
                           const ITweakNameResolver* name) const;

protected:
    // fields
    IObjectSP     topObject;
private:
    class CollectName;
    class GetQualifier;
    class CollectObject;
    SensMgrConst(const SensMgrConst &rhs);
    SensMgrConst& operator=(const SensMgrConst& rhs);
};

/** Implements generic 'external' tweaking of instrument and/or market data.
    The 'lite' version offers no support for restoring objects. This can
    be useful when a number of shifts are to be applied (eg scenarios) and it
    is cheaper to copy the original before shifting */
class RISKMGR_DLL SensMgr: public SensMgrConst{
public:
    /** takes a reference to the given object to allow tweaking etc */
    SensMgr(IObject*  topObject);

    /** takes a reference to the given object to allow tweaking etc */
    SensMgr(IObjectSP topObject);

    IObjectSP theFirst_mutable(CClassConstSP type,
                               const ITweakNameResolver* name) const;

    /** Shifts data inside object as indicated by sens. Shifted object
        is modified version of original. The pointer returned will always
        be the same as that used in the constructor */
    virtual IObjectSP shift(ITweakID* sens);

    /** same as above shift method but name for overrideName parameter is
        used rather than a call to ITweakNameResolver. If overrideName
        is null then that is treated the same nameResolver returning null */
    virtual IObjectSP shift(ITweakID* sens, OutputNameConstSP overrideName);
                                
    /** Indicates whether the object supplied was actually changed by the
        shift request */
    bool getShiftStatus() const;

    virtual ~SensMgr();

protected:
    class Shift;
    // fields
    bool          shifted;
private:
    class Find;
    SensMgr(const SensMgr &rhs);
    SensMgr& operator=(const SensMgr& rhs);
};

/** Same functionality as SensMgr but offers ability to restore an object
    after shifting it */
class RISKMGR_DLL SensMgrOpt: public SensMgr{
public:
    /** takes a reference to the given object to allow tweaking etc. The
     optimalShift flag drives how the algorithm for restoring the original
     object. If false, a clone of the original object is taken, if true,
     use is made of any restore methods and cloning is kept to a minimum */
    SensMgrOpt(IObject* topObject, bool optimalShift);

    /** Same as above but takes smart pointer */
    SensMgrOpt(IObjectSP topObject, bool optimalShift);

    /** Equivalent to SensMgrOpt(topObject, true) */
    SensMgrOpt(IObject*  topObject);

    /** Equivalent to SensMgrOpt(topObject, true) */
    SensMgrOpt(IObjectSP  topObject);

    /** Shifts data inside object as indicated by sens. Shifted object may
        be copy of original or may be modified version of original. Use
        restore to return object to its original state. If optimalShift
        was set in constructor then sens parameter will be dynamically cast
        to ITweakOptID which will be used if it succeeds. Otherwise it will be
        like sens is ITweakOptID but the restorableShift method always returns
        false. */
    virtual IObjectSP shift(ITweakID*    sens);

    /** Same as above but avoids the internal cast */
    IObjectSP shift(ITweakOptID* sens);
                                
    /** same as above shift method for ITweakID but name for
        overrideName parameter is used rather than a call to
        ITweakNameResolver. If overrideName is null then that is
        treated the same nameResolver returning null */
    virtual IObjectSP shift(ITweakID* sens, OutputNameConstSP overrideName);

    /** Same as above but avoids the internal cast */
    IObjectSP shift(ITweakOptID* sens, OutputNameConstSP overrideName);

    /** restores an object after a tweak */
    void restore();

    ~SensMgrOpt();
private:
    class OptimalRestoreShift;
    class OptimalShift;
    IObjectSP optimalShift(ITweakID* sens, ITweakOptID* optSens, 
                           OutputNameConstSP overrideName);
    IObjectSP copyAndShift(ITweakID*    sens, OutputNameConstSP overrideName);
    void optimalShiftRestore();
    
    // fields
    bool          doOptimalShift;
    OptimalShift* optShift; // used when using OptimalShift class for shifting

    // not implemented
    SensMgrOpt(const SensMgrOpt &rhs);
    SensMgrOpt& operator=(const SensMgrOpt& rhs);
};

FORWARD_DECLARE(SensMgrOpt)

DRLIB_END_NAMESPACE
#endif
