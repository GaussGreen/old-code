//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : JumpRateParallel.cpp
//
//   Description : for JumpRate parallel tweak
//
//   Date        : 25 Feb 2002
//
//----------------------------------------------------------------------------

#ifndef JUMP_RATE_PARALLEL_H
#define JUMP_RATE_PARALLEL_H
#include "edginc/ScalarShift.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for jump rate parallel */
class RISKMGR_DLL JumpRateParallel: public ScalarShift,
                    public Additive {
public:
    friend class JumpRateParallelHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;


    /** What an object must implement to be tweakable for JUMP_RATE_PARALLEL */
    class RISKMGR_DLL Shift{
    public:
        friend class JumpRateParallelHelper;
        static CClassConstSP const TYPE;
        virtual ~Shift();
        /** Returns the name of the vol - used to determine whether to tweak
            the object */
        virtual string sensName(JumpRateParallel* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(JumpRateParallel* shift) = 0;
    };
    typedef Shift IShift;

    /** What an object must implement to be able to perform a restorable
        tweak for JUMP_RATE_PARALLEL. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL RestorableShift: public virtual Shift{
    public:
        friend class JumpRateParallelHelper;
        static CClassConstSP const TYPE;
        virtual ~RestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(JumpRateParallel* shift) = 0;
    };
    typedef RestorableShift IRestorableShift;

    /** constructor with explicit shift size */
    JumpRateParallel(double     shiftSize);

    /** Once used to make a shift, this reports the appropriate divisor
        for this sensitivity */
    double divisor() const;

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    CClassConstSP shiftInterface() const;

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents which is also
        restorable */
    CClassConstSP restorableShiftInterface() const;


    // methods to shift given objects

    /** Returns true if the supplied object matches the supplied name
        for this sensitivity.  The object must implement the
        JumpRateParallel.Shift interface */
    virtual bool nameMatches(const OutputName&         name,
                             IObjectConstSP          obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP          obj);

    /**
     * @param obj The object to shift. The object must implement the
     JumpRateParallel.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /**
     * @param param1 The object to shift. The
     object must implement the JumpRateParallel.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

private:
    /** for reflection */
    JumpRateParallel();
    JumpRateParallel(const JumpRateParallel &rhs);
    JumpRateParallel& operator=(const JumpRateParallel& rhs);
};

typedef JumpRateParallel CJumpRateParallel;
typedef smartConstPtr<JumpRateParallel> JumpRateParallelConstSP;
typedef smartPtr<JumpRateParallel> JumpRateParallelSP;
typedef JumpRateParallelConstSP CJumpRateParallelConstSP;
typedef JumpRateParallelSP CJumpRateParallelSP;


DRLIB_END_NAMESPACE

#endif
