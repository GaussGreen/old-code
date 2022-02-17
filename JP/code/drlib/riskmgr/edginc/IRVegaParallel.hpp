//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRVegaParallel.hpp
//
//   Description : Controls calculation of parallel IR vega
//
//   Author      : Andrew J Swain
//
//   Date        : 22 February 2002
//
//
//----------------------------------------------------------------------------

#ifndef _IRVEGAPARALLEL_HPP
#define _IRVEGAPARALLEL_HPP
#include "edginc/ScalarShift.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for ir vega parallel */
class RISKMGR_DLL IRVegaParallel: public ScalarShift,
                      public Additive {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** What an object must implement to be tweakable for IRVEGA_PARALLEL */
    class RISKMGR_DLL Shift{
    public:
        friend class IRVegaParallelHelper;
        static CClassConstSP const TYPE;
        virtual ~Shift();
        /** Returns the name of the vol - used to determine whether to tweak
            the object */
        virtual string sensName(IRVegaParallel* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(IRVegaParallel* shift) = 0;
    };
    typedef Shift IShift;

    /** What an object must implement to be able to perform a restorable
        tweak for IRVEGA_PARALLEL. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL RestorableShift: public virtual Shift{
    public:
        friend class IRVegaParallelHelper;
        static CClassConstSP const TYPE;
        virtual ~RestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(IRVegaParallel* shift) = 0;
    };
    typedef RestorableShift IRestorableShift;

    /** constructor with explicit shift size */
    IRVegaParallel(double shiftSize);

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
        IRVegaParallel.Shift interface */
    virtual bool nameMatches(const OutputName&         name,
                             IObjectConstSP          obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP          obj);

    /**
     * @param obj The object to shift. The object must implement the
     IRVegaParallel.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /**
     * @param param1 The object to shift. The
     object must implement the IRVegaParallel.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

private:
    friend class IRVegaParallelHelper;
    /** for reflection */
    IRVegaParallel();
    IRVegaParallel(const IRVegaParallel &rhs);
    IRVegaParallel& operator=(const IRVegaParallel& rhs);
};

typedef smartConstPtr<IRVegaParallel> IRVegaParallelConstSP;
typedef smartPtr<IRVegaParallel> IRVegaParallelSP;

DRLIB_END_NAMESPACE

#endif
