//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ShortTermSqueezeTweak.hpp
//
//   Description : SHORT_TERM_SQUEEZE_TWEAK sensitivity
//
//   Author      : Oliver Brockhaus
//
//   Date        : 05 Mar 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDGSHORT_TERM_SQUEEZE_TWEAK_HPP
#define EDGSHORT_TERM_SQUEEZE_TWEAK_HPP
#include "edginc/ScalarShift.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Results.hpp"
#include "edginc/NotApplicable.hpp"
#include "edginc/Untweakable.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for SHORT_TERM_SQUEEZE_TWEAK - a scalar shift where the derivative is
    calculated via a one sided tweak operation 
    refers to correlation skew parameter used in XCBs */
class RISKMGR_DLL ShortTermSqueezeTweak:  public          ScalarShift,
                     virtual public  Additive 
{
public:
    friend class ShortTermSqueezeTweakHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** What an object must implement to be tweakable for SHORT_TERM_SQUEEZE_TWEAK */
    class RISKMGR_DLL Shift{
    public:
        friend class ShortTermSqueezeTweakHelper;
        static CClassConstSP const TYPE;
        virtual ~Shift();

        /** Returns the name of the stock/asset - used to determine
            whether to tweak the object */
        virtual string sensName(ShortTermSqueezeTweak* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(ShortTermSqueezeTweak* shift) = 0;
    };
    typedef Shift IShift;

    /** What an object must implement to be able to perform a restorable
        tweak for SHORT_TERM_SQUEEZE_TWEAK. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL RestorableShift: public virtual Shift{
    public:
        friend class ShortTermSqueezeTweakHelper;
        static CClassConstSP const TYPE;
        virtual ~RestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(ShortTermSqueezeTweak* shift) = 0;
    };
    typedef RestorableShift IRestorableShift;

    /** constructor with explicit shift size */
    ShortTermSqueezeTweak(double shiftSize);
    ShortTermSqueezeTweak(double shiftSize, bool useSign);

    /** Once used to make a shift, this reports the appropriate divisor
        for this sensitivity */
    double divisor() const;

    /** Validation */
    void validatePop2Object();

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    CClassConstSP shiftInterface() const;

    /** returns the interface identifying what an object has to do in order
        to support the tweak that this object represents which is also
        restorable */
    CClassConstSP restorableShiftInterface() const;


    // methods to shift given objects

    /** Returns true if the supplied object matches the supplied name
        for this sensitivity. The object must implement this sensitivity's
        Shift interface */
    virtual bool nameMatches(const OutputName&         name,
                             IObjectConstSP            obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list. The object must
        implement this sensitivity's Shift interface */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP            obj);

    /**
     * @param obj The object to shift. The object must implement the
     ShortTermSqueezeTweak.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    

    /**
     *
     * @param obj The object to shift. The
     object must implement the ShortTermSqueezeTweak.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

    /** Returns true if implementations of Shift should override the sign of
    the shift size when making the tweak */
    bool overrideShiftSizeSign() const; 

private:
    bool useShiftSizeSign;
    /** for reflection */
    ShortTermSqueezeTweak();
    ShortTermSqueezeTweak(const ShortTermSqueezeTweak &rhs);
    ShortTermSqueezeTweak& operator=(const ShortTermSqueezeTweak& rhs);
};

typedef smartConstPtr<ShortTermSqueezeTweak> ShortTermSqueezeTweakConstSP;
typedef smartPtr<ShortTermSqueezeTweak> ShortTermSqueezeTweakSP;

DRLIB_END_NAMESPACE

#endif
