//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaSurface.hpp
//
//   Description : Surface delta - ATM vol rides surface, smile remains constant
//                 for fund proxies
//
//   Author      : Andrew J Swain
//
//   Date        : 7 March 2003
//
//
//----------------------------------------------------------------------------

#ifndef DELTA_SURFACE_PROXY_HPP
#define DELTA_SURFACE_PROXY_HPP
#include <map>
#include "edginc/ScalarShift.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"
#include "edginc/TwoSidedDeriv.hpp"

DRLIB_BEGIN_NAMESPACE

class DeltaSurface;

/** Sens Control for proxy DeltaSurface - a scalar shift where the derivative is
    calculated via a two sided tweak operation */
class RISKMGR_DLL DeltaSurfaceProxy: public ScalarShift,
                         virtual public Additive,
                         virtual public ITwoSidedDeriv {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static string SECOND_ORDER_NAME;
    const static double DEFAULT_SHIFT;
    const static double MINIMUM_SHIFT;

    /** What an object must implement to be tweakable for DELTA_SURFACE_PROXY */
    class RISKMGR_DLL Shift{
    public:
        friend class DeltaSurfaceHelper;
        static CClassConstSP const TYPE;
        virtual ~Shift();

        /** Returns true if this object matches the supplied name with 
            respect to the sensitivity */
        virtual bool sensNameMatches(DeltaSurfaceProxy* shift,
                                     const OutputName&  name) const = 0;

        /** Appends the name(s) of this object with respect to
            the sensitivity to the supplied list. */
        virtual void sensAppendName(DeltaSurfaceProxy* shift,
                                    OutputNameArray&   namesList) const = 0;
      
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(DeltaSurfaceProxy* shift) = 0;
    };
    typedef Shift IShift;

    /** What an object must implement to be able to perform a restorable
        tweak for DELTA_SURFACE_PROXY. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL RestorableShift: public virtual Shift{
    public:
        friend class DeltaSurfaceHelper;
        static CClassConstSP const TYPE;
        virtual ~RestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(DeltaSurfaceProxy* shift) = 0;
    };
    typedef RestorableShift IRestorableShift;

    /** constructor with explicit shift size */
    DeltaSurfaceProxy(double shiftSize);

    /** calculates DeltaSurface via 2 sided tweak */
    virtual void calculate(TweakGroup*      tweakGroup,
                           CResults*        results);

    /** Scales DeltaSurface and gamma numbers in results object */
    virtual void scaleResult(Results*     results,     // (M)
                             double       scaleFactor) const;

    /** Adds DeltaSurface and gamma numbers in results object */
    virtual void addResult(Results*           results,     // (M)
                           const Results*     resultsToAdd,
                           double             scaleFactor) const;

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
     DeltaSurface.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    

    /**
     *
     * @param param1 The object to shift. The
     object must implement the DeltaSurface.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

    /** Returns the shift which has been made for the current pricing
        call */
    virtual ScalarShiftArray getComponentShifts() const;

    /** make a surface delta shift using own shift */
    DeltaSurface* makeDeltaSurfaceShift() const;

private:
    friend class DeltaSurfaceProxyHelper;

    /** for reflection */
    DeltaSurfaceProxy();
    DeltaSurfaceProxy(const DeltaSurface &rhs);
    DeltaSurfaceProxy& operator=(const DeltaSurface& rhs);
};

typedef smartConstPtr<DeltaSurfaceProxy> DeltaSurfaceProxyConstSP;
typedef smartPtr<DeltaSurfaceProxy> DeltaSurfaceProxySP;

DRLIB_END_NAMESPACE

#endif
