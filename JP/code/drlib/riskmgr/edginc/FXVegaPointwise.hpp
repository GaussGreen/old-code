//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FXVegaPointwise.hpp
//
//   Description : Controls calculation of FX Vega pointwise
//
//   Author      : Andrew J Swain
//
//   Date        : 3 February 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDG_FX_VEGA_POINTWISE_HPP
#define EDG_FX_VEGA_POINTWISE_HPP
#include "edginc/VectorShift.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for FX vega pointwise */
class RISKMGR_DLL FXVegaPointwise: public VectorShift,
                       public virtual Additive {
public:
    friend class FXVegaPointwiseHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;


    /** What an object must implement to be tweakable for FX_VEGA_POINTWISE */
    class RISKMGR_DLL IShift{
    public:
        friend class FXVegaPointwiseHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the vol - used to determine whether to tweak
            the object */
        virtual string sensName(FXVegaPointwise* shift) const = 0;

        /** Return the array of expiries (ie maturities/benchmark dates) that
            need to be tweaked for this vol */
        virtual ExpiryArrayConstSP sensExpiries(FXVegaPointwise* shift) const =0;

        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(FXVegaPointwise* shift) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        tweak for FX_VEGA_POINTWISE. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        friend class FXVegaPointwiseHelper;
        static CClassConstSP const TYPE;
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(FXVegaPointwise* shift) = 0;
    };

    /** constructor with explicit shift size */
    FXVegaPointwise(double     shiftSize);

    /** constructor with explicit shift size and name override (override
        allows a VEGA_POINTWISE calculation to be stored under, eg,
        VEGA_MATRIX) */
    FXVegaPointwise(const string& overrideName,
                    double     shiftSize);

    /** constructor with explicit shift size and explicit expiry to
        tweak. Also need to specify all the expiries. Building
        FXVegaPointwise this way is useful as it allows individal points to
        be tweaked. Note that the supplied expiries are overwritten by the
        vector shift calculate method */
    FXVegaPointwise(double                     shiftSize,
                    const ExpiryConstSP&       expiry,
                    const ExpiryArrayConstSP&  allExpiries);

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
                             IObjectConstSP          obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list. The object must
        implement this sensitivity's Shift interface */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP          obj);

     /**
     * @param obj The object to shift. The object must implement the
     FXVegaPointwise.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /** The supplied object is queried for the expiries array needed
        for doing vega pointwise and this array is returned. The supplied
        object must implement the FXVegaPointwise.Shift interface */
    virtual IObjectConstSP qualifier(IObjectConstSP obj);

    /**
     * @param param1 The object to shift. The
     object must implement the FXVegaPointwise.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

private:
    /** for reflection */
    FXVegaPointwise();
    FXVegaPointwise(const FXVegaPointwise &rhs);
    FXVegaPointwise& operator=(const FXVegaPointwise& rhs);
};

typedef smartConstPtr<FXVegaPointwise> FXVegaPointwiseConstSP;
typedef smartPtr<FXVegaPointwise> FXVegaPointwiseSP;


DRLIB_END_NAMESPACE

#endif
