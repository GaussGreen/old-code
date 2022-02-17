//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AssetAssetVegaPointwise.cpp
//
//   Description : Controls calculation of AssetVega pointwise
//
//   Author      : Andre Segger
//
//   Date        : 24 September 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDG_ASSET_VEGA_POINTWISE_HPP
#define EDG_ASSET_VEGA_POINTWISE_HPP
#include "edginc/VectorShift.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for vega pointwise */
class RISKMGR_DLL AssetVegaPointwise: public VectorShift,
                          public virtual Additive {
public:
    friend class AssetVegaPointwiseHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;


    /** What an object must implement to be tweakable for VEGA_POINTWISE */
    class RISKMGR_DLL IShift{
    public:
        friend class AssetVegaPointwiseHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the vol - used to determine whether to tweak
            the object */
        virtual string sensName(AssetVegaPointwise* shift) const = 0;

        /** Return the array of expiries (ie maturities/benchmark dates) that
            need to be tweaked for this vol */
        virtual ExpiryArrayConstSP sensExpiries(AssetVegaPointwise* shift) const =0;

        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(AssetVegaPointwise* shift) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        tweak for VEGA_POINTWISE. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        friend class AssetVegaPointwiseHelper;
        static CClassConstSP const TYPE;
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(AssetVegaPointwise* shift) = 0;
    };

    /** constructor with explicit shift size */
    AssetVegaPointwise(double     shiftSize);

    /** constructor with explicit shift size and name override (override
        allows a VEGA_POINTWISE calculation to be stored under, eg,
        VEGA_MATRIX) */
    AssetVegaPointwise(const string& overrideName,
                  double     shiftSize);

    /** constructor with explicit shift size and explicit expiry to
        tweak. Also need to specify all the expiries. Building
        AssetVegaPointwise this way is useful as it allows individal points to
        be tweaked. Note that the supplied expiries are overwritten by the
        vector shift calculate method */
    AssetVegaPointwise(double                     shiftSize,
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
     AssetVegaPointwise.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /** The supplied object is queried for the expiries array needed
        for doing vega pointwise and this array is returned. The supplied
        object must implement the AssetVegaPointwise.Shift interface */
    virtual IObjectConstSP qualifier(IObjectConstSP obj);

    /**
     * @param param1 The object to shift. The
     object must implement the AssetVegaPointwise.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

    /** implements a one sided vector derivative for each instance of the
        market data which is sensitive to this SensControl */
    virtual void calculate(TweakGroup*  tweakGroup,
                           CResults*    results);


private:
    /** for reflection */
    AssetVegaPointwise();
    AssetVegaPointwise(const AssetVegaPointwise &rhs);
    AssetVegaPointwise& operator=(const AssetVegaPointwise& rhs);
};

typedef smartConstPtr<AssetVegaPointwise> AssetVegaPointwiseConstSP;
typedef smartPtr<AssetVegaPointwise> AssetVegaPointwiseSP;


DRLIB_END_NAMESPACE

#endif
