
#ifndef EDG_ASSET_VEGA_PARALLEL_H
#define EDG_ASSET_VEGA_PARALLEL_H
#include "edginc/ScalarShift.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"
#include "edginc/TweakGroup.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for vega parallel */
class RISKMGR_DLL AssetVegaParallel: public ScalarShift,
                         public virtual Additive {
public:
    friend class AssetVegaParallelHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;


    /** What an object must implement to be tweakable for VEGA_PARALLEL */
    class RISKMGR_DLL Shift{
    public:
        friend class AssetVegaParallelHelper;
        static CClassConstSP const TYPE;
        virtual ~Shift();
        /** Returns the name of the vol - used to determine whether to tweak
            the object */
        virtual string sensName(AssetVegaParallel* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(AssetVegaParallel* shift) = 0;
    };
    typedef Shift IShift;

    /** What an object must implement to be able to perform a restorable
        tweak for VEGA_PARALLEL. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL RestorableShift: public virtual Shift{
    public:
        friend class AssetVegaParallelHelper;
        static CClassConstSP const TYPE;
        virtual ~RestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(AssetVegaParallel* shift) = 0;
    };
    typedef RestorableShift IRestorableShift;

    /** constructor with explicit shift size */
    AssetVegaParallel(double     shiftSize);

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
        VegaParallel.Shift interface */
    virtual bool nameMatches(const OutputName&         name,
                             IObjectConstSP          obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP          obj);

    /**
     * @param obj The object to shift. The object must implement the
     VegaParallel.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /**
     * @param param1 The object to shift. The
     object must implement the VegaParallel.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

    /** implements a one sided scalar derivative for each instance of the
        market data which is sensitive to this SensControl */
    virtual void calculate(TweakGroup*  tweakGroup, CResults*    results);

    void setEquityConversionFactor(double equityConvFactor);

    double getEquityConversionFactor() const;

    void setEquityName(const string& newEquityName);

    string getEquityName() const; 


private:
    /** for reflection */
    AssetVegaParallel();
    AssetVegaParallel(const AssetVegaParallel &rhs);
    AssetVegaParallel& operator=(const AssetVegaParallel& rhs);

    /** private fields */
    double equityConversionFactor; // $unregistered
    string equityName; // $unregistered
};

typedef AssetVegaParallel CAssetVegaParallel;
typedef smartConstPtr<AssetVegaParallel> AssetVegaParallelConstSP;
typedef smartPtr<AssetVegaParallel> AssetVegaParallelSP;
typedef AssetVegaParallelConstSP CAssetVegaParallelConstSP;
typedef AssetVegaParallelSP CAssetVegaParallelSP;

DRLIB_END_NAMESPACE

#endif
