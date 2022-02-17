//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VegaSkewParallel.hpp
//
//   Description : Controls calculation of Vega Skew Parallel
//
//   Author      : Mark A Robson
//
//   Date        : 7 March 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_VEGA_SKEW_PARALLEL_HPP
#define EDG_VEGA_SKEW_PARALLEL_HPP
#include "edginc/ScalarShift.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for vega skew parallel */
class RISKMGR_DLL VegaSkewParallel: public ScalarShift,
                        public Additive {
public:
    friend class VegaSkewParallelHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;


    /** What an object must implement to be tweakable for VEGA_SKEW_PARALLEL */
    class RISKMGR_DLL IShift{
    public:
        friend class VegaSkewParallelHelper;
        static CClassConstSP const TYPE;

        /** skew is calculated via a tweak of 
            -shift * log(strike/spot)/log(SKEW_NORMALISATION) */
        static const double SKEW_NORMALISATION;
        /** holds value of log(SKEW_NORMALISATION) */
        static const double LOG_SKEW_NORMALISATION;

        IShift();
        virtual ~IShift();
        /** Returns the name of the vol - used to determine whether to tweak
            the object */
        virtual string sensName(VegaSkewParallel* shift) const = 0;

        /** Name matching mechanism with default single name matching on sensName */
        virtual bool sensNameMatches(VegaSkewParallel* shift, const OutputName& name) const;

        /** Name adding mechanism with default single name adding of sensName */
        virtual void sensAppendName(VegaSkewParallel* shift, OutputNameArray& namesList) const;

        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(VegaSkewParallel* shift) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        tweak for VEGA_PARALLEL. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        friend class VegaSkewParallelHelper;
        static CClassConstSP const TYPE;
        IRestorableShift();
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(VegaSkewParallel* shift) = 0;
    };

    /** constructor with explicit shift size */
    VegaSkewParallel(double     shiftSize);

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
     VegaSkewParallel.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /**
     * @param param1 The object to shift. The
     object must implement the VegaSkewParallel.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

    /** Returns spot value of associated asset. Fails if spot has not been
        set */
    double getSpot() const;

    /** Stores the spot value of the associated asset. Allows vol to 
        implement skew tweak */
    void setSpot(double spotValue);

    /** return skew shift for given strike i.e.
        +shift * log(strike/spot)/log(SKEW_NORMALISATION) */
    double skewShift(double strike) const;

    /** static version of skewShift so can reuse in other skew greeks */
    static double skewShift(double shift, double strike, double spot);

private:
    /** for reflection */
    VegaSkewParallel();
    VegaSkewParallel(const VegaSkewParallel &rhs);
    VegaSkewParallel& operator=(const VegaSkewParallel& rhs);
    //// fields ////
    double spot;  // note: populated by asset during tweak
    bool   spotSet; // indicates if spot value has been set
};

typedef smartConstPtr<VegaSkewParallel> VegaSkewParallelConstSP;
typedef smartPtr<VegaSkewParallel> VegaSkewParallelSP;

DRLIB_END_NAMESPACE

#endif
