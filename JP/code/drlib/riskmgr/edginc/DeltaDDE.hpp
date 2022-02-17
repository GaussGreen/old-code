//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaDDE.hpp
//
//   Description : Delta sensitivity for DDE
//
//   Author      : Qing Hou
//
//   Date        : 16 Mar 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDG_DELTA_DDE_H
#define EDG_DELTA_DDE_H

#include "edginc/Delta.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for delta - a scalar shift where the derivative is
    calculated via a two sided tweak operation. Specifically for AssetDDE */
class RISKMGR_DLL DeltaDDE: public ScalarShift,
             virtual public Additive,
             virtual public IAssetSpotGreek,
             virtual public ITwoSidedDeriv {
public:
    friend class DeltaDDEHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static string SECOND_ORDER_NAME;
    const static double DEFAULT_SHIFT;
    const static double MINIMUM_SHIFT;

    string outputName() const { return NAME; }
    string output2ndOrderName() const { return SECOND_ORDER_NAME; }

    /** What an object must implement to be tweakable for DeltaDDE */
    class RISKMGR_DLL Shift{
    public:
        friend class DeltaDDEHelper;
        static CClassConstSP const TYPE;
        virtual ~Shift();
        /** Returns the name of the stock/asset - used to determine
            whether to tweak the object */
        virtual string sensName(DeltaDDE* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(DeltaDDE* shift) = 0;
    };
    typedef Shift IShift;

    /** What an object must implement to be able to perform a restorable
        tweak for DeltaDDE. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL RestorableShift: public virtual Shift{
    public:
        friend class DeltaDDEHelper;
        static CClassConstSP const TYPE;
        virtual ~RestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(DeltaDDE* shift) = 0;
    };
    typedef RestorableShift IRestorableShift;

    /** constructor with explicit shift size */
    DeltaDDE(double shiftSize);

    // own special calculate
    virtual void calculate(TweakGroup*      tweakGroup,
                           CResults*        results);

        /** Scales delta and gamma numbers in results object */
    virtual void scaleResult(Results*     results,     // (M)
                             double       scaleFactor) const;

    /** Adds delta and gamma numbers in results object */
    virtual void addResult(Results*           results,     // (M)
                           const Results*     resultsToAdd,
                           double             scaleFactor) const;

    /** Once used to make a shift, this reports the appropriate divisor
        for this sensitivity */
    double divisor() const;

    bool nameMatches(const OutputName&         name,
                     IObjectConstSP            obj);

    void appendName(OutputNameArray&          namesList,
                    IObjectConstSP            obj);

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    CClassConstSP shiftInterface() const;

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents which is also
        restorable */
    CClassConstSP restorableShiftInterface() const;


    // methods to shift given objects

    /**
     * @param obj The object to shift. The object must implement the
     VegaParallel.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    

    /**
     *
     * @param param1 The object to shift. The
     object must implement the DeltaDDE.RestorableShift interface
    */
    virtual void restore(IObjectSP obj);

    /** Returns the shift which has been made for the current pricing
        call */
    virtual ScalarShiftArray getComponentShifts() const;
    
private:
    /** for reflection */
    DeltaDDE();
};

typedef smartConstPtr<DeltaDDE> DeltaDDEConstSP;
typedef smartPtr<DeltaDDE> DeltaDDESP;


DRLIB_END_NAMESPACE

#endif
