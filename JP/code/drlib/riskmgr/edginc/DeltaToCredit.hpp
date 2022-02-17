//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaToCredit.hpp
//
//   Description : DeltaToCredit sensitivity
//
//   Author      : André Segger
//
//   Date        : 15 Sep 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDG_DELTA_TO_CREDIT_HPP
#define EDG_DELTA_TO_CREDIT_HPP
#include "edginc/ScalarShift.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"
#include "edginc/TwoSidedDeriv.hpp"
#include "edginc/AssetSpotGreek.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for DeltaToCredit - a scalar shift where the derivative is
    calculated via a two sided tweak operation */
class RISKMGR_DLL DeltaToCredit: public ScalarShift,
                     virtual public Additive,
                     virtual public IAssetSpotGreek,
                     virtual public ITwoSidedDeriv {
public:
    friend class DeltaToCreditHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static string SECOND_ORDER_NAME;
    const static double DEFAULT_SHIFT;
    const static double MINIMUM_SHIFT;

    /** What an object must implement to be tweakable for DELTA_TO_CREDIT */
    class RISKMGR_DLL Shift{
    public:
        friend class DeltaToCreditHelper;
        static CClassConstSP const TYPE;
        virtual ~Shift();
        /** Returns the name of the stock/asset - used to determine
            whether to tweak the object */
        virtual string sensName(DeltaToCredit* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(DeltaToCredit* shift) = 0;
    };
    typedef Shift IShift;

    /** What an object must implement to be able to perform a restorable
        tweak for DELTA. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL RestorableShift: public virtual Shift{
    public:
        friend class DeltaToCreditHelper;
        static CClassConstSP const TYPE;
        virtual ~RestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(DeltaToCredit* shift) = 0;
    };
    typedef RestorableShift IRestorableShift;

    /** constructor with explicit shift size */
    DeltaToCredit(double     shiftSize);

    /** Constructor with explicit name and shift size, useful for things
        like DeltaNextDay, e.g. */
    DeltaToCredit(double     shiftSize,
                  string     newName);

    /** constructor with explicit shift size, model and control pointers */
    DeltaToCredit(double     shiftSize,
                  IModel*    model,
                  Control*   control);

    /** calculates deltaToCredit via 2 sided tweak */
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
     VegaParallel.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /** @param param1 The object to shift. The
        object must implement the DeltaToCredit.RestorableShift interface */
    virtual void restore(IObjectSP obj);

    /** Returns the shift which has been made for the current pricing
        call */
    virtual ScalarShiftArray getComponentShifts() const;

protected:
    /** Questionable whether we should derive from this class - CrossGamma
        is currently */
    DeltaToCredit(const CClassConstSP& clazz,
                  const string&        outputName,
                  const double&        shiftSize);
    
private:
    /** for reflection */
    DeltaToCredit();
    DeltaToCredit(const DeltaToCredit& rhs);
    DeltaToCredit& operator=(const DeltaToCredit& rhs);
};

typedef DeltaToCredit CDeltaToCredit;
typedef smartConstPtr<DeltaToCredit> DeltaToCreditConstSP;
typedef smartPtr<DeltaToCredit> DeltaToCreditSP;
typedef DeltaToCreditConstSP CDeltaToCreditConstSP;
typedef DeltaToCreditSP CDeltaToCreditSP;

DRLIB_END_NAMESPACE

#endif
