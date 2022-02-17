//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaSurface.hpp
//
//   Description : Surface delta - ATM vol rides surface, smile remains constant
//
//   Author      : Andrew J Swain
//
//   Date        : 24 February 2003
//
//
//----------------------------------------------------------------------------

#ifndef DELTA_SURFACE_HPP
#define DELTA_SURFACE_HPP
#include <map>
#include "edginc/ScalarShift.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"
#include "edginc/TwoSidedDeriv.hpp"

DRLIB_BEGIN_NAMESPACE


/** Sens Control for DeltaSurface - a scalar shift where the derivative is
    calculated via a two sided tweak operation */
class RISKMGR_DLL DeltaSurface: public ScalarShift,
                    virtual public Additive,
                    virtual public ITwoSidedDeriv {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static string SECOND_ORDER_NAME;
    const static double DEFAULT_SHIFT;
    const static double MINIMUM_SHIFT;

    /** What an object must implement to be tweakable for DELTA_SURFACE */
    class RISKMGR_DLL Shift{
    public:
        friend class DeltaSurfaceHelper;
        static CClassConstSP const TYPE;
        Shift();
        virtual ~Shift();
        /** Returns the name of the stock/asset - used to determine
            whether to tweak the object */
        virtual string sensName(DeltaSurface* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(DeltaSurface* shift) = 0;
    };
    typedef Shift IShift;

    /** What an object must implement to be able to perform a restorable
        tweak for DELTA_SURFACE. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL RestorableShift: public virtual Shift{
    public:
        friend class DeltaSurfaceHelper;
        static CClassConstSP const TYPE;
        virtual ~RestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(DeltaSurface* shift) = 0;
    };
    typedef RestorableShift IRestorableShift;

    /** constructor with explicit shift size */
    DeltaSurface(double shiftSize);

    /** constructor with explicit shift size, model and control pointers */
    DeltaSurface(double   shiftSize,
                 IModel*  model,
                 Control* control);

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
                             IObjectConstSP          obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list. The object must
        implement this sensitivity's Shift interface */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP          obj);

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

    /** Returns spot value of associated asset. Fails if spot has not been
        set */
    double getSpot() const;

    /** Stores the spot value of the associated asset. Allows vol to 
        implement skew tweak */
    void setSpot(double spotValue);

    void setSpot(double        spotValue,   // asset spot
                 const string& eqName,      // asset name
                 const string& volName);    // asset's vol name
   
private:
    friend class DeltaSurfaceHelper;
    /** This is a work around for solaris assembler - it falls over if
        symbol names become too long (eg map<string, string>). Work around
        is to inherit of string (NB string is a template class) */
    class RISKMGR_DLL MyString: public string{
    public:
        MyString(const string& s);
        MyString();
    };
    //// Ideally typedef map<string, string> TweakNameMap;
    typedef map<MyString, MyString> TweakNameMap;

    string getEqName(const string& volName);

    /** for reflection */
    DeltaSurface();
    DeltaSurface(const DeltaSurface &rhs);
    DeltaSurface& operator=(const DeltaSurface& rhs);

    //// fields ////
    double       spot;     // note: populated by asset during tweak
    bool         spotSet;  // indicates if spot value has been set

    // lookup table of equity vs vol names
    TweakNameMap volNameMap; // $unregistered
};

typedef smartConstPtr<DeltaSurface> DeltaSurfaceConstSP;
typedef smartPtr<DeltaSurface> DeltaSurfaceSP;

DRLIB_END_NAMESPACE

#endif
