//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CEVPowerPointwise.cpp
//
//   Description : for CEVPower pointwise
//
//   Date        : 25 Feb 2002
//
//
//----------------------------------------------------------------------------

#ifndef CEVPOWER_POINTWISE_HPP
#define CEVPOWER_POINTWISE_HPP
#include "edginc/VectorShift.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for CEVPower pointwise */
class RISKMGR_DLL CEVPowerPointwise: public VectorShift,
                     public Additive {
public:
    friend class CEVPowerPointwiseHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;


    /** What an object must implement to be tweakable for CEVPower */
    class RISKMGR_DLL IShift{
    public:
        friend class CEVPowerPointwiseHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the vol - used to determine whether to tweak
            the object */
        virtual string sensName(CEVPowerPointwise* shift) const = 0;

        /** Return the array of expiries (ie maturities/benchmark dates) that
            need to be tweaked for this vol */
        virtual ExpiryArrayConstSP sensExpiries(CEVPowerPointwise* shift) const =0;

        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(CEVPowerPointwise* shift) = 0;
    };

    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        friend class CEVPowerPointwiseHelper;
        static CClassConstSP const TYPE;
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(CEVPowerPointwise* shift) = 0;
    };

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents which is also
        restorable */
    CClassConstSP restorableShiftInterface() const;

    virtual void restore(IObjectSP obj);

    /** constructor with explicit shift size */
    CEVPowerPointwise(double     shiftSize);

    /** constructor with explicit shift size and name override */
    CEVPowerPointwise(const string& overrideName,
                  double     shiftSize);

    /** Once used to make a shift, this reports the appropriate divisor
        for this sensitivity */
    double divisor() const;

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    CClassConstSP shiftInterface() const;

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
     CEVPowerPointwise.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /** The supplied object is queried for the expiries array needed
        for doing CEVPower pointwise and this array is returned. The supplied
        object must implement the CEVPowerPointwise.Shift interface */
    virtual IObjectConstSP qualifier(IObjectConstSP obj);

private:
    /** for reflection */
    CEVPowerPointwise();
    CEVPowerPointwise(const CEVPowerPointwise &rhs);
    CEVPowerPointwise& operator=(const CEVPowerPointwise& rhs);
};

typedef smartConstPtr<CEVPowerPointwise> CEVPowerPointwiseConstSP;
typedef smartPtr<CEVPowerPointwise> CEVPowerPointwiseSP;


DRLIB_END_NAMESPACE

#endif
