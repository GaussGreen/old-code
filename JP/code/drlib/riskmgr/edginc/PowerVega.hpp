//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PowerVega.hpp
//
//   Description : vol shift scenario - shift every point (across strikes) 
//                 of each maturity by:
//                 Shift * T^-Power where T is the time (in years) between now 
//                 and the relevant maturity, all of them at the same time.
//
//   Author      : Andrew J Swain
//
//   Date        : 28 August 2001
//
//
//----------------------------------------------------------------------------

#ifndef POWERVEGA__HPP
#define POWERVEGA__HPP
#include "edginc/ScalarPerturbation.hpp"

DRLIB_BEGIN_NAMESPACE

/** vol shift scenario - shift every point (across strikes) 
    of each maturity by:
    Shift * T^-Power where T is the time (in years) between now 
    and the relevant maturity, all of them at the same time. */
class RISKMGR_DLL PowerVega: public ScalarPerturbation {
public:
    static CClassConstSP const TYPE;

    /** What an object must implement to support PowerVega */
    class RISKMGR_DLL Shift{
    public:
        friend class PowerVegaHelper;
        static CClassConstSP const TYPE;
        Shift();
        virtual ~Shift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(PowerVega* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(PowerVega* shift) = 0;
    };
    typedef Shift IShift;

    /** constructor with explicit shift */
    PowerVega(double power, double shift);

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
     PowerVega.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /** return the power shift (shift * T^-p) */
    double powerShift(double years) const;

private:
    friend class PowerVegaHelper;
    /** for reflection */
    PowerVega();
    PowerVega(const PowerVega &rhs);
    PowerVega& operator=(const PowerVega& rhs);
    double power;
};


typedef smartConstPtr<PowerVega> PowerVegaConstSP;
typedef smartPtr<PowerVega> PowerVegaSP;

DRLIB_END_NAMESPACE

#endif
