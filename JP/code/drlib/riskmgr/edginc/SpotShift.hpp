//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SpotShift.hpp
//
//   Description : spot shift scenario - shift spot by supplied % value
//
//   Author      : Mark A Robson
//
//   Date        : 16 Jul 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SPOTSHIFT__HPP
#define EDR_SPOTSHIFT__HPP
#include "edginc/ScalarPerturbation.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for spot shift scenario - shift spot by supplied % value */
class RISKMGR_DLL SpotShift: public ScalarPerturbation {
public:
    static CClassConstSP const TYPE;

    /** What an object must implement to support SpotShift */
    class RISKMGR_DLL Shift{
    public:
        friend class SpotShiftHelper;
        static CClassConstSP const TYPE;
        virtual ~Shift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(SpotShift* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(SpotShift* shift) = 0;
    };
    typedef Shift IShift;

    /** constructor with explicit % spot shift amount */
    SpotShift(double spot);

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
     SpotShift.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
private:
    friend class SpotShiftHelper;
    /** for reflection */
    SpotShift();
    SpotShift(const SpotShift &rhs);
    SpotShift& operator=(const SpotShift& rhs);
};


typedef smartConstPtr<SpotShift> SpotShiftConstSP;
typedef smartPtr<SpotShift> SpotShiftSP;

DRLIB_END_NAMESPACE

#endif
