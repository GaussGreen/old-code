//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolParallelShift.hpp
//
//   Description : vol shift scenario - add parallel shift to vol
//
//   Author      : Andrew J Swain
//
//   Date        : 5 June 2001
//
//
//----------------------------------------------------------------------------

#ifndef VOLPARALLELSHIFT__HPP
#define VOLPARALLELSHIFT__HPP
#include "edginc/ScalarPerturbation.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for vol shift scenario - add parallel shift to vol */
class RISKMGR_DLL VolParallelShift: public ScalarPerturbation {
public:
    static CClassConstSP const TYPE;
    const static double MIN_SPOT_VOL;
    const static double MIN_FWD_VOL;

    /** What an object must implement to support VolParallelShift */
    class RISKMGR_DLL Shift{
    public:
        friend class VolParallelShiftHelper;
        static CClassConstSP const TYPE;
        Shift();
        virtual ~Shift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(VolParallelShift* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(VolParallelShift* shift) = 0;
    };
    typedef Shift IShift;

    /** constructor with explicit shift */
    VolParallelShift(double shift);

    //// overridden to shift fund proxies if flag is set
    virtual bool findAndShift(IObjectSP         objectToShift, 
                              OutputNameConstSP name);
 
    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    CClassConstSP shiftInterface() const;

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
     VolParallelShift.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
private:
    friend class VolParallelShiftHelper;
    /** for reflection */
    VolParallelShift();
    VolParallelShift(const VolParallelShift &rhs);
    VolParallelShift& operator=(const VolParallelShift& rhs);

    bool shiftFundProxies;
};


typedef smartConstPtr<VolParallelShift> VolParallelShiftConstSP;
typedef smartPtr<VolParallelShift> VolParallelShiftSP;

DRLIB_END_NAMESPACE

#endif
