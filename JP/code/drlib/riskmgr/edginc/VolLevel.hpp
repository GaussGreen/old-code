//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolLevel.hpp
//
//   Description : vol level scenario - set vol to supplied value
//
//   Author      : Andrew J Swain
//
//   Date        : 8 May 2001
//
//
//----------------------------------------------------------------------------

#ifndef VOLLEVEL__HPP
#define VOLLEVEL__HPP
#include "edginc/ScalarPerturbation.hpp"

DRLIB_BEGIN_NAMESPACE

/** Perturbation for vol level scenario - set vol to supplied value */
class RISKMGR_DLL VolLevel: public ScalarPerturbation{
public:
    static CClassConstSP const TYPE;

    /** What an object must implement to support VolLevel */
    class RISKMGR_DLL Shift{
    public:
        friend class VolLevelHelper;
        static CClassConstSP const TYPE;
        Shift();
        virtual ~Shift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(VolLevel* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(VolLevel* shift) = 0;
    };
    typedef Shift IShift;

    /** constructor with explicit vol level */
    VolLevel(double vol);

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
     VolLevel.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
private:
    friend class VolLevelHelper;
    /** for reflection */
    VolLevel();
    VolLevel(const VolLevel &rhs);
    VolLevel& operator=(const VolLevel& rhs);
};


typedef smartConstPtr<VolLevel> VolLevelConstSP;
typedef smartPtr<VolLevel> VolLevelSP;

DRLIB_END_NAMESPACE

#endif
