//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ParSpreadLevel.hpp
//
//   Description : CDS Par Spread level scenario - set spread to supplied value
//
//   Author      : Andrew McCleery
//
//   Date        : 16 March 2004
//
//
//----------------------------------------------------------------------------

#ifndef PAR_SPREAD_LEVEL_HPP
#define PAR_SPREAD_LEVEL_HPP
#include "edginc/ScalarPerturbation.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for scenario - set CDS Par Spread to supplied value */
class RISKMGR_DLL ParSpreadLevel: public ScalarPerturbation {
public:
    static CClassConstSP const TYPE;

    /** What an object must implement to support ParSpreadLevel */
    class RISKMGR_DLL IShift{
    public:
        friend class ParSpreadLevelHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(ParSpreadLevel* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(ParSpreadLevel* shift) = 0;
    };

    /** constructor with explicit cs level */
    ParSpreadLevel(double spread);

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    CClassConstSP shiftInterface() const;

    // methods to shift given objects

    /** Returns true if the supplied object matches the supplied name
        for this sensitivity. The object must implement this sensitivity's
        Shift interface */
    virtual bool nameMatches(const OutputName&         name,
                             IObjectConstSP            obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list. The object must
        implement this sensitivity's Shift interface */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP            obj);

    /**
     * @param obj The object to shift. The object must implement the
     ParSpreadLevel.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
private:
    friend class ParSpreadLevelHelper;
    /** for reflection */
    ParSpreadLevel();
    ParSpreadLevel(const ParSpreadLevel &rhs);
    ParSpreadLevel& operator=(const ParSpreadLevel& rhs);
};


typedef smartConstPtr<ParSpreadLevel> ParSpreadLevelConstSP;
typedef smartPtr<ParSpreadLevel> ParSpreadLevelSP;

DRLIB_END_NAMESPACE

#endif
