//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : CreditSpreadLevel.hpp
//
//   Description : credit spread level scenario - set cs to supplied value
//
//   Author      : Tycho von Rosenvinge and Andre Segger
//
//   Date        : 23 April 2002
//
//
//----------------------------------------------------------------------------

#ifndef CREDIT_SPREAD_LEVEL_HPP
#define CREDIT_SPREAD_LEVEL_HPP
#include "edginc/ScalarPerturbation.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for cs level scenario - set cs to supplied value */
class RISKMGR_DLL CreditSpreadLevel: public ScalarPerturbation {
public:
    static CClassConstSP const TYPE;

    /** What an object must implement to support CreditSpreadLevel */
    class RISKMGR_DLL Shift{
    public:
        friend class CreditSpreadLevelHelper;
        static CClassConstSP const TYPE;
        virtual ~Shift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(CreditSpreadLevel* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(CreditSpreadLevel* shift) = 0;
    };
    typedef Shift IShift;

    /** constructor with explicit cs level */
    CreditSpreadLevel(double spread);

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
     CreditSpreadLevel.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
   
private:
    double shiftSize; // poor name - legacy reasons $unregistered
    friend class CreditSpreadLevelHelper;
    /** for reflection */
    CreditSpreadLevel();
    CreditSpreadLevel(const CreditSpreadLevel &rhs);
    CreditSpreadLevel& operator=(const CreditSpreadLevel& rhs);
};

typedef smartConstPtr<CreditSpreadLevel> CreditSpreadLevelConstSP;
typedef smartPtr<CreditSpreadLevel> CreditSpreadLevelSP;

DRLIB_END_NAMESPACE

#endif
