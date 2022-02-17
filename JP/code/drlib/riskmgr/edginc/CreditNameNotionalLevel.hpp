//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : NotionalLevel.hpp
//
//   Description : Change in name notional scenario 
//
//   Author      : Linus Thand 
//
//   Date        : 20 June 2006
//
//----------------------------------------------------------------------------

#ifndef CREDIT_NAME_NOTIONAL_LEVEL__HPP
#define CREDIT_NAME_NOTIONAL_LEVEL__HPP
#include "edginc/ScalarPerturbation.hpp"

DRLIB_BEGIN_NAMESPACE

/** A scenario which replaces a notional in a CDO portfolio. This is 
 * used by managed trades where clients can request a change in notionals in 
 * a porfolio, and where another variable can be solved for for PV neutrality 
 * (e.g a parallel shift in strikes). **/
 
class CreditNameNotionalLevel : public ScalarPerturbation  {
 public:
    static CClassConstSP const TYPE;

    /** What an object must implement to support NotionalLevel */
    class Shift{
     public:
        static CClassConstSP const TYPE;
        virtual string sensName(CreditNameNotionalLevel* shift) const = 0;
        virtual bool   sensShift(CreditNameNotionalLevel* shift) = 0;
    };

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    CClassConstSP shiftInterface() const;

    /** Returns true if the supplied object matches the supplied name
        for this sensitivity. The object must implement this sensitivity's
        Shift interface */
    virtual bool nameMatches(const OutputName& name,
                             IObjectConstSP    obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list. The object must
        implement this sensitivity's Shift interface */
    virtual void appendName(OutputNameArray& namesList,
                            IObjectConstSP   obj);

    /**
     * @param obj The object to shift. The object must implement the
     NotionalLevel.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);

 private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    CreditNameNotionalLevel();
    ~CreditNameNotionalLevel();
    CreditNameNotionalLevel(const CreditNameNotionalLevel& rhs); 
    CreditNameNotionalLevel& operator=(const CreditNameNotionalLevel& rhs);
};

FORWARD_DECLARE(CreditNameNotionalLevel);

DRLIB_END_NAMESPACE

#endif
