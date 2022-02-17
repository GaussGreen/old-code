//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BorrowLevel.hpp
//
//   Description : borrow level scenario - overwrites the borrow curve with a flat curve
//
//   Author      : André Segger
//
//   Date        : 01 July 2002
//
//
//----------------------------------------------------------------------------

#ifndef BORROWLEVEL_HPP
#define BORROWLEVEL_HPP
#include "edginc/ScalarPerturbation.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for borrow curve shift scenario - adds a parallel shift */
class RISKMGR_DLL BorrowLevel: public ScalarPerturbation {
public:
    static CClassConstSP const TYPE;

    /** What an object must implement to support BorrowLevel */
    class RISKMGR_DLL Shift{
    public:
        friend class BorrowLevelHelper;
        static CClassConstSP const TYPE;
        virtual ~Shift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(BorrowLevel* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(BorrowLevel* shift) = 0;
    };
    typedef Shift IShift;

    /** constructor with explicit shift */
    BorrowLevel(double shift);

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
     BorrowLevel.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
private:
    friend class BorrowLevelHelper;
    /** for reflection */
    BorrowLevel();
    BorrowLevel(const BorrowLevel &rhs);
    BorrowLevel& operator=(const BorrowLevel& rhs);
};

typedef smartConstPtr<BorrowLevel> BorrowLevelConstSP;
typedef smartPtr<BorrowLevel> BorrowLevelSP;

DRLIB_END_NAMESPACE

#endif
