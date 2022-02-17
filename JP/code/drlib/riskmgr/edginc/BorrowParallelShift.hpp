//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BorrowParallelShift.hpp
//
//   Description : borrow curve shift scenario - adds a parallel shift
//
//   Author      : Andrew J Swain
//
//   Date        : 10 April 2002
//
//
//----------------------------------------------------------------------------

#ifndef BORROWPARALLELSHIFT__HPP
#define BORROWPARALLELSHIFT__HPP
#include "edginc/ScalarPerturbation.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for borrow curve shift scenario - adds a parallel shift */
class RISKMGR_DLL BorrowParallelShift: public ScalarPerturbation {
public:
    static CClassConstSP const TYPE;

    /** What an object must implement to support BorrowParallelShift */
    class RISKMGR_DLL Shift{
    public:
        friend class BorrowParallelShiftHelper;
        static CClassConstSP const TYPE;
        virtual ~Shift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(BorrowParallelShift* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(BorrowParallelShift* shift) = 0;
    };
    typedef Shift IShift;

    /** constructor with explicit shift */
    BorrowParallelShift(double shift);

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
     BorrowParallelShift.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
private:
    friend class BorrowParallelShiftHelper;
    /** for reflection */
    BorrowParallelShift();
    BorrowParallelShift(const BorrowParallelShift &rhs);
    BorrowParallelShift& operator=(const BorrowParallelShift& rhs);
};

typedef smartConstPtr<BorrowParallelShift> BorrowParallelShiftConstSP;
typedef smartPtr<BorrowParallelShift> BorrowParallelShiftSP;

DRLIB_END_NAMESPACE

#endif
