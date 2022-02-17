//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RhoBorrowPointwise.hpp
//
//   Description : Controls calculation of Rho borrow pointwise
//
//   Author      : Stephen Hope
//
//   Date        : 9 March 2001
//
//
//----------------------------------------------------------------------------

#ifndef RHOBORROWPOINTWISE_HPP
#define RHOBORROWPOINTWISE_HPP
#include "edginc/VectorShift.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for rho pointwise */
class RISKMGR_DLL RhoBorrowPointwise: public VectorShift,
                          public Additive {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** What an object must implement to be tweakable for RHO_BORROW_POINTWISE */
    class RISKMGR_DLL IShift{
    public:
        friend class RhoBorrowPointwiseHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the borrow curve - used to determine whether 
            to tweak the object */
        virtual string sensName(RhoBorrowPointwise* shift) const = 0;

        /** Return the array of expiries (ie maturities/benchmark dates) that
            need to be tweaked for this  borrow curve */
        virtual ExpiryArrayConstSP sensExpiries(RhoBorrowPointwise* shift) const =0;

        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(RhoBorrowPointwise* shift) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        tweak for RHO_BORROW_POINTWISE. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        friend class RhoBorrowPointwiseHelper;
        static CClassConstSP const TYPE;
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(RhoBorrowPointwise* shift) = 0;
    };

    /** constructor with explicit shift size */
    RhoBorrowPointwise(double shiftSize);

    /** Once used to make a shift, this reports the appropriate divisor
        for this sensitivity */
    double divisor() const;

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    CClassConstSP shiftInterface() const;

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents which is also
        restorable */
    CClassConstSP restorableShiftInterface() const;


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
     RhoBorrowPointwise.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /** The supplied object is queried for the expiries array needed
        for doing rho pointwise and this array is returned. The supplied
        object must implement the RhoBorrowPointwise.Shift interface */
    virtual IObjectConstSP qualifier(IObjectConstSP obj);

    /**
     * @param obj The object to query
     * @return Whether the given object implements 
     RhoBorrowPointwise.RestorableShift
     */
    virtual bool restorableShift(IObjectConstSP& obj) const;

    /**
     * @param param1 The object to shift. The
     object must implement the RhoBorrowPointwise.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

private:
    friend class RhoBorrowPointwiseHelper;
    /** for reflection */
    RhoBorrowPointwise();
    RhoBorrowPointwise(const RhoBorrowPointwise &rhs);
    RhoBorrowPointwise& operator=(const RhoBorrowPointwise& rhs);
};

typedef smartConstPtr<RhoBorrowPointwise> RhoBorrowPointwiseConstSP;
typedef smartPtr<RhoBorrowPointwise> RhoBorrowPointwiseSP;


DRLIB_END_NAMESPACE

#endif
