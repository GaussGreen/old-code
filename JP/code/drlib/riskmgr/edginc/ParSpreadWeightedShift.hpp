//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ParSpreadWeightedShift.hpp
//
//   Description : Scenario shift: Par Spread curve benchmarks can be shifted by
//                 different amounts. If the user defined shiftSizes lie between
//                 Par Spread curve benchmarks then we interpolate. 
//                 Extrapolate flat off the ends.
//
//   Author      : Jose Hilera
//
//   Date        : 15-June-2005
//
//
//----------------------------------------------------------------------------

#ifndef _ParSpreadWEIGHTSHIFT__HPP
#define _ParSpreadWEIGHTSHIFT__HPP

#include "edginc/MultiExpiryShift.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

/** Scenario shift: Par Spread curve benchmarks can be shifted by
    different amounts. If the user defined shiftSizes lie between
    Par Spread benchmarks then we interpolate. Extrapolate flat off
    the ends. */
class RISKMGR_DLL ParSpreadWeightedShift: public MultiExpiryShift {
public:    
    static CClassConstSP const TYPE;

    /** What an object must implement to support ParSpreadWeightedShift */
    class RISKMGR_DLL IShift{
    public:
        friend class ParSpreadWeightedShiftHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(ParSpreadWeightedShift* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(ParSpreadWeightedShift* shift) = 0;
    };

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
     PowerVega.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
 
    /**
     * Abstract method to shift the DoubleArray passed as an argument.
     * The implementation in derived classes will invoke the shiftArray 
     * method with 3 parameters (see below) indicating whether or not 
     * the shift is additive (as opposed to multiplicative) */
    virtual void shiftArray(ExpiryArraySP rate_expiries, 
                            DoubleArray* rates,
                            DateTime& today) = 0;

    /**
     * Shifts the DoubleArray passed as 3rd argument in the dates
     * indicated by the 2nd argument. The first argument determines
     * whether the shift is additive or multiplicative.
     * Note this method is not virtual: subclasses are expected to call
     * back to it with the right arguments, rather than overriding it */
    void shiftArray(bool isAdditiveShift, 
                    ExpiryArraySP rate_expiries, 
                    DoubleArray* rates,
                    DateTime& today);

 
    virtual void validatePop2Object();

    // what's the shift for a given date ? */
    virtual double shiftSize(const DateTime& today,
                             const DateTime& shiftDate) const;

protected:
    ParSpreadWeightedShift(CClassConstSP clazz);

private:
    friend class ParSpreadWeightedShiftHelper;
    // for reflection
    ParSpreadWeightedShift();
};

DRLIB_END_NAMESPACE

#endif
