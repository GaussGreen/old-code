//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolBenchmarkShift.hpp
//
//   Description : vol shift scenario - add benchmark shift to vol
//
//   Author      : Andrew J Swain
//
//   Date        : 5 June 2001
//
//
//----------------------------------------------------------------------------

#ifndef VOLBENCHMARKSHIFT__HPP
#define VOLBENCHMARKSHIFT__HPP
#include "edginc/ScalarPerturbation.hpp"
#include "edginc/Expiry.hpp"

DRLIB_BEGIN_NAMESPACE

class VegaPointwise;

/** Sens Control for vol shift scenario - add benchmark shift to vol */
class RISKMGR_DLL VolBenchmarkShift: public ScalarPerturbation {
public:
    static CClassConstSP const TYPE;

    /** What an object must implement to support VolBenchmarkShift */
    class RISKMGR_DLL Shift{
    public:
        friend class VolBenchmarkShiftHelper;
        static CClassConstSP const TYPE;
        Shift();
        virtual ~Shift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(VolBenchmarkShift* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(VolBenchmarkShift* shift) = 0;
    };
    typedef Shift IShift;

    /** constructor with explicit shift */
    VolBenchmarkShift(ExpirySP expiry, double shift);

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
     VolBenchmarkShift.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /** does this match the given expiry ? */
    bool expiryEquals(const Expiry* expiry) const;
   
    /** is this the last shift for given market data ? */
    bool lastShift() const;

    /** next two methods useful for param vol scenarios */
    /** store all expiries on vol to be shifted */
    void cacheExpiries(const ExpiryArrayConstSP& expiries);

    ExpiryConstSP getExpiry() const;

    ExpiryArrayConstSP getAllExpiries() const;

private:
    friend class VolBenchmarkShiftHelper;
    /** for reflection */
    VolBenchmarkShift();
    VolBenchmarkShift(const VolBenchmarkShift &rhs);
    VolBenchmarkShift& operator=(const VolBenchmarkShift& rhs);

    ExpiryConstSP expiry;
    bool          isLastShift;
    ExpiryArrayConstSP allExpiries; // transient
};


typedef smartConstPtr<VolBenchmarkShift> VolBenchmarkShiftConstSP;
typedef smartPtr<VolBenchmarkShift> VolBenchmarkShiftSP;

DRLIB_END_NAMESPACE

#endif
