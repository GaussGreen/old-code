//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolRelativeShift.hpp
//
//   Description : Scenario shift where a relative shift (e.g. a 10% shift
//                 to 20% vol gives 22% vol) is applied for each benchmark.
//                 Shifts are defined for ranges of benchmarks (e.g. <= 1Y,
//                 1Y -> 5Y, 5Y ->  30Y etc).
//
//   Author      : Andrew J Swain
//
//   Date        : 24 April 2003
//
//
//----------------------------------------------------------------------------

#ifndef _VOLRELATIVESHIFT__HPP
#define _VOLRELATIVESHIFT__HPP
#include "edginc/MultiExpiryStepShift.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolProcessedBS.hpp"

DRLIB_BEGIN_NAMESPACE
class TweakGroup;
class Results;

/** Scenario shift where a relative shift (e.g. a 10% shift
to 20% vol gives 22% vol) is applied for each benchmark.
Shifts are defined for ranges of benchmarks (e.g. <= 1Y,
1Y -> 5Y, 5Y ->  30Y etc).
*/

class MARKET_DLL VolRelativeShift: public MultiExpiryStepShift {
public:    
    static CClassConstSP const TYPE;

    // what's the shift for a given date ? */
    virtual double shiftSize(const DateTime& today,     // to anchor expiries
                             const DateTime& shiftDate) const;

    /** What an object must implement to support VolRelativeShift */
    class MARKET_DLL IShift{
    public:
        friend class VolRelativeShiftHelper;
        static CClassConstSP const TYPE;
        IShift();
        virtual ~IShift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(VolRelativeShift* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(VolRelativeShift* shift) = 0;
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
 

    virtual void validatePop2Object();

    /** Override clone method to copy non registered fields over */
    virtual IObject* clone() const;

    /** Stores the spot value of the associated asset */
    void setSpot(const DateTime& today, const Asset* asset);

private:
    friend class VolRelativeShiftHelper;
    // for reflection
    VolRelativeShift();

    // fields
    double moneyness;
    bool   gotVol; // indicates if moneyness vol has been set
    CVolProcessedBSSP vol;  // moneyness vol $unregistered
};


typedef smartConstPtr<VolRelativeShift> VolRelativeShiftConstSP;
typedef smartPtr<VolRelativeShift> VolRelativeShiftSP;

DRLIB_END_NAMESPACE

#endif
