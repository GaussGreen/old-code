//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CtsDivOverride.hpp
//
//   Description : [Scenario] Override dividends with a continuous
//                 dividend (for converts, surprise, surprise)
//                 default case is a single cts div d
//                 can also specify cts div d which grows by g% every year until some
//                 end date when it stops growing
//
//   Author      : Mark A Robson
//
//   Date        : 8 August 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDR_CTSDIVOVERRIDE__HPP
#define EDR_CTSDIVOVERRIDE__HPP

#include "edginc/ScalarPerturbation.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/DividendList.hpp"


DRLIB_BEGIN_NAMESPACE

/** Perturbation for cts div override scenario - set cts div to supplied value*/
class MARKET_DLL CtsDivOverride: public ScalarPerturbation {
public:
    static CClassConstSP const TYPE;

    /** What an object must implement to support CtsDivOverride */
    class MARKET_DLL IShift{
    public:
        friend class CtsDivOverrideHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(CtsDivOverride* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(CtsDivOverride* shift) = 0;
    };

    /** builds a list of dividends given the cts dividend required */
    DividendListSP dividends(const DateTime& today) const;

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
     CtsDivOverride.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
private:
    friend class CtsDivOverrideHelper;
    /** for reflection */
    CtsDivOverride();
    CtsDivOverride(const CtsDivOverride &rhs);
    CtsDivOverride& operator=(const CtsDivOverride& rhs);

    virtual void validatePop2Object();

    double     growthRate;
    ExpirySP   growthEnd;
    bool       isGrowing;

};


typedef smartConstPtr<CtsDivOverride> CtsDivOverrideConstSP;
typedef smartPtr<CtsDivOverride> CtsDivOverrideSP;

DRLIB_END_NAMESPACE

#endif
