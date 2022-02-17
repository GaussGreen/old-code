//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : ParSpreadLevelAtGivenDate.hpp
//
//   Description : CDS Par Spread level scenario
//                 1. Get spread S corresponding to given date (doing linear interpolation)
//                 2. Flatten spread curve to this value S
//
//   Author      : Antoine Gregoire
//
//   Date        : June 2005
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_PAR_SPREAD_LEVEL_AT_GIVEN_DATE_H
#define QLIB_PAR_SPREAD_LEVEL_AT_GIVEN_DATE_H

#include "edginc/Perturbation.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL ParSpreadLevelAtGivenDate : public Perturbation {
public:

    static CClassConstSP const TYPE;

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    virtual CClassConstSP shiftInterface() const;

    /** Returns true if the supplied object matches the supplied name
        for this sensitivity.  The object must implement the
        VegaParallel.Shift interface */
    virtual bool nameMatches(const OutputName&         name,
                             IObjectConstSP            obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP            obj);

    /** Shifts the object (which supports being tweaked
        by this type of sens control) using given shift. The return value
        indicates whether or not components of this object need to be
        tweaked ie true: infrastructure should continue to recurse through
        components tweaking them; false: the infrastructure shouldn't
        touch any components within this object */
    virtual bool shift(IObjectSP obj);

    /** What an object must implement to be tweakable for ParSpreadLevelAtGivenDate */
    class RISKMGR_DLL IShift{
    public:
        static CClassConstSP const TYPE;
        virtual ~IShift();
        virtual bool sensShift(ParSpreadLevelAtGivenDate* shift) = 0;
        virtual string sensName(ParSpreadLevelAtGivenDate* shift) const = 0;
    };

    /** Returns the level date */
    DateTime getLevelDate() const;

private:
    // Reflection
    ParSpreadLevelAtGivenDate();
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);                               
                           
    // FIELDS                           
    DateTime levelDate;
};

DRLIB_END_NAMESPACE

#endif //QLIB_PAR_SPREAD_LEVEL_AT_GIVEN_DATE_H
