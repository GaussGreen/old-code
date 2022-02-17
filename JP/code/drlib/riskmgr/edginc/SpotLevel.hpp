//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SpotLevel.hpp
//
//   Description : spot level scenario - set spot to supplied value
//
//   Author      : Andrew J Swain
//
//   Date        : 25 April 2001
//
//
//----------------------------------------------------------------------------

#ifndef SPOTLEVEL__HPP
#define SPOTLEVEL__HPP
#include "edginc/ScalarPerturbation.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for spot level scenario - set spot to supplied value */
class RISKMGR_DLL SpotLevel: public ScalarPerturbation {
public:
    static CClassConstSP const TYPE;

    /** What an object must implement to support SpotLevel */
    class RISKMGR_DLL Shift{
    public:
        friend class SpotLevelHelper;
        static CClassConstSP const TYPE;
        virtual ~Shift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(SpotLevel* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(SpotLevel* shift) = 0;
    };
    typedef Shift IShift;

    /** constructor with explicit spot level */
    SpotLevel(double spot);

    /** constructor with explicit spot level and floorNegFwdPrice*/
    SpotLevel(double spot, bool floorNegFwdPrice);

    void validatePop2Object();

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
     SpotLevel.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);

    //// given today's date, when is spot defined
    virtual DateTime spotDate(const DateTime& today) const;

    bool zeroNegFwd() const;
protected:
    SpotLevel(CClassConstSP clazz, double spot);   
private:
    friend class SpotLevelHelper;
    /** for reflection */
    SpotLevel();
    SpotLevel(const SpotLevel &rhs);
    SpotLevel& operator=(const SpotLevel& rhs);

    bool floorNegFwdPrice; // what to use for setWillZeroNegFwd
};


typedef smartConstPtr<SpotLevel> SpotLevelConstSP;
typedef smartPtr<SpotLevel> SpotLevelSP;

DRLIB_END_NAMESPACE

#endif
