//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SpotLevelProbability.hpp
//
//   Description : spot level scenario - for probability supplied, use the
//                 vol of the asset to determine the spot level that encloses
//                 this confidence interval
//
//   Author      : Andrew J Swain
//
//   Date        : 16 June 2002
//
//
//----------------------------------------------------------------------------

#ifndef SPOTLEVELPROBABILITY__HPP
#define SPOTLEVELPROBABILITY__HPP
#include "edginc/Perturbation.hpp"
#include "edginc/Asset.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for spot level scenario where spot depends on 
    confidence interval & vol */
class MARKET_DLL SpotLevelProbability: public Perturbation {
public:
    static CClassConstSP const TYPE;

    /** what the asset should set its spot to */
    double spotLevel(const DateTime& today, const Asset* asset) const;

    /** validation code to be called after object construction */
    virtual void validatePop2Object();

    /** What an object must implement to support SpotLevelProbability */
    class MARKET_DLL Shift{
    public:
        friend class SpotLevelProbabilityHelper;
        static CClassConstSP const TYPE;
        virtual ~Shift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(SpotLevelProbability* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(SpotLevelProbability* shift) = 0;
    };
    typedef Shift IShift;

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
     SpotLevelProbability.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
private:
    friend class SpotLevelProbabilityHelper;
    /** for reflection */
    SpotLevelProbability();
    SpotLevelProbability(const SpotLevelProbability &rhs);
    SpotLevelProbability& operator=(const SpotLevelProbability& rhs);

    // fields
    double probability; // defines confidence interval
    double moneyness;   // given spot, defines strike for vol interp
    string maturity;    // with time, defines date for vol interp
    string time;
    OutputNameArraySP toTweak; // to be removed
};


typedef smartConstPtr<SpotLevelProbability> SpotLevelProbabilityConstSP;
typedef smartPtr<SpotLevelProbability> SpotLevelProbabilitySP;

DRLIB_END_NAMESPACE

#endif
