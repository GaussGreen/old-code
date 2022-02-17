//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : BCStrikeMappingTweakBase.hpp
//
//   Description : Tweak to shift/set the Base Correlation Strike Mapping
//
//   Author      : Jose Hilera
//
//   Date        : 11 October 2005
//
//----------------------------------------------------------------------------


#ifndef QLIB_BCSTRIKEMAPPINGTWEAKBASE_HPP
#define QLIB_BCSTRIKEMAPPINGTWEAKBASE_HPP

#include "edginc/ScalarPerturbation.hpp"

DRLIB_BEGIN_NAMESPACE

/** Perturbation for strike mapping scenarios */
class RISKMGR_DLL BCStrikeMappingTweakBase: public ScalarPerturbation {
public:
    static CClassConstSP const TYPE;

    /** What an object must implement to be tweakable for types derived from
        BCStrikeMappingTweakBase. */
    class RISKMGR_DLL IShift{
    public:
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name - used to determine
            whether to tweak the object */
        virtual string sensName(BCStrikeMappingTweakBase* shift) const = 0;
        
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface, false otherwise */
        virtual bool sensShift(BCStrikeMappingTweakBase* shift) = 0;
    };

    /** Shifts the original StrikeMapping and returns the new one */
    virtual double applyShift(double& unadjStrikeMapping) = 0;

    /** returns the interface identifying what an object has to do in order
        to support the tweak that this object represents */
    CClassConstSP shiftInterface() const;

    /** Returns true if the supplied object matches the supplied name
        for this sensitivity.  The object must implement this sensitivity's
        Shift interface */
    virtual bool nameMatches(const OutputName& name,
                             IObjectConstSP    obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list */
    virtual void appendName(OutputNameArray&  namesList,
                            IObjectConstSP    obj);
     
    /**
     * @param obj The object to shift. The object must implement this 
     sensitivity's Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
protected:
    BCStrikeMappingTweakBase(CClassConstSP clazz = TYPE);
    
private:
    /** for reflection */
    static void load(CClassSP& clazz);
    BCStrikeMappingTweakBase(const BCStrikeMappingTweakBase&);
    BCStrikeMappingTweakBase& operator=(const BCStrikeMappingTweakBase&);
};


DRLIB_END_NAMESPACE

#endif
