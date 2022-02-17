//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Perturbation.hpp
//
//   Description : Defines the ability to shift a named piece of data
//
//   Author      : Mark A Robson/Andrew J Swain
//
//   Date        : 19 Jan 2005
//
//
//----------------------------------------------------------------------------


#ifndef EDR_PERTURBATION_HPP
#define EDR_PERTURBATION_HPP
#include "edginc/OutputName.hpp"
#include "edginc/TweakID.hpp"
#include "edginc/TweakNameResolver.hpp"

DRLIB_BEGIN_NAMESPACE

/**  Defines the ability to shift a named piece of data */
class RISKMGR_DLL IPerturbation: virtual public IObject {
public:
    /* Define some strings that can be used by multiple scenarios to configure
     * the type of tweak to be performed.
     * Added here to enforce a consistent use of theses strings, and to avoid
     * creating them on the fly wherever needed. */
    static const string ABSOLUTE; // (A) ie, value += shift
    static const string RELATIVE; // (R) ie, value *= (1 + shift)
    static const string SET;      // (S) ie, value  = shift

    static CClassConstSP const TYPE;
    IPerturbation();
    
    ~IPerturbation();

    /** apply this Perturbation to the object with the specified name contained
        within the supplied object. This is, for clarity, a non-restorable
        shift ie there is no mechanism for undoing the shift. The return
        value indicates if anything was actually shifted (true => yes) */
    virtual bool findAndShift(IObjectSP         objectToShift, 
                              OutputNameConstSP name) = 0;

    /** whether the Perturbation needs to be applied before the
        "getMarket()" phase */
    virtual bool applyBeforeGetMarket() const;
  
private:
    static void load(CClassSP& clazz);
};

// typedef for smart pointers to Perturbations
typedef smartConstPtr<IPerturbation> IPerturbationConstSP;
typedef smartPtr<IPerturbation> IPerturbationSP;

/** Base class which most implementations can derive from. To do:
    add comment here saying what methods you have to implement */
class RISKMGR_DLL Perturbation: public CObject,
                    public virtual IPerturbation,
                    public virtual ITweakID,
                    public virtual ITweakNameResolver{
public:
    static CClassConstSP const TYPE;
    ~Perturbation();

    /** Returns this */
    virtual ITweakNameResolver* nameResolver();

    /** returns the name identifying the market data to be shifted. Returns
        null if not set */
    virtual OutputNameConstSP getMarketDataName() const;

    /** Does nothing */
    virtual void reset();

    /** apply this Perturbation to the object with the specified name contained
        within the supplied object. Returns true if something was shifted */
    virtual bool findAndShift(IObjectSP         objectToShift, 
                              OutputNameConstSP name);

protected:
    Perturbation(CClassConstSP clazz);

    OutputNameConstSP marketDataName; /* somewhat bogus - ideally remove
                                         getMarketDataName from 
                                         ITweakNameResolver */
private:
    static void load(CClassSP& clazz);
};


DRLIB_END_NAMESPACE
#endif
