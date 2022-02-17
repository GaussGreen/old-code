//----------------------------------------------------------------------------
//
//   Filename    : ISDAConvention.hpp
//
//   Description : Classes for handling ISDA holiday adjustment 
//
//   Author      : Ian Stares   
//
//   Date        : June 30 2006
//

//
//----------------------------------------------------------------------------

#ifndef ISDACONVN_HPP
#define ISDACONVN_HPP

#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/SamplingConvention.hpp"

DRLIB_BEGIN_NAMESPACE

// class for single factor ISDA adjustment
class MARKET_DLL ISDAConvention1Factor :   public CObject,
                              public virtual SamplingConvention  {
public:
    static CClassConstSP const TYPE;

    static const string OMIT;

    // gets an observation date given holidays and a sample date
    // typically used once we're down in an asset history
    // Returns false if obs is to be omitted
    virtual bool observationDate(const DateTime& sampleDate, 
                                 const Holiday*  hols,
                                 DateTime*       obsDate) const;

    virtual bool rollAssetsTogether() const;

    virtual bool scheduleMoves() const;

    virtual bool isUnadjusted() const;

    virtual bool isOmit() const;

    virtual void validatePop2Object();

protected:
    // for ISDAConventionNFactor inheritance
    ISDAConvention1Factor(CClassConstSP const clazz);
    bool                omit; // are we using the OMIT adjustment method?

private:
    ISDAConvention1Factor();
    static void load(CClassSP& clazz);
    static IObject* defaultISDA1F();

    string              isdaBDA; // string for bad day adjustment (may be OMIT)
    bool                moveSchedule; // does schedule move with the observations
    //transient fields
    BadDayConventionSP  bda; // bad day adjustment object
};

// class for multi factor ISDA adjustment
class MARKET_DLL ISDAConventionNFactor :  public ISDAConvention1Factor {
public:
    static CClassConstSP const TYPE;

    virtual void validatePop2Object();

    virtual bool rollAssetsTogether() const;

private:
    ISDAConventionNFactor();
    static void load(CClassSP& clazz);
    static IObject* defaultISDANF();

    // do the assets adjust for hoidays on their own?
    bool rollSeparately;
};
typedef smartPtr<ISDAConventionNFactor> ISDAConventionNFactorSP;

DRLIB_END_NAMESPACE

#endif




