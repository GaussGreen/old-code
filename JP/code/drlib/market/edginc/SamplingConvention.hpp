//----------------------------------------------------------------------------
//
//   Filename    : SamplingConvention.hpp
//
//   Description : Classes for handling holiday adjustment 
//                 for use with MarketObservable and ObservableHistory
//                 see AssetHistory QLib Design Proposal - Ian Stares Jan 2006
//
//   Author      : Ian Stares   
//
//   Date        : February 2 2006
//
//
//----------------------------------------------------------------------------

#ifndef SAMPCONVN_HPP
#define SAMPCONVN_HPP

#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Holiday.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL SamplingConvention :  public virtual IObject {
public:
    static CClassConstSP const TYPE;

    // gets an observation date given holidays and a sample date
    // typically used once we're down in an asset history
    // Returns false if obs is to be omitted
    virtual bool observationDate(const DateTime& sampleDate, 
                                 const Holiday*  hols,
                                 DateTime*       obsDate) const = 0;

    virtual bool rollAssetsTogether() const = 0;

    virtual bool scheduleMoves() const = 0;

    virtual bool isUnadjusted() const = 0;

    virtual bool isOmit() const = 0;

private:
    static void load(CClassSP& clazz);
};
typedef smartPtr<SamplingConvention> SamplingConventionSP;

class MARKET_DLL UnadjustedConvention :  public CObject,
                              public virtual SamplingConvention {
public:
    static CClassConstSP const TYPE;
    
    // gets an observation date given holidays and a sample date
    // typically used once we're down in an asset history
    virtual bool observationDate(const DateTime& sampleDate, 
                                 const Holiday*  hols,
                                 DateTime*       obsDate) const;

    virtual bool rollAssetsTogether() const;

    virtual bool scheduleMoves() const;

    virtual bool isUnadjusted() const;

    virtual bool isOmit() const;

    UnadjustedConvention();
private:

    static void load(CClassSP& clazz);

    static IObject* defaultUnadjustedConvention();
};

DRLIB_END_NAMESPACE

#endif




