//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyVolBase.hpp
//
//   Description : Energy volatility object created with an underlying energy 
//                 2 factor vol surface. Based on drcommodityvolcurve.h
//                
//
//   Author      : Sean Chen
//
//   Date        : June 01, 2005
//
//----------------------------------------------------------------------------
#ifndef ENERGY_VOL_BASE_HPP
#define ENERGY_VOL_BASE_HPP

#include "edginc/DateTime.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/DataDictionary.hpp"
#include "edginc/Maths.hpp"


#include <string>
using namespace std;

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL EnergyVolBase: public MarketObject, 
                     virtual public IGetMarket
{

public:

    friend class EnergyVolBaseHelper;

    static CClassConstSP const TYPE;

    virtual ~EnergyVolBase();

    void validatePop2Object();

    virtual double getATMVol(const DateTime& expiry, const DateTime& futureExpiry) const = 0; 

    virtual double getSmileVolByStrike(const DateTime& expiry, double strike) const = 0;
    virtual double getSmileVolByDelta(const DateTime& expiry, double delta) const = 0;

    DateTime getBaseDate() const;

    string getName() const;

    virtual void getMarket(const IModel* model, const MarketData* market);

protected:
    
    EnergyVolBase(const CClassConstSP& clazz);

    string name;
    DateTime baseDate;
    
};

typedef smartConstPtr<EnergyVolBase> EnergyVolBaseConstSP;
typedef smartPtr<EnergyVolBase> EnergyVolBaseSP;

// support for wrapper class
typedef MarketWrapper<EnergyVolBase> EnergyVolBaseWrapper;
    
DRLIB_END_NAMESPACE

#endif
