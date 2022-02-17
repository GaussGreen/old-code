//   Filename    : PhysicalAndCashSettlement.hpp
//
//   Description : written for DailyAccumulator.cpp
//
//   Author      : Keith Law
//
//   Date        : 8 Dec 2006
//
//
//   $Log:
//----------------------------------------------------------------------------

#include "edginc/Settlement.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE


class MARKET_DLL PhysicalAndCashSettlement : public CObject{
public:
    static CClassConstSP const TYPE;
    friend class PhysicalAndCashSettlementHelper;

    PhysicalAndCashSettlement(); //default constructor
  
    bool                    isPhysicalAndCash; // physical and cash settle
    double                  physicalWeight;
    double                  cashWeight;
    bool                    phyDeliveredAtSpot;
};

typedef smartPtr<PhysicalAndCashSettlement> PhysicalAndCashSettlementSP;

DRLIB_END_NAMESPACE

