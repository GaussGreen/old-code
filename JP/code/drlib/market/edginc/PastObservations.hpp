//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : PastObservations.hpp
//
//   Description : Records historic (overriding) values for an asset
//                 including source/observation type for centralised sampling
//
//   Author      : Ian Stares
//
//   Date        : 15 May 2006
//

//
//----------------------------------------------------------------------------

#ifndef EDR_PASTOBS_HPP
#define EDR_PASTOBS_HPP

#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/ObservationType.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL PastObservations : public CObject{
public:
    static CClassConstSP const TYPE;

    // fields - public for ease of building PastValues
    CashFlowArray       samples;
    string              source;
    ObservationTypeSP   obsType;
    string              fxSource;
    ObservationTypeSP   fxObsType;

    const CashFlowArray& getSamples();

    PastObservations(CashFlowArray& smpls);

    static IObjectSP fromCashFlowArray(const IObjectSP& object, 
                                       CClassConstSP    requiredType);

private:
    PastObservations();
    static IObject* defaultPastObs();
    static void load(CClassSP& clazz);
};

typedef smartConstPtr<PastObservations> IPastObservationsConstSP;
typedef smartPtr<PastObservations> PastObservationsSP;

// support for arrays
typedef array<PastObservationsSP, PastObservations> PastObservationsArray;

DRLIB_END_NAMESPACE

#endif

