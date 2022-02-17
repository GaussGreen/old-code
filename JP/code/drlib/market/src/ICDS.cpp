//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : ICDS.cpp
//
//   Description : Interface for a CDS object
//
//   Date        : 18 January 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ICDS.hpp"
#include "edginc/IForwardRatePricer.hpp"
#include "edginc/IDecretionCurve.hpp"
#include "edginc/IBadDayAdjuster.hpp"

DRLIB_BEGIN_NAMESPACE

/*=============================================================================
 * ICDS
 *===========================================================================*/
ICDS::~ICDS(){}

void ICDS::load(CClassSP& clazz){
    REGISTER_INTERFACE(ICDS, clazz);
    EXTENDS(ICreditVanillaInstrument);
    EXTENDS(IHasCreditFeeLeg);
    EXTENDS(IHasCreditContingentLeg);
}

CClassConstSP const ICDS::TYPE = CClass::registerInterfaceLoadMethod(
    "ICDS", typeid(ICDS), load);

/**This has been overridden to set getPV = getFeeLegPV + getContingentLegPV.
     * If you override this, you should still ensure that this is the case!*/
double ICDS::getPV(
    const DateTime&              valuationDate,
    const DateTime&              settlementDate,
    const IDiscountCurveRisky&   crv,
    const IDecretionCurveConstSP prepay,
    IForwardRatePricerSP         model,
    IBadDayAdjusterConstSP       bda
    ) const
{
    ICreditFeeLegSP feeLeg = this->getFeeLeg();
    ICreditContingentLegSP ctgLeg = this->getContingentLeg();

    //legs may be optional in the parent structure
    double feePv = 0.0;
    if (!!feeLeg)
    {
        DateTime earliestRiskyDate = (!ctgLeg) ? valuationDate
                                             : ctgLeg->firstObservationStartDate();
        DateTime latestRiskyDate = (!ctgLeg) ? feeLeg->getLastPayDate()
                                             : ctgLeg->lastObservationEndDate();

        feePv = feeLeg->getFeeLegPV(valuationDate, settlementDate,
                                    earliestRiskyDate, latestRiskyDate,
                                    crv, crv, prepay,
                                    true, getAccrualDcc(), false, model);
    }

    double ctgPv = 0.0;
    if (!!ctgLeg)
    {
        ctgPv = ctgLeg->getContingentLegPV(valuationDate, settlementDate, crv, bda);
    }

    return feePv + ctgPv;
}

/**This has been overridden to set getPV = getFeeLegPV + getContingentLegPV.
 * If you override this, you should still ensure that this is the case!*/
double ICDS::getPV(
    const DateTime&              valuationDate,
    const IDiscountCurveRisky&   crv,
    const IDecretionCurveConstSP prepay,
    IForwardRatePricerSP         model,
    IBadDayAdjusterConstSP       bda
    ) const
{
    ICreditFeeLegSP feeLeg = this->getFeeLeg();
    ICreditContingentLegSP ctgLeg = this->getContingentLeg();

    //legs may be optional in the parent structure
    double feePv = 0.0;
    if (!!feeLeg)
    {
        DateTime earliestRiskyDate = (!ctgLeg) ? valuationDate
                                             : ctgLeg->firstObservationStartDate();
        DateTime latestRiskyDate = (!ctgLeg) ? feeLeg->getLastPayDate()
                                             : ctgLeg->lastObservationEndDate();

        feePv = feeLeg->getFeeLegPV(valuationDate,
                                    earliestRiskyDate, latestRiskyDate,
                                    crv, crv, prepay,
                                    true, getAccrualDcc(), model);
    }

    double ctgPv = 0.0;
    if (!!ctgLeg)
    {
        ctgPv = ctgLeg->getContingentLegPV(valuationDate, crv, bda);
    }

    return feePv + ctgPv;
}

/*=============================================================================
 * ICDSConvention
 *===========================================================================*/
ICDSConvention::~ICDSConvention(){}

void ICDSConvention::load(CClassSP& clazz){
    REGISTER_INTERFACE(ICDSConvention, clazz);
    EXTENDS(IObject);
}

CClassConstSP const ICDSConvention::TYPE = CClass::registerInterfaceLoadMethod(
    "ICDSConvention", typeid(ICDSConvention), load);

DRLIB_END_NAMESPACE
