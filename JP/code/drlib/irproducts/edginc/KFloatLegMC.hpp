#ifndef _KFLOATLEG_MC_HPP
#define _KFLOATLEG_MC_HPP

#include "edginc/MonteCarlo.hpp"
#include "edginc/SVGenIProdCreator.hpp"
//#include "edginc/KFloatLeg.hpp"
#include "edginc/SVGenKFloat.hpp"

DRLIB_BEGIN_NAMESPACE

// class KFloatLeg;

class KFloatLegMC : public MCProductClient {
private:
    SVGenKFloatSP floatSVGen;
    SVKFloatSP floatSV;
public:
//     KFloatLegMC(const KFloatLeg* inst, SimSeriesSP simSeries, InstrumentSettlementSP instSettle);
    // we don't want to have dependency on the instrument, so instrument passes all we need.
    KFloatLegMC( //const KFloatLeg* inst, 
                SimSeriesSP simSeries,
                InstrumentSettlementSP instSettle,
                const DateTime&   valueDate, // = inst->getValueDate()
                CouponSchedDatesSP sched, // = inst->sched
                YieldCurveConstSP       yc, //  inst->discount.getSP()
                IProdCreatorSP    _index, // inst->index
                DoubleArraySP     notionals, // inst->notionals
                DoubleArraySP     weights, // inst->weights
                DoubleArraySP     spreads,
                DoubleArraySP     dcfs,
                const RateType::Enum& rateType,
                const DateTimeArray&  principalDates, //inst->principalDates
                const DoubleArray     principalPayments);

    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen);
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const;
    virtual void payoff(const IPathGenerator* pathGen, IMCPrices& prices);
};

DRLIB_END_NAMESPACE

#endif
