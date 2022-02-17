#include "edginc/config.hpp"
#include "edginc/KFloatLegMC.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/IndexSpecIR.hpp"

DRLIB_BEGIN_NAMESPACE

KFloatLegMC::KFloatLegMC( //const KFloatLeg* inst, 
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
                          const DoubleArray     principalPayments
                        ) :
        MCProductClient( IMultiMarketFactors::asMulti( yc ).get(),
                         valueDate,
                         yc.get(),
                         IRefLevelSP( IRefLevel::Util::makeZero( valueDate ) ),  // fix!
                         simSeries,  // fix!
                         IPastValuesSP( IPastValues::Util::makeTrivial( valueDate, 0.0 ) ),  // fix
                         instSettle.get(),
                         sched->pay.back() )
        //MCProductClient(inst->getValueDate(), inst->discount.get(), simSeries)
{
    SVGenDiscFactorSP dfResetGen( new SVGenDiscFactor( valueDate, yc, sched->resetEff.getDates() ) );
    /*For now assume that the underlying index in KFloatLeg is of type IndexSpecIR.  In general,
    it's required to be only an IProdCreatorSP.  Have to FIX THIS */
    IndexSpecIRSP index( IndexSpecIRSP::dynamicCast( _index ) );
    SVGenDiscFactorSP dfPayGen( new SVGenDiscFactor( valueDate, yc, sched->pay.getDates() ) );
    SVGenDiscFactorSP dfOffsetGen( new SVGenDiscFactor( valueDate, yc, sched->pay.getDates() ) );
    SVGenIndexSpecIRSP indexSpecIRGen( new SVGenIndexSpecIR( dfResetGen, dfOffsetGen, index->tenor ) );
    floatSVGen = SVGenKFloatSP( new SVGenKFloat( sched, 
                                                 notionals, 
                                                 weights, 
                                                 spreads,
                                                 dcfs, 
                                                 rateType, 
                                                 principalDates, 
                                                 principalPayments,
                                                 indexSpecIRGen,
                                                 dfResetGen, 
                                                 dfPayGen ) );
}

void KFloatLegMC::pathGenUpdated( IStateVariableGen::IStateGen* newPathGen )
{
    //floatSV = SVKFloatSP::dynamicCast(floatSVGen->getSVIProdCreator(newPathGen));
    floatSV = DYNAMIC_POINTER_CAST<SVKFloat>( floatSVGen->getSVIProdCreator( newPathGen ) );
}

void KFloatLegMC::collectStateVars( IStateVariableCollectorSP svCollector ) const
{
    svCollector->append( floatSVGen.get() );
}

void KFloatLegMC::payoff( const IPathGenerator* pathGen, IMCPrices& prices )
{
    prices.add( floatSV->elem( 0 ) );
}

DRLIB_END_NAMESPACE
