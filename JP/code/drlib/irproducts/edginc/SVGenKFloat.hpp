#ifndef SVGEN_KFLOAT_HPP
#define SVGEN_KFLOAT_HPP

#include "edginc/SVGenKComponent.hpp"
#include "edginc/SVGenIndexSpecIR.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CouponSched.hpp"

DRLIB_BEGIN_NAMESPACE

class  IRPRODUCTS_DLL SVGenKFloat : public SVGenKComponent {
private :
    SVGenIndexSpecIRSP irSVGen;

    /*The following variables copied from KFloatLeg.hpp */
    CouponSchedDatesSP sched;

    DoubleArraySP   notionals;  // notional
    DoubleArraySP   weights;    // coupon is weights[]*index + spreads[]
    DoubleArraySP   spreads;    
    DoubleArraySP   dcfs;       // If dcfs are provided, dcc is ignored
    // for each coupon - end

    //bool isRateSimple;    // simple or continuous compounding rate
    RateType::Enum rateType;
    DayCountConventionSP dcc; // day count convention

    /****************** transient fields ************/
    // principal payment fields calculated in the setup function
    DateTimeArray principalDates;
    DoubleArray principalPayments;

    SVGenDiscFactorSP dfResetGen, dfPayGen;

public:
    SVGenKFloat(CouponSchedDatesSP sched, 
        DoubleArraySP notionals, 
        DoubleArraySP weights, 
        DoubleArraySP spreads, 
        DoubleArraySP dcfs, 
        //bool isRateSimple, 
        const RateType::Enum& rateType,
        DateTimeArray principalDates, 
        DoubleArray principalPayments,
        SVGenIndexSpecIRSP irSVGen,
        SVGenDiscFactorSP dfResetGen,
        SVGenDiscFactorSP dfPayGen) :
    irSVGen(irSVGen), sched(sched), notionals(notionals), weights(weights), 
        spreads(spreads), dcfs(dcfs), rateType(rateType), principalDates(principalDates),
        principalPayments(principalPayments),dfResetGen(dfResetGen), dfPayGen(dfPayGen) {}

    virtual SVIProdCreatorSP getSVIProdCreator(IStateVariableGen::IStateGen* pathGen) const;
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const;
};
DECLARE_REF_COUNT(SVGenKFloat);

class  IRPRODUCTS_DLL SVKFloat : public virtual SVKComponent {
public:
    virtual double elem(size_t index);
    virtual bool doingPast() const;
    SVKFloat(CouponSchedDatesSP sched, 
        DoubleArraySP notionals, 
        DoubleArraySP weights, 
        DoubleArraySP spreads, 
        DoubleArraySP dcfs, 
        //bool isRateSimple, 
        const RateType::Enum& rateType,
        DateTimeArray principalDates, 
        DoubleArray principalPayments,
        SVIndexSpecIRSP irSV,
        SVDiscFactorSP dfReset,
        SVDiscFactorSP dfPay) :
    irSV(irSV), sched(sched), notionals(notionals), weights(weights), 
        spreads(spreads), dcfs(dcfs), rateType(rateType), principalDates(principalDates),
        principalPayments(principalPayments), dfReset(dfReset), dfPay(dfPay) {} 

private:
    SVIndexSpecIRSP irSV;
    CouponSchedDatesSP sched;

    DoubleArraySP   notionals;  // notional
    DoubleArraySP   weights;    // coupon is weights[]*index + spreads[]
    DoubleArraySP   spreads;    
    DoubleArraySP   dcfs;       // If dcfs are provided, dcc is ignored
    // for each coupon - end

    //bool isRateSimple;    // simple or continuous compounding rate
    RateType::Enum rateType;

    /****************** transient fields ************/
    // principal payment fields calculated in the setup function
    DateTimeArray principalDates;
    DoubleArray principalPayments;

    SVDiscFactorSP dfReset, dfPay;
};

DECLARE_REF_COUNT(SVKFloat);

DRLIB_END_NAMESPACE
#endif
