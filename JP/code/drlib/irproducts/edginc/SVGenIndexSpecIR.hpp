#ifndef SVGEN_INDEX_SPEC_HPP
#define SVGEN_INDEX_SPEC_HPP

#include "edginc/SVGenIProdCreator.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/SVGenDiscFactor.hpp"

DRLIB_BEGIN_NAMESPACE

class  IRPRODUCTS_DLL SVGenIndexSpec : public SVGenIProdCreator {
};
DECLARE_REF_COUNT(SVGenIndexSpec);

class  IRPRODUCTS_DLL SVGenIndexSpecIR : public SVGenIndexSpec {
private:
    SVGenDiscFactorSP dfResetGen, dfOffsetGen;
    SVDiscFactorSP dfReset, dfOffset;
    MaturityPeriodSP tenor;
public:
    SVGenIndexSpecIR(SVGenDiscFactorSP dfResetGen, SVGenDiscFactorSP dfOffsetGen, MaturityPeriodSP tenor);
    virtual SVIProdCreatorSP getSVIProdCreator(IStateVariableGen::IStateGen* pathGen) const;
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const;
};
DECLARE_REF_COUNT(SVGenIndexSpecIR);

class  SVIndexSpecIR : public SVIProdCreator {
public:
    virtual double elem(size_t index);
    virtual bool doingPast() const;
    SVIndexSpecIR(SVDiscFactorSP dfReset, SVDiscFactorSP dfOffset, MaturityPeriodSP tenor) : dfReset(dfReset), dfOffset(dfOffset), tenor(tenor) {}
private:
    SVDiscFactorSP dfReset, dfOffset;
    MaturityPeriodSP tenor;
};
DECLARE_REF_COUNT(SVIndexSpecIR);

DRLIB_END_NAMESPACE

#endif
