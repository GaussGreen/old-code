#include "edginc/config.hpp"
#include "edginc/SVGenIndexSpecIR.hpp"

DRLIB_BEGIN_NAMESPACE

SVGenIndexSpecIR::SVGenIndexSpecIR(SVGenDiscFactorSP dfGen, 
                                   SVGenDiscFactorSP dfOffsetGen, 
                                   MaturityPeriodSP tenor) :
                                dfResetGen(dfGen), dfOffsetGen(dfOffsetGen), tenor(tenor) {}

void SVGenIndexSpecIR::collectStateVars(IStateVariableCollectorSP svCollector) const {
    svCollector->append(dfResetGen.get());
    svCollector->append(dfOffsetGen.get());
}

SVIProdCreatorSP SVGenIndexSpecIR::getSVIProdCreator(IStateVariableGen::IStateGen* pathGen) const {
    SVDiscFactorSP dfReset(dfResetGen->getSVDiscFactor(pathGen));
    SVDiscFactorSP dfOffset(dfOffsetGen->getSVDiscFactor(pathGen));

    return SVIProdCreatorSP(new SVIndexSpecIR(dfReset, dfOffset, tenor));
}

double SVIndexSpecIR::elem(size_t index) {
    double dfResetVal = dfReset->getDF(index);
    double dfOffsetVal = dfOffset->getDF(index);
    double freq = tenor->annualFrequency();
    //double value = (dfResetVal/dfOffsetVal - 1)*tenor->annualFrequency();
    double value = (log(dfResetVal) - log(dfOffsetVal))*freq;
    return value;
    //return ((*dfReset)[index]/(*dfOffset)[index] - 1)/tenor->annualFrequency();
}

bool SVIndexSpecIR::doingPast() const {return dfReset->doingPast();}

DRLIB_END_NAMESPACE
