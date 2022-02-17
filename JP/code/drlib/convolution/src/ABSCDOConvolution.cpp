 //----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : ABSCDOConvolution.cpp
//
//   Description : Convolution Algorithm
//
//   Date        : Aug 2004
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CCMSkew.hpp"
#include "edginc/CCMCalibration.hpp"
#include "edginc/CCMConvolution.hpp"
#include "edginc/ConvolutionProduct.hpp"
#include "edginc/ABSCDOConvolution.hpp"
#include "edginc/TrancheLossCalculator.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE
#define TINY 1e-12
#define REALLY_TINY 1e-64
#define TOWER_PROBA_EPSI      1e-12
#define SKEW_FCUMUL_NB_POINT  100L

typedef ABSCDOConvolution::ConditionalCurve        ConditionalCurve;        // for ease
typedef ABSCDOConvolution::NameParam               NameParam;               // for ease
typedef ABSCDOConvolution::NameParamSP             NameParamSP;             // for ease
typedef ABSCDOConvolution::NameParamArray          NameParamArray;          // for ease
typedef ABSCDOConvolution::MarketScenario          MarketScenario;          // for ease
typedef ABSCDOConvolution::MarketScenarioSP        MarketScenarioSP;        // for ease
typedef ABSCDOConvolution::MarketScenarioArray     MarketScenarioArray;     // for ease

//// constructor that zeroes everything
NameParam::NameParam():
    survival(0.0), beta(0.0), qM(0.0),
    ntl(0.0), R(0.0)
{
}
 
double NameParam::u2cprobs(double s, double M, double b)
{
    if (fabs(s-1.) < 3e-16)
        s = 1.;
    if (fabs(b) < 3e-16)
        return s;
    double t = CCMSkew::fInv(s, b, qM, SKEW_FCUMUL_NB_POINT);
    return CCMSkew::condSurvivalProba(1.0, t, b, qM, M);
}

// create conditional curves for this name under a given scenario
void NameParam::calcConditionalCurve(const double lossMF, const double decMF, 
                                     const DateTimeArray& timeline, 
                                     ConditionalCurve& curve, 
                                     double lossMarketWeight,
                                     double decretionMarketWeight)
{
    try {
        // merge convolution timeline with decretion timeline
        DateTimeArray finalTimeline = 
            DateTime::merge(timeline, *(decretion->getStepDates()));
        int numDates = finalTimeline.size(); // exclude value date
        curve.resize(numDates);

        DateTime valueDate = finalTimeline[0];
        curve.addSlice(valueDate, 
                       defaultrate->calcDefaultPV(valueDate, valueDate),
                       decretion->pv(valueDate, valueDate));
        for (int i = 1; i < numDates; ++i) {
            DateTime currentDate = finalTimeline[i];
            double unConditionalSurvivialProb = 
                defaultrate->calcDefaultPV(valueDate, currentDate);
            double conditionalSurvivalProb = u2cprobs(unConditionalSurvivialProb, lossMF, beta);
            double unConditionalDecretion = decretion->pv(valueDate, currentDate);
            double conditionalDecretion = u2cprobs(unConditionalDecretion, 
                                                   lossMF * lossMarketWeight + decMF * decretionMarketWeight, 
                                                   decBeta);
            curve.addSlice(currentDate, conditionalSurvivalProb, conditionalDecretion);
        }
    } catch (exception& e) {
        throw ModelException(e, "NameParam::calcConditionalCurve");
    }
}

void NameParam::getDistrib(const DateTime currentDate, const ConditionalCurve& curve, int& tpt,
                           double* bottomLoss, double* bottomProb, double* topLoss, double* topProb)
{
    try {
        tpt = 0;
        double leftSurvivalProb = curve.getSlice(0).survivialProb;
        double rightSurvivalProb;
        double startDecretion = curve.getSlice(0).decretion;
        double rightDecretion;
        for (int i = 1; i < curve.getSize(); ++i) {
            const ConditionalCurve::Slice& slice = curve.getSlice(i);
            if (slice.date > currentDate) break;
            rightSurvivalProb = slice.survivialProb;
            rightDecretion = slice.decretion;
            bottomProb[tpt] = topProb[tpt] = 1 - (leftSurvivalProb - rightSurvivalProb);
            bottomLoss[tpt] = ntlLoss * rightDecretion;
            topLoss[tpt] = ntlRecovery * rightDecretion + ntl * (startDecretion - rightDecretion);
            ++tpt;
            leftSurvivalProb = rightSurvivalProb;
        }
        // handle the last segment
        bottomProb[tpt] = topProb[tpt] = 1 - rightSurvivalProb;
        bottomLoss[tpt] = 0;
        topLoss[tpt] = ntl * (startDecretion - rightDecretion);
        ++tpt;
    } catch (exception& e) {
        throw ModelException(e, "NameParam::getDistrib");
    }
}

/////////////////////////////////////////////////////////////////////////////////////

MarketScenario::MarketScenario(
                   const double lossMF, const double lossWeight,
                   const double decMF, const double decWeight)
  : lossMF(lossMF), lossWeight(lossWeight), decMF(decMF), decWeight(decWeight) 
{
}

void MarketScenario::createConditionalCurves(NameParamArray& basketInfo, 
                                             const DateTimeArray& timeline,
                                             double lossMarketWeight,
                                             double decretionMarketWeight)
{
    try {
        conditionalCurves.reserve(basketInfo.size());
        ConditionalCurve tmp; size_t i;
        for (i = 0; i < basketInfo.size(); ++i) {
            conditionalCurves.push_back(tmp);
        }
        for (i = 0; i < basketInfo.size(); ++i) {
            basketInfo[i]->calcConditionalCurve(lossMF, decMF, timeline, 
                                                conditionalCurves[i], 
                                                lossMarketWeight,
                                                decretionMarketWeight);
        }
    } catch (exception& e) {
        throw ModelException(e, "MarketScenario::createConditionalCurves");
    }
}

void MarketScenario::setDistribution(const int nbLoss)
{
    const char routine[] = "MarketScenario::setDistribution";
    try
    {
        //add the distribution
        bottomLossDistrib = LossDistributionSP(new LossDistribution(nbLoss, 0));
        topLossDistrib = LossDistributionSP(new LossDistribution(nbLoss, 0));
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

void MarketScenario::accumulateDistribution(DoubleArray& dist, LossDistribution::DistributionType type)
{
    const char routine[] = "MarketScenario::accumulateDistribution";
    try
    {
        LossDistributionSP distrib;
        switch (type) {
        case LossDistribution::BOTTOMLOSS:
            distrib = bottomLossDistrib; break;
        case LossDistribution::TOPLOSS:
            distrib = topLossDistrib; break;
        }
        const DoubleArray& lossDist = distrib->getLossDistribution();
        //double lossI;
        for (int i=0; i<lossDist.size(); i++)
        {
            dist[i] += lossDist[i] * lossWeight;
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/////////////////////////////////////////////////////////////////////////////////////

void GaussQuadIntegMethod::calcPoints(int numPoints) {
    if (numPoints > MAX_GAUSS_INTEG_POINTS)
        throw ModelException("num of points too big! max is 500");
    imsl_d_gauss_quad_rule(numPoints, weights, points, IMSL_HERMITE, 0);
    for (int i = 0; i < numPoints; ++i) {
        points[i] *= sqrt(2.);
        weights[i] /= sqrt(Maths::PI); 
    }
}


void ABSCDOConvolution::calcABSCDOLossDistribution(DateTime currentDate, int lossSize,
                                                   const MarketScenarioMatrix& scenarios,
                                                   const NameParamArray& basket)
{
    static const int MAXDECINTERVALS = 30 * 12;
    int tpt;
    double lossBottom[MAXDECINTERVALS]; double probBottom[MAXDECINTERVALS];
    double lossTop[MAXDECINTERVALS]; double probTop[MAXDECINTERVALS];

    try {
        int nbName = basket.size();
        // for each decretion market factor
        for (size_t i = 0; i < scenarios.size(); ++i) {
            const MarketScenarioArray& currentDecScenario = scenarios[i];
            // for each loss market factor
            for (size_t j = 0; j < currentDecScenario.size(); ++j) {
                const MarketScenarioSP& currentLossScenario = currentDecScenario[j];
                currentLossScenario->setDistribution(lossSize);
                for (int k = 0; k < nbName; ++k) {
                    basket[k]->getDistrib(currentDate, 
                                          currentLossScenario->getConditionalCurve(k),
                                          tpt, lossBottom, probBottom,
                                          lossTop, probTop);
                    currentLossScenario->getLossDistribution()->convoluteDistribution(
                                          tpt, probBottom, lossBottom);
                    currentLossScenario->getTopLossDistribution()->convoluteDistribution(
                                          tpt, probTop, lossTop);
                }
            }
        }
    } catch (exception& e) {
        throw ModelException(e, "ABSCDOConvolution::calcABSCDOLossDistribution");
    }
}

void ABSCDOConvolution::calcABSCDOTrancheLossFast(DateTime currentDate, double K1, double K2, 
                                              double lossUnit, double port,
                                              const MarketScenarioMatrix& scenarios,
                                              const NameParamArray& basket,
                                              double& loss, double& topLoss)
{
    loss = topLoss = 0;

    static const int MAXDECINTERVALS = 30 * 12;
    int tpt;
    double lossBottom[MAXDECINTERVALS]; double probBottom[MAXDECINTERVALS];
    double lossTop[MAXDECINTERVALS]; double probTop[MAXDECINTERVALS];

    try {
        // for each decretion market factor
        int nbName = basket.size();
        for (size_t i = 0; i < scenarios.size(); ++i) {
            const MarketScenarioArray& currentDecScenario = scenarios[i];
            double lossMFLoss = 0;
            double lossMFTopLoss = 0;
            // for each loss market factor
            for (size_t j = 0; j < currentDecScenario.size(); ++j) {
                const MarketScenarioSP& currentLossScenario = currentDecScenario[j];
                double nameMeanBottom, nameMeanTop, nameVarBottom, nameVarTop, nameMaxBottom, nameMaxTop;
                double portMeanBottom(0), portMeanTop(0), portVarBottom(0), portVarTop(0);
                double portMaxBottom(0), portMaxTop(0);
                for (int k = 0; k < nbName; ++k) {
                    basket[k]->getDistrib(currentDate, 
                                          currentLossScenario->getConditionalCurve(k),
                                          tpt, lossBottom, probBottom,
                                          lossTop, probTop);

                    // update mean and variance
                    nameMeanBottom = nameMeanTop = nameVarBottom = nameVarTop = 0;
                    nameMaxBottom = nameMaxTop = 0;
                    for (int j = 0; j < tpt; ++j) {
                        nameMeanBottom += lossBottom[j] * (1 - probBottom[j]);
                        nameMeanTop += lossTop[j] * (1 - probTop[j]);
                        nameVarBottom += lossBottom[j] * lossBottom[j] * (1 - probBottom[j]);
                        nameVarTop += lossTop[j] * lossTop[j] * (1 - probTop[j]);
                        if (nameMaxBottom < lossBottom[j]) nameMaxBottom = lossBottom[j];
                        if (nameMaxTop < lossTop[j]) nameMaxTop = lossTop[j];
                    }
                    nameVarBottom -= nameMeanBottom * nameMeanBottom;
                    nameVarTop -= nameMeanTop * nameMeanTop;
                    portMeanBottom += nameMeanBottom;
                    portMeanTop += nameMeanTop;
                    portVarBottom += nameVarBottom;
                    portVarTop += nameVarTop;
                    portMaxBottom += nameMaxBottom;
                    portMaxTop += nameMaxTop;
                }

                // normalize results
				if (portMaxBottom == 0) {
					portMeanBottom = portVarBottom = 0;
				} else {
					portMeanBottom /= portMaxBottom;
					portVarBottom /= portMaxBottom * portMaxBottom;
				}

				if (portMaxTop == 0) {
					portMeanTop = portVarTop = 0;
				} else {
					portMeanTop /= portMaxTop;
					portVarTop /= portMaxTop * portMaxTop;
				}

                // do beta approximation
                portMaxBottom *= lossUnit;
                portMaxTop *= lossUnit;
                double bottomLo = ccmBetaTrancheLossInteg(K1/portMaxBottom, portMeanBottom, portVarBottom);
                double bottomHi = ccmBetaTrancheLossInteg(K2/portMaxBottom, portMeanBottom, portVarBottom);
                double topLo = ccmBetaTrancheLossInteg((port - K2)/portMaxTop, portMeanTop, portVarTop);
                double topHi = ccmBetaTrancheLossInteg((port - K1)/portMaxTop, portMeanTop, portVarTop);
                double bottomEL = (bottomLo - bottomHi)*portMaxBottom;
                double topEL = (topLo - topHi)*portMaxTop;

                lossMFLoss += currentLossScenario->getLossWeight() * bottomEL;
                lossMFTopLoss += currentLossScenario->getLossWeight() * topEL;
            }
            loss += currentDecScenario[0]->getDecWeight() * lossMFLoss;
            topLoss += currentDecScenario[0]->getDecWeight() * lossMFTopLoss;
        }
    } catch (exception& e) {
        throw ModelException(e, "ABSCDOConvolution::calcABSCDOTrancheLossFast");
    }
}

void ABSCDOConvolution::calcABSCDOTrancheLossSaddle(
                                              DateTime currentDate, double k1, double k2, 
                                              double lossUnit, double port,
                                              int lossMarketFactors, const double* lossMF, const double* lossWeights,
                                              const MarketScenarioMatrix& scenarios,
                                              const NameParamArray& basket,
                                              double& loss, double& topLoss, int numSaddlePoints)
{
    loss = topLoss = 0;

    static const int MAXDECINTERVALS = 30 * 12;
    int tpt;
    double lossBottom[MAXDECINTERVALS]; double probBottom[MAXDECINTERVALS];
    double lossTop[MAXDECINTERVALS]; double probTop[MAXDECINTERVALS];

    try {
        int nbName = basket.size();
        // for each decretion market factor
        for (size_t i = 0; i < scenarios.size(); ++i) {
            SaddlePoint lossSaddle(lossMarketFactors, lossMF, lossWeights,
                                   numSaddlePoints, nbName, 
                                   port, lossUnit, k1, k2);
            SaddlePoint topLossSaddle(lossMarketFactors, lossMF, lossWeights,
                                      numSaddlePoints, nbName, 
                                      port, lossUnit, port-k2, port-k1);

            const MarketScenarioArray& currentDecScenario = scenarios[i];
            // for each loss market factor
            for (size_t j = 0; j < currentDecScenario.size(); ++j) {
                const MarketScenarioSP& currentLossScenario = currentDecScenario[j];
                for (int k = 0; k < nbName; ++k) {
                    basket[k]->getDistrib(currentDate, 
                                          currentLossScenario->getConditionalCurve(k),
                                          tpt, lossBottom, probBottom,
                                          lossTop, probTop);
                    lossSaddle.setParams(j, k, tpt, lossBottom, probBottom);
                    topLossSaddle.setParams(j, k, tpt, lossTop, probTop);
                }
            }
            loss += currentDecScenario[0]->getDecWeight() * lossSaddle.getExpLoss();
            topLoss += currentDecScenario[0]->getDecWeight() * topLossSaddle.getExpLoss();
        }
    } catch (exception& e) {
        throw ModelException(e, "ABSCDOConvolution::calcABSCDOTrancheLossSaddle");
    }
}


void ABSCDOConvolution::calcABSCDOLossDensities(
    DoubleArray&      density,
    LossDistribution::DistributionType type,
    const MarketScenarioMatrix& scenarios,
    int nbName,
    int lossSize)
{
    const char routine[] = "ABSCDOConvolution::calcABSCDOLossDensities";
    try
    {
        density.clear();

        // treat the case nbName == 0 separately
        if (nbName == 0)
        {
            density.resize(1);
            density[0] = 1.0;
        }    
        else
        {
            density.resize(lossSize);
            density[0] = 0;

            // for each decretion market factor
            for (size_t i = 0; i < scenarios.size(); ++i) {
                const MarketScenarioArray& currentDecScenario = scenarios[i];

                DoubleArray work;
                work.clear();
                work.resize(lossSize);

                for (size_t j = 0; j < currentDecScenario.size(); j++)
                {
                    currentDecScenario[j]->accumulateDistribution(work, type);
                }
            
                // update loss density and loss density cond cpty survival 
                for (int k = 0; k < lossSize; ++k)
                {
                    density[k] += currentDecScenario[0]->getDecWeight() * work[k];
                }
            }
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

////////////////////////////////////////////////////////////////////////////////////


DRLIB_END_NAMESPACE
