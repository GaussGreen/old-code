//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : ABSCDOConvolution.hpp
//
//   Description : Convolution Algorithm
//
//   Date        : Fed 2006
//
//----------------------------------------------------------------------------

#ifndef EDR_ABSCDOConvolution_HPP
#define EDR_ABSCDOConvolution_HPP

#include "edginc/DefaultRates.hpp"
#include "edginc/LossDistribution.hpp"
#include "edginc/IDecretionCurve.hpp"
#include "edginc/SaddlePoint.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Maths.hpp"
#include "edginc/imsl.h"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE
class ConvolutionProduct;

/* ---------------------------------------------------------------------------
 * Integ Method
 */
#define CCM_INTEG_LOW       -7.5
#define CCM_INTEG_UP         7.5
#define CCM_INTEG_NB         101
struct CONVOLUTION_DLL IGIntegMethod{
    double l;    /* lower bound */
    double u;    /* upper bound */
    int    n;    /* nb of steps */
    double step; /* step size   */
    IGIntegMethod():l(CCM_INTEG_LOW), u(CCM_INTEG_UP), n(CCM_INTEG_NB),
                    step((u-l)/(n-1)){}

    double weightCalc(double m){
        return step * N1Density(m);
    }
};

#define MAX_GAUSS_INTEG_POINTS  500
struct CONVOLUTION_DLL GaussQuadIntegMethod {
    double weights[MAX_GAUSS_INTEG_POINTS];
    double points[MAX_GAUSS_INTEG_POINTS];
    void calcPoints(int numPoints);
};

class CONVOLUTION_DLL ABSCDOConvolution
{
public:
    // per scenario conditional curve (one per name)
    class CONVOLUTION_DLL ConditionalCurve
    {
    public:
        struct CONVOLUTION_DLL Slice {
            DateTime date;
            double survivialProb; // survival probability at this date
            double decretion; // expected remaining notional at this date
            Slice(DateTime d, double prob, double dec)
                : date(d), survivialProb(prob), decretion(dec) {}
        };

        ConditionalCurve() {}

        void resize(int length)
        { curve.reserve(length); }

        void addSlice(DateTime d, double prob, double dec)
        { curve.push_back(Slice(d, prob, dec)); }

        int getSize() const
        { return curve.size(); }

        const Slice& getSlice(int i) const
        { return curve[i]; }
    private:
        std::vector<Slice> curve;
        friend class ABSCDOConvolution;
    };

    class CONVOLUTION_DLL NameParam {
    public:
        string nameId;        // curve name
        int    index;         // name index
        double survival;      // survival probability, to be replaced by default rate curve
        double beta;          // beta
        double decBeta;       // decretion beta
        double qM;            // skew parameter (not used yet)
        double ntl;           // name notional (expressed in lossunit)
        double R;             // name recovery
        double ntlLoss;       // N*(1-R)
        double ntlRecovery;   // N*R
        IDecretionCurveConstSP decretion;   // name decretion curve
        DefaultRatesSP         defaultrate; // default rate curve
        NameParam();

        // create conditional curve for a scenario
        void calcConditionalCurve(const double lossMF, const double decMF,
                                  const DateTimeArray& timeline,
                                  ConditionalCurve& curve,
                                  double lossMarketWeight,
                                  double decretionMarketWeight);

        // calculate distribution of top and borrom loss for a given name given conditional curve
        void getDistrib(
            const DateTime currentDate, const ConditionalCurve& curve, int& tpt,
            double* bottomLoss, double* bottomProb, double* topLoss, double* topProb);

    private:
        // calc conditional survivial prob from unconditional prob
        double u2cprobs(double s, double M, double b);
    };

    DECLARE_REF_COUNT(NameParam);

    class CONVOLUTION_DLL MarketScenario
    {
    public:
        MarketScenario() {};
        MarketScenario(const double lossMF, const double lossWeight,
            const double decMF, const double decWeight);

        void createConditionalCurves(NameParamArray& basketInfo,
                                    const DateTimeArray& timeline,
                                    double lossMarketWeight,
                                    double decretionMarketWeight);

        int  numOfCurves() const
        { return conditionalCurves.size(); }

        const ConditionalCurve& getConditionalCurve(int i) const
        { return conditionalCurves[i]; }

        double getLossWeight()
        { return lossWeight; }

        double getDecWeight()
        { return decWeight; }

        void setDistribution(const int nbLoss);

        LossDistributionSP getLossDistribution() { return bottomLossDistrib; }
        LossDistributionSP getTopLossDistribution() { return topLossDistrib; }

        void accumulateDistribution(DoubleArray& dist,
                                    LossDistribution::DistributionType type);

    private:
        double              lossMF;                  // loss market factor
        double              lossWeight;              // loss weight
        double              decMF;                   // decretion market factor
        double              decWeight;               // decretion weight

        std::vector<ConditionalCurve> conditionalCurves;

        LossDistributionSP  bottomLossDistrib;
        LossDistributionSP  topLossDistrib;
    };
    typedef refCountPtr<MarketScenario>            MarketScenarioSP;
    typedef vector<MarketScenarioSP>               MarketScenarioArray;
    typedef vector<MarketScenarioArray>            MarketScenarioMatrix;

    // static methods
    static void calcABSCDOLossDistribution(DateTime currentDate, int lossSize,
                                           const MarketScenarioMatrix& scenarios,
                                           const NameParamArray& basket);

    static void calcABSCDOLossDensities(DoubleArray&      density,
                                        LossDistribution::DistributionType type,
                                        const MarketScenarioMatrix& scenarios,
                                        int nbName,
                                        int lossSize);

    static void calcABSCDOTrancheLossFast(DateTime currentDate, double k1, double k2,
                                      double lossUnit, double port,
                                      const MarketScenarioMatrix& scenarios,
                                      const NameParamArray& basket,
                                      double& loss, double& lossTop);

    static void calcABSCDOTrancheLossSaddle(DateTime currentDate, double k1, double k2,
                                      double lossUnit, double port,
                                      int lossMarketFactors, const double* lossMF, const double* lossWeight,
                                      const MarketScenarioMatrix& scenarios,
                                      const NameParamArray& basket,
                                      double& loss, double& lossTop, int numSaddlePoints);
};

typedef refCountPtr<ABSCDOConvolution> ABSCDOConvolutionSP;

DRLIB_END_NAMESPACE
#endif
