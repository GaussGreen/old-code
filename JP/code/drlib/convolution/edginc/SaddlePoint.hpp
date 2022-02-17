//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : SaddlePoint.hpp
//
//   Description : SaddlePoint Algorithm
//
//   Date        : Fed 2006
//
//
//----------------------------------------------------------------------------
#ifndef EDR_SADDLEPOINT_HPP
#define EDR_SADDLEPOINT_HPP

#include "edginc/RootFinder.hpp"
#include <vector>

DRLIB_BEGIN_NAMESPACE

class SaddlePoint
{
public:
    SaddlePoint(size_t nbMarketFactors, const double* marketFactors, const double* mfWeights,
                int nbSaddlePoints, size_t nbNames, 
                double port, double lossUnit, double k1, double k2);
    virtual ~SaddlePoint() {};

    // set per name param for current marlet factor
    void setParams(int mfIndex, int nameIndex, int numPoints, double* loss, double* probs);

    // get tranche expected loss
    double getExpLoss();

private:
    class SolverHelper : public Func1D::WtDeriv {
    public:
        SolverHelper(int mfIndex, double k, SaddlePoint* driver) : 
          mfIndex(mfIndex), k(k), driver(driver) {};
        void operator()(double x, double& f, double& df) const {
            driver->DKANDDDK(mfIndex, x, f, df);
            f -= k;
        }
    private:
        int mfIndex;
        double k;
        SaddlePoint* driver;
    };

    double K(int mfIndex, double z);
    double DK(int mfIndex, double z);
    double DDK(int mfIndex, double z);
    void   DKANDDDK(int mfIndex, double z, double& k, double& dk) const;

    void interpolateSaddlePoints();
    void findSaddlePoints(int mfIndex);
    double findSaddlePoint(int mfIndex, double k);
    bool bracket(SolverHelper& calculator, double& left, double& right, bool adjustLeft);

    double getBaseLoss(int mfIndex, double saddle, double k);
    void calcExpLoss(int mfIndex);
    void calcTotalExpLoss(int mfIndex);

    size_t nbMarketFactors;
    size_t nbNames;
    int nbPoints; // nb points to explicitly calculate
    double port; // portfolio notional in $    
    double lossUnit;
    double k1, k2; // strikes in %

    struct Scenario {
        double theta;
        double a;
        double getDiscVal(double z) const { return theta*exp(a*z); }
    };

    struct MarketFactor {
        double factor;
        double weight;
        std::vector<std::vector<Scenario> > params;
        double totalExpLoss; // E[a|m]
        double saddlePoint1; // for k1
        double saddlePoint2;// for k2
        double I2;
        bool interp;
    };
    
    std::vector<MarketFactor> mf;
    
    friend class SolverHelper;
};

typedef refCountPtr<SaddlePoint> SaddlePointSP;

DRLIB_END_NAMESPACE
#endif
