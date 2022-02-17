//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : SaddlePoint.cpp
//
//   Description : SaddlePoint Algorithm
//
//   Date        : Feb 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Interpolator.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/SaddlePoint.hpp"
#include "edginc/Spline.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE
#define TINY 1e-12
#define REALLY_TINY 1e-64
#define REALLY_BIG 1e+64

////////////////////////////////////////////////////////////////////////////////

SaddlePoint::SaddlePoint(size_t nbMarketFactors, const double* marketFactors, const double* mfWeights,
                         int nbSaddlePoints, size_t nbNames, 
                         double port, double lossUnit, double k1, double k2)
  : nbMarketFactors(nbMarketFactors), nbNames(nbNames), nbPoints(nbPoints), port(port), 
    lossUnit(lossUnit), k1(k1/port), k2(k2/port)
{
    size_t i,j;

    mf.reserve(nbMarketFactors);
    for (i = 0; i < nbMarketFactors; ++i) {
        MarketFactor tmp;
        tmp.factor = marketFactors[i];
        tmp.weight = mfWeights[i];
        tmp.interp = true;
        tmp.params.reserve(nbNames);
        std::vector<Scenario> param;
        for (j = 0; j < nbNames; ++j)
            tmp.params.push_back(param);
        mf.push_back(tmp);
    }

    // calculate interpolation points
    // this uses a naive implementation of even spacing
    int stepSize = nbMarketFactors / nbSaddlePoints;
    i = 0;
    while (i < nbMarketFactors) {
        mf[i].interp = false;
        i += stepSize;
    }
    mf[nbMarketFactors-1].interp = false; // always explicitly calc the first and last points
}

double SaddlePoint::getExpLoss()
{
    const char routine[] = "SaddlePoint::getExpLoss";
    try {
        size_t i;
    
        // calculate saddle point for all interp points
        for (i = 0; i < nbMarketFactors; ++i) {
            calcTotalExpLoss(i);
            if (!mf[i].interp)
                findSaddlePoints(i);
        }

        // interpolate all other saddle points
        interpolateSaddlePoints();

        // calc tranche expected loss for all market factors
        for (i = 0; i < nbMarketFactors; ++i)
            calcExpLoss(i);

        // sum over market factors
        double res = 0;
        for (i = 0; i < nbMarketFactors; ++i)
            res += mf[i].weight*mf[i].I2;

#ifdef DEBUG
        // xxx, for debug
        std::vector<double> finp;
        for (i = 0; i < nbMarketFactors; ++i)
            finp.push_back(mf[i].weight);
        finp.clear();
        for (i = 0; i < nbMarketFactors; ++i)
            finp.push_back(mf[i].I2);
#endif

        return res*port;
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

void SaddlePoint::setParams(int mfIndex, int nameIndex, int numPoints, double* loss, double* probs)
{
    ASSERT((size_t)mfIndex < nbMarketFactors && (size_t)nameIndex < nbNames);
    std::vector<Scenario>& param = mf[mfIndex].params[nameIndex];
    param.clear();
    param.reserve(numPoints);
    Scenario sce;
    for (int i = 0; i < numPoints; ++i) {
        sce.theta = 1 - probs[i];
        sce.a = loss[i]*lossUnit/port;
        param.push_back(sce);
    }
}

///////////////////////////////////////////////////////////////////////////////

double SaddlePoint::K(int mfIndex, double z) 
{
    register double res = 0; 
    register double nsum;
    for (size_t i = 0; i < nbNames; ++i) {
        std::vector<Scenario>& param = mf[mfIndex].params[i];
        nsum = 0;
        for (size_t j = 0; j < param.size(); ++j) {
            nsum += param[j].getDiscVal(z);
        }
        res += log(nsum);
    }
    return res/nbNames;
}

double SaddlePoint::DK(int mfIndex, double z) 
{
    register double res = 0;
    register double nsum;
    register double ndsum;
    register double discval;
    for (size_t i = 0; i < nbNames; ++i) {
        std::vector<Scenario>& param = mf[mfIndex].params[i];
        nsum = ndsum = 0;
        for (size_t j = 0; j < param.size(); ++j) {
            discval = param[j].getDiscVal(z);
            nsum += discval;
            ndsum += param[j].a * discval;
        }
        res += ndsum / nsum;
    }
    return res / nbNames;
}

double SaddlePoint::DDK(int mfIndex, double z) 
{
    register double res = 0;
    register double nsum;
    register double ndsum;
    register double nddsum;
    register double discval;
    for (size_t i = 0; i < nbNames; ++i) {
        std::vector<Scenario>& param = mf[mfIndex].params[i];
        nsum = ndsum = nddsum = 0;
        for (size_t j = 0; j < param.size(); ++j) {
            discval = param[j].getDiscVal(z);
            nsum += discval;
            ndsum += param[j].a * discval;
            nddsum += param[j].a * param[j].a * discval; 
        }
        res += (nddsum * nsum - ndsum * ndsum) / (nsum * nsum);
    }
    return res / nbNames;
}

void SaddlePoint::DKANDDDK(int mfIndex, double z, double& f, double& df) const
{
    f = 0;
    df = 0;
    register double nsum;
    register double ndsum;
    register double nddsum;
    register double discval;
    for (size_t i = 0; i < nbNames; ++i) {
        const std::vector<Scenario>& param = mf[mfIndex].params[i];
        nsum = ndsum = nddsum = 0;
        for (size_t j = 0; j < param.size(); ++j) {
            discval = param[j].getDiscVal(z);
            nsum += discval;
            ndsum += param[j].a * discval;
            nddsum += param[j].a * param[j].a * discval; 
        }
        f += ndsum/nsum;
        df += (nddsum * nsum - ndsum * ndsum) / (nsum * nsum);
    }
    f /= nbNames;
    df /= nbNames;
}

void SaddlePoint::calcTotalExpLoss(int mfIndex)
{
    mf[mfIndex].totalExpLoss = 0; 
    for (size_t i = 0; i < nbNames; ++i) {
        std::vector<Scenario>& param = mf[mfIndex].params[i];
        for (size_t j = 0; j < param.size(); ++j) {
             mf[mfIndex].totalExpLoss += param[j].a * param[j].theta;
        }
    }
}

bool SaddlePoint::bracket(SolverHelper& calculator, double& left, double& right, bool adjustLeft)
{
    static const int MAXSCALE = 10; // max scale 2^10
    try {
        double ddf; double lval; double rval;
        int iter = 0;
        do {
            calculator(left, lval, ddf);
            calculator(right, rval, ddf);
            if (!Maths::finite(lval) || !Maths::finite(rval))
                return false;
            if (lval * rval > 0) {
                if (adjustLeft)
                    left *= 2;
                else
                    right *= 2;
            } else {
                return true;
            }
        } while (++iter <= MAXSCALE);
        return false;
    } catch (exception& e) {
        return false;
    }
}

double SaddlePoint::findSaddlePoint(int mfIndex, double k)
{
    const char routine[] = "SaddlePoint::findSaddlePoint";

    try {
        double totalExpLoss = mf[mfIndex].totalExpLoss;

        // locate saddle points using root solver
        if (fabs(k - totalExpLoss) < TINY)
            return 0;

        SolverHelper helper(mfIndex, k/nbNames, this);

        double left = k < totalExpLoss ? -80000 : -0.5;
        double right = k < totalExpLoss ? 0.5 : 80000;
        if (!bracket(helper, left, right, k < totalExpLoss))
            return REALLY_BIG;

        double result;
        RtSafe_solve(helper, 
                     left, 
                     right, 
                     1e-4, // accuracy 
                     true, // throw on Error
                     result, 
                     RtSafe::default_MAXIT); // max iteration

        if (result == 0)
            throw ModelException("Saddle point is singular!");
        
        return result;
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

double SaddlePoint::getBaseLoss(int mfIndex, double saddle, double k)
{
    const char routine[] = "SaddlePoint::getBaseLoss";

    try {
        double totalExpLoss = mf[mfIndex].totalExpLoss;

        // adapted from PCM saddle point code
        double ahat = totalExpLoss - k;
        if (saddle == REALLY_BIG) {
            if (ahat < 0) // strike too high, [0,k] loss same as total loss
                return totalExpLoss;
            else // strike too low, [0,k] loss almost 0
                return 0;
        }

        double what = fabs(saddle) * sqrt(2 * saddle * k - 2 * nbNames * K(mfIndex, saddle))/saddle;
        double uhat = saddle * sqrt(nbNames * DDK(mfIndex, saddle));
        if (fabs(what) < REALLY_TINY || fabs(saddle*ahat*uhat) < REALLY_TINY) {
            if (ahat < 0) // strike too high, [0,k] loss same as total loss
                return totalExpLoss;
            else // strike too low, [0,k] loss almost 0
                return 0;
        }
		double norm_coeff = pow(what,-3) - 1/what + 1/(saddle*ahat*uhat);
          
        double rslt = totalExpLoss - ahat*(N1(-what) + norm_coeff*N1Density(what));
        return rslt;
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

void SaddlePoint::findSaddlePoints(int mfIndex)
{
    mf[mfIndex].saddlePoint1 = findSaddlePoint(mfIndex, k1);
    mf[mfIndex].saddlePoint2 = findSaddlePoint(mfIndex, k2);
}


void SaddlePoint::interpolateSaddlePoints()
{
    const char routine[] = "SaddlePoint::interpolateSaddlePoints";
    static const int SAFEGUARD = 100;
    try {
        size_t i, j;
        
        CubicSplineInterpECnd interpolator;
        std::vector<double> xinp;
        std::vector<double> finp;

        // collect input for k1
        j = 0;
        for (i = 0; i < nbMarketFactors; ++i) {
            if (!mf[i].interp) {
                ++j;
                xinp.push_back(mf[i].factor);
                finp.push_back(mf[i].saddlePoint1);
            }
        }
        Interpolator::InterpolantConstSP interpolant = 
            interpolator.computeInterp(&xinp[0], &finp[0], j);
        // fill interpolated value, force calc the smaller points
        for (i = 0; i < nbMarketFactors; ++i) {
            if (mf[i].interp) {
                double interpVal = interpolant->value(mf[i].factor);
                if (fabs(interpVal) < SAFEGUARD)
                    interpVal = findSaddlePoint(i, k1);
                mf[i].saddlePoint1 = interpVal;
            }
        }

#ifdef DEBUG
        // xxx, for debug
        finp.clear();
        for (i = 0; i < nbMarketFactors; ++i)
            finp.push_back(mf[i].saddlePoint1);
#endif
        // collect input for k2
        j = 0; xinp.clear(); finp.clear();
        for (i = 0; i < nbMarketFactors; ++i) {
            if (!mf[i].interp) {
                ++j;
                xinp.push_back(mf[i].factor);
                finp.push_back(mf[i].saddlePoint2);
            }
        }
        interpolant = interpolator.computeInterp(&xinp[0], &finp[0], j);
        // fill interpolated value, force calc the closest point to 0
        for (i = 0; i < nbMarketFactors; ++i) {
            if (mf[i].interp) {
                double interpVal = interpolant->value(mf[i].factor);
                if (fabs(interpVal) < SAFEGUARD)
                    interpVal = findSaddlePoint(i, k2);
                mf[i].saddlePoint2 = interpVal;
            }
        }

#ifdef DEBUG
        // xxx, for debug
        finp.clear();
        for (i = 0; i < nbMarketFactors; ++i)
            finp.push_back(mf[i].saddlePoint2);
#endif
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

void SaddlePoint::calcExpLoss(int mfIndex)
{
    const char routine[] = "SaddlePoint::calculateExpLoss";
    try {
        //xxx for debug
        if (mfIndex == 16) {
            double foo;
            foo = 2;
        }
        double totalExpLoss = mf[mfIndex].totalExpLoss;
        bool needCorrection = false;

        // tranche [0, k1]
        double bl1 = getBaseLoss(mfIndex, mf[mfIndex].saddlePoint1, k1);
        if (bl1 < 0 || bl1 > totalExpLoss) {
            needCorrection = true;
        }

        // tranche [0, k2]
        double bl2 = getBaseLoss(mfIndex, mf[mfIndex].saddlePoint2, k2);
        if (bl2 < 0 || bl2 > totalExpLoss) {
            needCorrection = true;
        }

        if (bl2 < bl1)
            needCorrection = true;
        
        if (needCorrection && mf[mfIndex].interp) {
            // try to correct interpolation error
            findSaddlePoints(mfIndex);
            bl1 = getBaseLoss(mfIndex, mf[mfIndex].saddlePoint1, k1);
            bl2 = getBaseLoss(mfIndex, mf[mfIndex].saddlePoint2, k2);
        }

        ASSERT(0 <= bl1 && bl1 <= totalExpLoss);
        ASSERT(0 <= bl2 && bl2 <= totalExpLoss);
        ASSERT(bl2 >= bl1);
        mf[mfIndex].I2 = bl2 - bl1; // result
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

DRLIB_END_NAMESPACE
