//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMRatesFactor.cpp
//
//   Description : Helper for SRM - used for holding intermediate data
//
//   Author      : Mark A Robson
//
//   Date        : 14 June 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMRatesFactor.hpp"
//#include "edginc/SRMRatesUtil.hpp"
//#include "edginc/VolProcessedBSIR.hpp"
//#include "edginc/DoubleMatrix.hpp"

DRLIB_BEGIN_NAMESPACE


void SRMRatesFactorModel::kFactor(
    double             rBarDateFrom,
    double             rBarDateTo,
    double             del_t,
    vector<double>&    k) const  // (0)
{
//    double denominator = srmUtil.rBar(dateFrom);
//    double ratio = srmUtil.rBar(dateTo) / denominator;
//    double del_t = SRMRatesUtil::yearFrac(dateFrom, dateTo);
//    if (Maths::isZero(denominator)){
//        throw ModelException("SRMRatesUtil::Model::kFactor", "Zero forward rate "
//                             " on "+dateFrom.toString());
//    }

    double ratio = rBarDateTo / rBarDateFrom;
    for (unsigned int i = 0; i < factors.size(); i++){
        k[i] = factors[i].kFactor(ratio, del_t);
    }
}

/** Returns 'mean reversion integral' as determined by model parameters */
void SRMRatesFactorModel::meanReversionIntegral(
    double delt_1, 
    double delt_2, 
    vector<double>& mr) const  //(0)
{
    mr.resize(factors.size());
    for (unsigned int i = 0; i < mr.size(); i++){
        mr[i] = factors[i].meanReversionIntegral(delt_1, delt_2);
    }
}


// implementing 1-factor variance/cov calculations

double SRMRatesFactorModel1F::variance(double delt_1, double delt_2) const
{
    return factors.front().variance(delt_1, delt_2);
}

/** qFwdSpotVolSq = q * Fwd * SpotVol^2 */
double SRMRatesFactorModel1F::covariance(
    double                delt_1,
    double                delt_2,
    double                qFwdSpotVolSq,
    const vector<double>& gFac) const
{
    double mr = factors.front().meanReversionIntegral(delt_1, delt_2);
    return (qFwdSpotVolSq * factors[0].alpha * gFac.front() * mr);
}

void SRMRatesFactorModel1F::xiFactor(
    double&               xi, 
    double&               der, 
    double&               maxXi, 
    const vector<double>& kfac, 
    const vector<double>& xiByFactor, 
    double                rbar_ratio) const
{
    /* transfer values to output */
    xi =  xiByFactor[0]/factors[0].alpha; /* xi Factor equals xi1 in the
                                             1 factor case */
    der = kfac[0]; /* the derivative of xi with respect to tau is
                      the kfactor  */
    maxXi = kfac[0]; /* Deterministic part of f(t,tau) vol squared */
}

double SRMRatesFactorModel1F::volWeight() const
{
    return factors[0].alpha;
}


// implementing 2-factor variance/cov calculations
    
double SRMRatesFactorModel2F::variance(double delt_1, double delt_2) const
{
    double var = factors[0].variance(delt_1, delt_2) +
        factors[1].variance(delt_1, delt_2);
    var += 2.0 * factors[0].alpha * factors[1].alpha * rho[0] *
        SRMRatesFactor::VFAC(delt_1, delt_2, factors[0].beta + factors[1].beta);
    return var;
}

double SRMRatesFactorModel2F::covariance(
    double                delt_1,
    double                delt_2,
    double                qFwdSpotVolSq,
    const vector<double>& gFac) const
{
    double mr1 = factors[0].meanReversionIntegral(delt_1, delt_2);
    double mr2 = factors[1].meanReversionIntegral(delt_1, delt_2);
    double alpha1 = factors[0].alpha;
    double alpha2 = factors[1].alpha;
    return (qFwdSpotVolSq * (
                gFac[0] * alpha1 * mr1 + gFac[0] * alpha1 * mr2 * rho[0] +
                gFac[1] * alpha2 * mr2 + gFac[1] * alpha2 * mr1 * rho[0]));
}

void SRMRatesFactorModel2F::xiFactor(
    double&               xi, 
    double&               der, 
    double&               maxXi, 
    const vector<double>& kfac, 
    const vector<double>& xiByFactor, 
    double                rbar_ratio) const
{
    double alpha1 = factors[0].alpha;
    double alpha2 = factors[1].alpha;
    /* calculate the Xi factor   */
    xi = (kfac[0] * alpha1 * xiByFactor[0] + 
          rho[0] * kfac[0] * alpha1 *  xiByFactor[1]
          + kfac[1] * alpha2 * rho[0] * xiByFactor[0] + 
          kfac[1] * alpha2 *  xiByFactor[1]) / rbar_ratio;

    /* calculate the derivative of xi with respect to tau */
    der =(-factors[0].beta * kfac[0] * alpha1 * 
          (xiByFactor[0] + rho[0] * xiByFactor[1])
          -factors[1].beta * kfac[1] * alpha2 * 
          (xiByFactor[1] + rho[0] * xiByFactor[0])
          + kfac[0] * alpha1 * 
          (alpha1 * kfac[0] + rho[0] * alpha2 * kfac[1])
          + kfac[1] * alpha2 * 
          (alpha2 * kfac[1] + rho[0] * alpha1 * kfac[0])) / rbar_ratio;
    
    /* Deterministic part of f(t,tau) vol squared         */
    maxXi = (Maths::square(alpha1 * kfac[0]) + 
             Maths::square(alpha2 * kfac[1]) 
             + 2 * alpha1 * alpha2 *  rho[0] * kfac[0] * kfac[1]) /
        rbar_ratio;
}

double SRMRatesFactorModel2F::volWeight() const
{
    double alpha0 = factors[0].alpha;
    double alpha1 = factors[1].alpha;
    double weight =
        alpha0*alpha0 + alpha1*alpha1 + alpha0*alpha1 * 2.0 * rho[0];
    return sqrt(weight);
}


// implementing 3-factor variance/cov calculations

double SRMRatesFactorModel3F::variance(double delt_1, double delt_2) const
{
    double var = factors[0].variance(delt_1, delt_2) +
        factors[1].variance(delt_1, delt_2) +
        factors[2].variance(delt_1, delt_2);
    double alpha1 = factors[0].alpha;
    double alpha2 = factors[1].alpha;
    double alpha3 = factors[2].alpha;
    double beta1 = factors[0].beta;
    double beta2 = factors[1].beta;
    double beta3 = factors[2].beta;
    var += 
        2.0 * alpha1 * alpha2 * rho[0] * 
        SRMRatesFactor::VFAC(delt_1, delt_2, beta1 + beta2) + 
        2.0 * alpha1 * alpha3 * rho[1] *
        SRMRatesFactor::VFAC(delt_1, delt_2, beta1 + beta3) +
        2.0 * alpha2 * alpha3 * rho[2] *
        SRMRatesFactor::VFAC(delt_1, delt_2, beta2 + beta3);
    return var;
}

double SRMRatesFactorModel3F::covariance(
    double                delt_1,
    double                delt_2,
    double                qFwdSpotVolSq,
    const vector<double>& gFac) const
{
    double mr1 = factors[0].meanReversionIntegral(delt_1, delt_2);
    double mr2 = factors[1].meanReversionIntegral(delt_1, delt_2);
    double mr3 = factors[2].meanReversionIntegral(delt_1, delt_2);
    double alpha1 = factors[0].alpha;
    double alpha2 = factors[1].alpha;
    double alpha3 = factors[2].alpha;
    return (qFwdSpotVolSq * (
            alpha1 * gFac[0] * (mr1 + rho[0] * mr2 + rho[1] * mr3) +
            alpha2 * gFac[1] * (rho[0] * mr1 + mr2 + rho[2] * mr3) +
            alpha3 * gFac[2] * (rho[1] * mr1 + rho[2] * mr2 + mr3)));
}

void SRMRatesFactorModel3F::xiFactor(
    double&               xi, 
    double&               der, 
    double&               maxXi, 
    const vector<double>& kfac, 
    const vector<double>& xiByFactor, 
    double                rbar_ratio) const
{
    // for ease
    double alpha1 = factors[0].alpha;
    double alpha2 = factors[1].alpha;
    double alpha3 = factors[2].alpha;
    xi = ( kfac[0] * alpha1 * (xiByFactor[0] + rho[0] * xiByFactor[1] + 
                               rho[1] * xiByFactor[2])
           + kfac[1] * alpha2 * (rho[0] * xiByFactor[0] + xiByFactor[1] +
                                 rho[2] * xiByFactor[2])
           + kfac[2] * alpha3 * (rho[1] * xiByFactor[0] + 
                                 rho[2] * xiByFactor[1] + xiByFactor[2]) ) /
        rbar_ratio;
    
    /* calculate the derivative of xi with respect to tau */
    der =( -factors[0].beta * kfac[0] * alpha1 *
           (xiByFactor[0] + rho[0] * xiByFactor[1]+ rho[1] * xiByFactor[2])
           -factors[1].beta * kfac[1] * alpha2 * 
           (rho[0] * xiByFactor[0] + xiByFactor[1]+ rho[2] * xiByFactor[2])
           -factors[2].beta * kfac[2] * alpha3 * 
           (rho[1] * xiByFactor[0] +rho[2] * xiByFactor[1] + xiByFactor[2])
           + kfac[0] * alpha1 * ( kfac[0] * alpha1 + rho[0] * kfac[1] * 
                                  alpha2 + rho[1] * kfac[2] * alpha3)
           + kfac[1] * alpha2 * ( rho[0] * kfac[0] * alpha1 + kfac[1] *
                                  alpha2 + rho[2] * kfac[2] * alpha3)  
           + kfac[2] * alpha3 * ( rho[1] * kfac[0] * alpha1 + rho[2] *
                                  kfac[1] * alpha2 + kfac[2] * alpha3)) /
        rbar_ratio;	
    
    
    /* Deterministic part of f(t,tau) vol squared */		
    maxXi = (Maths::square(alpha1 * kfac[0]) + 
             Maths::square(alpha2 * kfac[1])+Maths::square(alpha3 * kfac[2])
             + 2 * alpha1 * alpha2 *  rho[0] * kfac[0] * kfac[1]
             + 2 * alpha1 * alpha3 *  rho[1] * kfac[0] * kfac[2]
             + 2 * alpha2 * alpha3 *  rho[2] * kfac[1] * kfac[2]) / 
        rbar_ratio;
}

double SRMRatesFactorModel3F::volWeight() const
{
    double alpha0 = factors[0].alpha;
    double alpha1 = factors[1].alpha;
    double alpha2 = factors[2].alpha;
    double weight = alpha0*alpha0 + alpha1*alpha1 + alpha2*alpha2 +
        alpha0*alpha1 * 2.0 * rho[0] + alpha0*alpha2 * 2.0 * rho[1] +
        alpha1*alpha2 * 2.0 * rho[2];
    return sqrt(weight);
}

DRLIB_END_NAMESPACE

