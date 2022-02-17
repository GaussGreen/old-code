//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMRatesFactor.hpp
//
//   Description : Helper for SRM - used for holding intermediate data
//
//   Author      : Mark A Robson
//
//   Date        : 14 June 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SRMRATESFACTOR_HPP
#define EDR_SRMRATESFACTOR_HPP

#include "edginc/Maths.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/SRMConstants.hpp"

DRLIB_BEGIN_NAMESPACE

class SRMRatesHJMUtil;


// let's assume a lot of in-lining here, so keeping body inside
struct SRMRatesFactor{
//    friend class SRMRatesUtil;
//    friend class SRMRatesUtil::Model;
//    friend class SRMRatesUtil::Model1F;
//    friend class SRMRatesUtil::Model2F;
//    friend class SRMRatesUtil::Model3F;

    SRMRatesFactor():alpha(0.0), beta(0.0) {}
    double kFactor(double rBarRatio, double yearFrac) const
    {
        return (rBarRatio * exp(-beta * yearFrac));
    }

    double variance(double delt_1, double delt_2) const
    {
        return (Maths::square(alpha) * VFAC(delt_1,delt_2, 2.0 * beta));
    }

    double meanReversionIntegral(double delt_1, double delt_2) const
    {
        return (alpha * VFAC(delt_1,delt_2, beta));
    }

    double covariance(double delt_1, double delt_2) const
    {
        /* mean-reversion term integral */
        double mr  =  alpha * VFAC(delt_1,delt_2, beta);
        return (alpha * mr);
    }

    /*****	ExpDecay  ********************************************/
    /*
    *   Exponential decay function.
    */
    double expDecay(double t) const
    {
        /* If beta=0 or t=0, ExpDecay(a,t)=1 */
        if (fabs (beta * t) < 1E-7) { // tidy up - to do
            return 1.0;
        }
        return ((1. - exp (- beta * t)) / beta / t);
    }/* ExpDecay */

    //// the integral of exp(-b*t) between t1 and t2
    static double VFAC(double t1, double t2, double b)
    {
        // not sure why the fabs is here?.
        // moreover why not SRMUtil::GFAC(t1, t2-t1, b)?
        return (b<=SRMConstants::SRM_TINY ?
                 fabs(t2-t1): ((exp(-b*t2) - exp(-b*t1))/b));
    }

    // For efficiency, leaving these public.
    double         alpha;
    double         beta;
    vector<double> GfactorArrayIR;
};


class SRMRatesFactorModel :     public virtual VirtualDestructorBase
{
    friend class SRMRatesHJMUtil;
public:

    SRMRatesFactorModel(int numFactors):
        factors(numFactors),
        rho((numFactors * (numFactors - 1))/2)
        {}

    virtual ~SRMRatesFactorModel(){}

    void kFactor(
        double             rBarDateFrom,
        double             rBarDateTo,
        double             del_t,
        vector<double>&    k) const;  // (0)

    /** Returns the number of factors */
    int numFactors() const{ return factors.size(); }

    /** Returns 'mean reversion integral' as determined by model parameters */
    void meanReversionIntegral(double delt_1, double delt_2, vector<double>& mr) const;

    /** Returns 'variance' as determined by model parameters */
    virtual double variance(double delt_1, double delt_2) const = 0;

    /** qFwdSpotVolSq = q * Fwd * SpotVol^2 */
    virtual double covariance(
        double                delt_1,
        double                delt_2,
        double                qFwdSpotVolSq,
        const vector<double>& gFac) const = 0;

    virtual void xiFactor(
        double&               xi,
        double&               der,
        double&               maxXi,
        const vector<double>& kfac,
        const vector<double>& xiByFactor,
        double                rbar_ratio) const = 0;

    //// need a comment
    virtual double volWeight() const = 0;

protected:
    vector<SRMRatesFactor> factors;
    vector<double> rho; // correlations
};
typedef smartPtr<SRMRatesFactorModel> SRMRatesFactorModelSP;

class SRMRatesFactorModel1F: public SRMRatesFactorModel {
public:
    SRMRatesFactorModel1F(): SRMRatesFactorModel(1){}

    virtual double variance(double delt_1, double delt_2) const;

    /** qFwdSpotVolSq = q * Fwd * SpotVol^2 */
    virtual double covariance(
        double                delt_1,
        double                delt_2,
        double                qFwdSpotVolSq,
        const vector<double>& gFac) const;

        virtual void xiFactor(
            double&               xi,
            double&               der,
            double&               maxXi,
            const vector<double>& kfac,
            const vector<double>& xiByFactor,
            double                rbar_ratio) const;

        virtual double volWeight() const;
};

class SRMRatesFactorModel2F: public SRMRatesFactorModel {
public:
    SRMRatesFactorModel2F(): SRMRatesFactorModel(2) {}

    virtual double variance(double delt_1, double delt_2) const;

    virtual double covariance(
        double                delt_1,
        double                delt_2,
        double                qFwdSpotVolSq,
        const vector<double>& gFac) const;

    virtual void xiFactor(
        double&               xi,
        double&               der,
        double&               maxXi,
        const vector<double>& kfac,
        const vector<double>& xiByFactor,
        double                rbar_ratio) const;

    virtual double volWeight() const;
};


class SRMRatesFactorModel3F: public SRMRatesFactorModel {
public:
    SRMRatesFactorModel3F(): SRMRatesFactorModel(3){}

    virtual double variance(double delt_1, double delt_2) const;

    virtual double covariance(
        double                delt_1,
        double                delt_2,
        double                qFwdSpotVolSq,
        const vector<double>& gFac) const;

    virtual void xiFactor(
        double&               xi,
        double&               der,
        double&               maxXi,
        const vector<double>& kfac,
        const vector<double>& xiByFactor,
        double                rbar_ratio) const;

    virtual double volWeight() const;
};

DRLIB_END_NAMESPACE
#endif
