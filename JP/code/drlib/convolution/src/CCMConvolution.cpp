//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CCMConvolution.cpp
//
//   Description : Convolution Algorithm
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CCMConvolution.hpp"
#include "edginc/TrancheLossCalculator.hpp"
#include "edginc/CCMCalibration.hpp"
#include "edginc/CCMSkew.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Addin.hpp"
#include "edginc/LossDistribution.hpp"
#include <algorithm>
#include "edginc/Hashtable.hpp" // for hash_string
#include "edginc/Nrfns.hpp"
//#include <iostream>
//#include <fstream>

DRLIB_BEGIN_NAMESPACE
#define TINY 1e-12
#define REALLY_TINY 1e-64
#define TOWER_PROBA_EPSI      1e-12
#define SKEW_FCUMUL_NB_POINT  100L

/* ---------------------------------------------------------------------------
 * Integ Method
 */
#define CCM_INTEG_LOW       -7.5
#define CCM_INTEG_UP         7.5
#define CCM_INTEG_NB         101

typedef CCMConvolution::LgdParam                LgdParam;                // for ease
typedef CCMConvolution::LgdParamSP              LgdParamSP;              // for ease
typedef CCMConvolution::LgdParamArray           LgdParamArray;           // for ease
typedef CCMConvolution::NameParam               NameParam;               // for ease
typedef CCMConvolution::NameParamSP             NameParamSP;             // for ease
typedef CCMConvolution::NameParamArray          NameParamArray;          // for ease
typedef CCMConvolution::RFLNameParam               RFLNameParam;               // for ease
typedef CCMConvolution::RFLNameParamSP             RFLNameParamSP;             // for ease
typedef CCMConvolution::RFLNameParamArray          RFLNameParamArray;          // for ease
typedef CCMConvolution::MarketFactor            MarketFactor;            // for ease
typedef CCMConvolution::MarketFactorSP          MarketFactorSP;          // for ease
typedef CCMConvolution::MarketFactorArray       MarketFactorArray;       // for ease

//CCMConvolution & inner class type registration
CClassConstSP const CCMConvolution::TYPE = CClass::registerClassLoadMethod(
    "CCMConvolution", typeid(CCMConvolution), load);

DEFINE_TEMPLATE_TYPE(CCMConvolutionArray);

CClassConstSP const CCMConvolution::LgdParam::TYPE = CClass::registerClassLoadMethod(
    "CCMConvolution::LgdParam", typeid(CCMConvolution::LgdParam), load);

DEFINE_TEMPLATE_TYPE_WITH_NAME("CCMConvolution::LgdParamArray", LgdParamArray);

CClassConstSP const CCMConvolution::NameParam::TYPE = CClass::registerClassLoadMethod(
    "CCMConvolution::NameParam", typeid(CCMConvolution::NameParam), load);

DEFINE_TEMPLATE_TYPE_WITH_NAME("CCMConvolution::NameParamArray", NameParamArray);

CClassConstSP const CCMConvolution::RFLNameParam::TYPE = CClass::registerClassLoadMethod(
    "CCMConvolution::RFLNameParam", typeid(CCMConvolution::RFLNameParam), load);

template<> CClassConstSP const RFLNameParamArray::TYPE = CClass::registerClassLoadMethod(
    "CCMConvolution::RFLNameParamArray", typeid(RFLNameParamArray), load);

CClassConstSP const CCMConvolution::MarketFactor::TYPE = CClass::registerClassLoadMethod(
    "CCMConvolution::MarketFactor", typeid(CCMConvolution::MarketFactor), load);

DEFINE_TEMPLATE_TYPE_WITH_NAME("CCMConvolution::MarketFactorArray", MarketFactorArray);

//Convolution addin function type registration
CClassConstSP const ConvolutionAddin::TYPE = CClass::registerClassLoadMethod(
    "ConvolutionAddin", typeid(ConvolutionAddin), load);

CClassConstSP const DeReConvoluteAddin::TYPE = CClass::registerClassLoadMethod(
    "DeReConvoluteAddin", typeid(DeReConvoluteAddin), load);

CClassConstSP const DensityAddin::TYPE = CClass::registerClassLoadMethod(
    "DensityAddin", typeid(DensityAddin), load);

CClassConstSP const DensityAddin::Results::TYPE = CClass::registerClassLoadMethod(
    "DensityAddin::Results", typeid(DensityAddin::Results), load);

//// constructor that zeroes everything
LgdParam::LgdParam():
    CObject(TYPE),
    LGD1(0), LGD2(0),
    LGD1d(0), LGD2d(0),
    T1(0), W1ind(0),
    isT1Calibrated(false), lc(0)
{
}
//// constructor that zeroes everything
NameParam::NameParam():
    CObject(TYPE),
    survival(0.0), beta(0.0), qM(0.0), pdep(0.0), indep(0.0), 
    cataRecFactor(0.0), lambda(0.0), betaRec(0.0), ntl(0.0), R(0.0),
    lgdNotional(0), lgdFloor(0), lgdCap(0)
{
}

//// constructor that zeroes everything
NameParam::NameParam(const CClassConstSP clazz):
    CObject(clazz),
    survival(0.0), beta(0.0), qM(0.0), pdep(0.0), indep(0.0), 
    cataRecFactor(0.0), lambda(0.0), betaRec(0.0), ntl(0.0), R(0.0),
    lgdNotional(0), lgdFloor(0), lgdCap(0)
{
}

//// constructor that zeroes everything
RFLNameParam::RFLNameParam():
    NameParam(TYPE), rflParam(0)
{
}

void NameParam::threshold(double pgauss, double &tgauss) const
{
    tgauss = CCMSkew::fInv(pgauss, beta, qM, SKEW_FCUMUL_NB_POINT);
}

void RFLNameParam::threshold(double pgauss, double &tgauss) const
{
    rflParam->threshold(pgauss, tgauss, beta);
}

//calculate probabilities
void NameParam::probas(double& pind, double& pgauss) const
{
    if (fabs(survival-1.) < 3e-16) {
        pind   = 1.;
        pgauss = 1.;
    } else {
        double ci = log(pdep)/log(survival); 
        pind   = pow(survival, (1.-ci)*indep);
        pgauss = survival / (pdep*pind);
    }
}

//calculate probabilities
void NameParam::probas(double& pind, double& pgauss, double& tgauss) const
{
    /* compute probas */
    probas(pind, pgauss);
    /* initialize threshold */
    threshold(pgauss, tgauss);
}

void NameParam::condSurvivalProba(  double pind,
                                    double tgauss,
                                    double M,
                                    double&      pm) const
{
    pm = CCMSkew::condSurvivalProba(pind,
                                    tgauss,
                                    beta,
                                    qM,
                                    M);
}

void RFLNameParam::condSurvivalProba(   double pind,
                                        double tgauss,
                                        double M,
                                        double&      pm) const
{
    rflParam->condSurvivalProba(    pind,
                                    tgauss,
                                    beta,
                                    M,
                                    pm);
}


//calculate probability conditional on a market factor
void NameParam::probCond(double      pind,
                         double      tgauss,
                         const LgdParam& recInfo,
                         double      M,
                         double&     pm,
                         double&     weight1) const
{
    condSurvivalProba(pind,
                      tgauss,
                      M,
                      pm);

    weight1 = CCMCalibration::recoveryWeightCalc(recInfo.W1ind,
                                                 recInfo.T1,
                                                 recInfo.isT1Calibrated,
                                                 indep,
                                                 betaRec,
                                                 qM,
                                                 M);
}

bool NameParam::operator()(const NameParamSP& e1, const NameParamSP& e2)
{
    return (e2->pdep < e1->pdep); // must act like a less than operation
}

bool NameParam::equals(const NameParamSP rhs) const
{
    if (rhs->getClass() != TYPE)
    {
        return false;
    }
    bool equal = true;

    //is this the same as rhs ?
    equal = equal && (nameId        == rhs->nameId);
    equal = equal && (survival      == rhs->survival);
    equal = equal && (beta          == rhs->beta);
    equal = equal && (qM            == rhs->qM);
    equal = equal && (pdep          == rhs->pdep);
    equal = equal && (indep         == rhs->indep);
    equal = equal && (cataRecFactor == rhs->cataRecFactor);
    equal = equal && (lambda        == rhs->lambda);
    equal = equal && (betaRec       == rhs->betaRec);
    equal = equal && (ntl           == rhs->ntl);
    equal = equal && (R             == rhs->R);
    equal = equal && (lgdNotional   == rhs->lgdNotional);
    equal = equal && (lgdFloor      == rhs->lgdFloor);
    equal = equal && (lgdCap        == rhs->lgdCap);

    return equal;
}

bool RFLNameParam::equals(const NameParamSP rhs) const
{
    if (rhs->getClass() != TYPE)
    {
        return false;
    }
    else
    {
        RFLNameParam* rflRhs = DYNAMIC_CAST(RFLNameParam,rhs.get());
        return NameParam::equals(rhs) && rflParam->equalTo(rflRhs->rflParam.get());
    }
}

int NameParam::hashCode() const
{
    int hCode  = 0;
    hCode ^= CDouble::hashCode(indep);
    hCode ^= CDouble::hashCode(beta);
    hCode ^= CDouble::hashCode(survival);
    hCode ^= CDouble::hashCode(pdep);
    hCode ^= CDouble::hashCode(betaRec);
    hCode ^= CDouble::hashCode(cataRecFactor);
    hCode ^= CDouble::hashCode(lambda);
    hCode ^= CDouble::hashCode(lgdNotional);
    hCode ^= CDouble::hashCode(lgdFloor);
    hCode ^= CDouble::hashCode(lgdCap);
    hCode ^= CDouble::hashCode(ntl);
    hCode ^= CDouble::hashCode(qM);
    hCode ^= CDouble::hashCode(R);
    hCode ^= hash_string(nameId);
    return hCode;
}

int RFLNameParam::hashCode() const
{
    int hCode  = NameParam::hashCode();
    //hCode ^= rflParam->hashCode();
    return hCode;
}



bool NameParam::equalsUntweakable(const NameParamSP rhs) const
{
    bool equal = true;

    //are the fields that we cant consider for quick greeks the same ?
    // - all but survival, beta and recovery
    equal = equal && (nameId        == rhs->nameId);
    equal = equal && (qM            == rhs->qM);
    equal = equal && (pdep          == rhs->pdep);
    equal = equal && (indep         == rhs->indep);
    equal = equal && (cataRecFactor == rhs->cataRecFactor);
    equal = equal && (lambda        == rhs->lambda);
    equal = equal && (betaRec       == rhs->betaRec);
    equal = equal && (ntl           == rhs->ntl);
    equal = equal && (lgdNotional   == rhs->lgdNotional);
    equal = equal && (lgdFloor      == rhs->lgdFloor);
    equal = equal && (lgdCap        == rhs->lgdCap);

    return equal;
}

bool RFLNameParam::equalsUntweakable(const NameParamSP rhs) const
{
    throw ModelException("RFLNameParam::equalsUntweakable","Not implemented"); 
}

//check for differences in the two states
//returns
//  -1 : more than 1 difference
//   0 : no difference
//   n : the index+1 of the NameParam meeting the condition
int NameParam::singleTweakDiff(const NameParamArray& state1,
                               const NameParamArray& state2)
{
    int nmIdx = 0;

    //invalid for different lengths of array
    if (state1.size() == state2.size())
    {
        //look for first difference
        for (int i=0; i<state1.size() && nmIdx==0; i++)
        {
            //are they equal ?
            if (!(state1[i]->equals(state2[i])))
            {
                //if not, do they differ only by
                //tweakable fields
                if (state1[i]->equalsUntweakable(state2[i]))
                {
                    //if so, keep the index 
                    nmIdx = i+1;
                }
                else
                {
                    //a difference we cannot deal with
                    nmIdx = -1;
                }
            }
        }

        //if we found a contender
        if (nmIdx > 0)
        {
            //check remaining names for no differences
            for (int i=nmIdx; i<state1.size() && nmIdx !=-1; i++)
            {
                //are they equal ?
                if (!(state1[i]->equals(state2[i])))
                {
                    //invalidates
                    nmIdx = -1;
                }
            }
        }
        else
        {
            //no contender - states are identical
        }
    }
    else
    {
        nmIdx = -1;
    }

    return nmIdx;
}

MarketFactor::MarketFactor() :
    CObject(TYPE),
    weight(0)
{}

MarketFactor::MarketFactor(const double m,
                           const double w) : CObject(TYPE)
{
    //store critical information
    M = m;
    weight = w;

    distribs.resize(0);
}

void MarketFactor::addDistribution(const int nbLoss, const int maxShortIdx)
{
    const char routine[] = "MarketFactor::addDistribution";
    try
    {
        //add storage for the next catastrophic scenario
        int numScen = distribs.size() + 1;

        //add to the arrays
        ////distribs.resize(numScen);

        if (numScen == 1)
        {
            //add the distribution
            LossDistributionSP ld = LossDistributionSP(new LossDistribution(nbLoss,maxShortIdx));
            distribs.push_back(ld);
        }
        else
        {
            //copy the previous distribution
            //to determine the starting state for the new case
            //relying on default copy constructor
            LossDistributionSP ld = LossDistributionSP(new LossDistribution(*(distribs[numScen-2].get())));
            distribs.push_back(ld);
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

LossDistributionSP MarketFactor::getLastDistribution()
{
    return distribs.back();
}

void MarketFactor::deconvolute(const int          nmIdx,
                               const NameParam&   nameInfo,
                               const LgdParam&    recInfo,
                               const double       name_pind,
                               const double       name_tgauss)
{
    //calculate params conditional on M
    double name_pm;
    double name_weight1;

    nameInfo.probCond(name_pind,
                      name_tgauss,
                      recInfo,
                      M,
                      name_pm,
                      name_weight1);

    int numScen = distribs.size();

    //deconvolution against all loss distributions
    for (int i=0; i<numScen; i++)
    {
        //only if this name has contributed to the loss distribution
        if (nmIdx+1 <= distribs[i]->getLastRunName())
        {
            distribs[i]->deconvolute(nmIdx,
                                     name_pm,
                                     name_weight1,
                                     recInfo.LGD1,
                                     recInfo.LGD2);
        }
    }
}

//convolute a new version of a name
void MarketFactor::convolute(const int          nmIdx,
                             const NameParam&   nameInfo,
                             const LgdParam&    recInfo,
                             const double       name_pind,
                             const double       name_tgauss)
{
    //calculate params conditional on M
    double name_pm;
    double name_weight1;

    nameInfo.probCond(name_pind,
                      name_tgauss,
                      recInfo,
                      M,
                      name_pm,
                      name_weight1);

    int numScen = distribs.size();

    //deconvolution against all loss distributions
    for (int i=0; i<numScen; i++)
    {
        //only if this name has contributed to the loss distribution
        if (nmIdx+1 <= distribs[i]->getLastRunName())
        {
            distribs[i]->convolute(nmIdx,
                                   name_pm,
                                   name_weight1,
                                   recInfo.LGD1,
                                   recInfo.LGD2);
        }
    }
}


void MarketFactor::accumulateDistribution(int scenario, DoubleArray& dist)
{
    const char routine[] = "MarketFactor::accumulateDistribution(dist)";
    try
    {
        //add the distribution for the required scenario
        const DoubleArray& lossDist = distribs[scenario]->getLossDistribution();

        for (int i=0; i<lossDist.size(); i++)
        {
            dist[i] += lossDist[i] * weight;
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

void MarketFactor::accumulateDistribution(int scenario, DoubleArray& dist, double weight2, DoubleArray& distCond)
{
    const char routine[] = "MarketFactor::accumulateDistribution(dist,distCond)";
    try
    {
        //add the distribution for the required scenario
        const DoubleArray& lossDist = distribs[scenario]->getLossDistribution();

        double lossI;

        for (int i=0; i<lossDist.size(); i++)
        {
            //save an array access
            lossI = lossDist[i];

            dist[i]     += lossI * weight;
            distCond[i] += lossI * weight2;
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

double MarketFactor::getFactor()
{
    return M;
}

double MarketFactor::getWeight()
{
    return weight;
}

double MarketFactor::getWeight2(const double pindCpty,
                                const double tgaussCpty,
                                const double cptyBeta,
                                const double cptyQM)
{
    double pmCpty = CCMSkew::condSurvivalProba(pindCpty,
                                               tgaussCpty,
                                               cptyBeta,
                                               cptyQM,
                                               M);
    double weight2 = pmCpty * weight;

    return weight2;
}

int MarketFactor::getLastRunName()
{
    LossDistributionSP distrib = getLastDistribution();
    return distrib->getLastRunName();
}

/* Helper function that calculates the two integer loss points that replace the 
 * non integer loss input in convolution
 * Uses probability of Default, not Survival
 * solves: { q1+q2 = q
 *         { q1*l1+q2*l2 = q*l
 */
static void solveSystem( 
    double    l,         /* (I) non int LGD value        */
    double    q,         /* (I) probability of default   */
    int&      l1,        /* (O) lower loss               */
    double&   q1,        /* (O) proba associated with l1 */
    int&      l2,        /* (O) upper loss               */
    double&   q2)        /* (O) proba associated with l2 */
{
    ASSERT(fabs(q-.5) <= .5);

    l1 = (long)floor(l);
    l2 = (long)ceil(l); 
    ASSERT(abs(l1-l2) <=1);
    if (fabs(l - l1) < TINY) {
        l2 = l1;
    } else if (fabs(l - l2) < TINY) {
        l1 = l2;
    }
    
    /* Note that we have l1=l2 if l is integer */
    if (l1 == l2) {
        q1 = q;
        q2 = 0.;
    } else {
        q1 = q * (l - l2)/(l1 - l2);
        q2 = q - (q1);
    }
    ASSERT(q1>=0. && q1<=1.);
    ASSERT(fabs(q1 * l1 + q2 * l2 - q*l)<1e-10);
}


/* ---------------------------------------------------------------------------
 * Integrate the loss density conditional on M
 */

/*
 * calculate loss density conditional for indep-gauss copula (i.e., 
 * conditional on a certain level of dependent default)
 * handles a cpty
 */
void CCMConvolution::gaussIndepLossDensityCalc(const int nbName,       /* (I) number of names (excl. cpty)                          */
                                               const DoubleArray& pind,
                                               const DoubleArray& tgauss,
                                               const LgdParamArray& recInfo)
{
    const char routine[] = "CCMConvolution::gaussIndepLossDensityCalc";
    try{
        /* 2. Integrate convolution of losses */
        int numInteg = integrands.size();

        for (int i=0; i<numInteg; ++i)
        {
            MarketFactorSP mf = integrands[i];

            /* convolution of m+1 to nth name onto convolution of m names*/
            // cf convolution from ccm2
            int lastRunName = mf->getLastRunName();

            //required to bypass MarketFactor::convolute
            LossDistributionSP distrib = mf->getLastDistribution();
            
#ifdef PRINT_CONDITIONAL_PROBABILITIES
            // WARNING - WARNING - WARNING - WARNING WARNING - WARNING
            // print debug information if this flag is defined
            // this is just for convenience
            FILE *fp = fopen("C:\\debug.txt","a");
            fprintf(fp, "%ld\n",i);
#endif
            for (int k = lastRunName; k<nbName; k++)
            {
                //bypass MarketFactor::convolute for efficiency
                //mf->convolute(k, recInfo[k]->LGD1, recInfo[k]->LGD2);
                int    l1      = recInfo[k]->LGD1;
                int    l2      = recInfo[k]->LGD2;
                double pm;
                double weight1;
                //get pm & weight1 conditional on M
                orderedBasket[k]->probCond(pind[k],tgauss[k],*recInfo[k],mf->getFactor(),pm,weight1);
#ifdef PRINT_CONDITIONAL_PROBABILITIES
                fprintf(fp, "%.10lf\t",pm);
#endif
                distrib->convolute(k, pm, weight1, l1, l2);
            }
#ifdef PRINT_CONDITIONAL_PROBABILITIES
            fprintf(fp, "\n");
            fclose(fp);
#endif
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}
    
/**
 * Convolution algorithm.
 * Takes the copula and the raw loss parameters as input
 * Calibrates the correlated recovery model parameters as a first step.
 * For each level of loss caused by the dependence copula, do the 
 * independence/gaussian convolution conditionnaly on each value of the 
 * market variable M.
 *
 * To handle recovery correlation with M, the convolution will use two 
 * loss levels for each default: LGD1[] and LGD2[], which are the same 
 * for any M, but with weights that depend on M. 
 * These weight are a combination of the independence copula weights 
 * and the gaussian copula weights, with the approximation that the repar-
 * tition of defaults between the two copulae is independent of M.
 */
CCMConvolution::CCMConvolution() : CObject(TYPE), nbName(0), lossSize(0), maxShortIdx(0), maxLongIdx(0)
{
}

CCMConvolution::CCMConvolution(
    const bool            useFastMethod,
    const NameParamArray& basketInfo,    /* (I) name loss amt info[nbName]              */
    int                   maxLongIndex,  /* (I) density long size                       */
    int                   maxShortIndex) /* (I) density shorts size                     */
    : CObject(TYPE), expLoss(0), expLossCond(0)
{
    const char routine[] = "CCMConvolution::CCMConvolution";

    useFast = useFastMethod;
    nbName    = basketInfo.size(); // for ease
    maxLongIdx = maxLongIndex;
    maxShortIdx = maxShortIndex;
    lossSize  = maxLongIdx + maxShortIdx + 1L;

//#define DUMP_INPUTS
#ifdef DUMP_INPUTS
    //create arrays to capture inputs
    //for creating unit test cases

    StringArray nameId(nbName);
    DoubleArray survival(nbName);
    DoubleArray beta(nbName);
    DoubleArray qM(nbName);
    DoubleArray pdep(nbName);
    DoubleArray indep(nbName);
    DoubleArray cataRecFactor(nbName);
    DoubleArray lambda(nbName);
    DoubleArray betaRec(nbName);
    DoubleArray ntl(nbName);
    DoubleArray R(nbName);
    DoubleArray lgdNotional(nbName);
    DoubleArray lgdFloor(nbName);
    DoubleArray lgdCap(nbName);

    for (int ii=0; ii<nbName; ii++)
    {
        nameId[ii]        = basketInfo[ii]->nameId;
        survival[ii]      = basketInfo[ii]->survival;
        beta[ii]          = basketInfo[ii]->beta;
        qM[ii]            = basketInfo[ii]->qM;
        pdep[ii]          = basketInfo[ii]->pdep;
        indep[ii]         = basketInfo[ii]->indep;
        cataRecFactor[ii] = basketInfo[ii]->cataRecFactor;
        lambda[ii]        = basketInfo[ii]->lambda;
        betaRec[ii]       = basketInfo[ii]->betaRec;
        ntl[ii]           = basketInfo[ii]->ntl;
        R[ii]             = basketInfo[ii]->R;
        lgdNotional[ii]   = basketInfo[ii]->lgdNotional;
        lgdFloor[ii]      = basketInfo[ii]->lgdFloor;
        lgdCap[ii]        = basketInfo[ii]->lgdCap;
    }

#endif

    try
    {
        /* 
         * 1. need to sort pdep, pskew, beta, lossAmt, qz by ascending pdep 
         */
        orderedBasket = basketInfo;
        sort(orderedBasket.begin(), orderedBasket.end(), NameParam());

        /* set up new market factor containers */
        GaussianIntegrationMethod integMethod = GaussianIntegrationMethod(CCM_INTEG_LOW,CCM_INTEG_UP,CCM_INTEG_NB);
        double M = integMethod.l;

        for (int i=0; i<integMethod.n; i++)
        {
            MarketFactorSP mf = MarketFactorSP(
                new MarketFactor(
                    M,
                    integMethod.weightCalc(M)));

            integrands.push_back(mf);
            M += integMethod.step;
        }

        if (useFast)
        {
            lossSize = 1;
        }
        else
        {
            lossSize = maxLongIdx + maxShortIdx + 1;
        }

        scenarioSureLoss.clear();
        scenarioWeight.clear();

        expLoss = 0.;
        expLossCond = 0.;
    }
    catch (exception& e)
    {
        throw ModelException(e, routine);
    }
}

void CCMConvolution::calcLossDensityRecM4()
{
    const char routine[] = "CCMConvolution::calcLossDensityRecM4";

    try {
        int j;

        /* treat the case nbName == 0 separately */
        if (nbName == 0) return;

        /*
         * 2. calculation of probas assigned to each of the 3 copula
         */
        DoubleArray pgauss(nbName); /* treatment is needed */
        DoubleArray pind(nbName);   /* for defaulted names  */
        DoubleArray tgauss(nbName);

        int i;

        for (i=0; i<nbName; ++i) {
            NameParamSP b = orderedBasket[i];
            if (b->ntl == 0) continue; /* skip defaulted name */
            b->probas(pind[i], pgauss[i], tgauss[i]);
        }

        /* check all proba in [0,1] */
        for (i=0; i<nbName; ++i)
        {
            const string& name = orderedBasket[i]->nameId;
            if(!(fabs(orderedBasket[i]->pdep   - .5) < .5 + 3e-16))
            {
                throw ModelException(routine, "pdep too high for name #"+name);
            }
            if(!(fabs(pgauss[i] - .5) < .5 + 3e-16)) 
            {
                throw ModelException(routine, "pgauss too high for name #"+name);
            }
            if (!(fabs(pind[i]   - .5) < .5 + 3e-16))
            {
                throw ModelException(routine, "pind too high for name #"+name);
            }
        }

        /*
         * 3. Calibration of Recovery: recInfo index is 
         *    relevant to orderedBasket
         */
        LgdParamArray recInfo;
        for (i=0; i<nbName; ++i) {

            LgdParamSP ri = LgdParamSP(new LgdParam());

            CCMCalibration::recoveryModelCalibrate(
                *orderedBasket[i],
                pind[i],
                pgauss[i],
                useFast ? false : true, // is loss amt discretized
                ri); 

            recInfo.push_back(ri);
        }

        /* 
         * 4. Convolution
         * do a convolution algorithm conditional on a dependent default scenario.
         * for n names, there will be upto n dependent default scenario (hopefully
         * less if many pdep are identical).
         * The loss is then shifted by the sure loss amount and multiplied by 
         * pdep. The sure loss amount is a double, which is discretized into
         * the 2 nearest integers shift1 and shift2.
         * We look a t n+1 states, whose proba is 
         * pdep[j-1] - pdep[j] where j in [0..n]
         *      pdep[j] < pdep[j-1] for all j
         *      pdep[-1] = 1, pdep[n] = 0 as extreme conditions
         */

        /* deal with case of a catastrophic default separately */
        hasCatastrophicDefault = false;
        double sureLoss = 0.;
        for (j=1; j<=nbName; ++j)
        {
            sureLoss += recInfo[j-1]->lc;
        }

        double sumW = 1.0 - orderedBasket[0]->pdep;

        // store the sure loss in this scenario
        scenarioSureLoss.push_back(sureLoss);
        //store the weight (used slightly differently to the other scenarios)
        scenarioWeight.push_back(sumW);

        if (sureLoss!=0 && sumW>=TOWER_PROBA_EPSI) /* keep even 1e-13 to keep total expected loss */
        {
            // record the catastrophic case
            hasCatastrophicDefault = true;
        }

        /* Iterate over the possible dependent default scenario */
        for (j=1; j<=nbName; ++j) {
            double weight = j==nbName ? orderedBasket[nbName-1]->pdep:
                orderedBasket[j-1]->pdep - orderedBasket[j]->pdep;

            /* calculate sure loss in this senario */
            sureLoss -= recInfo[j-1]->lc;

            if (fabs(weight) < TOWER_PROBA_EPSI) {
                continue;
            }
            sumW  += weight;
            ASSERT(weight>=0.);

            // store the sure loss in this scenario
            scenarioSureLoss.push_back(sureLoss);
            //store the weight
            scenarioWeight.push_back(weight);

            for (int k=0; k<integrands.size(); k++)
            {
                integrands[k]->addDistribution(lossSize, maxShortIdx);
            }

            // calc loss density conditional on this level of dependent default
            gaussIndepLossDensityCalc(j, pind, tgauss, recInfo);
        }
        
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

void CCMConvolution::getDensities(
    const NameParamSP counterparty,  /* (I) optional. Counterparty info is separate */
    DoubleArray&      density,
    DoubleArray&      densityCond)
{
    const char routine[] = "CCMConvolution::getDensities";
    try
    {
        if (useFast)
        {
            throw ModelException("CCMConvolution::getDensityCond","Fast convolution used : density not calculated");
        }

        density.clear();
        densityCond.clear();

        /* treat the case nbName == 0 separately */
        if (nbName == 0)
        {
            density.resize(1);
            densityCond.resize(1);
            density[0]     = 1.0;
            densityCond[0] = 1.0;
        }    
        else
        {
            double pgaussCpty = 1.; //counterparty gaussian probability
            double pindCpty   = 1.; //counterparty independent probability
            double tgaussCpty = 0.; //counterparty thresholded gaussian probability
            double sumW2      = 0.;
            DoubleArray weight2;    //per market factor weights

            if (!!counterparty) {
                if (fabs(counterparty->survival-1.) < 3e-16)
                {
                    pindCpty   = 1.;
                    pgaussCpty = 1.;
                }
                else
                {
                    double ci  = log(counterparty->pdep)/log(counterparty->survival); 
                    pindCpty   = pow(counterparty->survival, (1.-ci)*counterparty->indep);
                    pgaussCpty = counterparty->survival / (counterparty->pdep*pindCpty);
                }
                tgaussCpty = CCMSkew::fInv(pgaussCpty,
                                           counterparty->beta,
                                           counterparty->qM,
                                           SKEW_FCUMUL_NB_POINT);

                //now calculate per market factor cpty weight2's
                weight2.resize(integrands.size());
                for (int j=0; j<integrands.size(); j++)
                {
                    double w2 = integrands[j]->getWeight2(pindCpty,
                                                          tgaussCpty,
                                                          counterparty->beta,
                                                          counterparty->qM);
                    sumW2 += w2;

                    weight2[j] = w2;
                }
            }

            density.resize(lossSize);
            densityCond.resize(lossSize);

            unsigned int numScenario = scenarioSureLoss.size();
            int shift1, shift2;
            double  q1, q2;

            /* deal with case of a catastrophic default separately */
            if (hasCatastrophicDefault)
            {
                solveSystem(scenarioSureLoss[0], scenarioWeight[0], shift1, q1, shift2, q2);
                density[shift1+maxShortIdx]   = q1;
                density[shift2+maxShortIdx]   += q2;
                if (!!counterparty)
                {
                    densityCond[shift1+maxShortIdx] = q1;
                    densityCond[shift2+maxShortIdx] += q2;
                }
            }
            else
            {
                density[0] = scenarioWeight[0];
                if (!!counterparty) densityCond[0] = scenarioWeight[0];
            }

            DoubleArray work1;
            DoubleArray work2;

            //now deal with dependent defaults
            for (unsigned int k = 1; k<numScenario; k++)
            {
                work1.clear();
                work2.clear();
                work1.resize(lossSize);
                work2.resize(lossSize);

                solveSystem(scenarioSureLoss[k], 1., shift1, q1, shift2, q2);

                //get the weighted distributions across all market factors for this scenario
                double q1ByW2 = 1.;
                double q2ByW2 = 1.;

                for (int j=0; j<integrands.size(); j++)
                {
                    if (!!counterparty)
                    {
                        integrands[j]->accumulateDistribution(k-1, work1, weight2[j], work2);
                    }
                    else
                    {
                        integrands[j]->accumulateDistribution(k-1, work1);
                    }
                }

                if (!!counterparty)
                {
                    q1ByW2 = q1/sumW2;
                    q2ByW2 = q2/sumW2;
                }

                /* update loss density and loss density cond cpty survival */
                for (int i=-maxShortIdx; i<=maxLongIdx; ++i)
                {
                    if (i+shift1 >= -maxShortIdx && i+shift1 <=  maxLongIdx)
                    {
                        density[i+shift1+maxShortIdx] += 
                            scenarioWeight[k] * q1 * work1[i+maxShortIdx];
                        if (!!counterparty)
                            densityCond[i+shift1+maxShortIdx] +=
                                scenarioWeight[k] * q1ByW2 * work2[i+maxShortIdx];
                        ASSERT(fabs(density[i+shift1+maxShortIdx]-.5) <= .5);
                    }

                    if (i+shift2 >= -maxShortIdx && i+shift2 <=  maxLongIdx)
                    {
                        density[i+shift2+maxShortIdx] += 
                            scenarioWeight[k] * q2 * work1[i+maxShortIdx];
                        if (!!counterparty)
                            densityCond[i+shift2+maxShortIdx] +=
                                scenarioWeight[k] * q2ByW2 * work2[i+maxShortIdx];
                        ASSERT(fabs(density[i+shift2+maxShortIdx]-.5) <= .5);
                    }
                }
            } //scenario loop

            if (!counterparty) densityCond = density;

        } //nbName = 0
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

void CCMConvolution::deconvolute(const NameParam& name)
{
    const char routine[] = "CCMConvolution::deconvolute";

    try
    {
        //establish the index of this name in the ordered basket
        int nmIdx=0;
        while ((nmIdx < nbName) && (orderedBasket[nmIdx]->nameId != name.nameId)) nmIdx++;

        if (nmIdx >= nbName)
        {
            throw ModelException(routine, "Failed to find " + name.nameId + " in existing portfolio.");
        }

        //calculate new copula probabilities
        double pind = 0.;
        double pgauss = 0.;
        double tgauss = 0.;

        name.probas(pind, pgauss, tgauss);

        //calculate recovery information for this name
        LgdParamSP ri = LgdParamSP(new LgdParam());

        CCMCalibration::recoveryModelCalibrate(
            name,
            pind,
            pgauss,
            useFast ? false : true, // is loss amt discretized
            ri);

        //deconvolute this name for each market factor
        for (int i=0; i<integrands.size(); i++)
        {
            integrands[i]->deconvolute(nmIdx,
                                       name,
                                       *ri,
                                       pind,
                                       tgauss);
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

void CCMConvolution::reconvolute(const NameParam& oldName,
                                 const NameParam& newName)
{
    const char routine[] = "CCMConvolution::reconvolute";

    try
    {
        deconvolute(oldName);
        convolute(newName);
    }
    catch (exception& e)
    {
        throw ModelException(e, routine);
    }
}

void CCMConvolution::convolute(const NameParam& name)
{
    const char routine[] = "CCMConvolution::convolute";

    try
    {
        //establish the index of this name in the ordered basket
        int nmIdx=0;
        while ((nmIdx < nbName) && (orderedBasket[nmIdx]->nameId != name.nameId)) nmIdx++;

        if (nmIdx >= nbName)
        {
            throw ModelException(routine, "Failed to find " + name.nameId + " in existing portfolio.");
        }

        //calculate new copula probabilities
        double pind = 0.;
        double pgauss = 0.;
        double tgauss = 0.;

        name.probas(pind, pgauss, tgauss);

        //calculate recovery information for this name
        LgdParamSP ri = LgdParamSP(new LgdParam());

        CCMCalibration::recoveryModelCalibrate(
            name,
            pind,
            pgauss,
            useFast ? false : true, // is loss amt discretized
            ri);

        //reconvolute this name for each market factor
        for (int i=0; i<integrands.size(); i++)
        {
            integrands[i]->convolute(nmIdx,
                                     name,
                                     *ri,
                                     pind,
                                     tgauss);
        }

    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

double& CCMConvolution::getExpLoss()
{
    if (!useFast)
    {
        throw ModelException("CCMConvolution::getDensityCond","Usual convolution used : expLoss not calculated");
    }
    return expLoss;
}

double& CCMConvolution::getExpLossCond()
{
    if (!useFast)
    {
        throw ModelException("CCMConvolution::getDensityCond","Usual convolution used : expLoss not calculated");
    }
    return expLossCond;
}

/**
 * Fast convolution algorithm
 */

#undef MAX
#define MAX(a,b) ((a)>(b) ? (a) : (b))

#undef MIN
#define MIN(a,b) ((a)<(b) ? (a) : (b))

#undef MAX3
#define MAX3(a,b,c) MAX(MAX(a,b),(c))

/* ----------------------------------------------------------------------
 * Beta incomplete function and closed form integral
 */
/*
 *  A better approx of a binomial can be done with an incomplete
 *  beta function. Utilities for incomplete beta function NRC p 227
 * 
 *  \todo use cephes as this seems to have troubles with numerical
 *        values like mu=1e-9 s=1e-7... the cephes implementation
 *        looks better.
 */
/* incomplete beta function continued fraction expansion */
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
static double betacf(double a, double b, double x)
{
    static char routine[] = "betacf";

    int m,m2;
    double aa,c,d,del,h,qab,qam,qap;

    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    c=1.0;
    d=1.0-qab*x/qap;
    if (fabs(d) < FPMIN) d=FPMIN;
    d=1.0/d;
    h=d;
    for (m=1;m<=MAXIT;m++) {
        m2=2*m;
        aa=m*(b-m)*x/((qam+m2)*(a+m2));
        d=1.0+aa*d;
        if (fabs(d) < FPMIN) d=FPMIN;
        c=1.0+aa/c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        h *= d*c;
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d=1.0+aa*d;
        if (fabs(d) < FPMIN) d=FPMIN;
        c=1.0+aa/c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        del=d*c;
        h *= del;
        if (fabs(del-1.0) < EPS) break;
    }
    if (m > MAXIT) 
        throw ModelException(routine, "a or b too big, or MAXIT too small");

    return h;
}
#undef MAXIT
#undef EPS
#undef FPMIN

/* incomplete beta function */
static double betai(double a, double b, double x)
{
    static char routine[] = "betai";

    double bt,y;
    double out;

    if (x < 0.0 || x > 1.0) 
        throw ModelException (routine, "Bad x=" + Format::toString(x));

    if (x == 0.0 || x == 1.0) 
        bt=0.0;
    else
        bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));

    if (x < (a+1.0)/(a+b+2.0)) 
    {
        y = betacf(a,b,x);
        out = bt * y / a;
    }
    else
    {
        y = betacf(b,a,1.0-x);
        out = 1.0 - bt * y / b;
    }
    return out;
}


/**
 * Calculate E((L-K)^+)
 * Assuming that L follows a beta distribution
 * with a given expected loss and variance.
 *
 * we use an "incomplete beta" integ close formula 
 */
double ccmBetaTrancheLossInteg(
    double K,      /* (I) lower strike   */
    double mean,   /* (I) expected loss  */
    double var)    /* (I) loss variance  */
{
    double a,b,x, y1,y2;                          
    double e;
    /* strike is >= 1 */
    if (1.-K<1e-10)
        return 0.0;
    /* strike is <= 0 */
    if (K<1e-10)
        return mean;
    /* variance is 0 */
    if (var<1e-10)
        return MAX(mean-K,0.);
    x = mean*(1.-mean)/var - 1.; 
    ASSERT(mean>=0. && mean<=1.);
    /* variance is max variance */
    ASSERT(x>=0.); 
    if (x<1e-12)
        return (1.-K)*mean;
    /* general case */
    a = mean*x;
    b = (1.-mean)*x;
    y1 = betai(b,a+1.,1.-K);
    y2 = betai(b,a   ,1.-K);

    e = a/(a+b) * y1 - K * y2;
    ASSERT(e>=-1.e-13 && e <= 1.-K+ 1.e-13);
    e = MAX(0.,e);
    e = MIN(1.-K, e);
    return e;
}

/* ----------------------------------------------------------------------
 * "Inspired" from convrecm4.c
 */


/* 
   Convolution routine as used at each value of M, returns tranche expected loss 
   for a set of independent default events.
*/
static double fastConvolution(
    const vector<double> &survivalProba, /* (I) name survival probas    */
    const LgdParamArray  p,             /* (I) loss parameters         */
    const vector<double> &weight1,       /* (I) weight on L1            */
    double               K1,            /* (I) lower strike            */
    double               K2,            /* (I) upper strike            */
    double               sureLoss)      /* (I) sure loss               */
{
    long   n = survivalProba.size();
    double out;
    long   i;
    double s, q1,q2, l1,l2;
    double expected, variance, allPos, allNeg, L1, L2, sureTrancheLoss;
    
    double offset = MIN(0.,K1)-MIN(0,K2);

    /* allPos =   positive notional */
    /* allNeg = - negative notional */
    
    allPos = 0.;
    allNeg = 0.;

    for(i=0; i<n; ++i )
    {
        allPos   += MAX3(MAX(p[i]->LGD1d,0.),MAX(p[i]->LGD2d,0.),MAX(p[i]->LGD1d-p[i]->LGD2d,0.));
        allNeg   += -MAX3(MAX(-p[i]->LGD1d,0.),MAX(-p[i]->LGD2d,0.),MAX(-p[i]->LGD1d+p[i]->LGD2d,0.));
    }

    sureTrancheLoss = MIN(MAX(sureLoss-K1,0.), K2-K1);

    /* transform strikes for past losses */
    K1 -= sureLoss;
    K2 -= sureLoss;
    if (fabs(allPos - allNeg)<TINY) 
        return sureTrancheLoss+offset;

    /* E(L) = sum_i{p_i L_i}
       V(L) = sum_i{p_i (1-p_i) L_i^2} */
    variance = 0.;
    expected = 0.;
    for(i=0; i<n; ++i )
    {
        s  = survivalProba[i];
        q1 = (1-s) * weight1[i];
        q2 = (1-s) * (1-weight1[i]);
        
        l1 = p[i]->LGD1d;
        l2 = p[i]->LGD2d;
        
        expected += l1*q1 + l2*q2;
        variance += l1*l1*q1*(1-q1) + l2*l2*q2*(1-q2) - 2 * l1*l2*q1*q2;
    }
    
    /*
     * renormalisation to get loss distribution on [0,1]
     * can handle short names
     */
    expected = (expected - allNeg) / (allPos - allNeg);
    K1       = (K1 - allNeg) / (allPos - allNeg);
    K2       = (K2 - allNeg) / (allPos - allNeg);
    variance /= (allPos - allNeg)*(allPos - allNeg);
    

    L1 = ccmBetaTrancheLossInteg(K1, expected, variance);
    L2 = ccmBetaTrancheLossInteg(K2, expected, variance);


    out = sureTrancheLoss+(L1-L2)*(allPos-allNeg)+offset;
    return out;
}

   
/*
 * calculate loss density conditional for indep-gauss copula (i.e., 
 * conditional on a certain level of dependent default)
 * handles a cpty
 */
static void gaussIndepLossDensityCalcFastConv(
    long                   nbName,      /* (I) number to process            */
    const LgdParamArray    recModel,    /* (I) name loss amt info           */
    const NameParamArray   basketinfo,  /* (I) name parameters              */
    const NameParamSP      cpty,         /* (I) counterparty info is separate*/
    const DoubleArray&     pind,         /* (I) indep survival proba         */
    double                 pindCpty,     /* (I) cpty indep survival proba    */
    const DoubleArray&     pgauss,       /* (I) gauss survival proba         */
    double                 pgaussCpty,   /* (I) cpty gauss survival proba    */
    double                 K1,           /* (I) lower strike                 */
    double                 K2,           /* (I) upper strike                 */
    double                 sureLoss,     /* (I) sure loss                    */
    double                &L,            /* (O) tranche loss                 */
    double                &Lcond)        /* (O) same for cpty conditional
                                            density */
{
    const char routine[] = "gaussIndepLossDensityCalcFastConv";

    try
    {
        int     i, j;
        double  M, sumW1, sumW2, w1, w2, tgaussCpty, LM;
        vector<double> pm(nbName);
        vector<double> tgauss(nbName);
        vector<double> weight1(nbName);
        GaussianIntegrationMethod integMethod = GaussianIntegrationMethod(CCM_INTEG_LOW,CCM_INTEG_UP,CCM_INTEG_NB);
        
        /* initialize threshold */
        for (i=0; i<nbName; ++i)
        {
            (basketinfo[i])->threshold(pgauss[i], tgauss[i]);
        }
        if (!!cpty)
        {
            cpty->threshold(pgaussCpty, tgaussCpty);
        }
        else
        {
            tgaussCpty = 0.;
        }
            
        
        /* 2. Integrate convolution of losses */
        L     = 0.;
        Lcond = 0.;
        for (i=0, M=integMethod.l, sumW1=sumW2=0.; i<integMethod.n; ++i)
        {
            for (j=0; j<nbName; ++j)
            {
                basketinfo[j]->condSurvivalProba(  pind[j],
                                                   tgauss[j],
                                                   M,
                                                   pm[j]);
                
                weight1[j] = CCMCalibration::recoveryWeightCalc(
                    recModel[j]->W1ind,
                    recModel[j]->T1,
                    recModel[j]->isT1Calibrated,
                    basketinfo[j]->indep,
                    basketinfo[j]->betaRec,
                    basketinfo[j]->qM,
                    M);
            }
#ifdef DEBUG_CCM_TXT
            if (fabs(M+3.15)<TINY)
            {
                FILE *f = fopen("ccmfast.txt","a");
                if (f == NULL) goto done;
                fprintf(f,"########### single name proba #################\n");
                for (j=0; j<nbName; ++j)
                    fprintf(f,"%s, %.16lf\n", basketinfo[j]->nameId.c_str(),
                            pm[j]);
                fclose(f);
            }
#endif
            
            /* convolution of n names */
            LM = fastConvolution(pm, recModel, weight1,
                                 K1, K2, sureLoss);
            
#ifdef DEBUG_CCM_TXT
            if (fabs(M+3.15)<TINY)
            {
                FILE *f = fopen("ccmfast.txt","a");
                fprintf(f,"########### tranche loss #################\n");
                fprintf(f, "K1=%lf, K2=%lf, sureLoss=%lf, ELM=%lf\n",
                        K1, K2, sureLoss, LM);
                fclose(f);
            }
#endif
            
            w1 = integMethod.weightCalc(M);
            sumW1 += w1;
            L += w1 * LM;
            
            if (!!cpty)
            {
                //double pmCpty = CCMSkew::condSurvivalProba(
                //    pindCpty, tgaussCpty, cpty->beta, cpty->qM, M);
                double pmCpty;
                cpty->condSurvivalProba(pindCpty, tgaussCpty, M, pmCpty);
                w2 = w1 * pmCpty;
                sumW2 += w2;
                Lcond += w2 * LM;
            }
            
            M += integMethod.step;
        }

        if (!cpty)
        {
            Lcond = L;
        }
        else
        {
            Lcond /=sumW2;
        }

    } catch(exception &e) {
        throw ModelException(e,routine);
    }
}

/**
 * Convolution algorithm.
 * Takes the copula and the raw loss parameters as input
 * Calibrates the correlated recovery model parameters as a first step.
 * For each level of loss caused by the dependence copula, do the 
 * independence/gaussian convolution conditionnaly on each value of the 
 * market variable M.
 *
 * To handle recovery correlation with M, the convolution will use two 
 * loss levels for each default: LGD1[] and LGD2[], which are the same 
 * for any M, but with weights that depend on M. 
 * These weight are a combination of the independence copula weights 
 * and the gaussian copula weights, with the parroximation that the repar-
 * tition of defaults between the two copulae is indepndent of M.
 */
void CCMConvolution::calcFastTrancheLoss(
    const NameParamSP counterparty,  /* (I) optional. Counterparty info is separate */
    double            histLoss,   /* (I) historical loss L(0,t) in $    */
    double            K1,         /* (I) lower strike                 */
    double            K2)         /* (I) upper strike                 */
{
    const char routine[] = "CCMConvolution::calcFastTrancheLoss";

    try
    {
        long    i, j;
        double  sureLoss, sumW;
        double  pindCpty = 1.;
        double  pgaussCpty = 1.;

        /* treat the case nbName == 0 separately */
        if (nbName == 0)
        {
            expLoss     = Maths::creditCollar(histLoss,K1,K2);//MIN(MAX(histLoss-K1,0),K2-K1);
            expLossCond = expLoss;
            return;
        }    

        /*
         * 2. calculation of probas assigned to each of the 3 copula
         */
        DoubleArray pgauss(nbName); /* treatment is needed */
        DoubleArray pind(nbName); /* for defaulted names  */

        for (i=0; i<nbName; ++i) {
            NameParamSP b = orderedBasket[i];
            if (b->ntl == 0) continue; /* skip defaulted name */
            b->probas(pind[i], pgauss[i]);
        }

        /* check all proba in [0,1] */
        for (i=0; i<nbName; ++i)
        {
            const string& name = orderedBasket[i]->nameId;
            if(!(fabs(orderedBasket[i]->pdep   - .5) < .5 + 3e-16))
            {
                throw ModelException(routine, "pdep too high for name #"+name);
            }
            if(!(fabs(pgauss[i] - .5) < .5 + 3e-16)) 
            {
                throw ModelException(routine, "pgauss too high for name #"+name);
            }
            if (!(fabs(pind[i]   - .5) < .5 + 3e-16))
            {
                throw ModelException(routine, "pind too high for name #"+name);
            }
        }

        /*
         * 3. Calibration of Recovery: recInfo index is 
         *    relevant to orderedBasket
         */
        LgdParamArray recInfo;
        recInfo.reserve(nbName);
        for (i=0; i<nbName; ++i) {

            LgdParamSP ri = LgdParamSP(new LgdParam());

            CCMCalibration::recoveryModelCalibrate(
                *orderedBasket[i],
                pind[i],
                pgauss[i],
                false, // is loss amt discretized
                ri); 

            recInfo.push_back(ri);
        }

        if (!!counterparty)
        {
            counterparty->probas(pindCpty, pgaussCpty);
        }

        /* 
         * 4. Convolution
         * do a convolution algorithm conditional on a dependent default scenario.
         * for n names, there will be upto n dependent default scenario (hopefully
         * less if many pdep are identical).
         * The loss is then shifted by the sure loss amount and multiplied by 
         * pdep. The sure loss amount is a double, which is discretized into
         * the 2 nearest integers shift1 and shift2.
         * We look a t n+1 states, whose proba is 
         * pdep[j-1] - pdep[j] where j in [0..n]
         *      pdep[j] < pdep[j-1] for all j
         *      pdep[-1] = 1, pdep[n] = 0 as extreme conditions
         */
        expLoss = expLossCond = 0.;
        
        /* deal with case of a catastrophic default separately */
        sureLoss = histLoss;
        for (j=1; j<=nbName; ++j) 
            sureLoss += recInfo[j-1]->lc;
        sumW = 1.0 - orderedBasket[0]->pdep;
        
        expLoss = sumW * Maths::creditCollar(sureLoss,K1,K2);//MIN(MAX(sureLoss-K1,0),K2-K1);
        expLossCond = expLoss;
        
        /* Iterate over the possible dependent default scenario */
        for (j=1; j<=nbName; ++j) 
        {
            double Lgauss, LgaussND;
            double weight = j==nbName ? orderedBasket[nbName-1]->pdep 
                : orderedBasket[j-1]->pdep - orderedBasket[j]->pdep;
            
            if (fabs(weight) < TOWER_PROBA_EPSI) 
                continue;
            sumW  += weight;
            ASSERT(weight>=0.);
            
            /* calculate sure loss in this senario */
            sureLoss = histLoss;
            for (i=j; i<nbName; ++i)
                sureLoss += recInfo[i]->lc;
            
            /* calc loss density conditional on this level of dependent default */
            gaussIndepLossDensityCalcFastConv(
                j, 
                recInfo,
                orderedBasket,                
                counterparty,
                pind,
                pindCpty,
                pgauss,
                pgaussCpty,
                K1,
                K2,
                sureLoss,
                Lgauss,
                LgaussND);
            
            expLoss     += weight * Lgauss;
            expLossCond += weight * LgaussND;
        }
        ASSERT(fabs(sumW-1.) < TOWER_PROBA_EPSI * (nbName+1));
    } catch(exception &e) {
        throw ModelException(e, routine);
    }
}

void CCMConvolution::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(CCMConvolution, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCCMConvolution);

    //fields
    FIELD( useFast,       "false=>usual convolution");
    FIELD( nbName,        "the number of names in the basket");
    FIELD( lossSize,      "the size of the loss distributions");
    FIELD( maxShortIdx,   "the short bound of the distribution");
    FIELD( maxLongIdx,    "the long bound of the distribution");
    FIELD( orderedBasket, "the basket, ordered by pdep");

    //transient fields
    FIELD( integrands,               "the loss distribution samples");
    FIELD( hasCatastrophicDefault,   "can all names be in default simultaneously ?");
    FIELD( scenarioSureLoss,         "store per-scenario sure loss amounts");
    FIELD( scenarioWeight,           "store per-scenario weights");
    FIELD( expLoss,                  "the expected loss");
    FIELD( expLossCond,              "and that conditional upon the (optional) counterparty's survival");
    //FIELD_MAKE_TRANSIENT(integrands);
    //FIELD_MAKE_TRANSIENT(hasCatastrophicDefault);
    //FIELD_MAKE_TRANSIENT(scenarioSureLoss);
    //FIELD_MAKE_TRANSIENT(scenarioWeight);
    //FIELD_MAKE_TRANSIENT(expLoss);
    //FIELD_MAKE_TRANSIENT(expLossCond);
}

IObject* CCMConvolution::defaultCCMConvolution()
{
    return new CCMConvolution();
}

void LgdParam::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(CCMConvolution::LgdParam, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultLgdParam);

    //transient fields
    FIELD( LGD1,           "LGD1");
    FIELD( LGD2,           "LGD2");
    FIELD( LGD1d,          "LGD1d");
    FIELD( LGD2d,          "LGD2d");
    FIELD( T1,             "T1");
    FIELD( W1ind,          "W1ind");
    FIELD( isT1Calibrated, "isT1Calibrated");
    FIELD( lc,             "lc");
}

IObject* LgdParam::defaultLgdParam()
{
    return new LgdParam();
}

void NameParam::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(CCMConvolution::NameParam, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultNameParam);

    //transient fields
    FIELD( nameId,        "curve id");
    FIELD( survival,      "survival probability");
    FIELD( beta,          "beta");
    FIELD( qM,            "skew parameter (not used yet)");
    FIELD( pdep,          "survival proba assigned to dep copula = p^c");
    FIELD( indep,         "determines survival proba assigned to ind copula = p^(1-c)*a");
    FIELD( cataRecFactor, "cata recovery factor, such that recovery_cata = rc*recovery");
    FIELD( lambda,        "recovery dispersion factor");
    FIELD( betaRec,       "beta_recovery");
    FIELD( ntl,           "name notional (expressed in lossunit)");
    FIELD( R,             "name recovery");
    FIELD( lgdNotional,   "lgdNotional");
    FIELD( lgdFloor,      "lgdFloor");
    FIELD( lgdCap,        "lgdCap");
}

void RFLNameParam::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(CCMConvolution::RFLNameParam, clazz);
    SUPERCLASS(NameParam);
    EMPTY_SHELL_METHOD(defaultRFLNameParam);

    //transient fields
    FIELD(rflParam,         "RFL Parameters");
}


IObject* NameParam::defaultNameParam()
{
    return new NameParam();
}

IObject* RFLNameParam::defaultRFLNameParam()
{
    return new RFLNameParam();
}

void CCMConvolution::MarketFactor::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(CCMConvolution::MarketFactor, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultMarketFactor);

    //transient fields
    FIELD( M,        "the market factor itself");
    FIELD( weight,   "the weight associated with M");
    FIELD( distribs, "loss distributions");
}

IObject* CCMConvolution::MarketFactor::defaultMarketFactor()
{
    return new MarketFactor();
}

void DensityAddin::Results::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(DensityAddin::Results, clazz);
    SUPERCLASS(CObject);

    //fields
    FIELD(density, "loss density");
    FIELD(densityCond, "loss density conditional upon cpty");
    EMPTY_SHELL_METHOD(defaultDensityAddinResults);
}

IObject* DensityAddin::Results::defaultDensityAddinResults()
{
    return new DensityAddin::Results();
}

void ConvolutionAddin::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(ConvolutionAddin, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConvolutionAddin);

    //fields
    FIELD (useFast,       "false=usual convolution; true=fast convolution");
    FIELD (nameId,        "name identifier");
    FIELD (survival,      "name probability of survival");
    FIELD (beta,          "name beta");
    FIELD (qM,            "qM");
    FIELD (pdep,          "pdep");
    FIELD (indep,         "indep");
    FIELD (cataRecFactor, "cataRecFactor");
    FIELD (lambda,        "lambda");
    FIELD (betaRec,       "betaRec");
    FIELD (ntl,           "name notional in loss units");
    FIELD (R,             "name recovery");
    FIELD (lgdNotional,   "lgdNotional");
    FIELD (lgdFloor,      "lgdFloor");
    FIELD (lgdCap,        "lgdCap");
    FIELD (maxLongIdx,    "maxLongIdx");
    FIELD (maxShortIdx,   "maxShortIdx");

    Addin::registerObjectMethod("CONVOLUTION",
                                Addin::RISK,
                                "Returns convolution results",
                                true, //require a handle name
                                Addin::returnHandle,
                                &ConvolutionAddin::run);
}

void DeReConvoluteAddin::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(DeReConvoluteAddin, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultDeReConvoluteAddin);

    //fields
    FIELD (action,  "-1 remove old; 0 remove old then add new; 1 add new");
    FIELD        (convo,   "the original convolution results");
    FIELD        (oldName, "the name to deconvolute");
    FIELD        (newName, "the name to convolute");

    Addin::registerObjectMethod("DERECONVOLUTE",
                                Addin::RISK,
                                "Returns convolution results after applying action for name",
                                true, //require a handle name
                                Addin::returnHandle,
                                &DeReConvoluteAddin::run);
}

void DensityAddin::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(DensityAddin, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultDensityAddin);

    //fields
    FIELD (convo, "convolution input");
    FIELD (cpty, "optional counterparty");
    FIELD_MAKE_OPTIONAL(cpty);

    Addin::registerObjectMethod("CONVOLUTION_DENSITY",
                                Addin::RISK,
                                "Returns loss density and conditional loss density",
                                false,
                                Addin::expandMulti, //expand out the results immediately
                                &DensityAddin::run);
}

IObject* ConvolutionAddin::defaultConvolutionAddin()
{
    return new ConvolutionAddin();
}

IObject* DeReConvoluteAddin::defaultDeReConvoluteAddin()
{
    return new DeReConvoluteAddin();
}

IObject* DensityAddin::defaultDensityAddin()
{
    return new DensityAddin();
}

ConvolutionAddin::ConvolutionAddin(): CObject(TYPE) {}

DeReConvoluteAddin::DeReConvoluteAddin(): CObject(TYPE) {}

DensityAddin::DensityAddin(): CObject(TYPE) {}

IObjectSP ConvolutionAddin::run()
{
    static string method = "ConvolutionAddin::run";

    int nbName = nameId.size();

    if (useFast)
    {
        throw ModelException(method, "fast convolution not yet implemented");
    }

    //now construct portfolio
    NameParamArray names;
    for (int i=0; i<nbName; i++)
    {
        NameParamSP nm = NameParamSP(new NameParam());

        nm->nameId        = nameId[i];
        nm->survival      = survival[i];
        nm->beta          = beta[i];
        nm->qM            = qM[i];
        nm->pdep          = pdep[i];
        nm->indep         = indep[i];
        nm->cataRecFactor = cataRecFactor[i];
        nm->lambda        = lambda[i];
        nm->betaRec       = betaRec[i];
        nm->ntl           = ntl[i];
        nm->R             = R[i];
        nm->lgdNotional   = lgdNotional[i];
        nm->lgdFloor      = lgdFloor[i];
        nm->lgdCap        = lgdCap[i];

        names.push_back(nm);
    }

    //construct the convolutor
    CCMConvolutionSP convoOut = CCMConvolutionSP(
        new CCMConvolution(
            false,
            names,
            maxLongIdx,
            maxShortIdx));

    //and run it
    convoOut->calcLossDensityRecM4();

    //and return
    return convoOut;
}

IObjectSP DeReConvoluteAddin::run()
{
    static string method = "DeReConvoluteAddin::run";

    //deep copy the input convolution
    CCMConvolutionSP convoOut = CCMConvolutionSP(convo.clone());

    //apply the action
    if (action == -1) convoOut->deconvolute(*oldName);
    if (action == 0)  convoOut->reconvolute(*oldName,*newName);
    if (action == 1)  convoOut->convolute(*newName);

    //and return
    return convoOut;
}

IObjectSP DensityAddin::run()
{
    static string method = "DensityAddin::run";

    //set up the density arrays - will be resized in the convolution
    ResultsSP out = ResultsSP(new Results());

    //then get the densities
    convo->getDensities(cpty, out->density,out->densityCond);

    //and return
    return out;
}

DRLIB_END_NAMESPACE
