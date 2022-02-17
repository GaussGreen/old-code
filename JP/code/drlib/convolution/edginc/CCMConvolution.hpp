//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CCMConvolution.hpp
//
//   Description : Convolution Algorithm
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------
#ifndef EDR_CCMCONVOLUTION_HPP
#define EDR_CCMCONVOLUTION_HPP

#include "edginc/LossDistribution.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/RflOnlyParameters.hpp"

DRLIB_BEGIN_NAMESPACE

class CONVOLUTION_DLL CCMConvolution : public CObject
{
public:
    static CClassConstSP const TYPE;

    /* Structure that puts together results of recovery calibration */
    class CONVOLUTION_DLL LgdParam : public CObject
    {
    public:
        int    LGD1;
        int    LGD2;
        double LGD1d;
        double LGD2d;
        double T1;
        double W1ind;
        bool   isT1Calibrated;
        double lc;

        //// constructor that zeroes everything
        LgdParam();

        //----------------
        // CObject methods
        //----------------
        static CClassConstSP const TYPE;
        static void load(CClassSP& clazz);
        static IObject* defaultLgdParam();
    };

    typedef smartPtr<LgdParam>          LgdParamSP;
    typedef array<LgdParamSP, LgdParam> LgdParamArray; //note array of smart pointers

    /**
     * Note that the loss given default is calculated through a payoff formula
     * for a given name i, we have
     * lgdNotional[i]=lgdPayoff[i][0],lgdFloor[i]=lgdPayoff[i][1],lgdCap[i]=lgdPayoff[i][2]
     * lgd[i]=min(max(lgdNotional-R[i](M), lgdFloor), lgdCap)
     *
     * note that lambda and beta_r are really model parameter, we stored them 
     * here because we do not need the cpty recovery details.
     */

    //forward declaration for operator()
    class NameParam;
    typedef smartPtr<NameParam> NameParamSP;
    typedef smartConstPtr<NameParam> NameParamConstSP;
    typedef array<NameParamSP, NameParam> NameParamArray; //note array of smart pointers
    typedef smartPtr<NameParamArray> NameParamArraySP;

    class CONVOLUTION_DLL NameParam : public CObject
    {
    public:
        /* identifier that is useful for error messages */
        string nameId;        /* curve id                           */
        /* probability of survival */
        double survival;      /* survival probability */
        /* model params */
        double beta;          /* beta */
        double qM;            /* skew parameter (not used yet) */
        double pdep;          /* survival proba assigned to dep copula = p^c */
        double indep;         /* determines survival proba assigned 
                                 to ind copula = p^(1-c)*a */
        double cataRecFactor; /* cata recovery factor, such that
                                 recovery_cata = rc*recovery  */
        double lambda;        /* recovery dispersion factor */
        double betaRec;       /* beta_recovery */
        /* portfolio definition */
        double ntl;           /* name notional (expressed in lossunit) */
        double R;             /* name recovery                         */
        double lgdNotional;
        double lgdFloor;
        double lgdCap;

        //// constructor that zeroes everything
        NameParam();
        NameParam(CClassConstSP clazz);
        virtual ~NameParam(){};
        //calculate probabilities
        void probas(double& pind, double& pgauss, double& tgauss) const;
        void probas(double& pind, double& pgauss) const;

        // compute the threshold tgauss from pgauss, i.e Fx-1(pgauss)
        virtual void threshold(double pgauss, double &tgauss) const;

        //calculate probability conditional on a market factor
        /*virtual*/ void probCond(double pind,
                              double tgauss,
                              const LgdParam& recInfo,
                              double M,
                              double& pm,
                              double& weight1) const;

        // compute the conditional survival probability pm
        virtual void condSurvivalProba(  double pind,
                                         double tgauss,
                                         double M,
                                         double&      pm) const;

        //to allow sorting
        bool operator()(const NameParamSP& e1, const NameParamSP& e2);

        //to allow comparison, for optimised convolution
        virtual bool equals(const NameParamSP rhs) const;
        virtual bool equalsUntweakable(const NameParamSP rhs) const;
        virtual int  hashCode() const;
        static int singleTweakDiff(const NameParamArray& state1, const NameParamArray& state2);

        //----------------
        // CObject methods
        //----------------
        static CClassConstSP const TYPE;
        static void load(CClassSP& clazz);
        static IObject* defaultNameParam();
    };

    //forward declaration for operator()
    class RFLNameParam;
    typedef smartPtr<RFLNameParam> RFLNameParamSP;
    typedef smartConstPtr<RFLNameParam> RFLNameParamConstSP;
    typedef array<RFLNameParamSP, RFLNameParam> RFLNameParamArray; //note array of smart pointers

    class CONVOLUTION_DLL RFLNameParam : public NameParam
    {
    public:

        //// constructor that zeroes everything
        RFLNameParam();

        virtual ~RFLNameParam(){};
        // compute the threshold tgauss from pgauss, i.e Fx-1(pgauss)
        virtual void threshold(double pgauss, double &tgauss) const;

        // compute the conditional survival probability pm
        virtual void condSurvivalProba(  double pind,
                                    double tgauss,
                                    double M,
                                    double&      pm) const;

        //to allow comparison, for optimised convolution
        virtual bool equals(const NameParamSP rhs) const;
        virtual int  hashCode() const;
        virtual bool equalsUntweakable(const NameParamSP rhs) const;

        //----------------
        // CObject methods
        //----------------
        static CClassConstSP const TYPE;
        static void load(CClassSP& clazz);
        static IObject* defaultRFLNameParam();

        /* model params */
        RflOnlyParametersConstSP rflParam;
    };

    //inner class to represent the sampling interval used in the numerical integration
    //of loss distributions
    class CONVOLUTION_DLL MarketFactor : public CObject
    {
    public:
        //default constructor to satisfy stl
        MarketFactor();

        //constructor
        MarketFactor(const double            m,
                     const double            w);

        void addDistribution(const int nbLoss,
                             const int maxShortIdx);

        LossDistributionSP getLastDistribution();

        //deconvolute a name
        void deconvolute(const int          nmIdx,
                         const NameParam&   nameInfo,
                         const LgdParam&    recInfo,
                         const double       name_pind,
                         const double       name_tgauss);

        //convolute a name
        void convolute(const int          nmIdx,
                       const NameParam&   nameInfo,
                       const LgdParam&    recInfo,
                       const double       name_pind,
                       const double       name_tgauss);

        void accumulateDistribution(int scenario, DoubleArray& dist);
        void accumulateDistribution(int scenario, DoubleArray& dist,
                                    double weight2, DoubleArray& distCond);

        double getFactor();
        double getWeight();
        //double getPm(int nmIdx) const;
        //double getWeight1(int nmIdx) const;
        double getWeight2(const double pindCpty,
                          const double tgaussCpty,
                          const double cptyBeta,
                          const double cptyQM);

        int getLastRunName();

        //----------------
        // CObject methods
        //----------------
        static CClassConstSP const TYPE;
        static void load(CClassSP& clazz);
        static IObject* defaultMarketFactor();

    private:

        //fields
        double                    M;          //the market factor itself
        double                    weight;     //the weight associated with M

        //DoubleArray               pm;         //per name probability of survival conditional on M
        //DoubleArray               weight1;    //per name weight


        LossDistributionArray     distribs;   //loss distributions
                                              //   - elements represent the loss units, from max short to max long loss
                                              //   - traditionally the convolution algorithm would accumulate this by
                                              //     looping over each default scenario in turn, and within that over each market
                                              //     factor in turn
                                              //   - it turns out that for a given market factor, the convolution of the
                                              //     nth name scenario requires first convoluting n-1 names, and then convoluting the
                                              //     nth name
                                              //   - so we can reuse the unweighted loss distribution from the previous name scenario
                                              //     to reduce the overall calculation time
                                              //   - therefore during the pricing run, this array will be "work in progress" until
                                              //     all name scenarios have been run
                                              //   - for more efficient tweaking, we require that we retain the final n-name scenario
                                              //     which is the state that this array will be left in
    };

    //cant use auto_ptr in an stl container....
    typedef smartPtr<MarketFactor>              MarketFactorSP;
    typedef array<MarketFactorSP, MarketFactor> MarketFactorArray; //note array of smart pointers

    //constructors
    CCMConvolution();
    CCMConvolution(const bool              useFastMethod, /* (I) false=>"usual" convolution */
                   const NameParamArray&   basketInfo,    /* (I) name loss amt info[nbName] */
                   const int               maxLongIdx,    /* (I) density long size */
                   const int               maxShortIdx);  /* (I) density shorts size */

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
    void calcLossDensityRecM4();

    /**
     * Fast Convolution algorithm.
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
    void calcFastTrancheLoss(
        const NameParamSP counterparty,  /* (I) optional. Counterparty info is separate */
        double histLoss,   /* (I) historical loss L(0,t) in $	*/
        double K1,         /* (I) lower strike                 */
        double K2);        /* (I) upper strike                 */

    void getDensities(
        const NameParamSP counterparty,  /* (I) optional. Counterparty info is separate */
        DoubleArray& density,
        DoubleArray& densityCond);

    //methods to apply to an existing convolution
    void deconvolute(
        const NameParam& name);
    void reconvolute(
        const NameParam& oldName,
        const NameParam& newName);
    void convolute(
        const NameParam& name);

    double&       getExpLoss();
    double&       getExpLossCond();

    //----------------
    // CObject methods
    //----------------
    static void load(CClassSP& clazz);
    static IObject* defaultCCMConvolution();

private:

    //private methods
    void gaussIndepLossDensityCalc(const int nbName,
                                   const DoubleArray& pind,
                                   const DoubleArray& tgauss,
                                   const LgdParamArray& recInfo);

    //fields
    bool                useFast;                //false=>"usual" convolution
    int                 nbName;                 //the number of names in the basket
    int                 lossSize;               //the size of the loss distributions
    int                 maxShortIdx;            //the short bound of the distribution
    int                 maxLongIdx;             //the long bound of the distribution

    NameParamArray      orderedBasket;          //the basket, ordered by pdep
    //DoubleArray         pgauss;                 //per name gaussian probability
    //DoubleArray         pind;                   //per name independent probability
    //DoubleArray         tgauss;                 //per name thresholded gaussian probability

    //LgdParamArray       recInfo;                //per name calibrated recovery information
    MarketFactorArray   integrands;             //the loss distribution samples

    bool                hasCatastrophicDefault; //can all names be in default simultaneously ?
    DoubleArray         scenarioSureLoss;       //store per-scenario sure loss amounts
    DoubleArray         scenarioWeight;         //store per-scenario weights

    double              expLoss;                //the expected loss
    double              expLossCond;            //and that conditional upon the (optional) counterparty's survival
};

typedef smartPtr<CCMConvolution> CCMConvolutionSP;
typedef array<CCMConvolutionSP, CCMConvolution> CCMConvolutionArray; //note array of SP's

//for unit testing the convolution code base
class CONVOLUTION_DLL ConvolutionAddin : public CObject
{
public:
    static CClassConstSP const TYPE;

    IObjectSP run();

private:

    //----------------
    // CObject methods
    //----------------

    static void load(CClassSP& clazz);

    //-------------------------
    // ConvolutionAddin methods
    //-------------------------
    static IObject* defaultConvolutionAddin();

    //for reflection
    ConvolutionAddin();

    //fields
    bool        useFast;    //TRUE => run "fast" convolution method
                            //FALSE=> run "usual" convolution method

    //to form individual CCMConvolution::NameParam objects
    StringArray nameId;
    DoubleArray survival;
    DoubleArray beta;
    DoubleArray qM;
    DoubleArray pdep;
    DoubleArray indep;
    DoubleArray cataRecFactor;
    DoubleArray lambda;
    DoubleArray betaRec;
    DoubleArray ntl;
    DoubleArray R;
    DoubleArray lgdNotional;
    DoubleArray lgdFloor;
    DoubleArray lgdCap;

    //density bounds
    int         maxLongIdx;
    int         maxShortIdx;
};

//to unit test the removal/addition of names from/to the baseline convolution
class CONVOLUTION_DLL DeReConvoluteAddin : public CObject
{
public:
    static CClassConstSP const TYPE;

    IObjectSP run();

private:

    //----------------
    // CObject methods
    //----------------

    static void load(CClassSP& clazz);

    //-------------------------
    // DeconvoluteAddin methods
    //-------------------------
    static IObject* defaultDeReConvoluteAddin();

    //for reflection
    DeReConvoluteAddin();

    //fields
    int                         action;  //-1 remove; 0 no action; 1 add
    CCMConvolutionSP            convo;   //the original convolution results
    CCMConvolution::NameParamSP oldName; //the name to deconvolute
    CCMConvolution::NameParamSP newName; //the name to convolute
};

class CONVOLUTION_DLL DensityAddin : public CObject
{
public:
    static CClassConstSP const TYPE;

    class CONVOLUTION_DLL Results : public CObject
    {
    public:
        //--------
        // CObject
        //--------
        static CClassConstSP const TYPE;
        static void load(CClassSP& clazz);

        //-------------------------
        // ConvolutionAddin methods
        //-------------------------
        static IObject* defaultDensityAddinResults();

        Results() : CObject(TYPE) {};

        //fields
        DoubleArray density;
        DoubleArray densityCond;
    };

    typedef smartPtr<Results> ResultsSP;

    IObjectSP run();

private:

    //----------------
    // CObject methods
    //----------------

    static void load(CClassSP& clazz);

    //-------------------------
    // DensityAddin methods
    //-------------------------
    static IObject* defaultDensityAddin();

    //for reflection
    DensityAddin();

    //fields
    CCMConvolutionSP            convo;   //the convolution results
    CCMConvolution::NameParamSP cpty;    //the optional counterparty information

};

double ccmBetaTrancheLossInteg(
    double K,      /* (I) lower strike   */
    double mean,   /* (I) expected loss  */
    double var);    /* (I) loss variance  */

DRLIB_END_NAMESPACE
#endif
