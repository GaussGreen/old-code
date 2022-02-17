//----------------------------------------------------------------------------
//
//   Group       : QR Equities 
//
//   Filename    : DependenceGauss.hpp
//
//   Description : Holds Basic Dependence (Maker) and POISSON and CLAYTON
//
//----------------------------------------------------------------------------
#ifndef DEPENDENCE_HPP
#define DEPENDENCE_HPP

#include <map>
#include "edginc/Class.hpp"
#include "edginc/Function.hpp"
#include "edginc/MultiFactors.hpp"
#include "edginc/MarketDataFetcher.hpp"

DRLIB_BEGIN_NAMESPACE

/***************************/
/* BASIC DEPENDENCE OBJECT */
/***************************/
class MCARLO_DLL Dependence : public CObject {
public:
    static CClassConstSP const TYPE;

    /* validation */
    virtual void validatePop2Object() = 0;

    /* correlate a matrix of random numbers */
    virtual void correlateSeries( DoubleMatrix& noise, int pathIdx) = 0;

    /** Returns the [forward] correlations between dates[index] and
        dates[index+1] that are implied by this dependence. (ie it
        is essentially E(Zi.Zj) where Zi and Zj are two correlated series).
        The index must be a valid index into the array of dates with which
        this object was created */
    virtual DoubleMatrix getCorrelations(int index) = 0;
    virtual double getCorrelationIJ(int i, int j, int index) = 0;
	
	virtual int getNumAssets() = 0;

private:
    /* Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

public:
	/* constructor */
    Dependence();
    Dependence(const CClassConstSP& clazz);
};

typedef smartConstPtr<Dependence> DependenceConstSP;
typedef smartPtr<Dependence> DependenceSP;
#ifndef QLIB_DEPENDENCE_CPP
EXTERN_TEMPLATE(class MCARLO_DLL_SP smartPtr<Dependence>);
#else
INSTANTIATE_TEMPLATE(class MCARLO_DLL smartPtr<Dependence>);
#endif

/**********************/
/* POISSON DEPENDENCE */
/**********************/
class MCARLO_DLL Poisson : public Dependence {
public:

    Poisson();
    ~Poisson();

	/* validation */
    virtual void validatePop2Object();

    // constructor
    Poisson( double crashRate );

    /* correlate a matrix of random numbers */
    virtual void correlateSeries( DoubleMatrix& noise, int pathIdx);

    /* see description in pure virtual declaration */
    virtual DoubleMatrix getCorrelations(int index);
	virtual double getCorrelationIJ(int i, int j, int index);

	virtual int getNumAssets() {return 0;}

    // calculates the cumulative Poisson distribution with intensity crashRate*time
    void numJumpsPreprocess( DoubleArray& timeInYears, int maxJumps );

    // calculates the number of jumps between dates given as datesInYears
    void numJumps( DoubleMatrix& noise );
    void numJumps( double* noise, int nbRows );

private:
    double          crashRate;
    DoubleMatrix    cumulative; // matrix P[X(j)<=i],i=0,1,..., 
                                // X(j) with intensity crashRate*t[j]
};

typedef smartConstPtr<Poisson> PoissonConstSP;
typedef smartPtr<Poisson> PoissonSP;
#ifndef QLIB_DEPENDENCE_CPP
EXTERN_TEMPLATE(class MCARLO_DLL_SP smartPtr<Poisson>);
#else
INSTANTIATE_TEMPLATE(class MCARLO_DLL smartPtr<Poisson>);
#endif

/**********************/
/* CLAYTON DEPENDENCE */
/**********************/
class MCARLO_DLL Clayton : public Dependence {
public:

    Clayton();
    ~Clayton();

    /* validation */
    void validatePop2Object();

    /** constructor */
    Clayton(
        const double        theta,
        const DoubleMatrix& alpha );

    /* correlate a matrix of random numbers */
    virtual void correlateSeries( DoubleMatrix& noise, int pathIdx);

    /* see description in pure virtual declaration */
    virtual DoubleMatrix getCorrelations(int index);
	virtual double getCorrelationIJ(int i, int j, int index);
	
	virtual int getNumAssets() {return nbAssets;}

private:
    int             nbAssets;
    int             nbSteps;
    double          theta; 
    DoubleMatrix    alpha;

    double          thetaUsed; 
    DoubleArray     tempA;
    DoubleArray     tempU;
    DoubleArray     tempX;
    double*         gammaRandom;
};

/**************************/
/* BASIC DEPENDENCE MAKER */
/**************************/

class MCARLO_DLL DependenceMaker : public CObject {
public:
    static CClassConstSP const TYPE;

    /** modifies MDF if necessary, default implementation empty */
    virtual void modifyMarketDataFetcher(MarketDataFetcherSP mdf); 

    class MCARLO_DLL ISupport {
    public:
        virtual ~ISupport() {}
    };
    
    /** interface that the ? must implement */
    virtual DependenceSP createDependence(const DependenceMaker::ISupport* support) const = 0;

    // Invoked when Class is 'loaded'
    static void load(CClassSP& clazz);            

    DependenceMaker(const CClassConstSP& clazz);
};

typedef smartConstPtr<DependenceMaker> DependenceMakerConstSP;
typedef smartPtr<DependenceMaker> DependenceMakerSP;
#ifndef QLIB_DEPENDENCE_CPP
EXTERN_TEMPLATE(class MCARLO_DLL_SP smartPtr<DependenceMaker>);
#else
INSTANTIATE_TEMPLATE(class MCARLO_DLL smartPtr<DependenceMaker>);
#endif

/********************/
/* BASIC SKEW MAKER */
/********************/

class MCARLO_DLL SkewMaker : public DependenceMaker {
public:
    static CClassConstSP const TYPE; 

    SkewMaker(const CClassConstSP& clazz);
    
    static void load(CClassSP& clazz);                
};

typedef smartPtr<SkewMaker> SkewMakerSP;

#ifndef QLIB_DEPENDENCE_CPP
EXTERN_TEMPLATE(class MCARLO_DLL_SP smartPtr<SkewMaker>);
#else
INSTANTIATE_TEMPLATE(class MCARLO_DLL smartPtr<SkewMaker>);
#endif

/***********/
/* for IMS */
/***********/

class DependenceMakerGauss;
typedef smartPtr<DependenceMakerGauss> DependenceMakerGaussSP;
class DependenceMakerGaussSrm;
typedef smartPtr<DependenceMakerGaussSrm> DependenceMakerGaussSrmSP;
class DependenceMakerGaussTerm;
typedef smartPtr<DependenceMakerGaussTerm> DependenceMakerGaussTermSP;
class DependenceMakerGaussTermSrm;
typedef smartPtr<DependenceMakerGaussTermSrm> DependenceMakerGaussTermSrmSP;
class DependenceMakerLocalCorr;
typedef smartPtr<DependenceMakerLocalCorr> DependenceMakerLocalCorrSP;

#define DEPENDENCEMAKER_TYPE_GAUSS          "GAUSS"
#define DEPENDENCEMAKER_TYPE_GAUSS_SRM      "GAUSSSRM"
#define DEPENDENCEMAKER_TYPE_GAUSS_TERM     "GAUSSTERM"
#define DEPENDENCEMAKER_TYPE_GAUSS_TERM_SRM "GAUSSTERMSRM"
#define DEPENDENCEMAKER_TYPE_LOCAL_CORR     "LOCALCORR"

class MCARLO_DLL DependenceMakerWrapper : public CObject,
							   virtual public ITypeConvert {
public: 
    ~DependenceMakerWrapper();
	string                          dependenceMakerType; 
	DependenceMakerGaussSP          dependenceMakerGauss;
	DependenceMakerGaussSrmSP       dependenceMakerGaussSrm;
	DependenceMakerGaussTermSP      dependenceMakerGaussTerm;
	DependenceMakerGaussTermSrmSP   dependenceMakerGaussTermSrm;
    DependenceMakerLocalCorrSP      dependenceMakerLocalCorr;

private:
    DependenceMakerWrapper(const DependenceMakerWrapper& rhs);
    DependenceMakerWrapper& operator=(const DependenceMakerWrapper& rhs);

	DependenceMakerSP             realDependenceMaker;
    
public:
	static CClassConstSP const TYPE;

	/** validatePop2Object is called first */
	void validatePop2Object();

	/** convert is called second */
	virtual void convert(IObjectSP&    object,
					     CClassConstSP requiredType) const;

	/** Invoked when Class is 'loaded' */
	static void load(CClassSP& clazz);

	// for reflection
	DependenceMakerWrapper();

	static IObject* defaultDependenceMakerWrapper();

	DependenceMakerWrapper(DependenceMakerGaussSP dependenceMaker); 
	DependenceMakerWrapper(DependenceMakerGaussTermSP dependenceMaker); 
    DependenceMakerWrapper(DependenceMakerLocalCorrSP dependenceMaker);
};

typedef smartPtr<DependenceMakerWrapper> DependenceMakerWrapperSP;

DRLIB_END_NAMESPACE
#endif // DEPENDENCE_HPP
