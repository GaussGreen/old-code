//----------------------------------------------------------------------------
//
//   Group       : QR Equities 
//
//   Filename    : DependenceLocalCorr.hpp
//
//   Description : Holds LocalCorr Dependence (Maker)
//
//----------------------------------------------------------------------------
#ifndef DEPENDENCE_LOCAL_CORR_HPP
#define DEPENDENCE_LOCAL_CORR_HPP

#include <map>
#include "edginc/Class.hpp"
#include "edginc/Dependence.hpp"
#include "edginc/MCPathConfig.hpp"
#include "edginc/LocalCorrSqueeze.hpp"

DRLIB_BEGIN_NAMESPACE

/*************************/
/* LOCAL CORR DEPENDENCE */
/*************************/
class IMCRandom;
typedef refCountPtr<IMCRandom> IMCRandomSP;

class MCARLO_DLL LocalCorr : public Dependence {
public:    
    LocalCorr();
    ~LocalCorr();    

    /** proper constructor */
    LocalCorr(DateTimeSP            valueDate,
              int                   nbIF,              
              IMCRandomSP           randomGenIF,
              DateTimeArraySP       idiosynFactorTimeline,
              int                   nbMF,
              IMCRandomSP           randomGenMF, 
              DateTimeArraySP       marketFactorTimeline,
              IntArraySP            fwdCorrIndexArray,
              DoubleMatrixArraySP   fwdCorrBetaFactorArray,
              IntArraySP            indexArrayMF, 
              LocalCorrSqueezeArray localCorrSqueezeArray,
              IntArray              localCorrSqueezeIndexArray);
    /** validation */
    virtual void validatePop2Object();

    /** correlate a matrix of random numbers */
    virtual void correlateSeries( DoubleMatrix& noise, int pathIdx );

    /**  see description in pure virtual declaration */
    virtual DoubleMatrix getCorrelations(int index);
	virtual double getCorrelationIJ(int i, int j, int index);
	virtual int getNumAssets() {return nbIF;}

protected:
    /** fields */
    int             nbIF;
    int             nbTimeStepsIF; 
    IMCRandomSP     randomGenIF;
    
    int             nbMF;
    int             nbTimeStepsMF; 
    IMCRandomSP     randomGenMF;
    
    DoubleMatrixArraySP storedIntegralsMF;
    //CDoubleMatrixSP     testStoredIntegralsMF;
    CDoubleMatrixSP     storedMF;

    int             nbRegions;

    DoubleArraySP   deltaTimeIF;
    DoubleArraySP   sqrtDeltaTimeIF;
    DoubleArraySP   deltaTimeMF;
    DoubleArraySP   sqrtDeltaTimeMF;

    LocalCorrSqueezeArray   localCorrSqueezeArray;
    IntArray                localCorrSqueezeIndexArray;
    
    IntArraySP             fwdCorrIndexArray;
    DoubleMatrixArraySP    fwdCorrBetaFactorArray;
    IntArraySP             indexArrayMF;          

    CDoubleMatrixSP storeCorrRN;
};

/*******************************/
/* LOCAL CORR DEPENDENCE MAKER */
/*******************************/

class IRandom;
typedef smartPtr<IRandom> IRandomSP;

class RandomCache;
typedef smartPtr<RandomCache> RandomCacheSP;

class MCARLO_DLL DependenceMakerLocalCorr : public SkewMaker {                                            
public:
    static CClassConstSP const TYPE;

    static const string DEFAULT_MARKET_FACTOR_FREQ;  
    static const int DEFAULT_NB_FACTORS; 
    static const double DEFAULT_PRECISION; 
    static const int DEFAULT_MAX_ITER; 

    virtual void validatePop2Object();    

    void modifyMarketDataFetcher(MarketDataFetcherSP mdf); 

    class MCARLO_DLL Support : virtual public DependenceMaker::ISupport {
    public:        
        virtual DateTimeArray getSimDates() const = 0;
        virtual DoubleMatrix getFwdVarAtDates(bool interpolateAtmFwd) const = 0;
        virtual const IMultiFactors* getMultiFactors() const = 0;
        virtual DependenceMakerSP getDependenceMaker() const = 0;
        virtual int getNbPastDates() const = 0;
        virtual const IRandomSP getRandomGenerator() const = 0;
        virtual const MCPathConfig::RandomCacheSP getIdiosynFactorRandomCache() const = 0;
        virtual const MCPathConfig::RandomCacheSP getMarketFactorRandomCache() const = 0;
        virtual bool carefulRandoms() const = 0;        
    };

    /** interface that the instrument must implement */
    virtual DependenceSP createDependence(const DependenceMaker::ISupport* support) const;     

    /** invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);    

    /** for reflection */
    DependenceMakerLocalCorr();
    DependenceMakerLocalCorr(const CClassConstSP& clazz);
    
    static IObject* defaultDependenceMakerLocalCorr();    

protected:
    /** registered fields */
    string                  samplingFrequency;      
    int                     nbFactors;
    double                  precision;
    int                     maxIter;        
};

typedef smartConstPtr<DependenceMakerLocalCorr> DependenceMakerLocalCorrConstSP;
typedef smartPtr<DependenceMakerLocalCorr> DependenceMakerLocalCorrSP;


DRLIB_END_NAMESPACE
#endif // DEPENDENCE_LOCAL_CORR_HPP
