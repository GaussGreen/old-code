//----------------------------------------------------------------------------
//
//   Group       : QR Equities 
//
//   Filename    : DependenceGauss.hpp
//
//   Description : Holds Gauss (SRM) Dependence (Maker)
//
//----------------------------------------------------------------------------
#ifndef DEPENDENCE_GAUSS_HPP
#define DEPENDENCE_GAUSS_HPP

#include <map>
#include "edginc/Class.hpp"
#include "edginc/Dependence.hpp"
#include "edginc/Function.hpp"
#include "edginc/SparseMatrix.hpp"
#include "edginc/MultiFactors.hpp"
#include "edginc/RevertTypeConvert.hpp"

DRLIB_BEGIN_NAMESPACE

/************************************************************/
/* MODES FOR FWDCORRELATION						            */
/* FwdCorrMode -> FwdCorrModeStandard & FwdCorrModeSimdates */
/************************************************************/

class MCARLO_DLL FwdCorrMode : public CObject {
public:
    static CClassConstSP const TYPE;
    
    /** invoked when class is loaded */
    static void load(CClassSP& clazz);
    virtual DoubleArray getBarriers() const = 0;
    virtual void checkInputs() const = 0;
    virtual bool doInterpolateAtmFwd() const = 0;

    /** for reflection */
    FwdCorrMode();
    FwdCorrMode(const CClassConstSP& type);
};
typedef smartPtr<FwdCorrMode> FwdCorrModeSP;

class MCARLO_DLL FwdCorrModeStandard : public FwdCorrMode {
public:
    static CClassConstSP const TYPE;

    virtual void checkInputs() const;

    virtual DoubleArray getBarriers() const;

    /** invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** for reflection */
    FwdCorrModeStandard();

    /** default constructor */
    FwdCorrModeStandard(const CClassConstSP& clazz);

    static IObject* defaultFwdCorrModeStandard();

    string getFwdCorrInterval() const; 
    int getNbStubDays() const; 
    virtual bool doInterpolateAtmFwd() const;

protected:
    string              fwdCorrInterval;
    int                 nbStubDays;
    double              fwdMaxSqError;
    bool                interpolateAtmFwd; // vol interp atm fwd or atm spot
};
typedef smartPtr<FwdCorrModeStandard> FwdCorrModeStandardSP;

class MCARLO_DLL FwdCorrModeSimdates : public FwdCorrMode {
public:
    static CClassConstSP const TYPE;

    virtual void checkInputs() const;

    virtual DoubleArray getBarriers() const;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** for reflection */
    FwdCorrModeSimdates();
    
    /** default constructor */
    FwdCorrModeSimdates(const CClassConstSP& clazz);

    /** constructor for backwards compatibility */
    FwdCorrModeSimdates(double barrierLow, double barrierHigh);

    static IObject* defaultFwdCorrModeSimdates();

    virtual bool doInterpolateAtmFwd() const;

protected:
    double              fwdMaxSqErrorLow;
    double              fwdMaxSqErrorHigh;
    bool                interpolateAtmFwd; // vol interp atm fwd or atm spot
};
typedef smartPtr<FwdCorrModeSimdates> FwdCorrModeSimdatesSP;

class MCARLO_DLL FwdCorrModeWrapper : public CObject {
public: // how can I have this protected or private?
    string                  fwdCorrModeType;
    FwdCorrModeStandardSP   fwdCorrModeStandard;
    FwdCorrModeSimdatesSP   fwdCorrModeSimdates;
    
private:
    FwdCorrModeSP           realFwdCorrMode;

public:
    static CClassConstSP const TYPE;

    /** validatePop2Object is called first */
    void validatePop2Object();

    /** additional helper method for constructor & validatePop2Object */
    void init();

    FwdCorrModeSP getFwdCorrMode() const;

    /** invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
     
    /** for reflection */
    FwdCorrModeWrapper();
    FwdCorrModeWrapper(const CClassConstSP& clazz);

    static IObject* defaultFwdCorrModeWrapper();
};
typedef smartPtr<FwdCorrModeWrapper> FwdCorrModeWrapperSP;

/***************************************************************/
/* GAUSS DEPENDENCE                                            */
/* independent steps, assets with step independent correlation */
/***************************************************************/
class MCARLO_DLL Gauss : public Dependence {
public:
    Gauss();
    ~Gauss();

    // validation
    virtual void validatePop2Object();

    /**** constructor ****/
    Gauss(const DoubleMatrix& correlation); 

    // correlate a matrix of random numbers
    virtual void correlateSeries( DoubleMatrix& noise, int pathIdx );    

    //// see description in pure virtual declaration
    virtual DoubleMatrix getCorrelations(int index);
	virtual double getCorrelationIJ(int i, int j, int index);

	virtual int getNumAssets() {return nbAssets;}

    // interfaced part of Beta-correlations setup available for ExtendedSparseGauss
    // will see if it is needed, maybe not
   // virtual void setBetaCorrelations(const SparseDoubleMatrix& betas, const DoubleArray& maxcorr) {}

protected:
    DoubleMatrix    correlation;
    int             nbAssets;
    DoubleMatrix    correlCoeffs;
    bool            isTrivial; // corr matrix is idenity => thus skip correlate ... 
};

/**************************/
/* GAUSS DEPENDENCE MAKER */
/**************************/

class DependenceMakerGauss; 
typedef smartConstPtr<DependenceMakerGauss> DependenceMakerGaussConstSP;
typedef smartPtr<DependenceMakerGauss> DependenceMakerGaussSP;

class MCARLO_DLL DependenceMakerGauss : public DependenceMaker,
                             virtual public IRevertTypeConvert {
public:
    static CClassConstSP const TYPE;

    class MCARLO_DLL Support : virtual public DependenceMaker::ISupport {
    public:
        /** retrieve correlation matrix */
        virtual CDoubleMatrixConstSP getGaussData() const = 0; 		
    };
    
    /** interface that the instrument must implement */
    virtual DependenceSP createDependence(const DependenceMaker::ISupport* support) const;     

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** constructor */
    DependenceMakerGauss();    
    static IObject* defaultDependenceMakerGauss();

    /** for the IMS interface (see IRevertTypeConvert) */
    virtual IObjectSP revert(const string& interfaceType) const;
};

/******************************/
/* GAUSS SRM DEPENDENCE MAKER */
/******************************/

class MCARLO_DLL DependenceMakerGaussSrm : public DependenceMaker {
public:
    static CClassConstSP const TYPE;

    void validatePop2Object();
    void init();
    void checkInputs() const;

    /** modifies MDF if necessary */
    virtual void modifyMarketDataFetcher(MarketDataFetcherSP mdf); 

    void setCorrTermStructureMode(bool doCorrTermStructure);
    bool getCorrTermStructureMode() const;

    class MCARLO_DLL Support : virtual public DependenceMaker::ISupport {
    public:        
        /** general helper method */
        virtual int nbEqEqAssets() const = 0;

        /** no correlation term structure, no correlation mapping */
        virtual vector<SparseDoubleMatrixSP> createSparseGaussMatrixArray(
            const DependenceMakerGaussSrm* dependenceMaker) const = 0;
    
        /** no correlation term structure, yes correlation mapping */
        virtual DateTimeArray getSimDates() const = 0;
        virtual vector<SparseDoubleMatrixSP> createSparseGaussMatrixArray(
            const DependenceMakerGaussSrm*  dependenceMaker,
            const IntArray&                 fwdCorrIndexArray,
            const DateTimeArray&            fwdCorrDatesArray) const = 0;
        
        /** yes correlation term structure and yes/no correlation mapping */
        virtual DoubleMatrix getFwdVarAtDates() const = 0;
        virtual void getCorrelationData(
            CorrelationCommonArray& corrObjArray,
            CorrelationTermArray&   corrTermObjArray) const = 0; 
        
        virtual vector<SparseDoubleMatrixSP> createSparseGaussTermMatrixArray(
            const DependenceMakerGaussSrm*      dependenceMaker,
            DoubleMatrixArraySP                 fwdCorrelations,
            const IntArray&                     fwdCorrIndexArray,
            const DateTimeArray&                fwdCorrDatesArray) const = 0;

        /** further methods */
        virtual SparseDoubleMatrixSP getBetaCorrelations(
            const DependenceMakerGaussSrm* dependenceMaker) const {return SparseDoubleMatrixSP();}
        virtual double getMaxBetaCorr(const DependenceMakerGaussSrm* dependenceMaker) const {return 0.8;}
    };

    /** interface that the instrument must implement */
    virtual DependenceSP createDependence(const DependenceMaker::ISupport* support) const;     

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** constructor */
    DependenceMakerGaussSrm();
    DependenceMakerGaussSrm(CClassConstSP clazz);
    static IObject* defaultDependenceMakerGaussSrm();
    
    bool            doCorrMapping() const; 
    double          getEigenValueFloor() const; 
    double          getMaxSqError() const; 
    FwdCorrModeSP   getFwdCorrModeUsed() const;
    
protected:   
    /** registered fields */
    bool                    corrMapping;
    bool                    corrTermStructure;
    double                  eigenValueFloor;
    double                  maxSqError; 
    FwdCorrModeWrapperSP    fwdCorrModeDetails;
     /** transient fields */
    FwdCorrModeSP           fwdCorrModeUsed;
};

typedef smartConstPtr<DependenceMakerGaussSrm> DependenceMakerGaussSrmConstSP;
typedef smartPtr<DependenceMakerGaussSrm> DependenceMakerGaussSrmSP;

DRLIB_END_NAMESPACE
#endif // DEPENDENCE_GAUSS_HPP
