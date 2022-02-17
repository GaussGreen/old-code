//----------------------------------------------------------------------------
//
//   Group       : QR Equities 
//
//   Filename    : DependenceGaussTerm.hpp
//
//   Description : Holds GaussTerm (SRM) Dependence (Maker)
//
//----------------------------------------------------------------------------
#ifndef DEPENDENCE_GAUSS_TERM_HPP
#define DEPENDENCE_GAUSS_TERM_HPP

#include <map>
#include "edginc/Class.hpp"
#include "edginc/Dependence.hpp"
#include "edginc/Function.hpp"
#include "edginc/SparseMatrix.hpp"
#include "edginc/MultiFactors.hpp"
#include "edginc/RevertTypeConvert.hpp"

DRLIB_BEGIN_NAMESPACE

/*************************************************************/
/* GAUSS TERM DEPENDENCE                                     */
/* independent steps, assets with step dependent correlation */
/*************************************************************/
class MCARLO_DLL GaussTerm : public Dependence {
public:    
    GaussTerm();
    ~GaussTerm();

    // validation
    virtual void validatePop2Object();

    /**** constructor ****/
    GaussTerm(const DoubleMatrixArraySP& correlationTerm,
              const IntArray&            fwdCorrIndexArray);

    // correlate a matrix of random numbers
    virtual void correlateSeries( DoubleMatrix& noise, int pathIdx);

    //// see description in pure virtual declaration
    virtual DoubleMatrix getCorrelations(int index);
	virtual double getCorrelationIJ(int i, int j, int index);
	virtual int getNumAssets() {return nbAssets;}

protected:
    DoubleMatrixArraySP correlationTerm;
    IntArray            fwdCorrIndexArray;
    int                 nbAssets;
    int                 nbFwdCorrSteps;    
    vector<vector<vector<double> > > correlCoeffsTerm;

};

/*****************************************/
/* EXTENDED SPARSE GAUSS TERM DEPENDENCE */
/*****************************************/

/** Same as Gauss but can cope for case where the 'square root' of the
correlation matrix is supplied (and is not in lower triangular form) */
class MCARLO_DLL ExtendedSparseGaussTerm: public GaussTerm {
public:
	/** THIS IS HACKY -- NEEDS TO BE EXTENDED FOR SPARSEFWDCORRS */
	
	/** constructor */
	ExtendedSparseGaussTerm(const vector<SparseDoubleMatrixSP>& sqrtCorrs,
                            IntArray                            fwdCorrIndexArray);

	DoubleMatrix getCorrelations(int index); 

	double getCorrelationIJ(int i, int j, int index);

	void correlateSeries( DoubleMatrix& noise, int pathIdx); 

protected:
	vector<SparseDoubleMatrixSP>  correlCoeffsSparse;    
};

/*******************************/
/* GAUSS TERM DEPENDENCE MAKER */
/*******************************/

class MCARLO_DLL DependenceMakerGaussTerm : public DependenceMaker,
                                            virtual public IRevertTypeConvert {
public:
    static CClassConstSP const TYPE;

    virtual void validatePop2Object();

    /** modifies MDF if necessary */
    virtual void modifyMarketDataFetcher(MarketDataFetcherSP mdf); 

    void init();

    void checkInputs() const;

    static const string DEFAULT_FWD_CORR_FREQ;

    class MCARLO_DLL Support : virtual public DependenceMaker::ISupport {
    public:
        virtual DateTimeArray getSimDates() const = 0;
        virtual DoubleMatrix getFwdVarAtDates(bool interpolateAtmFwd) const = 0;
        virtual const IMultiFactors* getMultiFactors() const = 0;
    };

    /** helper methods */
    static void createFwdCorrIdxAndDatesArray(
        const DateTimeArray&           timeline, // valuedate + simdates
        const FwdCorrModeSP            fwdCorrMode,                                   
        const TimeMetricArray&         timeMetricArray,
        IntArray&                      fwdCorrIndexArray,     // output
        DateTimeArray&                 fwdCorrDatesArray);    // output                        

    static CDoubleMatrixSP fwdVarAtFwdCorrDates(
        const DateTimeArray&            timeline,
        const IntArray&                 fwdCorrIndexArray,
        const DoubleMatrix&             fwdVarAtSimDates);
    
    DoubleMatrixArraySP computeFwdCorrMatrixArray (
        IntArray&             fwdCorrIndexArray,      // output
        DateTimeArray&        fwdCorrDatesArray,      // output        
        CorrTermDataSP        corrTermData,
        const IMultiFactors*  mAsset,
        const DoubleMatrix&   fwdVarAtSimDates,
        const DateTimeArray&  timeline) const ;

    /** interface that the instrument must implement */
    virtual DependenceSP createDependence(const DependenceMaker::ISupport* support) const;     

    /** invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** for reflection */
    DependenceMakerGaussTerm();
    DependenceMakerGaussTerm(const CClassConstSP& clazz);

    static IObject* defaultDependenceMakerGaussTerm();

    // for the IMS interface (see IRevertTypeConvert)
    virtual IObjectSP revert(const string& interfaceType) const;

    FwdCorrModeSP getFwdCorrModeUsed() const;

protected:
    /** registered fields */
    double                  eigenValueFloor;
    double                  maxSqError;
    string                  fwdCorrMode; 
    FwdCorrModeWrapperSP    fwdCorrModeDetails;
    /** transient fields */
    FwdCorrModeSP           fwdCorrModeUsed;
};

typedef smartConstPtr<DependenceMakerGaussTerm> DependenceMakerGaussTermConstSP;
typedef smartPtr<DependenceMakerGaussTerm> DependenceMakerGaussTermSP;

/***********************************/
/* GAUSS TERM SRM DEPENDENCE MAKER */
/***********************************/
/** dummy class, since already available in IMS */
class MCARLO_DLL DependenceMakerGaussTermSrm : public DependenceMakerGaussSrm {
public:
    static CClassConstSP const TYPE;    

private:
    DependenceMakerGaussTermSrm();
    static void load(CClassSP& clazz);
    static IObject* defaultDependenceMakerGaussTermSrm();
};

DRLIB_END_NAMESPACE
#endif // DEPENDENCE_GAUSS_TERM_HPP
