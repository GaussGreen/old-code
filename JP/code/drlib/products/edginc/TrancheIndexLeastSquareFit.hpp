//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : TrancheIndexLeastSquareFit.hpp
//
//   Description : Objective function used to calibrate CDO models
//
//   Author      : Antoine Gregoire
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#ifndef TRANCHE_INDEX_LEAST_SQUARE_FIT_HPP
#define TRANCHE_INDEX_LEAST_SQUARE_FIT_HPP

#include "edginc/Calibrator.hpp"
#include "edginc/CDOQuotesBootstrapperTimeOnly.hpp"
#include "edginc/ICDOFineGrid.hpp"
#include "edginc/ICDOQuotesGenerator.hpp"
#include "edginc/ICDOQuotesWeights.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Objective function used to calibrate CDO models:
 * - Base correlations (uses bootstrapping)
 * - CID [not implemented]
 * 
 * 
 * This class is responsible for building "on the fly" the CDO instruments
 * given the CDOParSpreads market objects.
 * 
 * Additionaly, we can specify a "fine grid" [eg: a set of (strike, maturity)
 * points associated to a model] to extrapolate the existing tranche quotes
 * for all fine grid points. Then the extrapolated tranche quotes will be used
 * as input for the calibration.
 * 
 * */
// We derive from Calibrator::ObjFuncLeastSquare to allow more flexibility for
// future models (eg: CID).
// Base correlation calibration will just use N (=nb beta skews to calibrate)
// calls to calcValue(funcVals) where "funcVals" is a CDoubleArray of size 1.
// 
class PRODUCTS_DLL TrancheIndexLeastSquareFit :
	public Calibrator::ObjFuncLeastSquare,
    public Calibrator::ObjFunc::IGenericBootstrappable 
{
public:	
	/** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Destructor */
    virtual ~TrancheIndexLeastSquareFit();

    /**
     * Give a chance to do s'thing with the market
     * [Implements Calibrator::ObjFunc]
     * */
    virtual void getMarket(MarketData* market);
    
    /**
     * Returns an IObect that contains all the IAdjustable objects 
     * that the calibrator is to operate upon
     * [Implements Calibrator::ObjFunc]
     * */
	virtual IObjectSP getAdjustableGroup();
	
	/**
	 * Returns the number of functions
	 * [Implements Calibrator::ObjFuncLeastSquare]
	 * */
    virtual int getNbFuncs() const;

	/**
	 * Calculates the values of the objective functions
	 * [Implements Calibrator::ObjFuncLeastSquare]
	 * */
    virtual void calcValue(CDoubleArray& funcvals) const;

	/** override for calcValue() so that can put a floor on the objective function */
	double calcValue() const;


    /** main method used to retrieve a bootstrapper capable of bootstrapping */
    virtual IBootstrapperSP getBootstrapper(
        const Calibrator::InstanceIDArray & ids);

    

	/// Bootstrap types
	/** bootstrap in time and strike dimension */
	static const string TIME_STRIKE;

	/** bootstrap in time dimension only (global calibration in strike dim) */
	static const string TIME;

    /** set the current instruments valued in calcValue. Used for bootstrapping */
    void setCurrentInstruments(CInstrumentArraySP inst);

    /** access to bootstrapType */
    string getBootstrapType() const
    {
         return bootstrapType;
    }
    /** get cdoQuotesIterator */
     CDOQuotesBootstrapperSP  getQuotesIterator() const
     {
         return cdoQuotesIterator;
     }

  // Return model
  IModelSP Model() { return model; }  

  // Return current guess for parameters
  CreditEngineParametersConstSP  CurrentGuess() { return currentGuess; }

  // Return instruments
  CInstrumentArraySP GetInstruments() { return instruments; }

  // Return controls
  CControlArraySP GetControls() { return controls; }

private:
    /** types of objective function */
    // MTM^2
    static const string OBJ_FUNC_MTM;
    // (MTM / FeeLeg)^2
    static const string OBJ_FUNC_MTM_FEE_LEG;
	
    /** floor for objective function - sets maximum precision */
	static const double DEFAULT_OBJ_FUNC_FLOOR;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
    
    /** Private constructor (only build instances of that class using reflection) */
    TrancheIndexLeastSquareFit();
    
    /** Default constructor */
    static IObject* defaultConstructor();
	
	// ------
	// FIELDS
	// ------
	
	/** CDO model [Mandatory] */
	IModelSP model;
	
	/** Index tranche quotes [Mandatory] */
	ICDOQuotesGeneratorWrapper cdoQuotes;

	/** Index tranche weights [Optional] */
	ICDOQuotesWeightsWrapper cdoQuotesWeights;


	/** Defines the "fine grid" [Optional] */
	ICDOFineGridSP fineGrid;

	/** Defines the "fine grid" [Optional] */
  CreditEngineParametersSP initialGuess;

	/** iterator for CDO quotes */
	CDOQuotesBootstrapperSP cdoQuotesIterator;
	
    /** Instruments corresponding to CDO quotes [transient - but tweakable] */
	// these are used only for global calibration
    CInstrumentArraySP instruments;

    /**
     * Temporary field (for bootstrapping only):
     * contains status information on the current bootstrapping
     * for each "InstanceID"
     * [transient]
     * */
    IBootstrapperArraySP idBootstrappers;
    
    /** "control" used for pricing [transient] */
    CControlArraySP controls;
    
    /**
     * Flag to know if we are bootstrapping (=True) or not (=False).
     * Default is false.
     * [transient]
     * */
     bool doBootstrapping;

	 /**
	 * bootstrap type
	 */
	 string bootstrapType;

     /** flag to set how we value the objective function [optional] default = OBJ_FUNC_MTM*/
     string objFuncType;

	 /** floor for objective function */
	 double objFuncFloor;

	 /* Transient field */
	 mutable CDoubleArraySP funcs;

  // workspace 
  mutable CreditEngineParametersConstSP currentGuess;
};

DECLARE(TrancheIndexLeastSquareFit);

DRLIB_END_NAMESPACE

#endif /*TRANCHE_INDEX_LEAST_SQUARE_FIT_HPP*/

