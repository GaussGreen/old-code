//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IndexSkew.hpp
//
//   Description : 
//
//   Author      : Antoine Gregoire
//
//   Date        : March 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDR_INDEX_SKEW_HPP
#define EDR_INDEX_SKEW_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/SkewSurface.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/MappingFunction.hpp"
#include "edginc/BetaSkewParallel.hpp"
#include "edginc/TweakableWith.hpp"
#include "edginc/StrikeMappingTweak.hpp"
#include "edginc/BCStrikeMappingTweakBase.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL IndexSkew : public MarketObject,
                  public virtual BetaSkewParallel::IRestorableShift,
                  public virtual TweakableWith<StrikeMappingTwk>,
                  public virtual BCStrikeMappingTweakBase::IShift
{
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Destructor */
    ~IndexSkew();
    
    /** Called immediately after object constructed */
    virtual void validatePop2Object();
    
    /** Rebuild wide spreads */
    virtual void fieldsUpdated(const CFieldArray& fields);
    
    /** Populates the object with the market data that this object needs */
    virtual void getMarket(const IModel* model, const MarketData* market);
    
    /** 
     * Returns the name of this object
     * */
    virtual string getName() const;
    
    /** Calculates the default rates associated to this index */
    DefaultRatesSP defaultRates() const;
    
    /** Calculates implied par spread at given date */
    void impliedSpreadsAndDurations(
        const YieldCurveConstSP discount,
        const DateTimeArray& dates,
        DoubleArray& impliedSpreads, /* (Output) */
        DoubleArray& durations) const; /* (Output) */
    
    /** Returns the recovery associated to this index */
    double getRecovery() const;
    
    /**
     * Returns interpolated skew of given type (normal, fast...)
     * at point (strike, date)
     * */
    double getSkew(
        double strike, const DateTime& date, SkewSurface::SkewType type) const;

    /**
     * Returns interpolated skew of given type (normal, fast...)
     * at point (strike, date, spreadRatio)
     * */
    double getSkew(
        double strike, const DateTime& date, double spreadRatio, SkewSurface::SkewType type) const;

    /**
     * Returns "neighbours" of the given point P{strike, date}
     * i.e. points needed to interpolate at P
     * */
    BetaSkewGridPointSet getNeighbours(double strike, 
                                       const DateTime& date, 
                                       const int numberOfNeighbours) const;

    /**
     * Returns "neighbours" of the given point P{strike, date, spreadRatio}
     * for surface "surfaceName"
     * i.e. points needed to interpolate at P
     * */
    BetaSkewGridPointSet getNeighbours(double strike, 
                                       const DateTime& date, 
                                       string surfaceName, 
                                       double spreadRatio,
                                       const int numberOfNeighbours) const;

    /** Returns all points defining the surface */
    BetaSkewGridPointArrayConstSP getAllPoints(
        string skewSurfaceName, DoubleArrayConstSP spreadRatios) const;

    /** Returns all skew surface names */
    StringArraySP getSkewSurfaceNames() const;

    /** Returns the strike mapping parameter */
    double getStrikeMapping() const;
    
    /** Returns the squeeze function*/
    const IMappingFunctionSP getSqueeze() const;

    /** Returns the historical beta */
    double getHistoricalBeta() const;
    
    /** Returns 'true' when offSpreadBetaSkews and spreadRatios are populated */
    bool useOffSpreads() const;
    
    /** Strike mapping tweak support */
    virtual string sensName(StrikeMappingTwk* shift) const;

    /** Strike mapping tweak support */
    virtual bool sensShift(StrikeMappingTwk* shift);

    /** Base Correlation Strike mapping tweak support */
    virtual string sensName(BCStrikeMappingTweakBase* shift) const;

    /** Base Correlation Strike mapping tweak support */
    virtual bool sensShift(BCStrikeMappingTweakBase* shift);

    /**
     * Creates a bootstrapper corresponding to this IInstanceIDBootstrappable
     * [Implements IInstanceIDBootstrappable]
     * */
    virtual IInstanceIDBootstrapperSP createBootstrapperFromSkewSurface(
		CFieldConstSP fieldToCalibrate, 
        string bootstrapType) const;

	// -----------------------------------
	// METHODS FOR BETA SKEW SENSITIVITIES
	// -----------------------------------

	virtual string sensName(BetaSkewParallel* shift) const;
    virtual bool   sensShift(BetaSkewParallel* shift);
    virtual void   sensRestore(BetaSkewParallel* shift);

private:

    /**
     * Returns index of skew surface with name surfaceName in
     * offSpreadBetaSkews (or -1 if not found)
     * */
    int getSurfaceIndex(string surfaceName) const;
    
    /** Check if the surface with index "surfaceIdx" is in the neighbourhood of "spreadRatio" */
    bool isNeighbour(double spreadRatio, int surfaceIdx) const;
    
    /** Builds homogeneous wide spreads surfaces */
    void buildWideSpreads();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
    
    /** Only build instances of that class using reflection */
    IndexSkew();
    
    /** Default constructor */
    static IObject* defaultConstructor();
    
    // ----------------
    // MANDATORY FIELDS
    // ----------------

    /** 
     *  Beta skews (normal and fast) matrix for that index,
     *  implicitly equivalent of a beta skew surface for a 
     *  spread ratio equal to 1 */
    SkewSurfaceSP betaSkews;
    
    /** Historical beta */
    double historicalBeta;
    
    /** Mapping rule parameter (also called "q"), range [0,1] */
    double strikeMapping;

    /** Squeeze function : relative tweak apllied to the skew surface */
    IMappingFunctionSP squeeze;
    
    /** CDS par spreads curve associated to the index */
    ICDSParSpreadsWrapper cdsParSpreads;

    // ---------------
    // OPTIONAL FIELDS
    // ---------------

    /** name of the index */
    string name;
    
    /**
     * Array of skew surfaces for different spread ratios
     * equivalent of betaSkews for different average spreads 
     * */
    SkewSurfaceArraySP offSpreadBetaSkews;
    
    /** Spread Ratio corresponding to each skew surface */
    DoubleArraySP spreadRatios;
    
    // --------------------------
    //  LOCAL FIELD (NOT EXPOSED)
    // --------------------------
    
    // internal flag set to true when offSpreadBetaSkews and spreadRatios are populated
    bool useOffSpreadsFlag;
    // internal array of spreadRatios including 1
    DoubleArraySP allSpreadRatios;
    // internal array of SkewSurface including the original surface
    SkewSurfaceArraySP allSpreadBetaSkews;
};

// Support for smart pointers
typedef smartPtr<IndexSkew> IndexSkewSP;
typedef smartConstPtr<IndexSkew> IndexSkewConstSP;
#ifndef QLIB_INDEXSKEW_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<IndexSkew>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IndexSkew>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<IndexSkew>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IndexSkew>);
#endif

// Support for arrays
typedef array<IndexSkewSP, IndexSkew> IndexSkewArray;
typedef smartPtr<IndexSkewArray> IndexSkewArraySP;
typedef smartConstPtr<IndexSkewArray> IndexSkewArrayConstSP;
#ifndef QLIB_INDEXSKEW_CPP
EXTERN_TEMPLATE(class MARKET_DLL array<IndexSkewSP _COMMA_ IndexSkew>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<IndexSkewArray>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IndexSkewArray>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL array<IndexSkewSP _COMMA_ IndexSkew>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<IndexSkewArray>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IndexSkewArray>);
#endif

// Support for wrapper class
typedef MarketWrapper<IndexSkew> IndexSkewWrapper;
#ifndef QLIB_INDEXSKEW_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<IndexSkew>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<IndexSkew>);
#endif

/** Specialisations of arrayObjectCast */
template <> class MARKET_DLL arrayObjectCast<IndexSkewWrapper>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const IndexSkewWrapper& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(IndexSkewWrapper& value);

    /** Turns the IObjectSP into an IndexSkewWrapper */
    static IndexSkewWrapper fromIObject(IObjectSP& value);
};

// Support for wrapper arrays
typedef array<IndexSkewWrapper,IndexSkewWrapper> IndexSkewWrapperArray;
typedef smartPtr<IndexSkewWrapperArray> IndexSkewWrapperArraySP;
#ifndef QLIB_INDEXSKEW_CPP
EXTERN_TEMPLATE(class MARKET_DLL array<IndexSkewWrapper _COMMA_ IndexSkewWrapper>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IndexSkewWrapperArray>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL array<IndexSkewWrapper _COMMA_ IndexSkewWrapper>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IndexSkewWrapperArray>);
#endif

DRLIB_END_NAMESPACE

#endif
