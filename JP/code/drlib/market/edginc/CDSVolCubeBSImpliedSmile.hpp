//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : CDSVolCubeBSImpliedSmile.hpp
//
//   Description : CDS Option CDSVolCube with a Black-Scholes implied-volatility
//                 smile
//
//   Author      : Charles Morcom
//
//   Date        : 21 December 2005
//
//----------------------------------------------------------------------------

#ifndef QR_CDSVOLCUBEBSIMPLIEDSMILE_HPP
#define QR_CDSVOLCUBEBSIMPLIEDSMILE_HPP

#include "edginc/Object.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/CDSVolCube.hpp"
#include "edginc/CDSVolATMMatrix.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Array.hpp"
#include "edginc/Model.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/CRVolParallel.hpp"
#include "edginc/CRVolPointwise.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_WRAPPER(CDSVolCubeBSImpliedSmile)
FORWARD_DECLARE(MultiQDistribution)
FORWARD_DECLARE(MarketData)
FORWARD_DECLARE(TimeMetric)
FORWARD_DECLARE(CDSVolCubeMultiQSmile)
FORWARD_DECLARE(IVolProcessed)
FORWARD_DECLARE(CVolRequest)
FORWARD_DECLARE_WRAPPER(ICDSParSpreads)
//FORWARD_DECLARE_WRAPPER(CDSVolATMMatrix)
FORWARD_DECLARE(CDSVolCubeBSLocalSmile)
FORWARD_DECLARE(Expiry)
FORWARD_DECLARE(ExpiryPair)

/**Class defining CDS option volatility cube for an implied vol BS smile. */
class MARKET_DLL CDSVolCubeBSImpliedSmile : public MarketObject,
                                 virtual public ITweakableWithRespectTo<CRVolParallel>,
                                 virtual public ITweakableWithRespectTo<CRVolPointwise>,
                                 virtual public ICDSVolCube {
public:
    static CClassConstSP const TYPE;

    /** Returns name of vol */
    virtual string getName() const;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    virtual void validatePop2Object();

    virtual double interpolateSingleSmile(
        double strike,
        double atmForward,
        double t,
        double atmSigma,
        const CDSVolCubeBSLocalSmile& smile) const;

    /**Creates a weighted composite smile from a vector of nearby smiles and weights.*/
    virtual void combineSmiles(CDSVolCubeBSLocalSmile* compositeSmile, 
        const std::vector<CDSVolCubeBSLocalSmile*>&    closestSmiles, 
        const std::vector<double>&                     smileWeights)  const;

    //virtual TimeMetricSP getTimeMetric() const;
    //const DateTime& VolCube::getBaseDate() const;

    /**Convert to a multi-q representation with the same ATM vols, if possible to calibrate.
       Could write this as an implicit cast, instead, but maybe this is better.*/
    //CDSVolCubeMultiQSmileSP convertToCDSVolCubeMultiQ() const;

    /*=========================================================================
     * ICDSVol interface methods
     *=======================================================================*/
    virtual IVolProcessed* getProcessedVol(const CVolRequest*    volRequest,
                                           const ICDSParSpreads* cds) const;

    // Parallel credit vega tweak (sensitivity) support
    virtual string sensName(const CRVolParallel*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<CRVolParallel>&);
    
    // Pointwise credit vega tweak (sensitivity) support
    virtual string sensName(const CRVolPointwise*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<CRVolPointwise>& shift);
    virtual ExpiryPairArrayConstSP sensQualifiers(const CRVolPointwise*) const;

    virtual ~CDSVolCubeBSImpliedSmile();
protected:
    CDSVolCubeBSImpliedSmile();
    CDSVolCubeBSImpliedSmile(const CClassConstSP& clazz);
private:
    /*=========================================================================
     * DATA FIELDS
     *=======================================================================*/
    /**Name of object*/
    string name;
    /**This also inherits name, base date, etc. from ATM vol matrix */
    CDSVolATMMatrixWrapper atmMatrix;
    /**How should the strike numbers in the smile be interpolated? Delta? Ratio? Absolute?
       Should have value one of CDSVolCube::STRIKE_TYPE_ + ABSOLUTE, RATIO, or DELTA.*/
    string smileStrikeType;
    /**How should the volatility numbers in the smile be interpreted? Offset from ATM? 
       Multiplcative or additive? Absolute? Should have value one of
       CDSVolCube::VOL_TYPE_ + ABSOLUTE, ADDITIVE_OFFSET, or MULTIPLICATICE_OFFSET.*/
    string smileVolType;
    /**Array of smile points*/
    CDSVolCubeBSLocalSmileArray localSmiles;

    //mutable CDSVolCubeBSLocalSmileArray smileArray;
    CDSVolCubeBSLocalSmileArray smileArray; // $unregistered
    // number of option expiries
    mutable int opN; // $unregistered
    // number of ul expiries
    mutable int ulN; // $unregistered

    friend class CDSVolCubeBSImpliedSmileHelper;
};

/**Describes a smile at a particular ATM vol point, described by indices into 
   the ATM matrix vol expiry and underlying maturity arrays*/
class MARKET_DLL CDSVolCubeBSLocalSmile : public CObject {
public:
    static CClassConstSP const TYPE;
    virtual ~CDSVolCubeBSLocalSmile();
private:
    CDSVolCubeBSLocalSmile();
    static IObject* defaultCDSVolCubeBSLocalSmile();
    /*=========================================================================
     * DATA MEMBERS
     *=======================================================================*/
    /**Which expiry in the CDSVolATMMatrix optionExpiries is this smile attached to? 
       If not an exact match, then error.*/
    ExpirySP optionExpiry;
    /**Which expiry in the CDSVolATMMatrix ulExpiries is this smile attached to? 
       If not an exact match, then error.*/
    ExpirySP underlyingExpiry;
    /**Strikes interpreted according to the smileStrikeType in the vol cube 
       that contains this - must be the same number as volNumbers, and must be
       strictly ascending.*/
    DoubleArray strikeNumbers; 
    /**Vols interpreted according to the smileVolType in the vol cube that 
       contains this. Must be the same number as strikeNumbers + a few other
       checks conditional on the smileVolType */
    DoubleArray volNumbers;
    static void load(CClassSP& clazz);

    friend class CDSVolCubeBSImpliedSmile;
};

DRLIB_END_NAMESPACE
#endif
