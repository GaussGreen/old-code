//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : CDSVolCubeMultiQSmile.hpp
//
//   Description : CDS Option CDSVolCube with a Multi-Q Smile with
//                 an arbitrary number of qs.
//
//   Author      : Charles Morcom
//
//   Date        : 21 December 2005
//
//
//----------------------------------------------------------------------------

#ifndef QR_CDSVOLCUBEMULTIQSMILE_HPP
#define QR_CDSVOLCUBEMULTIQSMILE_HPP




#include "edginc/Object.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/CDSVolCube.hpp"
#include "edginc/CDSVolATMMatrix.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/Model.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_WRAPPER(CDSVolCubeMultiQSmile)
FORWARD_DECLARE(MultiQDistribution)
FORWARD_DECLARE(MarketData)
FORWARD_DECLARE(IVolProcessed)
FORWARD_DECLARE(CVolRequest)
FORWARD_DECLARE(Expiry)
FORWARD_DECLARE_WRAPPER(ICDSParSpreads)
//FORWARD_DECLARE_WRAPPER(CDSVolATMMatrix)
FORWARD_DECLARE(CDSVolCubeMultiQLocalSmile)

 

/**Class defining CDS option volatility cube for an implied vol BS smile. */
class MARKET_DLL CDSVolCubeMultiQSmile : public MarketObject,
                              virtual public ICDSVolCube {
public:
    static CClassConstSP const TYPE;

    /** Returns name of vol */
    virtual string getName() const;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Check validity */
    virtual void validatePop2Object();

    //virtual TimeMetricSP getTimeMetric() const;
    //const DateTime& VolCube::getBaseDate() const;

    /**Convert to an equivalent BS implied representation with the same ATM vols.*/
    //CDSVolCubeBSImpliedSmileSP convertToCDSVolCubeBSImplied() const;

    /*=========================================================================
     * ICDSVol interface methods
     *=======================================================================*/
    virtual IVolProcessed* getProcessedVol(const CVolRequest*    volRequest,
                                           const ICDSParSpreads* cds) const;

    /**Creates a weighted composite smile from a vector of nearby smiles and weights.*/
    virtual void combineSmiles(CDSVolCubeMultiQLocalSmile* compositeSmile, 
        const std::vector<CDSVolCubeMultiQLocalSmile*>&    closestSmiles, 
        const std::vector<double>&                     smileWeights)  const;

    virtual ~CDSVolCubeMultiQSmile();
protected:
    CDSVolCubeMultiQSmile();
    CDSVolCubeMultiQSmile(const CClassConstSP& clazz);

private:
    /*=========================================================================
     * DATA FIELDS
     *=======================================================================*/
    /**Name of object*/
    string name;
    /**This also inherits name, base date, etc. from ATM vol matrix */
    CDSVolATMMatrixWrapper atmMatrix;
    /**Array of smiles*/
    CDSVolCubeMultiQLocalSmileArray localSmiles;
  double meanReversion;

    //mutable CDSVolCubeBSLocalSmileArray smileArray;
    mutable CDSVolCubeMultiQLocalSmileArray smileArray; // $unregistered
    // number of option expiries
    mutable int opN; // $unregistered
    // number of ul expiries
    mutable int ulN; // $unregistered

    static void load(CClassSP& clazz);

    friend class CDSVolCubeMultiQSmileHelper;
};


/**Represents a parameterized multi-q smile for a given expiry/maturity tenor 
   point in a vol matrix.*/
class MARKET_DLL CDSVolCubeMultiQLocalSmile : public CObject {
public:
    virtual ~CDSVolCubeMultiQLocalSmile();
    static CClassConstSP const TYPE;
protected:
    CDSVolCubeMultiQLocalSmile();
    static IObject* defaultCDSVolCubeMultiQLocalSmile();
private:
    /*=========================================================================
     * DATA MEMBERS
     *=======================================================================*/
    /**Which expiry in the CDSVolATMMatrix optionExpiries is this smile attached to? 
       If not an exact match, then error.*/
    ExpirySP optionExpiry;
    /**Which expiry in the CDSVolATMMatrix ulExpiries is this smile attached to? 
       If not an exact match, then error.*/
    ExpirySP underlyingExpiry;
    /**Smile pesudo-deltas - N(K/F) - that define points where q changes.
       These must be strictly increasing and 0<delta<=1
       The corresponding q applies for deltas <= the q value. The last delta value
       is ignored: a safe, meaningful value to use is 1.*/   
    DoubleArray deltaAry;
    /**Qs on regions between pseudo-deltas. 
       These are "trader" qs rather than "QR" qs: 0 is lognormal; 1 is normal.*/
    DoubleArray qAry;

    // NON-INTERFACE FIELDS ///////////////////////////////////////////////////
    

    static void load(CClassSP& clazz);

    friend class CDSVolCubeMultiQSmile;
};

DRLIB_END_NAMESPACE
#endif
