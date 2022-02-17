//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : CreditIndexBasis.hpp
//
//   Description : Holds the information about the basis to be applied to a
//                 credit index
//
//   Author      : Jose Hilera
//
//   Date        : August 30, 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_CREDITINDEXBASIS_HPP
#define QLIB_CREDITINDEXBASIS_HPP

#include "edginc/Object.hpp"
#include "edginc/ExpiryResult.hpp"
#include "edginc/CDSParSpreadsAdjustment.hpp"
#include "edginc/CreditIndexBasisParallel.hpp"
#include "edginc/CreditIndexBasisPointwise.hpp"


DRLIB_BEGIN_NAMESPACE


/** Holds the index basis information required to transform the spreads of the 
 * names in an index into the index spreads.
 * Note that CreditIndexBasis is immutable, ie, once created it cannot be 
 * changed
 *
 ********************   DO NOT ADD ANY NON CONST METHODS  ******************
 *
 * (except for getMarket and validatePop2Object which are essentially 
 * part of the construction of the object) 
 */
class MARKET_DLL CreditIndexBasis : public CDSParSpreadsAdjustment,
                         virtual public ITweakableWithRespectTo<CreditIndexBasisParallel>,
                         virtual public ITweakableWithRespectTo<CreditIndexBasisPointwise> {

public:
    static CClassConstSP const TYPE;

    //values for controlling how index basis may apply
    //used in CreditMetricsModel and in MarketDataFetcherCDS
    //so defined here as a reasonable common source
    static const string APPLY_TO_NONE;              //for completeness, may be useful....
    static const string APPLY_TO_MAPPED_NAMES_ONLY; //non constituents will not be adjusted
    static const string APPLY_TO_ALL_NAMES;         //more than 1 index in the map => error)

    // SPs are copied into the structure, no clone is made
    CreditIndexBasis(const string&       name,
                     ExpiryArrayConstSP  indexExpiries, 
                     const DateTime&     valueDate);

    // Constructor to allow this CreditIndexBasis to be just a dummy adjustment
    CreditIndexBasis(bool useIndexBasis);

    // Construct an invalid index basis. Value of the parameter is ignored
    CreditIndexBasis(ExpiryArrayConstSP  indexExpiries,
                     const DateTime&     valueDate,
                     const bool          invalidate);

    virtual ~CreditIndexBasis();

    void getMarket(const IModel* model, const MarketData* market);
    void validatePop2Object();

    // Returns whether basis is being used or not
    bool used() const;

    // Returns the specified value date
    DateTime getValueDate() const;

    // Returns the name of this object
    virtual string getName() const;

    // Returns the array of expiries where the basis is defined
    ExpiryArrayConstSP getBasisExpiries() const;

    // Returns the whole array of index basis
    DoubleArrayConstSP getBasis() const;

    // Returns information about tenor mapping
    IntArraySP getMapping() const;

    //Utility method for mapping an arbitrary name to an
    //arbitrary index
    IntArraySP getMapping(ExpiryArrayConstSP nameExpiries) const;

    // Returns the basis on a specific expiry, interpolating if required
    // CAUTION: This method assumes the basis has been calculated already.
    double getBasisOnExpiry(const DateTime& expiry,
                            const string&   forName) const;

    // Set the basis weighting associated to a name
    // has no effect if the name is not found (requires prior initialisation/population)
    void setCalculatedForWeights(string name, double weight);

    // Return the weighting of the basis associated to the specified name
    // returns 1.0 if the name is not found
    double getCalculatedForWeight(const string& name) const;

    // Parallel tweak (sensitivity) support
    virtual string sensName(const CreditIndexBasisParallel*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<CreditIndexBasisParallel>&);

    // Pointwise tweak (sensitivity) support
    virtual string sensName(const CreditIndexBasisPointwise*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<CreditIndexBasisPointwise>& shift);
    virtual ExpiryWindowArrayConstSP sensQualifiers(const CreditIndexBasisPointwise*) const;

    /** Validation */
    virtual bool expiriesEqual(const ICDSBootstrappableWrapper& cdsToAdjust) const;

    //--------------------------------
    // CDSParSpreadsAdjustment methods
    //--------------------------------

    /** Validation */
    virtual bool expiriesValid(const ICDSBootstrappableWrapper& cdsToAdjust) const;

    /** Returns the expiries themselves */
    virtual ExpiryArrayConstSP getAdjustmentExpiries() const;

    /** Returns the adjusted recovery */
    virtual double getAdjustedRecovery (const ICDSBootstrappableWrapper& cdsToAdjust) const;

    /** Returns the adjusted par spread curve */
    virtual DoubleArrayConstSP getAdjustedParSpreads(const ICDSBootstrappableWrapper& cdsToAdjust) const;

    /** Returns a single adjusted rate, corresponding to expiry */
    virtual double getAdjustedParSpread(const ICDSBootstrappable* cdsToAdjust,
                                        ExpiryConstSP expiry) const; 

    /** Tweaking support */
    virtual const DoubleArrayConstSP adjustMaxTweakSizes(const ICDSBootstrappableWrapper& cdsToAdjust,
                                                         DoubleArrayConstSP maximumAllowed,
                                                         IExpiryRiskPropertyConstSP wrt) const;

protected:
    // Mechanism to map an arbitary date onto the index tenors
    // for interpolation
    bool mapDate(const DateTime& theDate,                // (I) the date to map
                 int&            index,                  // (O) the mapped index
                 bool&           extrapolateFlat) const; // (O) whether that requires flat extrapolation

    // Once mapped, we can interpolate
    double interpolate(const DateTime& date,
                       const int       index,
                       const bool      extrapolateFlat) const;

    // Set up the names for which this basis has been calculated
    // and initialise the weightings to 1.0
    void initialiseCalculatedForData();

    DoubleArraySP   basis;        // The actual basis. It is "protected" because 
                                  // derived classes need access to it

    DateTimeArraySP expiryDates;  // The expiries converted to dates as of value date
                                  // do once for more efficient interpolation

    IntArraySP      tenorMapping; // How the single name tenors were mapped to the
                                  // index tenors. Useful for debugging.
                                  // Retrieveable via addin method BASIS_GET_MAPPING

private:
    // For Reflection
    CreditIndexBasis();
    static void load (CClassSP& clazz);
    static IObject* defaultCreditIndexBasis();

    // Fields
    string               name;         // Name of this market object
    DateTime             valueDate;
    ExpiryArrayConstSP   expiries;     // The tenors the basis applies to
                                       // NB No offset or bdc adjustments

    StringArraySP        calculatedForNames;    //we reweight the basis for adjusted names 
    DoubleArraySP        calculatedForWeights;  //that fail to bootstrap so that the basis
                                                //can applied successfully after the event
                                                //Therefore we retain this weighting information
                                                //along with the basis

    // Transient fields
    bool validForUse;        /* In some cases, we may want to provide an
                              * instance that should not be used
                              * Accessor methods should check this field
                              * and throw if false */

    bool useIndexBasis;      /* If false, this object is a "dummy" basis adjustment
                              * in the sense that it is equivalent to making no 
                              * adjustment at all, ie, as if the basis was zero */
};


typedef smartConstPtr<CreditIndexBasis> CreditIndexBasisConstSP;
typedef smartPtr<CreditIndexBasis>      CreditIndexBasisSP;
typedef MarketWrapper<CreditIndexBasis> CreditIndexBasisWrapper;

DRLIB_END_NAMESPACE

#endif
