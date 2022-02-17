//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids 
//
//   Description : InstrumentCollection with weights and caching
//
//   Author      : Linus Thand 
//
//   Date        : 10 July 2006
//
//----------------------------------------------------------------------------


#ifndef EDG_WEIGHTED_INSTRUMENT_COLLECTION_H
#define EDG_WEIGHTED_INSTRUMENT_COLLECTION_H

#include "edginc/config.hpp"
#include "edginc/Model.hpp"
#include "edginc/CInstrumentCollection.hpp"
#include "edginc/WeightedInstrumentTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/** An IInstrumentCollection similar to ArrayInstrumentCollection, but with
 *  caches on each of the instruments, and weights associated with each cache.
 *  The weights are necessary when solving for reinvestment (e.g. a managed trade)
 *  through ImpliedScalarShiftMulti, since this collection implements 
 *  WeightedInstrumentShift. The caches are here so that each instrument in the
 *  collection doesn't have to be repriced for each iteration of the weight
 *  solving, which isn't necessary since the weight is only a scaling factor 
 *  that is applied on top  of the result.
 */


class RISKMGR_DLL WeightedInstrumentPriceCache : public CObject {
 public:
    static CClassConstSP const TYPE;
    void fieldsUpdated(const CFieldArray& fields);
    CResultsSP Price(IModelSP model, CControlSP c, const double weight);
    CInstrumentSP getInstrument(void);
    CInstrumentConstSP getInstrument(void) const;
    WeightedInstrumentPriceCache(CInstrumentSP inst);
    WeightedInstrumentPriceCache();
    ~WeightedInstrumentPriceCache();
 private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    WeightedInstrumentPriceCache(const WeightedInstrumentPriceCache& rhs); 
    WeightedInstrumentPriceCache& operator=(const WeightedInstrumentPriceCache& rhs);
    void markDirty(void);

    /** Fields **/
    CInstrumentSP inst;
    bool cacheDirty;
    CResultsSP cachedResults;

    CControlSP localCtrl;
    CControlSP originalCtrl;
    IModelSP localMdl;
    IModelSP originalMdl;
};

FORWARD_DECLARE(WeightedInstrumentPriceCache)


class RISKMGR_DLL WeightedInstrumentCollection : 
    public CInstrumentCollection,
    virtual public IRestorableWithRespectTo<WeightedInstrumentTweak>                             
{
 public:

    static CClassConstSP const TYPE;
    WeightedInstrumentCollection(); 
    ~WeightedInstrumentCollection();
    virtual void validatePop2Object();

    /**
    * Number of instruments in the collection
    *
    * May be 0.
    */

    int size() const;

    /**
     * An individual instrument from the collection
     */

    InstrumentSP operator [](int i);

    /**
     * An individual instrument from the collection
     */

    InstrumentConstSP operator [](int i) const;

    //@}

    /**
     * Housekeeping
     */

    //@{

    /**
     * Throw an exception if one of the instruments isn't in a valid state
     */

    virtual void Validate();

    /**
     * Call IInstrument::GetMarket() on all the instruments
     */

    void GetMarket(const IModel*, const CMarketDataSP);

    /**
     * Today
     */

    DateTime getValueDate() const;

    /**
     * Curve to use for discounting payments/values
     *
     * Taken from the first instrument in the underlying array; throws an
     * exception if it's empty
     */

    string discountYieldCurveName() const;

    /**
     * Max of IInstrument::endDate() over all instruments
     */

    DateTime endDate(const Sensitivity*) const;

    /**
     * Scale results of a pricing run
     */
    void scaleOutputs(CControlSP control, CResultsArraySP results);

    // ITweakableWithRespectTo<InstrumentWeight>
    virtual string sensName(const WeightedInstrumentTweak* tag) const;
    virtual TweakOutcome sensShift(const PropertyTweak<WeightedInstrumentTweak>& shift);
    // IRestorableWithRespectTo<InstrumentWeight>
    virtual void sensRestore(const PropertyTweak<WeightedInstrumentTweak>& shift);

    void Price(IModel *model, CControl *control, CResultsArraySP resultss);
 private:
    WeightedInstrumentCollection(const WeightedInstrumentCollection& rhs); 
    WeightedInstrumentCollection& operator=(const WeightedInstrumentCollection& rhs);
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
   
    /* Fields */

    // Caches for each of the instruments
    WeightedInstrumentPriceCacheArraySP priceCaches; 
    CInstrumentArraySP instruments;
    DoubleArray weights; 
};

FORWARD_DECLARE(WeightedInstrumentCollection);

DRLIB_END_NAMESPACE
#endif

