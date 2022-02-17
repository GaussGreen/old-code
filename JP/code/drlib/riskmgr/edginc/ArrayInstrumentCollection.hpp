/**
 * @file ArrayInstrumentCollection.hpp
 */

#ifndef DRLIB_ArrayInstrumentCollection_H
#define DRLIB_ArrayInstrumentCollection_H

#include "edginc/Instrument.hpp"
#include "edginc/CInstrumentCollection.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ArrayInstrumentCollection)

/**
 * Concrete array implementation of IInstrumentCollection
 *
 * QLib's "vectorised pricing" facility allows pricing/greeking of multiple
 * instruments in one pass.  The list of instruments to be priced is defined by
 * the IInstrumentCollection argument to MultiRiskMgrInterface --- see
 * IInstrumentCollection class documentation for an overview.
 *
 * ArrayInstrumentCollection is an implementation backed by an explicit array
 * of individual instruments.  Specialised implementations (such as
 * VanillaGridInstrumentCollection) can be defined to support optimised
 * pricing of structured sets of instruments.
 */

class RISKMGR_DLL ArrayInstrumentCollection: public CInstrumentCollection,
                                 public virtual ISensitiveStrikes {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP TYPE;

private:

    CInstrumentArraySP instruments;

    ArrayInstrumentCollection();
    static IObject* defaultArrayInstrumentCollection();

public:

    /**
     * A collection of instruments defined by an underlying array
     *
     * @param instruments   May not be empty.  It isn't copied, so don't
     *                      go changing it afterwards.
     */

    ArrayInstrumentCollection(CInstrumentArraySP instruments);

    ~ArrayInstrumentCollection();

    /**
     * Array access methods
     *
     * These are just passed to the access methods on the underlying array
     */

    //@{

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

    void Validate();

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
     * Union of IInstrument::getSensitiveStrikes() over all instruments
     */

    DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                      const IModel* model);

    /**
     * True unless one of the instruments is ISensitiveStrikes and
     * its ISensitiveStrikes::avoidVegaMatrix() is false.
     */

    bool avoidVegaMatrix(const IModel* model);

    //@}

    /**
     * Pricing
     */

    //@{

    /**
     * Price all the instruments in the array, putting the results for each in
     * the corresponding entry in @a resultss
     *
     * Calls IModel::PriceMulti() which by default simply uses IModel::Price()
     * to price each instrument separately.
     */

    void Price(IModel *model, CControl *control, CResultsArraySP resultss);

    /**
     * Scale results of a pricing run
     */

    void scaleOutputs(CControlSP control, CResultsArraySP results);

    //@}
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
