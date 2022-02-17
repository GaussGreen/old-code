/**
 * @file VanillaGridInstrumentCollection.hpp
 */

#ifndef QLIB_VanillaGridInstrumentCollection_H
#define QLIB_VanillaGridInstrumentCollection_H

#include "edginc/CInstrumentCollection.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Vanilla.hpp" 

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(VanillaGrid)
FORWARD_DECLARE(VanillaGridInstrumentCollection)

/**
 * An IInstrumentCollection comprising a matrix of Vanilla instruments with
 * given strikes and maturities
 *
 * Enables a whole lot of vanillas to be priced/greeked at once using
 * QLib's "vectorised pricing" facility (see IInstrumentCollection for
 * an overview).
 *
 * It's actually a wrapper round VanillaGrid, essentially just marshalling
 * the results of pricing the VanillaGrid "instrument" into an array of
 * normal-looking Results rather than a single Results with a matrix
 * output request.
 */

class PRODUCTS_DLL VanillaGridInstrumentCollection: public CInstrumentCollection,
                                       public virtual IMCIntoProduct /* FIXME so that we can call memoryUsagePerPath */ {

    VanillaGridInstrumentCollection();
    static IObject* defaultVanillaGridInstrumentCollection();
    static void load(CClassSP& clazz);

public:

    static CClassConstSP TYPE;

private:

    VanillaGridSP grid;

    mutable vector<pair<int, int> > _indices; // $unregistered
    const vector<pair<int, int> >& indices() const;
    mutable DoubleArrayConstSP _weights; // $unregistered

public:
    
    /**
     * A collection of Vanilla instruments with given strikes and maturities
     */

    VanillaGridInstrumentCollection(VanillaGridSP grid);

    /**
     * The grid
     */

    VanillaGridConstSP getGrid() const;

    /**
     * Array access methods
     */

    //@{

    /**
     * Number of vanillas in the grid
     */

    int size() const;

    /**
     * An individual Vanilla from the grid
     */

    InstrumentSP operator [](int i);

    /**
     * Weights assigned to each Vanilla in the grid
     */

    DoubleArrayConstSP weights() const;

    //@}

    /**
     * Housekeeping
     */

    //@{

    /**
     * Validate after construction from xml etc.
     */

    void validatePop2Object();

    /**
     * Throw an exception if the underlying VanillaGrid is in an invalid state
     */

    void Validate();

    /**
     * Call IInstrument::GetMarket() on the underlying VanillaGrid
     */

    void GetMarket(const IModel* model, const CMarketDataSP market);

    /**
     * Today
     */

    DateTime getValueDate() const;

    /**
     * Curve for discounting payoff
     */

    string discountYieldCurveName() const;

    /**
     * Date after which none of the Vanillas are sensitive to market data (with
     * respect to the supplied sensitivity)
     */

    DateTime endDate(const Sensitivity* sensControl) const;

    /**
     * IMCIntoProduct implementation
     *
     * Purpose of this is to give MonteCarlo::RunMulti() access to
     * VanillaGridMCSV::storagePerPath() --- which actually just returns 0.
     */

    IMCProduct* createProduct(const MonteCarlo* model) const;

    //@}

    /**
     * Pricing
     */

    //@{

    /**
     * Price all the Vanilla's in the grid
     *
     * Calls IModel::Price() on the VanillaGrid "instrument", and marshals the
     * OPTION_PRICE output request into the corresponding @a resultss entries.
     */

    void Price(IModel *model,
               CControl *control,
               CResultsArraySP resultss);

    /**
     * Scale results of a pricing run
     *
     * No-op.
     */

    void scaleOutputs(CControlSP control, CResultsArraySP results);

    //@}
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
