/**
 * @file InstrumentCollection.hpp
 */

#ifndef DRLIB_InstrumentCollection_H
#define DRLIB_InstrumentCollection_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Instrument.hpp"

DRLIB_BEGIN_NAMESPACE

class IModel;
class SensControl;
class Sensitivity;

FORWARD_DECLARE(IInstrumentCollection)

/**
 * A virtual array of instruments which can all be priced together
 *
 * This interface is the heart of QLib's "vectorized pricing" facility, which
 * allows sets of instruments to be priced/Greeked in one pass, exploiting
 * opportunities for optimisation such as
 *
 *    -  reusing MonteCarlo paths to price them all, rather than regenerating
 *       them for each;
 *
 *    -  potentially, reading off whole term structures of Vanilla prices
 *       from a tree (Tree1f) rather than pricing them individually with
 *       separate trees of different depths.
 *
 *
 * <H3>Overview of vectorized pricing</H3>
 *
 * The client entry point is MultiRiskMgrInterface, which accepts an "array"
 * of instruments to price defined via an IInstrumentCollection.
 *
 * The key internal entry point is IInstrumentCollection::Price(), which is
 * similar to IInstrument::Price() except that it accepts an array of Results
 * (one for each contained instrument) rather than just one.
 *
 * The boring implementation ArrayInstrumentCollection is a thin wrapper round
 * a concrete array of instruments; its ArrayInstrumentCollection::Price()
 * method simply calls Model::Price() for each element.  This (obviously)
 * doesn't give any significant performance gain over pricing each instrument
 * in a separate run, but is very general.
 *
 * More interesting implementations, like VanillaGridInstrumentCollection, are
 * not concrete arrays of instruments, but parameterised pseudo-arrays which
 * are never explicitly unpacked.  Their Price() implementations can take
 * advantage of the known structure of the collection to speed things up.
 *
 * Most of the entry points in the QLib framework are now defined in terms
 * of IInstrumentCollection rather than IInstrument: at the top level,
 *
 *    -  RiskMgrInterface is a wrapper round MultiRiskMgrInterface,
 *       calling it with an ArrayInstrumentCollection::singleton()
 * 
 * Further down, similarly:
 *
 *    -  RiskMgr::calculate() calls RiskMgr::calculateMulti()
 *    -  RiskMgr::run() calls RiskMgr::runMulti()
 *    -  Control::calculate() calls Control::calculateMulti()
 *    -  CModel::Run() calls CModel::RunMulti()
 *    -  CModel::getInstrumentAndModelMarket() calls
 *       CModel::getInstrumentsAndModelMarket(), plural
 *
 * At the bottom level, though, we revert by default to one-at-a time pricing,
 * because in the general case that's the only possibility:
 *
 *    -  CModel::PriceMulti() calls CModel::Price()
 * 
 * To implement fast special cases, you can override either
 * IInstrumentCollection::Price() or CModel::PriceMulti().
 *
 *
 * <H3>Using this interface</H3>
 *
 * Apart from Price(), and the obvious operator[]() and size(),
 * IInstrumentCollection defines some housekeeping methods similar to
 * IInstrument's, to support market data fetching and some other more obscure
 * QLib features.  See ArrayInstrumentCollection for obvious implementations
 * in terms of arbitrary instruments.
 *
 * You should inherit from CInstrumentCollection if you can, rather than
 * directly from IInstrumentCollection, because then you get emptyResults() for
 * free.
 */

class RISKMGR_DLL IInstrumentCollection: public virtual IObject {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP TYPE;

protected:

    IInstrumentCollection();

public:

    virtual ~IInstrumentCollection();

    /**
     * "Array" access methods
     *
     * An IInstrumentCollection can behave like an array of instruments,
     * although it doesn't have to be implemented as a concrete array
     */

    //@{

    /**
     * Number of instruments in the collection
     */

    virtual int size() const = 0;

    /**
     * An individual instrument from the collection
     *
     * Note that it doesn't have to exist explicitly (see e.g. implementation
     * of VanillaGridInstrumentCollection::operator[]()), so don't call this
     * unless you really want to.
     */

    virtual InstrumentSP operator [](int i) = 0;

    //@}

    /**
     * Housekeeping
     */

    //@{

    /**
     * Throw an exception if the object isn't in a valid state
     */

    virtual void Validate() = 0;

    /**
     * Equivalent to calling IInstrument::GetMarket() on the underlying
     * instruments
     */

    virtual void GetMarket(const IModel*, const CMarketDataSP) = 0;

    /**
     * Today
     */

    virtual DateTime getValueDate() const = 0;

    /**
     * Curve for discounting payments---must be common to all instruments
     */

    virtual string discountYieldCurveName() const = 0;

    /**
     * Date after which no instruments in the collection are sensitive to
     * market data (with respect to the supplied sensitivity)
     */

    virtual DateTime endDate(const Sensitivity*) const = 0;

    //@}

    /**
     * Pricing
     */

    //@{

    /**
     * Price all the instruments in the collection using the supplied IModel
     * and CControl
     *
     * Vector equivalent of IInstrument::Price: equivalent to pricing each
     * instrument #i separately and putting the results in the i'th entry of @a
     * resultss.  But you can use whatever tricks you can to make it go faster
     * in your particular implementation.
     *
     * @param resultss    On entry, assumed to have length == this->size()
     */

    virtual void Price(IModel *model, CControl *control, CResultsArraySP resultss) = 0;

    /**
     * Scale results of a pricing run
     *
     * Equivalent having each instrument in the collection scale the
     * corresponding entry in @a resultss
     */

    virtual void scaleOutputs(CControlSP control, CResultsArraySP resultss) = 0;

    //@}

    /**
     * Utilities
     */

    //@{

    /**
     * An array of empty Results objects, with length equal to our size()
     *
     * This is implemented as CInstrumentCollection::emptyResults(), so you
     * generally don't have to define it
     */

    virtual CResultsArraySP emptyResults() const = 0;

    /**
     * An IInstrumentCollection containing a single instrument
     */

    static IInstrumentCollectionSP singleton(InstrumentSP instrument);

    /**
     * An IInstrumentCollection containing a single instrument
     */

    static IInstrumentCollectionSP singleton(Instrument* instrument);

    //@}
};
#ifndef QLIB_IINSTRUMENTCOLLECTION_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<IInstrumentCollection>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<IInstrumentCollection>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<IInstrumentCollection>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<IInstrumentCollection>);
#endif

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
