/**
 * @file IResultsIdentifier.hpp
 */

#ifndef DRLIB_IResultsIdentifier_H
#define DRLIB_IResultsIdentifier_H

#include "edginc/Atomic.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Results_forward.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(Void)
FORWARD_DECLARE(Untweakable)
FORWARD_DECLARE(OutputName)
FORWARD_DECLARE(Expiry)
FORWARD_DECLARE(ExpiryWindow)
FORWARD_DECLARE(ExpiryPair)
FORWARD_DECLARE(ExpiryAndStrike)
FORWARD_DECLARE(BoxedInt)
FORWARD_DECLARE(IResultsIdentifier)

/**
 * A "generalized name" for an entry in a Results dictionary
 *
 * Storing an entry in a Results dictionary is not always just a matter of
 * knowing its category label (like "DELTA") and the market data name to which
 * it refers ("IBM") and calling a single Results method.  For instance,
 * entries which also need to be qualified by expiry (like "IBM 3M vega") have
 * to find their way into an ExpiryResultArray which is itself stored in the
 * Results under the appropriate packet/name.
 *
 * This abstract class provides an interface for Results storage which extends
 * beyond the scalar case in a uniform way.  You just use the appropriate SP()
 * method to create an IResultsIdentifier which you can then use as a name
 * for your result.
 *
 * IResultsIdentifier's are mostly created in
 * PerNameRiskPropertySensitivity::nameRiskQuantities(), which makes an
 * overloaded call to IResultsIdentifier::SP() to create an identifier for
 * each NamedRiskQuantity it returns.
 */

class RISKMGR_DLL IResultsIdentifier: public virtual IObject {
public:

    static CClassConstSP const TYPE;

    IResultsIdentifier();
    ~IResultsIdentifier();

public:

    /**
     * Whether a number has been stored against this name in @a results.
     */

    virtual bool exists(ResultsConstSP results) const = 0;

    /**
     * Store a number against this name in @a results.
     */

    virtual void store(ResultsSP results, double value) const = 0;

     /**
     * Store an IObjectConst against this name in @a results.
     */
    
    virtual void store(ResultsSP results, IObjectSP value) const = 0;

    /**
     * Store "untweakable" against this name in @a results.
     */

    virtual void storeUntweakable(CResultsSP results,
                                  UntweakableConstSP oops) const = 0;

    /**
     * Store "not applicable" against this name in @a results.
     */

    virtual void storeNotApplicableToName(CResultsSP results) const = 0;

    /**
     * Store "not applicable" against the whole packet in @a results
     */

    virtual void storeNotApplicableToInstrument(CResultsSP results) const = 0;

    /**
     * The packet name in which the result will be stored.
     *
     * It would be nice not to have this; it's used for filtering BasketDelta
     * results in Delta.cpp, please don't use it elsewhere.
     */

    virtual string packet() const = 0;

    /**
     * The name in which the result will be stored in the packet().
     *
     * It would be nice not to have this; it's used for filtering BasketDelta
     * results in Delta.cpp, please don't use it elsewhere.
     */

    virtual OutputNameConstSP entry() const = 0;

    /**
     * Constructors
     */

    //@{

    /**
     * A name for a scalar Results entry
     *
     * This a pretty thin wrapper round Results::storeScalarGreek().
     */

    static IResultsIdentifierSP SP(string packet, const string& entry);

    /**
     * A name for a scalar Results entry
     *
     * This a pretty thin wrapper round Results::storeScalarGreek().
     */

    static IResultsIdentifierSP SP(string packet, OutputNameConstSP entry,
                                   VoidConstSP = VoidConstSP(),
                                   VoidConstSP = VoidConstSP());

    /**
     * A name for an entry in a DoubleArray in a Results dictionary
     */

    static IResultsIdentifierSP SP(string packet, OutputNameConstSP entry,
                                   int index);

    /**
     * A name for an entry in a DoubleArray in a Results dictionary
     */

    static IResultsIdentifierSP SP(string packet, OutputNameConstSP entry,
                                   BoxedIntConstSP index);

    /**
     * A name for an entry in an ExpiryResultArray in a Results dictionary
     */

    static IResultsIdentifierSP SP(string packet, OutputNameConstSP entry,
                                   ExpiryConstSP expiry);

    /**
     * A name for an entry in an ExpiryResultArray in a Results dictionary
     *
     * Keys by ExpiryWindow::expiry and ignores ExpiryWindow::previous,
     * ExpiryWindow::next.
     */

    static IResultsIdentifierSP SP(string packet, OutputNameConstSP entry,
                                   ExpiryWindowConstSP expiry);

    /**
     * A name for an entry in an ExpiryPairResultArray in a Results dictionary
     */

    static IResultsIdentifierSP SP(string packet, OutputNameConstSP entry,
                                   ExpiryPairConstSP expiry);

    /**
     * A name for an entry in a MatrixResult in a Results dictionary
     */

    static IResultsIdentifierSP SP(string packet, OutputNameConstSP entry,
                                   ExpiryAndStrikeConstSP expiryAndStrike);

    /**
     * A name for an entry in an ExpiryPairResultArray in a Results dictionary
     *
     * Keys by ExpiryWindow::expiry and ignores ExpiryWindow::previous,
     * ExpiryWindow::next.
     */

    static IResultsIdentifierSP SP(string packet, OutputNameConstSP entry,
                                   ExpiryWindowConstSP expiry0,
                                   ExpiryWindowConstSP expiry1);
};

DECLARE(IResultsIdentifier)

DRLIB_END_NAMESPACE

#endif
