/**
 * @file IResultsFunction.hpp
 */

#ifndef DRLIB_IResultsFunction_H
#define DRLIB_IResultsFunction_H

#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

class Results;
FORWARD_DECLARE(OutputRequest)
FORWARD_DECLARE(IResultsFunction)

/**
 * A function which extracts a value from a Results dictionary.
 *
 * This class is used to say what quantity a greek is calculated with respect
 * to---typically price() of course.  See the @a derivand argument to
 * RiskPropertySensitivity::Deriv::Deriv().
 *
 * At a lower level it says what the value of a HypotheticalQuantity is in a
 * given state of the world.  See HypotheticalQuantity::quantity().
 */

class RISKMGR_DLL IResultsFunction: public virtual IObject {
public:

    static CClassConstSP const TYPE;

    /**
     * Given a Results, return a number from it
     */

    virtual double operator()(const Results *) const = 0;

    /**
     * Any special output requests required in the Results
     */

    virtual OutputRequestArrayConstSP outputRequests() const = 0;

    /**
     * A function which calls Results::retrievePrice() on its argument
     */

    static IResultsFunctionConstSP price();

    /**
     * A function returning the value of a certain OutputRequest in its Results
     * argument
     */

    static IResultsFunctionConstSP outputRequest(OutputRequestConstSP request);

    /**
     * A function returning the value of a certain OutputRequest in its Results
     * argument
     */

    static IResultsFunctionConstSP outputRequest(const string& request);

    /**
     * A boring function which just returns zero and ignores its argument
     */

    static IResultsFunctionConstSP zero();
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
