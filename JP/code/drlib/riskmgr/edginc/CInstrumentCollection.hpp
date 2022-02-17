/**
 * @file CInstrumentCollection.hpp
 */

#ifndef DRLIB_CInstrumentCollection_H
#define DRLIB_CInstrumentCollection_H

#include "edginc/IInstrumentCollection.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Abstract base class for IInstrumentCollection implementations
 *
 * Use this if you can for your IInstrumentCollection implementations,
 * so that you get the emptyResults() utility method defined for free.
 */

class RISKMGR_DLL CInstrumentCollection: public CObject,
                             public virtual IInstrumentCollection {

public:

    static CClassConstSP TYPE;

private:

    static void load(CClassSP& clazz);

protected:

    CInstrumentCollection(const CClassConstSP& clazz);

public:

    /**
     * An array of empty Results objects, with length equal to our size()
     */

    CResultsArraySP emptyResults() const;

    virtual ~CInstrumentCollection();
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
