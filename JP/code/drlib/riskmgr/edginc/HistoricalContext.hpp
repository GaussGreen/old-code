//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : HistoricalContext.hoo
//
//   Description : HistoricalContext interface
//
//   Date        : July 2005
//
//
//----------------------------------------------------------------------------


#ifndef EDR_HISTORICALCONTEXT_HPP
#define EDR_HISTORICALCONTEXT_HPP

#include "edginc/smartPtr.hpp"

DRLIB_BEGIN_NAMESPACE

/** Products use IHistoricalContext to store any path dependent information. 
    This allows products to work with path and slice based pricers.*/
class MCARLO_DLL IHistoricalContext {
public:

    virtual ~IHistoricalContext() {};

    /** A deep copy is needed so that values referred to via pointers
        are not shared across paths.  A copy of the object is made
        into the destination */
    virtual void deepCopyTo( IHistoricalContext* destination ) const = 0;
};

typedef refCountPtr<IHistoricalContext> IHistoricalContextSP;

class MCARLO_DLL IHistoricalContextGenerator {
public:
    // create a new historical context
    virtual IHistoricalContextSP createHistoricalContext() = 0;
    
    // get initial context
    virtual IHistoricalContextSP getInitialHistoricalContext() = 0;
    virtual ~IHistoricalContextGenerator() {}
};


DRLIB_END_NAMESPACE


#endif

