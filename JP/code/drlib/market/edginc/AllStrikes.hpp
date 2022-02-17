//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//
//----------------------------------------------------------------------------

#ifndef EDG_ALL_STRIKES_HPP
#define EDG_ALL_STRIKES_HPP

#include "edginc/config.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL IAllStrikes: virtual public IObject {
public:

    /** returns a double list of all strikes on the vol surface */
    virtual DoubleArraySP getAllStrikes() const = 0;

    static CClassConstSP const TYPE; 
};

// typedef for smart pointers to ISensitiveStrikes
typedef smartConstPtr<IAllStrikes> IAllStrikesConstSP;
typedef smartPtr<IAllStrikes> IAllStrikesSP;

DRLIB_END_NAMESPACE

#endif // ALL_STRIKES_HPP
