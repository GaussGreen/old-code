//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//
//----------------------------------------------------------------------------

#ifndef EDG_NEXT_STRIKE_HPP
#define EDG_NEXT_STRIKE_HPP

#include "edginc/config.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL INextStrike: virtual public IObject{
public:
    // DON'T USE THIS CLASS - IT SHOULD HAVE DIED YEARS AGO
    // IF YOU REALLY WANT TO SEE STRIKES ON A SURFACE USE IAllStrikes
    // INSTEAD

    /** given a current spot level, get the next strike on the vol surface where 
        the slope is non-differentiable */
    virtual double getNextStrike(const double& strike,
                                 bool          isUp,
                                 bool&         offSurface) const = 0;

    static CClassConstSP const TYPE; // in Object.cpp

};

// typedef for smart pointers to ISensitiveStrikes
typedef smartConstPtr<INextStrike> INextStrikeConstSP;
typedef smartPtr<INextStrike> INextStrikeSP;

DRLIB_END_NAMESPACE

#endif // NEXT_STRIKE_HPP
