#ifndef SVGEN_KCOMPONENT_HPP
#define SVGEN_KCOMPONENT_HPP

/*#include "edginc/Holiday.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CouponSched.hpp"*/

#include "edginc/SVGenIProdCreator.hpp"

DRLIB_BEGIN_NAMESPACE

class  IRPRODUCTS_DLL SVGenKComponent : public SVGenIProdCreator {
   // public : 
};
DECLARE_REF_COUNT(SVGenKComponent);


/*Now the definitions of state variables start*/
class  IRPRODUCTS_DLL SVKComponent : public SVIProdCreator {
};
DECLARE_REF_COUNT(SVKComponent);

DRLIB_END_NAMESPACE
#endif
