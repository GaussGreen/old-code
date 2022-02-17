#include "edginc/config.hpp"
#include "edginc/ZeroPair.hpp"
#include "edginc/ZeroCurve.hpp"

DRLIB_BEGIN_NAMESPACE

// default constructor
ZeroPair::ZeroPair() {}
    
// useful constructor
ZeroPair::ZeroPair(ZeroCurveSP discZC, ZeroCurveSP growZC) : 
    discZC(discZC), growZC(growZC) {}
    
// is there anybody in there?
bool ZeroPair::empty() const
{
    return (discZC.get() == 0);
}

const ZeroCurve* ZeroPair::get(bool useProjectionCurve) const
{
    return (useProjectionCurve ? growZC.get() : discZC.get());
}

DRLIB_END_NAMESPACE
