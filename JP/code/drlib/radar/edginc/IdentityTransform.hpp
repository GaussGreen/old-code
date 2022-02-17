#ifndef IDENTITY_TRANSFORM_H
#define IDENTITY_TRANSFORM_H

#include "edginc/IFittingVarTransform.hpp"
#include "edginc/RadarRepUtil.hpp"

DRLIB_BEGIN_NAMESPACE

class  RADAR_DLL IdentityTransform : public virtual IFittingVarTransform {
public:
    IdentityTransform(void) {}
    virtual TransformedArray operator () (const FittingArray& raw) const {
        TransformedArray tarr(raw.begin(), raw.end());
        return tarr;
    }
    virtual ~IdentityTransform() {}
};


DRLIB_END_NAMESPACE
#endif
