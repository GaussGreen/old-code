#ifndef IFITTING_VAR_TRANSFORM_H
#define IFITTING_VAR_TRANSFORM_H 

#include "edginc/DECLARE.hpp"
//#include "edginc/Format.hpp"
#include "edginc/RadarRepUtil.hpp"

// #include <vector>
// #include <map>
// #include <numeric>

DRLIB_BEGIN_NAMESPACE

/// FIXME: do we start with "market variables" and transform them into the fitting ones or start directly with the "fitting variables"?
class  RADAR_DLL IFittingVarTransform {
public:
	// FittingVarTransform(vector<FittingVariable*> transformations); 
	virtual TransformedArray operator () (const FittingArray& raw) const = 0;
	virtual ~IFittingVarTransform() {}
};

DECLARE_REF_COUNT(IFittingVarTransform);

DRLIB_END_NAMESPACE
#endif
