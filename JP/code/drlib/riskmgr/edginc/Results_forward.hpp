/**
 * @file Results_forward.hpp
 */

#ifndef DRLIB_Results_forward_H
#define DRLIB_Results_forward_H

#include "edginc/Array.hpp"

DRLIB_BEGIN_NAMESPACE

class Results;
typedef Results CResults;
typedef smartPtr<Results> ResultsSP;
typedef smartConstPtr<Results> ResultsConstSP;
#ifndef QLIB_RESULTS_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<Results>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<Results>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<Results>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<Results>);
#endif
typedef ResultsSP CResultsSP;
typedef ResultsConstSP CResultsConstSP;
typedef array<CResultsSP, CResults> CResultsArray;
typedef smartPtr<CResultsArray> CResultsArraySP;
typedef smartConstPtr<CResultsArray> CResultsArrayConstSP;
#ifndef QLIB_RESULTS_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL array<CResultsSP _COMMA_ CResults>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<CResultsArray>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<CResultsArray>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL array<CResultsSP _COMMA_ CResults>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<CResultsArray>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<CResultsArray>);
#endif

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
