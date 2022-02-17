/**
 * @file Control_forward.hpp
 */

#ifndef DRLIB_Control_forward_H
#define DRLIB_Control_forward_H

#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

class Control;
typedef Control CControl;
typedef smartConstPtr<CControl> CControlConstSP;
typedef smartPtr<CControl> CControlSP;
typedef array<CControlSP, CControl> CControlArray;
typedef smartPtr<CControlArray> CControlArraySP;
typedef smartConstPtr<CControlArray> CControlArrayConstSP;
#ifndef QLIB_CONTROL_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<CControl>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<CControl>);
EXTERN_TEMPLATE(class RISKMGR_DLL array<CControlSP _COMMA_ CControl>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<CControlArray>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<CControlArray>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<CControl>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<CControl>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL array<CControlSP _COMMA_ CControl>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<CControlArray>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<CControlArray>);
#endif

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
