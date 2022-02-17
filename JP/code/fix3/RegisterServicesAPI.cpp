// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2000 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// $Header$
//

#pragma hdrstop

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include "fix123API.h"

 
namespace IR {

LG_REGISTER_CLASS_TAG(_T_CURVE, _T_CURVEType)
LG_REGISTER_CLASS_TAG(_TREE_DATA, _TREE_DATAType)
LG_REGISTER_CLASS_TAG(_MKTVOL_DATA, _MKTVOL_DATAType)

LG_REGISTER_CLASS(MARKET_DATA)
LG_REGISTER_CLASS(MODEL_DATA)
LG_REGISTER_CLASS(BaseVolData)
LG_REGISTER_CLASS(SwapVolData)
LG_REGISTER_CLASS(ModelParametersData)
LG_REGISTER_CLASS(OPT_OUT_DATA_API)
LG_REGISTER_CLASS(POPULATE_MKTVOL_DATA)

LG_REGISTER_LIBRARY("IRLib")

} // IR

