/**********************************************************
*
*	kconstant.h
*
***********************************************************/

#ifndef _KCONSTANT_H
#define _KCONSTANT_H

#include "kstring.h"
#include "kplatdep.h"
#define	DLL_EXPORT  __declspec(dllexport)

const	KString __ContYieldDiv = "CONTINUOUS_DIV";
const	KString __FixAmountDiv = "FIXED_AMOUNT_DIV";
const	KString __FixedYieldDiv = "FIXED_YIELD_DIV";
const	KString __NonDiv		= "NON_DIVIDEND";

//often used Date types
const KString	__divDate("DIVIDEND_DATE");
const KString	__exerDate("EXERCISE_DATE");
const KString	__couponSetDate("COUPON_SET_DATE");
const KString	__couponPayDate("COUPON_PAY_DATE");
const KString	__principalSetDate("PRINCIPAL_SET_DATE");
const KString	__principalPayDate("PRINCIPAL_PAY_DATE");
const KString	__maturityDate("MATURITY_DATE");
const KString	__lastTreeDate("LAST_TREE_DATE");

//often used critical values
const KString	__coupon("COUPON_AMOUNT");
const KString	__principal("PRINCIPAL_AMOUNT");
//_couponPayDate - __couponSetDate
const KString	__couponSetPayDiff("COUPON_SET_PAY_DIFF");
const KString	__couponFactor("COUPON_FACTOR");
//__principalPayDate - __principalSetDate
const KString	__principalSetPayDiff("PRINCIPAL_SET_PAY_DIFF");
const KString	__timeFromValueDate("TIME_FROM_VALUE_DATE");
const KString	__timeTillNextTp("TIME_TILL_NEXT_TP");
const KString	__zeroPrice("ZERO_PRICE");
const KString	__fwdRateTillNextTp("FWD_RATE_TILL_NEXT_TP");
//critical values for equity tree
const KString	__eqTreeSpotVol("EQ_TREE_SPOTVOL");
const KString	__eqTreeLimit("EQ_TREE_LIMIT");
const KString	__eqTreeMaxLimit("EQ_TREE_MAX_LIMIT");
const KString	__eqFwdPrice("EQ_FWD_PRICE");
const KString	__eqTreeMidNode("EQ_TREE_MID_NODE");
const KString	__eqDividend("EQ_DIVIDEND");
//critical values for interest rate tree
const KString	__irTreeSpotVol1F("IR_TREE_SPOTVOL_1F");
const KString	__irTreeSpotVol2F("IR_TREE_SPOTVOL_2F");
const KString	__irTreeLimit1F("IR_TREE_LIMIT_1F");
const KString	__irTreeLimit2F("IR_TREE_LIMIT_2F");
const KString	__irTreeMaxLimit1F("IR_TREE_MAX_LIMIT_1F");
const KString	__irTreeMaxLimit2F("IR_TREE_MAX_LIMIT_2F");
const KString	__irTreeDrift("IR_TREE_DRIFT");
//critical values for spread tree
const KString	__spTreeSpotVol("SP_TREE_SPOTVOL");
const KString	__spTreeLimit("SP_TREE_LIMIT");
const KString	__spTreeMaxLimit("SP_TREE_MAX_LIMIT");
const KString	__spTreeDrift("SP_TREE_DRIFT");
const KString	__fwdSpreadTillNextTp("FWD_SPREAD_TILL_NEXT_TP");
const KString	__riskyZeroPrice("RISKY_ZERO_PRICE");
const KString	__riskyEqFwdPrice("RISKY_EQ_FWD_PRICE");
const KString	__riskyEqTreeMidNode("RISKY_EQ_TREE_MID_NODE");


//critical values for general tree
const KString	__d1TreeSpotVol("D1_TREE_SPOTVOL");
const KString	__d1TreeSpotVolShift("D1_TREE_SPOTVOLSHIFT");		//HY3.4v
const KString	__d1TreeLimit("D1_TREE_LIMIT");
const KString	__d1TreeMaxLimit("D1_TREE_MAX_LIMIT");
const KString	__d1FwdPrice("D1_FWD_PRICE");
const KString	__d1TreeMidNode("D1_TREE_MID_NODE");
const KString	__d1Dividend("D1_DIVIDEND");



//options
const KString	__exerciseDate("EXERCISE_DATE");
const KString	__optionStrike("OPTION_STRIKE");

#endif

