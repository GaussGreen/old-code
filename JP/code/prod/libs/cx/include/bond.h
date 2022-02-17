/*
*************************************************************************
**FILE NAME: bond.c
**
** Bond testing functions
**
**
*************************************************************************
*/

#ifndef CX_BOND_H
#define CX_BOND_H

#include "cx.h"


int CxBondTest(TBond* bond,
	       double * output);

int CxCashFlowListPV(TCashFlowList * cfl,
		     TDate valueDate,
		     TCurve* discCurve,
		     CxTCreditCurve *spreadCurve,
		     double         *pv);


int CxBondFeeLegPV(TBond* bond,
		   TDate valueDate,
		   TDate startDate,
		   TCurve* discCurve,
		   CxTCreditCurve *spreadCurve,
		   CxTRecoveryCurve *recoveryCurve,
		   TBoolean isClean,
		   double         *pv);

CxTContingentLeg* CxBondContingentLegMake
(TDate     startDate,
 TBond*     bond,  
 double    notional, 
 long      delay,    
 TBoolean  protectStart);

int CxBondContingentLegPV
(TDate             today,
 TDate             valueDate,
 TDate startDate,
 TBond* bond,
 double            notional,
 long              delay,
 TCurve           *discCurve,
 CxTCreditCurve   *spreadCurve,
 CxTRecoveryCurve *recoveryCurve,
 TBoolean          protectStart,
 double           *pv);

int CxBondPrice
(TDate             today,
 TDate             valueDate,
 TDate startDate,
 TBond* bond,
 double            notional,
 long              delay,
 TCurve           *discCurve,
 CxTCreditCurve   *spreadCurve,
 CxTRecoveryCurve *recoveryCurve,
 TBoolean          protectStart,
 TBoolean              isClean,
 double           *pv);

int CxBondFeeLegAccrualsOnDefaultPV
(TCashFlowList* cfl,
 TDate valueDate,
 TCurve* discCurve,
 CxTCreditCurve *spreadCurve,
 CxTRecoveryCurve* recoveryCurve,
 double *bondAccrualOnDefault);

#endif
