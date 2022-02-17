/*
 // Module:      driface
 // File:        cdiopt.h
 // Function:    CDI option pricer
 // Author:      Julia Chislenko July 2001

$Header$

#ifndef DRI_CDI_OPT
#define DRI_CDI_OPT

#include "drlstd.h"

#include "cgeneral.h"
#include "bastypes.h"
#include "macros.h"             // MAX 

//f---------------------------------------------------------------------
 / Pricing routine for CDI options.
 / Returns SUCCESS/FAILURE.
 //

DLL_EXPORT(int)
DriCDIOpt(   
char        instrType,   	// (I) P/R/S/C/F 
double      notional,           // (I) notional 
TDate       startDate,          // (I) accrual start 
TDate       endDate,            // (I) accrual end = expiration for cap/floor 
double      strike,             // (I) strike rate 
TBoolean    cashSettle,         // (I) TRUE=cash settle FALSE=physical for swaptions 
long        dcDenom,            // (I) day count denom for CDI compounding 
				       (act if 360 or 365, bus if 30 or 252) 
long        numCDIFixs,         // (I) num CDI fixings betw start and value date 
double     *cdiFixings,         // (I) CDI fixings betw start and value date 
TCurve     *zc,                 // (I) Brazil zero curve 
long        zeroInterpType,     // (I) zero interp type 
TCurve     *baseVolCrv,         // (I) used for caps, can be NULL otherwise 
TSwaptionMatrix2D *swVolMtx,    // (I) used for swaptions, can be NULL otherwise 
long        numModlPars,        // (I) number of model params 
double     *modlParams,         // (I) MRs, weights, corrs 
char       *holidays,           // (I) for business day compounding 
double     *pv);		// (O) present value 

//f---------------------------------------------------------------------
 // Dr-Wrapper for {\tt DriCDIOpt}.
 The argument {\tt dataFnam} specifies the name of the file
  containing the data (if it is NULL, the default name
  "cdiopt_w.dat" is used.
  Returns SUCCESS/FAILURE.
 

DLL_EXPORT(int)
DriCDIOptW(char *dataFnam);

#endif*/
