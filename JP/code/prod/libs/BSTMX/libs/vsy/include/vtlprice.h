/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtlprice.h
 * Function:	
 * Author:	David Liu
 *****************************************************************/
#ifndef	_vtlprice_H
#define	_vtlprice_H

#include "vpbase.h"

#include "kstdinc.h"
#include "kutilios.h"
#include "kmktcrv.h"
				// VP Tools

//--------------------------------------------------------------
//

void
KVPRootPrice(
	KMap(String, double) &results,	// (O) Results
    const	SharedPointer<KVPInstr>	&vpRoot,// (I) Root product
  	KMarketCurves	&mktCurves,	// (I) Curves and curve types
  	KVolDiag	&irVolDiag,	// (I) IR volatility data.
  	KMrParam	&irMrParam,	// (I) IR mr data.
  	KSmileParam	&irSmileParam,	// (I) IR skew data.
  	KVolDiag	&bsVolDiag,	// (I) Basis volatility data.
  	KMrParam	&bsMrParam,	// (I) Basis mr data.
  	KSmileParam	&bsSmileParam,	// (I) Basis skew data.
  	double		irBsCorr,	// (I) IR basis correlation.
  	KResetBank	&resetBank,	// (I) Rate reset bank
  	int		debugLevel);	// (I) Debug level

void KVPRootPrice_TMX(     
    KMap(String, double) &results,
    const   SharedPointer<KVPInstr> &vpRoot,// (I) Root product

    KMarketCurves   &mktCurves, // (I) Curves and curve types
 
    KVolDiag    &irVolDiag,     // (I) IR volatility data.
    KMrParam    &irMrParam,     // (I) IR mr data.
    KSmileParam &irSmileParam,  // (I) IR skew data.
    KVolDiag    &bsVolDiag,     // (I) Basis volatility data.
    KMrParam    &bsMrParam,     // (I) Basis mr data.
    KSmileParam &bsSmileParam,  // (I) Basis skew data.
    double      irBsCorr,       // (I) IR basis correlation.
    KResetBank  &resetBank,     // (I) Rate reset bank
    int         debugLevel);    // (I) Debug level
                   
/*=============================================================================
 * Returns factor volatilities for the given instrument and market
 *===========================================================================*/
void
KVPRootFactorVols(
    KVector(TDate)& volDates,
    KVector(KVector(double))& volRates,
    const	SharedPointer<KVPInstr>	&vpRoot,// (I) Root product
  	KMarketCurves	&mktCurves,	// (I) Curves and curve types
  	KVolDiag	&irVolDiag,	// (I) IR volatility data.
  	KMrParam	&irMrParam,	// (I) IR mr data.
  	KSmileParam	&irSmileParam,	// (I) IR skew data.
  	KVolDiag	&bsVolDiag,	// (I) Basis volatility data.
  	KMrParam	&bsMrParam,	// (I) Basis mr data.
  	KSmileParam	&bsSmileParam,	// (I) Basis skew data.
  	double		irBsCorr,	// (I) IR basis correlation.
  	KResetBank	&resetBank,	// (I) Rate reset bank
  	int		debugLevel);	// (I) Debug level


#endif




