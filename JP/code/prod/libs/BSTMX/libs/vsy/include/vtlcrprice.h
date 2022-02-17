/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtlcrprice.h
 * Function:	
 * Author:	David Liu
 *****************************************************************/
#ifndef	_vtlcrprice_H
#define	_vtlcrprice_H

#include "vpbase.h"

#include "kstdinc.h"
#include "kutilios.h"
#include "kmktcrv.h"
				// VP Tools

//--------------------------------------------------------------
//

void
KVPRootCreditPrice(
    KMap(String, double) &results,  // (O) Results
    const   SharedPointer<KVPInstr> &vpRoot,// (I) Root product
    KMarketCurves   &mktCurves, // (I) Curves and curve types
    KVolDiag    &irVolDiag, // (I) IR volatility data.
    KMrParam    &irMrParam, // (I) IR mr data.
    KSmileParam &irSmileParam,  // (I) IR skew data.
    KVolDiag    &crVolDiag, // (I) Credit volatility data.
    KMrParam    &crMrParam, // (I) Credit mr data.
    KSmileParam &crSmileParam,  // (I) Credit skew data.
    double      irCrCorr,   // (I) IR Credit correlation.
    KResetBank  &resetBank, // (I) Rate reset bank
  	int		debugLevel);	// (I) Debug level


/*=============================================================================
 * Returns factor volatilities for the given instrument and market
 *===========================================================================*/
void
KVPRootCreditFactorVols(
    KVector(TDate)& factorVolDates,
    KVector(KVector(double))& factorVols,
    const   SharedPointer<KVPInstr> &vpRoot,// (I) Root product
    KMarketCurves   &mktCurves, // (I) Curves and curve types
    KVolDiag    &irVolDiag, // (I) IR volatility data.
    KMrParam    &irMrParam, // (I) IR mr data.
    KSmileParam &irSmileParam,  // (I) IR skew data.
    KVolDiag    &crVolDiag, // (I) Credit volatility data.
    KMrParam    &crMrParam, // (I) Credit mr data.
    KSmileParam &crSmileParam,  // (I) Credit skew data.
    double      irCrCorr,   // (I) IR Credit correlation.
    KResetBank  &resetBank, // (I) Rate reset bank
  	int		debugLevel);	// (I) Debug level

#endif




