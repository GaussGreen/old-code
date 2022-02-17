/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtlogtest.h
 * Function:	
 * Author:	David Liu
 *****************************************************************/
#ifndef	_vtlogtest_H
#define	_vtlogtest_H

#include "vtlbase.h"

#include "kmodpar.h"
#include "kvoldat.h"
#include "kzcurve.h"

//--------------------------------------------------------------
/**
 * This routine performs PSA basis test.  It consists in calculating
 * from the tree the following values and comparing them with theoretical
 * results.
 * 1. the IR zeros,
 * 2. the expected forward rates of Libor and PSA basis rate,
 * 3. the lognormal volatilities of Libor, basis spread, as well as PSA basis
 *    rate, 
 * 4. the realized correlations between libor and spread, and between Libor 
 *    and PSA basis rates.

 * The routine prints a text format output of the test on
 * the user specified stream "os".
 *
 * <br><br><b>Arguments Details</b><br><br><ul>
 * <li><b> vt: </b>
 *	A virtual tree already initialized (but not calibrated).
 * 	The tree should already contain the discount and index curves
 *	with the corresponding names
 *	(these will not be added by the routine).
 * <li><b> market:</b>
 * 	market environment, including zero curves and vol curves.
 * <li><b> dealParam:</b>
 * 	model parameters. 
 * <li><b> resetDates:</b>
 *	A vector of reset date for the forward/options.
 * <li> <b>payDates:</b>
 * 	A vector of payment date for the forward/options.
 *        Must be of the same length as resetDates.
 * <li> <b>strikes:</b>
 * 	A vector of requested strike offset (i.e. expressed
 *	not as a distance from the ATM forward spread that will be internally 
 *	calculated).
 * <li> <b>floatRate1:</b>
 * 	The definition for the first rate used in the test.
 * <li> <b>curveName1:</b>
 * 	The name of the zero curve used to calculate
 *	the index first rate.
 * <li> <b>zcCurve1:</b>
 *	The first index zero curve (used to calculated the 
 *	the deterministic forwards).
 * <li> <b>floatRate2:</b>
 * 	The definition for the second rate used in the test.
 * <li> <b>curveName2:</b>
 * 	The name of the zero curve used to calculate
 *	the index seconrd rate.
 * <li> <b>zcCurve1:</b>
 *	The second index zero curve (used to calculated the 
 *	the deterministic forwards).
 * <li> <b>discCurveName:</b>
 *	The name of the zero curve used to discount.
 * <li> <b>discZcCurve:</b>
 * 	The zero curve used for discounting.
 * <li> <b>os:</b>
 * 	The stream where to write the test results.
 * 
 */

void KVTreeTestBasisRatesVolAndCorr(
	KVTree *vt,			// (I) 

	KVolDiag& irVolDiag,		// (I) ir volatility data
	KVolDiag& bsVolDiag,		// (I) basis volatility data 
	KMrParam&  treeModelParam,	// (I) full model parameter

	KVector(TDate) &resetDates,	// (I)
	KVector(TDate) &payDates,	// (I)
	KVector(double) &strikes,	// (I)	offset 

	KRate &floatRate0,		// (I)
	const String &curveName0,	// (I)
	const KZCurve &zcCurve0,	// (I)

	KRate &floatRate1,		// (I)
	const String &curveName1,	// (I)
	const KZCurve &zcCurve1,	// (I)

	const String &discCurveName,	// (I)
	const KZCurve &discZcCurve,	// (I)
	ostream &os);			// (I)



#endif




