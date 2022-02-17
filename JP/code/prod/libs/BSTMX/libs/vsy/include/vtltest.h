/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtltest.h
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_vtltest_H
#define	_vtltest_H

#include "vtlbase.h"
#include "kzcurve.h"


//--------------------------------------------------------------
/**
 * This routine performs a test consisting on comparing 
 * forward and discount factors computed form the zero curves
 * and from the tree, together with a calculation of the implied
 * Black-Scholes volatility smile for caps on the rate.
 * The routine prints a text format output of the test on
 * the user specified stream "os".
 *
 * <br><br><b>Arguments Details</b><br><br><ul>
 * <li><b> vt: </b>
 *	A virtual tree already initialized (but not calibrated).
 * 	The tree should already contain the discount and index curves
 *	with the corresponding names
 *	(these will not be added by the routine).
 * <li><b> resetDates:</b>
 *	A vector of reset date for the forward/options.
 * <li> <b>payDates:</b>
 * 	A vector of payment date for the forward/options.
 *        Must be of the same length as resetDates.
 * <li> <b>strikes:</b>
 * 	A vector of requested strike offset (i.e. expressed
 *	not as a distance from the ATM forward rate that will be internally 
 *	calculated.
 * <li> <b>floatRate:</b>
 * 	The definition for the rate used in the test.
 * <li> <b>curveName:</b>
 * 	The name of the zero curve used to calculate
 *	the index rate.
 * <li> <b>zcCurve:</b>
 *	The index zero curve (used to calculated the 
 *	the deterministic forwards).
 * <li> <b>discCurveName:</b>
 *	The name of the zero curve used to discount.
 * <li> <b>discZcCurve:</b>
 * 	The zero curve used for discounting.
 * <li> <b>os:</b>
 * 	The stream where to write the test results.
 * 
 */

void KVTreeTestForwardRates1(
	KVTree &vt,			// (I) 
	KVector(TDate) &resetDates,	// (I) 
	KVector(TDate) &payDates,	// (I) 
	KVector(double) &strikes,	// (I) 
	KRate &floatRate,		// (I) 
	const String &curveName,	// (I) 
	const KZCurve &zcCurve,		// (I) 
	const String &discCurveName,	// (I) 
	const KZCurve &discZcCurve,	// (I) 
	ostream &os);			// (I) 



//--------------------------------------------------------------
/**
 * This routine performs a test consisting in  calculating
 * the realized correlation of two index rates in
 * the tree, together with a calculation of the implied
 * normal volatility smile for spread options between the rates.
 * The routine prints a text format output of the test on
 * the user specified stream "os".
 *
 * <br><br><b>Arguments Details</b><br><br><ul>
 * <li><b> vt: </b>
 *	A virtual tree already initialized (but not calibrated).
 * 	The tree should already contain the discount and index curves
 *	with the corresponding names
 *	(these will not be added by the routine).
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

void KVTreeTestCorrRates(
	KVTree &vt,			// (I) 
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




