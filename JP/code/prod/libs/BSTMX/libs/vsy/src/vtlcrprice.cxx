//**************************************************************
// 
//  Main Credit Tree Pricing Routine
//
//**************************************************************
				// DR C++ Library
#include "kstdinc.h"
#include "kutilios.h"
#include "kmktcrv.h"

				// VTree Tools
#include "vtlbase.h"
				// Credit Tree
#include "kcrxtree.h"

extern	"C" {
#include "date_sup.h"

#include "drlstr.h"		// Strtok
#include "drltime.h"
#include "drlio.h"
#include "drlmem.h"
#include "drlproc.h"
#include "drlts.h"
#include "drlstr.h"
#include "drlgetop.h"		// getopt()

};




//--------------------------------------------------------------
//

void
KVPRootCreditPrice(
	KMap(String, double) &results,	// (O) Results
    const	SharedPointer<KVPInstr>	&vpRoot,// (I) Root product
  	KMarketCurves	&mktCurves,	// (I) Curves and curve types
  	KVolDiag	&irVolDiag,	// (I) IR volatility data.
  	KMrParam	&irMrParam,	// (I) IR mr data.
  	KSmileParam	&irSmileParam,	// (I) IR skew data.
  	KVolDiag	&crVolDiag,	// (I) Credit volatility data.
  	KMrParam	&crMrParam,	// (I) Credit mr data.
  	KSmileParam	&crSmileParam,	// (I) Credit skew data.
  	double		irCrCorr,	// (I) IR Credit correlation.
  	KResetBank	&resetBank,	// (I) Rate reset bank
  	int		debugLevel)	// (I) Debug level
{
static	char			routine[] = "KVPRootCreditPrice";

	KCrxTree		vt;

	KVolDiag		copyIRVolDiag, copyCRVolDiag;	// A copy of IR,CR vol

	SharedPointer<KVPAtom>	vpAtomRoot;
	SharedPointer<KVPToolAtom> vpToolRoot;	// Root tool

    try {


 
	//----------------------------------------------
        // Root instrument
	//----------------------------------------------
	ASSERT_OR_THROW(vpRoot != NULL);


	//----------------------------------------------
        // Perform Pricing
	//----------------------------------------------

	//
	// Initialize tree
	//
	vt.Initialize(
		mktCurves,
		irVolDiag,
		irMrParam,
		irSmileParam,
		crVolDiag,
		crMrParam,
		crSmileParam,
		irCrCorr,
		resetBank);
 

	//----------------------------------------------
	// Build instrument tree pricing tools
	//----------------------------------------------
	SharedPointerConvertTo(vpRoot, vpAtomRoot);	
	vpToolRoot = NewToolRecursive(vpAtomRoot, vt);
	ASSERT_OR_THROW(vpToolRoot != NULL);


	//
	// Set tree time line
	//
	vt.SetUpTimeline();


	//----------------------------------------------
	// CET
	//----------------------------------------------

	//
	// Make a copy of mIRVolDiag, which is subject to modification
	// by CET iteration
	//
	copyIRVolDiag = irVolDiag;


	KPirTreeCet(
		vt,
		mktCurves,
		copyIRVolDiag,
		irMrParam,
		irSmileParam,
		resetBank);

	if (debugLevel > 0) {
		dppLog << "irVolDiag (AFTER CET):\n" << irVolDiag << endl;
	}

    /*=========================================================================
     * CREDIT CET NEXT (Added by Charles Morcom, October 2005)
     *=======================================================================*/
    copyCRVolDiag = crVolDiag;
    KCrxTreeCet(
        vt,
        mktCurves,
        copyIRVolDiag,
        irMrParam,
        irSmileParam,
        copyCRVolDiag,
        crMrParam,
        crSmileParam,
        irCrCorr,
        resetBank
        );
    if (debugLevel > 0) {
        dppLog << "crVolDiag (AFTER CREDIT CET):\n" << crVolDiag << endl;
    }

	//
	// Calibrate tree
	//
	vt.Calibrate();

	//
	// Initialize tool
	//
	vpToolRoot->Initialize();


	//
	// Rollback
	//
	for (int tpIdx=vt.TPNum(); tpIdx >= 0; tpIdx--) {
		vt.Update(tpIdx);
		vpToolRoot->Update();
	}


	//----------------------------------------------
        // Output results
	//----------------------------------------------
	results = vpToolRoot->GetResults();


    }
    catch (KFailure) {
	throw KFailure ("%s: failed.\n", routine);
    }
}

void
KVPRootCreditFactorVols(
    KVector(TDate)& volDates, 
    KVector(KVector(double))& volRates,
    const	SharedPointer<KVPInstr>	&vpRoot,// (I) Root product
  	KMarketCurves	&mktCurves,	// (I) Curves and curve types
  	KVolDiag	&irVolDiag,	// (I) IR volatility data.
  	KMrParam	&irMrParam,	// (I) IR mr data.
  	KSmileParam	&irSmileParam,	// (I) IR skew data.
  	KVolDiag	&crVolDiag,	// (I) Credit volatility data.
  	KMrParam	&crMrParam,	// (I) Credit mr data.
  	KSmileParam	&crSmileParam,	// (I) Credit skew data.
  	double		irCrCorr,	// (I) IR Credit correlation.
  	KResetBank	&resetBank,	// (I) Rate reset bank
  	int		debugLevel)	// (I) Debug level
{
static	char			routine[] = "KVPRootCreditPrice";

	KCrxTree		vt;

	KVolDiag		copyIRVolDiag, copyCRVolDiag;	// A copy of IR,CR vol

	SharedPointer<KVPAtom>	vpAtomRoot;
	SharedPointer<KVPToolAtom> vpToolRoot;	// Root tool

    try {


 
	//----------------------------------------------
        // Root instrument
	//----------------------------------------------
	ASSERT_OR_THROW(vpRoot != NULL);


	//----------------------------------------------
        // Perform Pricing
	//----------------------------------------------

	//
	// Initialize tree
	//
	vt.Initialize(
		mktCurves,
		irVolDiag,
		irMrParam,
		irSmileParam,
		crVolDiag,
		crMrParam,
		crSmileParam,
		irCrCorr,
		resetBank);
 

	//----------------------------------------------
	// Build instrument tree pricing tools
	//----------------------------------------------
	SharedPointerConvertTo(vpRoot, vpAtomRoot);	
	vpToolRoot = NewToolRecursive(vpAtomRoot, vt);
	ASSERT_OR_THROW(vpToolRoot != NULL);


	//
	// Set tree time line
	//
	vt.SetUpTimeline();


	//----------------------------------------------
	// CET
	//----------------------------------------------

	//
	// Make a copy of mIRVolDiag, which is subject to modification
	// by CET iteration
	//
	copyIRVolDiag = irVolDiag;


	KPirTreeCet(
		vt,
		mktCurves,
		copyIRVolDiag,
		irMrParam,
		irSmileParam,
		resetBank);

	if (debugLevel > 0) {
		dppLog << "irVolDiag (AFTER CET):\n" << irVolDiag << endl;
	}

    /*=========================================================================
     * CREDIT CET NEXT (Added by Charles Morcom, October 2005)
     *=======================================================================*/
    copyCRVolDiag = crVolDiag;
    KCrxTreeCet(
        vt,
        mktCurves,
        copyIRVolDiag,
        irMrParam,
        irSmileParam,
        copyCRVolDiag,
        crMrParam,
        crSmileParam,
        irCrCorr,
        resetBank
        );
    if (debugLevel > 0) {
        dppLog << "crVolDiag (AFTER CREDIT CET):\n" << crVolDiag << endl;
    }

	//
	// Calibrate tree
	//
	vt.Calibrate();

	//
	// Initialize tool
	//
	//vpToolRoot->Initialize();


	//
	// Rollback
	//
	//for (int tpIdx=vt.TPNum(); tpIdx >= 0; tpIdx--) {
	//	vt.Update(tpIdx);
	//	vpToolRoot->Update();
	//}


	//----------------------------------------------
        // Output results
	//----------------------------------------------
	//results = vpToolRoot->GetResults();
    volDates = vt.SpotVolDates();
    volRates = vt.SpotVolatilities();


    }
    catch (KFailure) {
	throw KFailure ("%s: failed.\n", routine);
    }
}

