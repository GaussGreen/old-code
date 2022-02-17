//**************************************************************
// 
//  Main Tree Pricing Routine
//
//**************************************************************
				// DR C++ Library
#include "kstdinc.h"
#include "kutilios.h"
#include "kmktcrv.h"

				// VTree Tools
#include "vtlbase.h"
				// Basis Tree
#include "kbirtree.h"

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

  	int		debugLevel)	// (I) Debug level
{
static	char			routine[] = "KVPRootPrice";

	KBirTree		vt;

	KVolDiag		copyIRVolDiag;	// A copy of IR vol

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
		bsVolDiag,
		bsMrParam,
		bsSmileParam,
		irBsCorr,
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
KVPRootPrice_TMX(       KMap(String, double) &results,
                        const   SharedPointer<KVPInstr> &vpRoot,
                                                    // (I) Root product
                        KMarketCurves   &mktCurves, // (I) Curves and curve types
                        KVolDiag    &irVolDiag,     // (I) IR volatility data.
                        KMrParam    &irMrParam,     // (I) IR mr data.
                        KSmileParam &irSmileParam,  // (I) IR skew data.
                        KVolDiag    &bsVolDiag,     // (I) Basis volatility data.
                        KMrParam    &bsMrParam,     // (I) Basis mr data.
                        KSmileParam &bsSmileParam,  // (I) Basis skew data.
                        double      irBsCorr,       // (I) IR basis correlation.
                        KResetBank  &resetBank,     // (I) Rate reset bank
                        int         debugLevel)     // (I) Debug levelconst KMrParam& MrParam)
{
    static char routine[] = "TmxPrice";

	KBirTree		vt;

	KVolDiag		copyIRVolDiag;	// A copy of IR vol

	SharedPointer<KVPAtom>	vpAtomRoot;
	SharedPointer<KVPToolAtom> vpToolRoot;	// Root tool
    

    try {
    

    vt.SetTmxFlag(TRUE);
 
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
		bsVolDiag,
		bsMrParam,
		bsSmileParam,
		irBsCorr,
		resetBank);
 
    //----------------------------------------------
    // Initialize TMX engine
	//----------------------------------------------
    vt.WrapperEnvPack( mktCurves,
                       irVolDiag,
                       irMrParam,
                       irSmileParam,
                       vpRoot->GetDiscName());
    

	//----------------------------------------------
	// Build instrument tree pricing tools
	//----------------------------------------------
	SharedPointerConvertTo(vpRoot, vpAtomRoot);	
	vpToolRoot = NewToolRecursive(vpAtomRoot, vt);
	ASSERT_OR_THROW(vpToolRoot != NULL);
 

    // 
    // check IR factor
    //
    ASSERT_OR_THROW(vt.irDim() == 1);
  
    //
	// Set tree time line
	//
    vt.WrapTreeTimeLine(irVolDiag, irMrParam);
    vt.PrintTreeTimeLine();

    // 
    // CET
    //
    vt.TMXirTreeCet(irMrParam);

    // reset mVolDates and mFactoVol according to NmrDate
    vt.ResetBsSpotVol(mktCurves,
                    irVolDiag, 
                    irMrParam,
                    irSmileParam,
                    bsVolDiag,
                    bsMrParam,
                    bsSmileParam,
                    irBsCorr,
                    resetBank);

    // unwrap tmx tree
    vt.Calibrate();

    vpToolRoot->Initialize();

    //
    // Roll back
    //
    for (int tpIdx = vt.TPNum(); tpIdx >= 0; tpIdx --)
    {
        //cout << tpIdx << " CurrDate[" << vt.TPIdxCurrent() << "]= " 
        //     << DateAlibToWrapper(vt.TPDateCurrent()) << endl;
        vt.SetNmrToCcy(tpIdx);
        vt.Update(tpIdx);
        vpToolRoot->Update();
        vt.SetCcyToNmr();
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



//--------------------------------------------------------------
//

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

  	int		debugLevel)	// (I) Debug level
{
static	char			routine[] = "KVPRootPrice";

	KBirTree		vt;

	KVolDiag		copyIRVolDiag;	// A copy of IR vol

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
		bsVolDiag,
		bsMrParam,
		bsSmileParam,
		irBsCorr,
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


	//
	// Calibrate tree
	//
	vt.Calibrate();

	//
	// Initialize tool
	//
//	vpToolRoot->Initialize();


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

