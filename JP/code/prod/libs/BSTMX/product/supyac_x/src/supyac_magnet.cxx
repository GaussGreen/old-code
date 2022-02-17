/***************************************************************
 * Module:		BasisTree
 * Submodule:	
 * File:		supyac_magnet.cxx
 * Function:		Magnet wrapper
 * Author:		David Liu,		Feb, 2002
 ***************************************************************/


#define FIDR_LITERAL(x) FIDR_LITERAL_(x)
#define FIDR_LITERAL_(x) #x
const char* _VERSION_ = FIDR_LITERAL(U_VERSION);

// MAGNET
#undef  _MAGNET_WITH_GENERIC		// need for using Array<TZeroCurve*>
#define _MAGNET_WITHOUT_GENERIC

#include "Magnet/Magnet.h"
#include "General/General.h"


MAGNET_PREFIX("SY_")
MAGNET_CATEGORY("SUPER YACCTION")
//using CM::MagnetVersion;
//using CM::XLOper;
//using CM::RuntimeError;
using namespace CM;


// Market environment
#include "kstlutil.h"
#include "kmktcrv.h"
#include "kvoldat.h"
#include "kmodpar.h"
#include "krstbank.h"

// Products
#include "vpbase.h"	
#include "vpbundle.h"
#include "vpcashfl.h"
#include "vpoption.h"
#include "vpfleg.h"
#include "vpfleg2idx.h"
#include "vpfleg3idx.h"
#include "vpkio.h"
#include "vpkio2idx.h"
#include "vpprotleg.h"
#include "vpdefprotect.h"

#include "vtlprice.h"
#include "vtlcrprice.h"

#include "vtlprice.h"

extern "C" {
#include "zcurveo.h"
#include "crxerror.h"
};

#include <cstdio>


//
// Global variables defined for debugging and market
// environments
//
extern	int			debugLevel;	// debug flag

static	Array<String>nilArray(1,"nil");


//
// Market environment objects
//
class KMarketEnv : public Object {

public:
    /** Default constructor */
    KMarketEnv() {mIRBSCorr = 0e0; mIsCredit = false;}

    /** Pure IR Constructor */
    KMarketEnv(
        Date    today,
        const Array<const TZeroCurve*> &zcCurves, 
        const Array<int>    &zcTypes,
        const Array<String> &zcNames,   
        const Array<Date>   &volExpDates, 
        const Array<Date>   &volMatDates, 
        const Array<int>    &volFreqs,   
        const Array<double> &volRates,  
        String              &volType,  
        const Array<double> &modParams,	  
        String              &tmxFlag,
        const Array<double> &smileParams, 
        const Array<Date>   &irBmkLiqDates,
        const Array<int>    &treeParams, 
        const Array<SharedPointer<KRate> > &rsRates,
        const Array<Date>   &rsDates, 
        const Array<double> &rsValues
    );

    /** Basis Constructor */
    KMarketEnv(
        Date    today,
        const Array<const TZeroCurve*> &zcCurves, 
        const Array<int>    &zcTypes,
        const Array<String> &zcNames,   
            // IR diffuse curve
        const Array<Date>   &irVolExpDates, 
        const Array<Date>   &irVolMatDates, 
        const Array<int>    &irVolFreqs,   
        const Array<double> &irVolRates,  
        String              &irVolType,  
        const Array<double> &irModParams,	   
        String              &tmxFlag,
        const Array<double> &irSmileParams, 
        const Array<Date>   &irBmkLiqDates,

                // Basis curve
        const Array<String> &bsInfo,
                // Basis spread 
        const Array<Date>   &bsVolExpDates, 
        const Array<double> &bsVolRates,  
        const Array<int>    &bsVolFreqs,   
        String              &bsVolType,  
        const Array<double> &bsModParams,	   
        const Array<double> &bsSmileParams, 

        double              corrIRBS, 
                // Tree parameters
        const Array<int>    &treeParams, 
                // Reset Bank
        const Array<SharedPointer<KRate> > &rsRates,
        const Array<Date>	&rsDates, 
        const Array<double>	&rsValues
    );

    /** Credit Constructor */
    KMarketEnv(
        Date    today,
        const Array<const TZeroCurve*>     &zcCurves,
        const Array<int>       &zcTypes,
        const Array<String>    &zcNames,
        const Array<int>       &zcInterps,
                // IR diffuse curve
        const Array<Date>      &irVolExpDates,
        const Array<Date>      &irVolMatDates,
        const Array<int>       &irVolFreqs,
        const Array<double>    &irVolRates,
        String                 &irVolType,
        const Array<double>    &irModParams,
        String                 &tmxFlag,
        const Array<double>    &irSmileParams,
        const Array<Date>      &irBmkLiqDates,

                // Credit curve
        String                 &crIRDiscName,
        double                 recovery,
                // Credit vol
        const Array<Date>      &crVolExpDates,
        const Array<Date>      &crVolMatDates,
        const Array<int>       &crVolFreqs,
        const Array<double>    &crVolRates,
        String                 &crVolType,
        const Array<double>    &crModParams,
        const Array<double>    &crSmileParams,

        double                 corrIRCR,

                // Tree parameters
        const Array<int>      &treeParams,

                // Reset Bank
        const Array<SharedPointer<KRate> > &rsRates,
        const Array<Date>     &rsDates,
        const Array<double>   &rsValues
    );

    /** Is credit or basis tree?
     *  Will call different tree pricer.
     */
virtual bool IsCredit() const {return mIsCredit;}
    /** Is TMX or FIX tree?
      * Will call different tree pricer
      */
virtual bool IsTmx() const {return mIsTmx;}

    /** Destructor */
    ~KMarketEnv() {}

    /** Write to a stram in Yacction format */
virtual ostream& YacctionWrite(ostream& os, int indent = FALSE);

    /** Write to a stream. */
friend  ostream& operator<<(ostream& os, const KMarketEnv& mkt);


    /** Members of the class */
    KMarketCurves   mMktCurves;	// curves and curve types
    KVolDiag        mIRVolDiag;	// IR volatility data.
    KMrParam        mIRMrParam;	// IR mr data.
    KSmileParam     mIRSmileParam;	// IR skew data.
                        // Basis
    KVolDiag        mBSVolDiag;	// Basis volatility data.
    KMrParam        mBSMrParam;	// Basis mr data.
    KSmileParam     mBSSmileParam;	// Basis skew data.
    double          mIRBSCorr;	// IR basis correlation.
                        // Credit
    KVolDiag        mCRVolDiag;       // Credit volatility data.
    KMrParam        mCRMrParam;       // Credit mr data.
    KSmileParam     mCRSmileParam;    // Credit skew data.
    double          mIRCRCorr;        // IR-CR correlation.
                        // Reset bank
    KResetBank      mVPResetBank;	// Reset bank

                        // Is credit or basis/IR tree?
    bool            mIsCredit;
                        // Is TMX or FIX tree
    bool            mIsTmx;

private:
    /** IR Constructor, used by basis as well */
void    SetIRMktEnv(
        Date    today,
        const Array<const TZeroCurve*> 	&zcCurves, 
        const Array<int>    &zcTypes,
        const Array<String> &zcNames,   
        const Array<int>        &zcInterps,
        const Array<Date>   &volExpDates, 
        const Array<Date>   &volMatDates, 
        const Array<int>    &volFreqs,   
        const Array<double> &volRates,  
        String              &volType,  
        const Array<double> &modParams,	   
        String              &tmxFlag,
        const Array<double> &smileParams, 
        const Array<Date>   &bmkLiqDates,
        const Array<int>    &treeParams, 
        const Array<SharedPointer<KRate> > &rsRates,
        const Array<Date>   &rsDates, 
        const Array<double> &rsValues
    );
};


//
// Pure IR Constructor
//
void
KMarketEnv::SetIRMktEnv(
                    // KMarketCurves
    Date                today,	      // (I) 1. today's date	
    const Array<const TZeroCurve*>  &zcCurves,  // (I) 2. ALIB zero curve objs
    const Array<int>    &zcTypes,     // (I) 3. zc types (KV_DIFF, ..)
    const Array<String> &zcNames,     // (I) 4. array of zc names
    const Array<int>    &zcInterps,   // (I) 5. array of zc interps
                    // KVolDiag
    const Array<Date>   &volExpDates, // (I) 6. option expiration dates
    const Array<Date>   &volMatDates, // (I) 7. underlying mat dates
    const Array<int>    &volFreqs,    // (I) 8. underlying rate freqs
    const Array<double> &volRates,    // (I) 9. volatilities
    String              &volType,     // (I) 10.percentage vol or bp vol
                    // KMrParam
    const Array<double> &modParams,   // (I) 11. n,beta,alpha,rho,bbone
                    // KSmileParam
    String              &tmxFlag,     // (I) 12. 'Y'/'N'
    const Array<double> &smileParams, // (I) 13. q1, q2, qfsh, iteration
    const Array<Date>   &bmkLiqDates, // (I) 14  benchmark liquid dates (TMX)
                    // Tree Parameter
    const Array<int>    &treeParams,  // (I) 15. ppy, smooth, numStd
                    // KResetBank
    const Array<SharedPointer<KRate> > &rsRates,// (I) 16. array of float rates 
    const Array<Date>   &rsDates,     // (I) 17. array of reset dates
    const Array<double> &rsValues)    // (I) 18. array of reset values
{
static	char routine[] = "KMarketEnv::SetIRMktEnv";

	int	idx;

	KRate	rsRateNS;	// Strip off the spread from floating reset rate 

 
	//--------------------------------------------------------------
	//	Market Curves
	//--------------------------------------------------------------

	mMktCurves.ReadMagnet(today,
			      zcCurves,
			      zcTypes,
			      zcNames,
			      zcInterps);

	//--------------------------------------------------------------
	//	Vol Curves
	//--------------------------------------------------------------
	
	// Check size consistency
	ASSERT_OR_THROW(volExpDates.size() == volMatDates.size());
	ASSERT_OR_THROW(volExpDates.size() == volFreqs.size());
	ASSERT_OR_THROW(volExpDates.size() == volRates.size());

	mIRVolDiag = KVolDiag(DppDateArrayToTDateVector(volExpDates),
			      DppDateArrayToTDateVector(volMatDates),
			      DppVectorFromCMLIBArray(volFreqs),
			      DppVectorFromCMLIBArray(volRates));

	if (toupper(volType[0]) == 'N')
		mIRVolDiag.mVolType = NORMVOL;
	else if (toupper(volType[0]) == 'L')
		mIRVolDiag.mVolType = LOGVOL;
	else
		throw KFailure("%s: invalid vol type (%s).\n", 
				routine, volType.c_str());


	//--------------------------------------------------------------
	//	Model Parameters
	//--------------------------------------------------------------
	mIRMrParam.ReadMagnet(DppVectorFromCMLIBArray(modParams),
			      DppVectorFromCMLIBArray(treeParams));

	//--------------------------------------------------------------
	//	TMX or FIX tree
	//--------------------------------------------------------------
	if (toupper(tmxFlag[0]) == 'Y')
		mIsTmx = true;
	else
        mIsTmx = false;

    //--------------------------------------------------------------
	//	Smile Parameters
	//--------------------------------------------------------------
    if (IsTmx())
    {
	    mIRSmileParam.ReadMagnet_tmx(DppVectorFromCMLIBArray(smileParams),
                                    DppDateArrayToTDateVector(bmkLiqDates),
                                    mIRVolDiag.mVolDates);
        mIRSmileParam.CheckLiqBmkDates(mIRVolDiag);
    }
    else
    {
	    mIRSmileParam.ReadMagnet(DppVectorFromCMLIBArray(smileParams));
    }

	
	// -------------------------------------------------------------
	//   Basis
	// -------------------------------------------------------------
	mBSMrParam = 0;
	mBSMrParam.mNumFact = 0;
	mBSSmileParam = 0;
	mBSVolDiag = mIRVolDiag;
	mBSVolDiag = 0e0;

	mIRBSCorr = 0e0;


	// -------------------------------------------------------------
	//   Reset Bank
	// -------------------------------------------------------------

	// Check size consistency
	ASSERT_OR_THROW(rsRates.size() == rsDates.size());
	ASSERT_OR_THROW(rsRates.size() == rsValues.size());

	// Add all reset rates
	for (idx = 0; idx < rsRates.size(); idx++)
	{
		// Strip off spread from reset rate
		rsRateNS = (*rsRates[idx]);
		rsRateNS.SetSpread(0e0);

		mVPResetBank.Insert(rsRateNS, 
				    rsDates[idx], 
				    rsValues[idx]);
	}


	 
}



//
// Pure IR Constructor
//
KMarketEnv::KMarketEnv(
                    // KMarketCurves
    Date                today,        // (I) 1. today's date	
    const Array<const TZeroCurve*>  &zcCurves,   // (I) 2. ALIB zero curve objs
    const Array<int>    &zcTypes,     // (I) 3. zc types (KV_DIFF, ..)
    const Array<String> &zcNames,     // (I) 4. array of zc names
                    // KVolDiag
    const Array<Date>   &volExpDates, // (I) 5. option expiration dates
    const Array<Date>   &volMatDates, // (I) 6. underlying mat dates
    const Array<int>    &volFreqs,    // (I) 7. underlying rate freqs
    const Array<double> &volRates,    // (I) 8. volatilities
    String              &volType,     // (I) 9. percentage vol or bp vol
                    // KMrParam
    const Array<double> &modParams,   // (I) 10. n,beta,alpha,rho,bbone
                    // KSmileParam
    String              &tmxFlag,     // (I) 11. 'Y'/'N'
    const Array<double> &smileParams, // (I) 12. q1, q2, qfsh, iteration
                                      //         or MultiQ smile
    const Array<Date>   &bmkLiqDates, // (I) 13. Liquid benchmark dates
                    // Tree Parameter
    const Array<int>    &treeParams,  // (I) 14. ppy, smooth, numStd
                    // KResetBank
    const Array<SharedPointer<KRate> > &rsRates,// (I) 15. array of float rates 
    const Array<Date>   &rsDates,     // (I) 16. array of reset dates
    const Array<double> &rsValues)    // (I) 17. array of reset values
{
static	char routine[] = "KMarketEnv::KMarketEnv";
	
        int             i;
        Array<int>      zcInterps;
	
        // Default IR curve interp type
        for (i=0; i<zcTypes.size(); i++)
            zcInterps.push_back(GTO_LINEAR_INTERP);

	//
	// Construct the IR environment
	//
    SetIRMktEnv(
        today,
        zcCurves,
        zcTypes,
        zcNames,
        zcInterps,
        volExpDates,
        volMatDates,
        volFreqs,
        volRates,
        volType,
        modParams,
        tmxFlag,
        smileParams,
        bmkLiqDates,
        treeParams,
        rsRates,
        rsDates,
        rsValues);


    mIsCredit = false;

}




//
// Basis Constructor
//
KMarketEnv::KMarketEnv(
                    // KMarketCurves
    Date                today,	      // (I) 1. today's date	
    const Array<const TZeroCurve*> &zcCurves,   // (I) 2. ALIB zero curve objs
    const Array<int>    &zcTypes,     // (I) 3. zc types (KV_DIFF, ..)
    const Array<String> &zcNames,     // (I) 4. array of zc names
                    // KVolDiag
    const Array<Date>   &irVolExpDates, // (I) 5. option exp dates
    const Array<Date>   &irVolMatDates, // (I) 6. underlying mat dates
    const Array<int>    &irVolFreqs,    // (I) 7. underlying rate freqs
    const Array<double> &irVolRates,    // (I) 8. volatilities
    String              &irVolType,     // (I) 9. percentage or bp vols
                    // KMrParam
    const Array<double> &irModParams,   // (I) 10. n,beta,alpha,rho,bbone
                    // KSmileParam
    String              &tmxFlag,       // (I) 11. 'Y'/'N'
    const Array<double> &irSmileParams, // (I) 12. q1, q2, qfsh, iteration
                                        //         or MultiQ smile
    const Array<Date>   &irbmkLiqDates, // (I) 13. Liquid benchmark dates

                    // Basis curve
    const Array<String> &bsInfo,        /* (I) 14. Extra param for basis
						 * [1] Ref basis Libor curve
						 * [2] Ref basis Disc curve
						 * [3] DCC for basis rate
						 * [4] DCC for Libor rate
						 * [5] Basis type (Spread,
						 *         Percentage) */
                    // Basis spread 
    const Array<Date>   &bsVolExpDates, // (I) 15. option exp dates
    const Array<double> &bsVolRates,    // (I) 16. volatilities
    const Array<int>    &bsVolFreqs,    // (I) 17. underlying rate freqs
    String              &bsVolType,     // (I) 18. percentage or bp vols
                    // KMrParam
    const Array<double> &bsModParams,	// (I) 19. model parameters
                    // KSmileParam
    const Array<double> &bsSmileParams, // (I) 20. q1, q2, qfsh, iter

    double              corrIRBS,       // (I) 21. IR-BS correlation
                    // Tree Parameter
    const Array<int>    &treeParams,    // (I) 22. ppy, smooth, numStd
                    // KResetBank
    const Array<SharedPointer<KRate> > &rsRates,// (I) 23. array of float rates 
    const Array<Date>   &rsDates,     // (I) 24. array of reset dates
    const Array<double> &rsValues)    // (I) 25. array of reset values
{
static	char routine[] = "KMarketEnv::KMarketEnv";

        int             i;
        Array<int>      zcInterps;
	
        // Default IR curve interp type
        for (i=0; i<zcTypes.size(); i++)
            zcInterps.push_back(GTO_LINEAR_INTERP);


	// 
	// Construct the IR environment
	//
	SetIRMktEnv(
        today,
        zcCurves, 
        zcTypes,
        zcNames,   
        zcInterps,   
        irVolExpDates, 
        irVolMatDates, 
        irVolFreqs,   
        irVolRates,  
        irVolType,  
        irModParams,
        tmxFlag,
        irSmileParams, 
        irbmkLiqDates,
        treeParams,
        rsRates,
        rsDates, 
        rsValues);


	//--------------------------------------------------------------
	//	Basis Curves
	//--------------------------------------------------------------

	if (mMktCurves.IsBasis())
	{
	    //
	    // Basis curve info
	    //
	    mMktCurves.mLiborCVName = bsInfo[0];
	    mMktCurves.mBSDiscCVName = bsInfo[1];
	    mMktCurves.mBasisDCC     = KDayCc(bsInfo[2].c_str());
	    mMktCurves.mLiborDCC     = KDayCc(bsInfo[3].c_str());

		switch (toupper((bsInfo[4])[0]))
		{
        case 'P':
			mMktCurves.mBSType = PER_SPREAD;
			break;
		case 'S':
			mMktCurves.mBSType = SUB_SPREAD;
			break;
	    case 'A':
			mMktCurves.mBSType = ADD_SPREAD;
			break;
		}

	    mMktCurves.mBSDelayShift = (bsInfo.size() > 5 ? 
					atof(bsInfo[5].c_str()) : 0e0);
 

	    //
	    // 	Basis vol
	    //
	    mBSVolDiag.BasisVolDiag(bsVolType.c_str(),
				    DppDateArrayToTDateVector(bsVolExpDates),
				    DppVectorFromCMLIBArray(bsVolRates),
				    bsVolFreqs[0],
				    mIRVolDiag);


	   //
	   //	Basis Model Parameters
	   //
	   mBSMrParam.ReadMagnet(DppVectorFromCMLIBArray(bsModParams),
			         DppVectorFromCMLIBArray(treeParams));

	   //
	   //	Basis Smile Parameters
	   //
 
	   mBSSmileParam.ReadMagnet(DppVectorFromCMLIBArray(bsSmileParams));

	
	   mIRBSCorr = corrIRBS;



	}  // end of basis

    mIsCredit = false;

	   
}



//
// Credit Constructor
//
KMarketEnv::KMarketEnv(
                    // KMarketCurves
    Date                   today,          // (I) 1. today's date
    const Array<const TZeroCurve*>     &zcCurves,      // (I) 2. ALIB zero curve objs
    const Array<int>       &zcTypes,       // (I) 3. zc types (KV_DIFF, ..)
    const Array<String>    &zcNames,       // (I) 4. array of zc names
    const Array<int>       &zcInterps,     // (I) 5. array of zc interps
                    // KVolDiag
    const Array<Date>      &irVolExpDates, // (I) 6. option exp dates
    const Array<Date>      &irVolMatDates, // (I) 7. underlying mat dates
    const Array<int>       &irVolFreqs,    // (I) 8. underlying rate freqs
    const Array<double>    &irVolRates,    // (I) 9. volatilities
    String                 &irVolType,     // (I) 10. percentage or bp vols
                    // KMrParam
    const Array<double>    &irModParams,   // (I) 11. n,beta,alpha,rho,bbone
                    // KSmileParam
    String                 &tmxFlag,       // (I) 12  'Y'/'N'
    const Array<double>    &irSmileParams, // (I) 13. q1, q2, qfsh, iter or 
                                           //         MultiQ Smile
    const Array<Date>      &irBmkLiqDates, // (I) 14  TMX liquid benchmark dates 

                    // Credit curve
    String                 &crIRDiscName,  // (I) 15. IR discount curve
    double                 recovery,       // (I) 16. recovery rate
                    // Credit vol
    const Array<Date>      &crVolExpDates, // (I) 17. option exp dates
    const Array<Date>      &crVolMatDates, // (I) 18. underlying mat dates
    const Array<int>       &crVolFreqs,    // (I) 19. underlying rate freqs
    const Array<double>    &crVolRates,    // (I) 20. volatilities
    String                 &crVolType,     // (I) 21. percentage or bp vols
                    // KMrParam
    const Array<double>    &crModParams,   // (I) 22. model parameters
                    // KSmileParam
    const Array<double>    &crSmileParams, // (I) 23. q1, q2, qfsh, iter

    double                 corrIRCR,       // (I) 24. IR-CR correlation
                    // Tree Parameter
    const Array<int>       &treeParams,    // (I) 25. ppy, smooth, numStd
                    // KResetBank
    const Array<SharedPointer<KRate> > &rsRates,// (I) 26. array of float rates
    const Array<Date>      &rsDates,       // (I) 27. array of reset dates
    const Array<double>    &rsValues)      // (I) 28. array of reset values
{
static    char routine[] = "KMarketEnv::KMarketEnv";


    //
    // Construct the IR environment
    //

    if (toupper(tmxFlag[0] == 'Y'))
        throw KFailure("%s: Credit .\n",routine);
    SetIRMktEnv(
        today,
        zcCurves,
        zcTypes,
        zcNames,
        zcInterps,
        irVolExpDates,
        irVolMatDates,
        irVolFreqs,
        irVolRates,
        irVolType,
        irModParams,
        tmxFlag,
        irSmileParams,
        irBmkLiqDates,
        treeParams,
        rsRates,
        rsDates,
        rsValues);


    //--------------------------------------------------------------
    //    Credit Curves
    //--------------------------------------------------------------

    if (mMktCurves.IsCredit())
    {
        //
        // Credit curve info
        //
        mMktCurves.mIRDiscCVName = crIRDiscName;
        mMktCurves.mRecovery     = recovery;

        //--------------------------------------------------------------
        //    Vol Curves
        //--------------------------------------------------------------

        // Check size consistency
        ASSERT_OR_THROW(crVolExpDates.size() == crVolMatDates.size());
        ASSERT_OR_THROW(crVolExpDates.size() == crVolFreqs.size());
        ASSERT_OR_THROW(crVolExpDates.size() == crVolRates.size());

        mCRVolDiag = KVolDiag(DppDateArrayToTDateVector(crVolExpDates),
                              DppDateArrayToTDateVector(crVolMatDates),
                              DppVectorFromCMLIBArray(crVolFreqs),
                              DppVectorFromCMLIBArray(crVolRates));

        if (toupper(crVolType[0]) == 'N')
            mCRVolDiag.mVolType = NORMVOL;
        else if (toupper(crVolType[0]) == 'L')
            mCRVolDiag.mVolType = LOGVOL;
        else
            throw KFailure("%s: invalid credit vol type (%s).\n",
                           routine, crVolType.c_str());


        //--------------------------------------------------------------
        //    Model Parameters
        //--------------------------------------------------------------
        mCRMrParam.ReadMagnet(DppVectorFromCMLIBArray(crModParams),
                              DppVectorFromCMLIBArray(treeParams));

        //--------------------------------------------------------------
        //    Smile Parameters
        //--------------------------------------------------------------

        mCRSmileParam.ReadMagnet(DppVectorFromCMLIBArray(crSmileParams));


        mIRCRCorr = corrIRCR;

        mIsCredit = true;

    }  // end of basis



}




//---------------------------------------------------------------
//
ostream& 
operator<<(ostream& os, const KMarketEnv& mkt)
{

	os << "===========================================================================" << endl;
	os << "MARKET AND MODEL DATA " << endl;
	os << "===========================================================================" << endl;

	os << "---------------------------------------------------------------------------" << endl;
	os << "DATA: mktCurves:\n" << mkt.mMktCurves << endl;
	os << "---------------------------------------------------------------------------" << endl;
	os << "DATA: irVolDiag:\n" << mkt.mIRVolDiag << endl;
	os << "---------------------------------------------------------------------------" << endl;
	os << "DATA: irMrParam:\n" << mkt.mIRMrParam << endl;
	os << "---------------------------------------------------------------------------" << endl;
	os << "DATA: irSmileParam:\n" << mkt.mIRSmileParam << endl;
	os << "---------------------------------------------------------------------------" << endl;

	if (mkt.mMktCurves.mIsBasis)
	{
	    os << "DATA: bsVolDiag:\n" << mkt.mBSVolDiag << endl;
	    os << "---------------------------------------------------------------------------" << endl;
	    os << "DATA: bsMrParam:\n" << mkt.mBSMrParam << endl;
	    os << "---------------------------------------------------------------------------" << endl;
	    os << "DATA: bsSmileParam:\n" << mkt.mBSSmileParam << endl;
	    os << "---------------------------------------------------------------------------" << endl;
	    os << "DATA: irBsCorr:\n" << mkt.mIRBSCorr << endl;
	}

	if (mkt.IsCredit())
	{
	    os << "DATA: crVolDiag:\n" << mkt.mCRVolDiag << endl;
	    os << "---------------------------------------------------------------------------" << endl;
	    os << "DATA: crMrParam:\n" << mkt.mCRMrParam << endl;
	    os << "---------------------------------------------------------------------------" << endl;
	    os << "DATA: crSmileParam:\n" << mkt.mCRSmileParam << endl;
	    os << "---------------------------------------------------------------------------" << endl;
	    os << "DATA: irCRCorr:\n" << mkt.mIRCRCorr << endl;
	}


	os << "===========================================================================" << endl;
	os << "DATA: ResetBank:\n" << mkt.mVPResetBank << endl;
	os << "===========================================================================" << endl;

	return os;
}


//---------------------------------------------------------------
//
ostream& 
KMarketEnv::YacctionWrite(ostream& os, int indent)
{

	//
	// 1. Print the zc table
	//
    if (IsCredit())
	    mMktCurves.MAWCreditYacctionWrite(os);    // MAW curve files
    else
	    mMktCurves.DrWYacctionWrite(os);          // DrW curve files


    if (!IsCredit())     // Basis/IR
    {
        //
        // 2. Print the vol calibration index
        //
        mIRVolDiag.YacctionWrite(os);
   
        //
        // 3. Smile parameters
        //
        mIRSmileParam.YacctionWrite(os);
   
        //
        // 4. Model parameters
        //
        if (!mMktCurves.mIsBasis)
            mIRMrParam.YacctionWrite(os);
        else    // Basis
        {
            KMrParam        treeMrParam;    // full tree mr parameters
            treeMrParam = Correlate(mIRMrParam, mBSMrParam, mIRBSCorr);
            treeMrParam.YacctionWrite(os);

            os << "# Basis Number of factors" << endl;
            os << mBSMrParam.mNumFact << endl;

            // Basis zero curve info
            mMktCurves.BasisYacctionWrite(os);

            // Basis model parameter (backbone only)
            //
            mBSMrParam.BasisYacctionWrite(os);

            // Basis smile
            //
            mBSSmileParam.BasisYacctionWrite(os);

            // Basis spread vol curves
            //
            mBSVolDiag.BasisYacctionWrite(os);
        }
    }
    else        // Credit MAW deal file
    {
            // MAW tree model param
            mCRMrParam.MAWYacctionWrite(os);

            // MAW IR CET iteration
            mIRSmileParam.MAWYacctionWrite(os);

            // Credit reference curve
            os << "# IR reference discount curve name for CDS" << endl;
            os << mMktCurves.mIRDiscCVName << endl;
    }

    return (os);
}




//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
//x                                                                    x*/
//x     Set error log file                                             x*/
//x                                                                    x*/
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

MAGNET_DESCRIPTION("Set error log")
MAGNET_PARAMETER("", flag, "", "On if flag >=1")
MAGNET_PARAMETER("", fileName, "C:/error.log", "Error log file")
MAGNET_RESULT("The error log flag")

MAGNET_X_FUNCTION1(
	int,		ERR_LOG,		// Wrapper type and name
	int	flag,		=0)			// (I) 1. error flag
{
static	char routine[] = "MAGNET::ERR_LOG";
	
    char *fileName = NULL;

	// Enable error logging, and define the error messege output file.
	//
	if (flag) 
    {
        fileName=GtoErrMsgGetFileName();
        if (fileName != NULL)
		    DppSetToUseGtoErrMsg(fileName, TRUE);
        else
		    DppSetToUseGtoErrMsg("C:/error.log", TRUE);
    }
	else 
		DppErrMsgSet(DPP_ERR_MSG_OFF);

	return flag;

}



//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
//x                                                                    x*/
//x     Addin version                                                  x*/
//x                                                                    x*/
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

MAGNET_DESCRIPTION("Addin version")
MAGNET_RESULT("The addin version")

MAGNET_FUNCTION0(
	String,	VERSION)		// Wrapper type and name
{
static	char routine[] = "MAGNET::VERSION";

	String version = "SY " + String(_VERSION_) 
            + " COMPILED on "  + __DATE__ + " " + __TIME__
			+ " - (Linking to " + GtoVersion() 
			+ " + Magnet Version " + MagnetVersion() + ")";
	
	return version;
}




//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
//x                                                                    x*/
//x     KVP classes wrappers                                           x*/
//x                                                                    x*/
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//----------------------------------------------------------------------
// Construct a fixed coupon rate object
//

MAGNET_DESCRIPTION("Construct a fixed coupon rate")
MAGNET_PARAMETER("", name, "fixRate", "Name of the object.")
MAGNET_PARAMETER("", coupon, "", "Coupon rate")
MAGNET_RESULT("A Fixed Rate Handle")
MAGNET_SEE_ALSO("FLOATRATE", "See the floating rate")

MAGNET_X_FUNCTION2(
	SharedPointer<KRate>, FIXEDRATE,	// Wrapper type and name

	String		name,	="",		  // (I) 1. Name of the object

	double		coupon, MAGNET_MANDATORY) // (I) 2. Coupon rate
{
static	char routine[] = "MAGNET:FIXEDRATE";

	KRate *fixedRate = NULL;

  try {
	fixedRate = new KRate(coupon);
	ASSERT_OR_THROW(fixedRate);

	// If name is empty, then use pointer address
	if (name.empty())
		name = "fixRate_" + String(format("%p", fixedRate));

	fixedRate->SetName(name.c_str());

	return Raw2SharedPointer(fixedRate);
  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nFIXEDRATE: failed.\n");
	return 0;
  }

}



//----------------------------------------------------------------------
// Construct a KRate object with 0 spot offset.
//

MAGNET_DESCRIPTION("Construct a floating rate")
MAGNET_PARAMETER("", name, "floatRate", "Name of the object.")
MAGNET_PARAMETER("", maturity, "", "The rate maturity, e.g. 3M, 10, etc.")
MAGNET_PARAMETER("", frequency, "", "The rate frequency, e.g. 3M")
MAGNET_PARAMETER("", dayCC, "", "The rate day count convention, ACT/360, ACT/365, 30/360, ...")
MAGNET_PARAMETER("", idxZcName, "", "The rate index curve name")
MAGNET_RESULT("A Floating Rate Handle")
MAGNET_SEE_ALSO("FIXEDRATE", "See the fixed rate")

MAGNET_X_FUNCTION5(
	SharedPointer<KRate>, 	FLOATRATE,	// Wrapper type and name

					// (I) 1. Name of the object
	String		name,		="", 

					// (I) 2. Maturity
	const TDateInterval&	maturity,	MAGNET_MANDATORY,	
					// (I) 3. Frequency
	const TDateInterval&	frequency,	MAGNET_MANDATORY,
					// (I) 4. DCC
	String&			dayCc,		MAGNET_MANDATORY,	
					// (I) 5. Forward curve
	String&			idxZcName,	MAGNET_MANDATORY)	
{
static	char routine[] = "MAGNET::FLOATRATE";

	KRate		*floatRate= NULL;
	KDateInterval	spotOffset(0e0);
	double		spread = 0e0;
	double		weight = 1e0;

  try {

	// Create object
	floatRate = new KRate (
		idxZcName.c_str(),
		KDateInterval(maturity),
		KDateInterval(frequency),
		KDayCc(dayCc.c_str()),
		spotOffset,
		spread,
		weight);
	ASSERT_OR_THROW(floatRate);


	// If name is empty, then use pointer address
	if (name.empty())
		name = "floatRate_" + String(format("%p", floatRate));

	floatRate->SetName(name.c_str());

	return Raw2SharedPointer(floatRate);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nFLOATRATE: failed.\n");
	return 0;
  }

}





//----------------------------------------------------------------------
// Construct a KVPWBundle object.
//
MAGNET_DESCRIPTION("Construct a weighted sum of instruments")
MAGNET_PARAMETER("", name, "sum", "Name of the object.")
MAGNET_PARAMETER("", instruments, "", "Array of instrument objects")
MAGNET_PARAMETER("", weights, "", "Array of instrument weights")
MAGNET_PARAMETER("", discZcName, "Discount curve of the 1st dependency", "The discount zero curve name")
MAGNET_RESULT("A Sum Handle")
MAGNET_X_FUNCTION4(
	SharedPointer<KVPWBundle>,  SUM,	// Wrapper type and name
		
					// (I) 1. Name of the object
	String		name,		="",

					// (I) 2. array of instruments	  
	const Array<SharedPointer<KVPInstr> >&	instruments, MAGNET_MANDATORY,
					// (I) 3. array of instr weights
	const Array<double>&		weights,	MAGNET_MANDATORY,
					// (I) 4. discount curve name
	String		     		discZcName, 	= "Default") 
{
static	char routine[] = "MAGNET::SUM";

	KVPWBundle	*wb = new KVPWBundle("Sum");


	size_t		idx = 0;

  try {

	// Check size consistency
	ASSERT_OR_THROW(instruments.size() > 0);
	ASSERT_OR_THROW(instruments.size() == weights.size());

	// Set instrument name

	// Add all instruments
	for (idx = 0; idx < instruments.size(); ++idx)
		wb->AddDepWeight(instruments[idx], weights[idx]);

	//
	// Set discount curve
	// This is optional, just to be consistent with vpbase.
	//
        if (discZcName != "Default")
                wb->SetDiscName(discZcName.c_str());
	else	// set the discount curve to be that of the 1st component.
		wb->SetDiscName((instruments[0]->GetDiscName()).c_str());

	// If name is empty, then use pointer address
	if (name.empty())
		name = "sum_" + String(format("%p", wb));

	wb->SetName(name.c_str());

	return Raw2SharedPointer(wb);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nSUM: failed.\n");
	return 0;
  }

}



//----------------------------------------------------------------------
// Construct a fixed cash flows list.
//
MAGNET_DESCRIPTION("Construct cash flows")
MAGNET_PARAMETER("", name, "cashflows", "Name of the object.")
MAGNET_PARAMETER("", dates, "", "Array of cash flow dates")
MAGNET_PARAMETER("", amounts, "", "Array of cash flow amounts")
MAGNET_PARAMETER("", discZcName, "", "The discount zero curve name")
MAGNET_RESULT("A Cash Flow Handle")

MAGNET_X_FUNCTION4(
	SharedPointer<KVPCashFlows>,	CASHFLOW,  // Wrapper type and name

				// (I) 1. Name of the object
	String	name,		="",	 

				// (I) 2. array of cash flow dates
	const Array<Date>&   dates,	MAGNET_MANDATORY, 
				// (I) 3. array of cash flow amounts
	const Array<double>& amounts, 	MAGNET_MANDATORY, 
				// (I) 4. discount curve name
	String&		     discZcName, MAGNET_MANDATORY) 
{
static	char routine[] = "MAGNET::CASHFLOW";

	KVPCashFlows	*cf = NULL;

  try {
	

	ASSERT_OR_THROW(dates.size() == amounts.size());

	//
	// Construct object 
	//
	cf = new KVPCashFlows(
		name.c_str(),
		DppDateArrayToTDateVector(dates),	// cash flow dates
		DppVectorFromCMLIBArray(amounts),	// cash flow amounts
		discZcName.c_str());
	ASSERT_OR_THROW(cf);

	// If name is empty, then use pointer address
	if (name.empty())
		name = "sum_" + String(format("%p", cf));

	cf->SetName(name.c_str());

	return Raw2SharedPointer(cf);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nCASHFLOW: failed.\n");
	return 0;
  }

}





//**********************************************************************
//
// Floating leg constructors
//
//**********************************************************************


//----------------------------------------------------------------------
// Construct a KVPFloatLeg object (simple)
//
MAGNET_DESCRIPTION("Construct a floating leg with simple schedule")
MAGNET_PARAMETER("", name, "floatLeg", "Name of the object.")
MAGNET_PARAMETER("", startDate, "", "Start date")
MAGNET_PARAMETER("", maturityDate,  "", "Maturity date")
MAGNET_PARAMETER("", frequency, "", "Reset and payment frequency, e.g. 3M")
MAGNET_PARAMETER("", dayCc,     "", "Payment day count convention, e.g. ACT/360, 30/360, ...")
MAGNET_PARAMETER("", stubConv, "", " \
<P>Stub convention. </P> \
<OL> \
<LI>Floating rate: </LI></OL> \
<UL> \
<LI>(N)ONE \373 Fully accrued payment on the next payment date.</LI> \
<LI>(B)OND \373 Fully accrued payment received on the next payment date minus a partially accrued payment from the last reset to current paid on the current date.</LI> \
<LI>(S)IMPLE - Partial accrue of past reset for the remaining period, paid at the next payment date.</LI> \
<LI>(P)AR - Par rate accrued for the remaining period, and paid at the next payment date, priced at par.</LI> \
<LI>(O)N_SPOT_SIMPLE - Spot reset of the same floating rate (e.g. 3m Libor) accrued for the remaining period, and paid at the next payment date.</LI></UL> \
<P>NONE, BOND, and SIMPLE stubs require past reset information.</P> \
<P> </P> \
<OL> \
<LI>Fixed rate:</LI></OL> \
<UL> \
<LI>(N)ONE - Full coupon payment on the next payment date.</LI> \
<LI>(B)OND - Full coupon received on the next payment date minus a partially accrued payment from the last reset to current on the current date.</LI> \
<LI>(S)IMPLE - Partial accrue of coupon for the remaining period, pay at the next payment date.</LI></UL> \
")
MAGNET_PARAMETER("", rate,     "", "Pay rate object")
MAGNET_PARAMETER("", formula,  "nil", "Payment formula")
MAGNET_PARAMETER("", discZcName, "", "The discount zero curve name")
MAGNET_RESULT("A Floating Leg Handle")
MAGNET_SEE_ALSO("FLOATLEG_GENERAL", "See the floating leg with arbitrary schedule")

MAGNET_X_FUNCTION9(
	SharedPointer<KVPFloatLeg>, FLOATLEG_SIMPLE, // Wrap type and name

	String	name,		="",       	   // (I) 1. Name of the object

	Date	startDate,	MAGNET_MANDATORY,  // (I) 2. start date
	Date	maturityDate,	MAGNET_MANDATORY,  // (I) 3. maturity date
  const TDateInterval& frequency, MAGNET_MANDATORY,  // (I) 4. freqency
	String& dayCc,		MAGNET_MANDATORY,  // (I) 5. DCC
	String& stubConv,	MAGNET_MANDATORY,  // (I) 6. stub convention
  	SharedPointer<KRate> rate, MAGNET_MANDATORY,  // (I) 7. underlying rate
	String	formula,	="nil",		   // (I) 8. payoff formula
	String& discZcName, 	MAGNET_MANDATORY)  // (I) 9. discount curve
{
static	char routine[] = "MAGNET::FLOATLEG_SIMPLE";

	TBoolean	stubAtEnd = FALSE;
	KVPFloatLeg	*floatLeg = NULL;

  try {
	
	// Create floatLeg
	floatLeg = new KVPFloatLeg(
		name.c_str(),
		startDate,
		maturityDate,
		KDateInterval(frequency),
		KDayCc(dayCc.c_str()),
		KStubConv(stubConv.c_str()),
		stubAtEnd,
		rate,
		formula.c_str(),
		discZcName.c_str());   // discount curve
	ASSERT_OR_THROW(floatLeg);
	
	// If name is empty, then use pointer address
	if (name.empty())
		name = "floatLeg_" + String(format("%p", floatLeg));

	floatLeg->SetName(name.c_str());

	return Raw2SharedPointer(floatLeg);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nFLOATLEG_SIMPLE: failed.\n");
	return 0;
  }

}




//----------------------------------------------------------------------
// Construct a KVPFloatLeg object (complex arbitrary resets)
// Reset effective dates are given explicitly.

MAGNET_DESCRIPTION("Construct a floating leg with arbitrary schedule")
MAGNET_PARAMETER("", name, "floatLeg", "Name of the object.")
MAGNET_PARAMETER("", resetDates,    "", "Array of reset dates")
MAGNET_PARAMETER("", resetEffDates, "", "Array of reset effective dates")
MAGNET_PARAMETER("", accStDates,    "", "Array of accrual start dates")
MAGNET_PARAMETER("", accEndDates,   "", "Array of accrual end dates")
MAGNET_PARAMETER("", payDates,      "", "Array of payment dates")
MAGNET_PARAMETER("", notionals,     "", "Array of notionals")
MAGNET_PARAMETER("", formula, 	    "nil", "(Array of payment formula or a single formula applied to all the rate payments in each reset period.")
MAGNET_PARAMETER("", rates, 	    "", "Array of pay rates for each reset period or a single rate definition applied to all the resets")
MAGNET_PARAMETER("", dayCc, 	    "", "Payment day count convention")
MAGNET_PARAMETER("", stubConv,      "", "Stub convention")
MAGNET_PARAMETER("", discZcName,    "", "The discount zero curve name")
MAGNET_RESULT("A Floating Leg Handle")
MAGNET_SEE_ALSO("FLOATLEG_SIMPLE", "See the floating leg with simple schedule")
//
MAGNET_X_FUNCTION12(
	SharedPointer<KVPFloatLeg>, FLOATLEG_GENERAL, // Wrap type and name

					// (I) 1. Name of the object	
	String			name,		="",

					// (I) 2. array of reset dates	
	const Array<Date>&	resetDates,	MAGNET_MANDATORY,
					// (I) 3. array of reset effective dates
	const Array<Date>&	resetEffDates,	MAGNET_MANDATORY,,
					// (I) 4. array of accrual start dates
	const Array<Date>&	accStDates,	MAGNET_MANDATORY,
					// (I) 5. array of accrual end dates
	const Array<Date>&	accEndDates,	MAGNET_MANDATORY,
					// (I) 6. array of pay dates
	const Array<Date>&	payDates,	MAGNET_MANDATORY,
					// (I) 7. array of notionals
	const Array<double>&	notionals,	MAGNET_MANDATORY,
					// (I) 8. payoff formula
	const Array<String>&	formula,  	=nilArray,
	                                // (I) 9. underlying rate index
	const	Array<SharedPointer<KRate> >& rates,	MAGNET_MANDATORY,
					// (I) 10. DCC
	String&			dayCc,		MAGNET_MANDATORY,
					// (I) 11. stub convention
	String&			stubConv,	MAGNET_MANDATORY,
					// (I) 12. discount curve
	String&			discZcName,	MAGNET_MANDATORY)
{
static	char routine[] = "MAGNET::FLOATLEG_GENERAL";

	TBoolean	stubAtEnd = FALSE;

	KVPFloatLeg	*floatLeg = NULL;

	KVector(SharedPointer<KRate>) vRates;
	KVector(String) vFormula;

  try {
	//
	// 1. Single rate and single pay formula
	//
	if (formula.size() == 1 &&
	    rates.size()   == 1)
		floatLeg = new KVPFloatLeg(
			name.c_str(),
			DppDateArrayToTDateVector(resetDates),
			DppDateArrayToTDateVector(resetEffDates),
			DppDateArrayToTDateVector(accStDates),
			DppDateArrayToTDateVector(accEndDates),
			DppDateArrayToTDateVector(payDates),
			DppVectorFromCMLIBArray(notionals),
			KDayCc(dayCc.c_str()),
			KStubConv(stubConv.c_str()),
			rates[0],
			formula[0].c_str(),
			discZcName.c_str());   // discount curve
	//
	// 2. Single rate and array of pay formula
	//
	else if (formula.size() > 1 &&
		 rates.size()  == 1)
		floatLeg = new KVPFloatLeg(
			name.c_str(),
			DppDateArrayToTDateVector(resetDates),
			DppDateArrayToTDateVector(resetEffDates),
			DppDateArrayToTDateVector(accStDates),
			DppDateArrayToTDateVector(accEndDates),
			DppDateArrayToTDateVector(payDates),
			DppVectorFromCMLIBArray(notionals),
			DppVectorFromCMLIBArray(formula),
			KDayCc(dayCc.c_str()),
			KStubConv(stubConv.c_str()),
			rates[0],
			discZcName.c_str());   // discount curve
	//
	// 3. Array of rate and single pay formula
	//
	else if (formula.size() == 1 &&
		 rates.size()   > 1)
	{
		for (int idx=0; idx < rates.size(); idx++)
		{
			vFormula.insert(vFormula.end(), formula[0]);
			vRates.insert(vRates.end(), rates[idx]);
		}

		floatLeg = new KVPFloatLeg(
			name.c_str(),
			DppDateArrayToTDateVector(resetDates),
			DppDateArrayToTDateVector(resetEffDates),
			DppDateArrayToTDateVector(accStDates),
			DppDateArrayToTDateVector(accEndDates),
			DppDateArrayToTDateVector(payDates),
			DppVectorFromCMLIBArray(notionals),
			vFormula,
			vRates,
			KDayCc(dayCc.c_str()),
			KStubConv(stubConv.c_str()),
			discZcName.c_str());   // discount curve
	}
	//
	// 4. Array of rates and array of pay formula
	//
	else if (formula.size() > 1 &&
		 rates.size()   > 1)
	{
		for (int idx=0; idx < rates.size(); idx++)
			vRates.insert(vRates.end(), rates[idx]);

		floatLeg = new KVPFloatLeg(
			name.c_str(),
			DppDateArrayToTDateVector(resetDates),
			DppDateArrayToTDateVector(resetEffDates),
			DppDateArrayToTDateVector(accStDates),
			DppDateArrayToTDateVector(accEndDates),
			DppDateArrayToTDateVector(payDates),
			DppVectorFromCMLIBArray(notionals),
			DppVectorFromCMLIBArray(formula),
			vRates,
			KDayCc(dayCc.c_str()),
			KStubConv(stubConv.c_str()),
			discZcName.c_str());   // discount curve
	}
	else
	{
		throw KFailure("%s: invalid inputs: following "
			"floatLeg constructions are currently supported:\n"
			"1: Single rate object and single formula;\n"
			"2. Single rate object and array of formula;\n"
			"3. Array of rate objects and single formula;\n"
			"4. Array of rate objects and array of formula;\n");
	}


	ASSERT_OR_THROW(floatLeg);
	
	// If name is empty, then use pointer address
	if (name.empty())
		name = "floatLeg_" + String(format("%p", floatLeg));

	floatLeg->SetName(name.c_str());

	return Raw2SharedPointer(floatLeg);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nFLOATLEG_GENERAL: failed.\n");
	return 0;
  }

}





//----------------------------------------------------------------------
// Construct a KVPFloatLeg2Idx object (simple)
//
MAGNET_DESCRIPTION("Construct a floating leg with 2 rate indices with simple schedule")
MAGNET_PARAMETER("", name, "floatLeg2Idx", "Name of the object.")
MAGNET_PARAMETER("", startDate, "", "Start date")
MAGNET_PARAMETER("", maturityDate,  "", "Maturity date")
MAGNET_PARAMETER("", frequency, "", "Reset and payment frequency, e.g. 3M")
MAGNET_PARAMETER("", dayCc,     "", "Payment day count convention, e.g. ACT/360, 30/360, ...")
MAGNET_PARAMETER("", stubConv, "", " \
<P>Stub convention. </P> \
<OL> \
<LI>Floating rate: </LI></OL> \
<UL> \
<LI>(N)ONE \373 Fully accrued payment on the next payment date.</LI> \
<LI>(B)OND \373 Fully accrued payment received on the next payment date minus a partially accrued payment from the last reset to current paid on the current date.</LI> \
<LI>(S)IMPLE - Partial accrue of past reset for the remaining period, paid at the next payment date.</LI> \
<LI>(P)AR - Par rate accrued for the remaining period, and paid at the next payment date, priced at par.</LI> \
<LI>(O)N_SPOT_SIMPLE - Spot reset of the same floating rate (e.g. 3m Libor) accrued for the remaining period, and paid at the next payment date.</LI></UL> \
<P>NONE, BOND, and SIMPLE stubs require past reset information.</P> \
<P> </P> \
<OL> \
<LI>Fixed rate:</LI></OL> \
<UL> \
<LI>(N)ONE - Full coupon payment on the next payment date.</LI> \
<LI>(B)OND - Full coupon received on the next payment date minus a partially accrued payment from the last reset to current on the current date.</LI> \
<LI>(S)IMPLE - Partial accrue of coupon for the remaining period, pay at the next payment date.</LI></UL> \
")
MAGNET_PARAMETER("", rate1,     "", "Pay rate 1 object")
MAGNET_PARAMETER("", rate2,     "", "Pay rate 2 object")
MAGNET_PARAMETER("", formula,  "nil", "Payment formula")
MAGNET_PARAMETER("", discZcName, "", "The discount zero curve name")
MAGNET_RESULT("A Floating Leg Handle")
MAGNET_SEE_ALSO("FLOATLEG2_GENERAL", "See the floating leg with arbitrary schedule")

MAGNET_X_FUNCTION10(
	SharedPointer<KVPFloatLeg2Idx>, FLOATLEG2_SIMPLE, // Wrap type and name

	String	name,		="",   		   // (I) 1. Name of the object

	Date	startDate,	MAGNET_MANDATORY,  // (I) 2. start date
	Date	maturityDate,	MAGNET_MANDATORY,  // (I) 3. maturity date
	const   TDateInterval& frequency, MAGNET_MANDATORY,  // (I) 4. freqency
	String& dayCc,		MAGNET_MANDATORY,  // (I) 5. DCC
	String& stubConv,	MAGNET_MANDATORY,  // (I) 6. stub convention
  	SharedPointer<KRate> rate1, MAGNET_MANDATORY, // (I) 7. rate 1
  	SharedPointer<KRate> rate2, MAGNET_MANDATORY, // (I) 8. rate 2
	String	formula,	="nil",		   // (I) 9. payoff formula
	String& discZcName, 	MAGNET_MANDATORY)  // (I) 10. discount curve
{
static	char routine[] = "MAGNET::FLOATLEG2_SIMPLE";

	TBoolean	stubAtEnd = FALSE;
	KVPFloatLeg2Idx	*floatLeg2Idx = NULL;

  try {
	
	// Create floatLeg
	floatLeg2Idx = new KVPFloatLeg2Idx(
		name.c_str(),
		startDate,
		maturityDate,
		KDateInterval(frequency),
		KDayCc(dayCc.c_str()),
		KStubConv(stubConv.c_str()),
		stubAtEnd,
		rate1,
		rate2,
		formula.c_str(),
		discZcName.c_str());   // discount curve
	ASSERT_OR_THROW(floatLeg2Idx);
	
	// If name is empty, then use pointer address
	if (name.empty())
		name = "floatLeg2Idx_" + String(format("%p", floatLeg2Idx));

	floatLeg2Idx->SetName(name.c_str());

	return Raw2SharedPointer(floatLeg2Idx);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nFLOATLEG2_SIMPLE: failed.\n");
	return 0;
  }

}




//----------------------------------------------------------------------
// Construct a KVPFloatLeg object (complex arbitrary resets)
// Reset effective dates are given explicitly.

MAGNET_DESCRIPTION("Construct a floating leg with 2 rate indices and arbitrary schedule")
MAGNET_PARAMETER("", name, "floatLeg", "Name of the object.")
MAGNET_PARAMETER("", resetDates,    "", "Array of reset dates")
MAGNET_PARAMETER("", resetEffDates, "", "Array of reset effective dates")
MAGNET_PARAMETER("", accStDates,    "", "Array of accrual start dates")
MAGNET_PARAMETER("", accEndDates,   "", "Array of accrual end dates")
MAGNET_PARAMETER("", payDates,      "", "Array of payment dates")
MAGNET_PARAMETER("", notionals,     "", "Array of notionals")
MAGNET_PARAMETER("", formula, 	    "nil", "(Array of payment formula or a single formula applied to all the rate payments in each reset period.")
MAGNET_PARAMETER("", rates1, 	    "", "Array of first pay rates for each reset period or a single rate definition applied to all the resets")
MAGNET_PARAMETER("", rates2, 	    "", "Array of second pay rates for each reset period or a single rate definition applied to all the resets")
MAGNET_PARAMETER("", dayCc, 	    "", "Payment day count convention")
MAGNET_PARAMETER("", stubConv,      "", "Stub convention")
MAGNET_PARAMETER("", discZcName,    "", "The discount zero curve name")
MAGNET_RESULT("A Floating Leg Handle")
MAGNET_SEE_ALSO("FLOATLEG2_SIMPLE", "See the floating leg with simple schedule")
//
MAGNET_X_FUNCTION13(
	SharedPointer<KVPFloatLeg2Idx>, FLOATLEG2_GENERAL, // Wrap type and name

					// (I) 1. Name of the object	
	String			name,		="",

					// (I) 2. array of reset dates	
	const Array<Date>&	resetDates,	MAGNET_MANDATORY,
					// (I) 3. array of reset effective dates
	const Array<Date>&	resetEffDates,	MAGNET_MANDATORY,,
					// (I) 4. array of accrual start dates
	const Array<Date>&	accStDates,	MAGNET_MANDATORY,
					// (I) 5. array of accrual end dates
	const Array<Date>&	accEndDates,	MAGNET_MANDATORY,
					// (I) 6. array of pay dates
	const Array<Date>&	payDates,	MAGNET_MANDATORY,
					// (I) 7. array of notionals
	const Array<double>&	notionals,	MAGNET_MANDATORY,
					// (I) 8. payoff formula
	const Array<String>&		formula,  	=nilArray,
	                                // (I) 9. underlying rate index
	const	Array<SharedPointer<KRate> >& rates1,	MAGNET_MANDATORY,
	                                // (I) 10. underlying rate index
	const	Array<SharedPointer<KRate> >& rates2,	MAGNET_MANDATORY,
					// (I) 11. DCC
	String&			dayCc,		MAGNET_MANDATORY,
					// (I) 12. stub convention
	String&			stubConv,	MAGNET_MANDATORY,
					// (I) 13. discount curve
	String&			discZcName,	MAGNET_MANDATORY)
{
static	char routine[] = "MAGNET::FLOATLEG2_GENERAL";

	TBoolean	stubAtEnd = FALSE;

	KVPFloatLeg2Idx	*floatLeg2Idx = NULL;

	KVector(SharedPointer<KRate>) vRates1, vRates2;

  try {

	//
	// 1. Single rate and single pay formula
	//
	if (formula.size() == 1 &&
	    rates1.size()  == 1 &&
	    rates2.size()  == 1)
		floatLeg2Idx = new KVPFloatLeg2Idx(
			name.c_str(),
			DppDateArrayToTDateVector(resetDates),
			DppDateArrayToTDateVector(resetEffDates),
			DppDateArrayToTDateVector(accStDates),
			DppDateArrayToTDateVector(accEndDates),
			DppDateArrayToTDateVector(payDates),
			DppVectorFromCMLIBArray(notionals),
			KDayCc(dayCc.c_str()),
			KStubConv(stubConv.c_str()),
			rates1[0],
			rates2[0],
			formula[0].c_str(),
			discZcName.c_str());   // discount curve
	//
	// 2. Single rate and array of pay formula
	//
	else if (formula.size() > 1 &&
		 rates1.size() == 1 &&
		 rates2.size() == 1)
		floatLeg2Idx = new KVPFloatLeg2Idx(
			name.c_str(),
			DppDateArrayToTDateVector(resetDates),
			DppDateArrayToTDateVector(resetEffDates),
			DppDateArrayToTDateVector(accStDates),
			DppDateArrayToTDateVector(accEndDates),
			DppDateArrayToTDateVector(payDates),
			DppVectorFromCMLIBArray(notionals),
			DppVectorFromCMLIBArray(formula),
			KDayCc(dayCc.c_str()),
			KStubConv(stubConv.c_str()),
			rates1[0],
			rates2[0],
			discZcName.c_str());   // discount curve
	//
	// 3. Array of rates and array of pay formula
	//
	else if (formula.size() > 1 &&
		 rates1.size()  > 1 &&
		 rates2.size()  > 1)
	{
		for (int idx=0; idx < rates1.size(); idx++)
		{
			vRates1.insert(vRates1.end(), rates1[idx]);
			vRates2.insert(vRates2.end(), rates2[idx]);
		}

		floatLeg2Idx = new KVPFloatLeg2Idx(
			name.c_str(),
			DppDateArrayToTDateVector(resetDates),
			DppDateArrayToTDateVector(resetEffDates),
			DppDateArrayToTDateVector(accStDates),
			DppDateArrayToTDateVector(accEndDates),
			DppDateArrayToTDateVector(payDates),
			DppVectorFromCMLIBArray(notionals),
			DppVectorFromCMLIBArray(formula),
			vRates1,
			vRates2,
			KDayCc(dayCc.c_str()),
			KStubConv(stubConv.c_str()),
			discZcName.c_str());   // discount curve
	}
	else
	{
		throw KFailure("%s: invalid inputs: following "
			"floatLeg constructions are currently supported:\n"
			"1: Single rate object and single formula;\n"
			"2. Single rate object and array of formula;\n"
			"3. Array of rate objects and array of formula;\n");
	}


	ASSERT_OR_THROW(floatLeg2Idx);
	
	// If name is empty, then use pointer address
	if (name.empty())
		name = "floatLeg2Idx_" + String(format("%p", floatLeg2Idx));

	floatLeg2Idx->SetName(name.c_str());

	return Raw2SharedPointer(floatLeg2Idx);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nFLOATLEG2_GENERAL: failed.\n");
	return 0;
  }

}


MAGNET_DESCRIPTION("Construct a floating leg with 3 rate indices and arbitrary schedule")
MAGNET_PARAMETER("", name, "floatLeg", "Name of the object.")
MAGNET_PARAMETER("", resetDates,    "", "Array of reset dates")
MAGNET_PARAMETER("", resetEffDates, "", "Array of reset effective dates")
MAGNET_PARAMETER("", accStDates,    "", "Array of accrual start dates")
MAGNET_PARAMETER("", accEndDates,   "", "Array of accrual end dates")
MAGNET_PARAMETER("", payDates,      "", "Array of payment dates")
MAGNET_PARAMETER("", notionals,     "", "Array of notionals")
MAGNET_PARAMETER("", formula, 	    "nil", "(Array of payment formula or a single formula applied to all the rate payments in each reset period.")
MAGNET_PARAMETER("", rates1, 	    "", "First pay rate for each reset period or a single rate definition applied to all the resets")
MAGNET_PARAMETER("", rates2, 	    "", "Second pay rate for each reset period or a single rate definition applied to all the resets")
MAGNET_PARAMETER("", rates3, 	    "", "Third pay rate for each reset period or a single rate definition applied to all the resets")
MAGNET_PARAMETER("", dayCc, 	    "", "Payment day count convention")
MAGNET_PARAMETER("", stubConv,      "", "Stub convention")
MAGNET_PARAMETER("", discZcName,    "", "The discount zero curve name")
MAGNET_RESULT("A Floating Leg Handle")
MAGNET_SEE_ALSO("FLOATLEG2_GENERAL", "See the floating leg with simple schedule")
//
MAGNET_X_FUNCTION14(
	SharedPointer<KVPFloatLeg3Idx>, FLOATLEG3_GENERAL, // Wrap type and name

					// (I) 1. Name of the object	
	String			name,		="",

					// (I) 2. array of reset dates	
	const Array<Date>&	resetDates,	MAGNET_MANDATORY,
					// (I) 3. array of reset effective dates
	const Array<Date>&	resetEffDates,	MAGNET_MANDATORY,,
					// (I) 4. array of accrual start dates
	const Array<Date>&	accStDates,	MAGNET_MANDATORY,
					// (I) 5. array of accrual end dates
	const Array<Date>&	accEndDates,	MAGNET_MANDATORY,
					// (I) 6. array of pay dates
	const Array<Date>&	payDates,	MAGNET_MANDATORY,
					// (I) 7. array of notionals
	const Array<double>&	notionals,	MAGNET_MANDATORY,
					// (I) 8. payoff formula
	const Array<String>&		formula,  	=nilArray,
	                                // (I) 9. underlying rate index
	const	Array<SharedPointer<KRate> >& rates1,	MAGNET_MANDATORY,
	                                // (I) 10. underlying rate index
	const	Array<SharedPointer<KRate> >& rates2,	MAGNET_MANDATORY,
	                                // (I) 11. underlying rate index
	const	Array<SharedPointer<KRate> >& rates3,	MAGNET_MANDATORY,
					// (I) 12. DCC
	String&			dayCc,		MAGNET_MANDATORY,
					// (I) 13. stub convention
	String&			stubConv,	MAGNET_MANDATORY,
					// (I) 14. discount curve
	String&			discZcName,	MAGNET_MANDATORY)
{
static	char routine[] = "MAGNET::FLOATLEG3_GENERAL";

	TBoolean	stubAtEnd = FALSE;

	KVPFloatLeg3Idx	*floatLeg3Idx = NULL;

	KVector(SharedPointer<KRate>) vRates1, vRates2, vRates3;

  try {

	//
	// Single rate and array of pay formula
	//
	if (formula.size() > 1 &&
	    rates1.size() == 1 &&
	    rates2.size() == 1 &&
	    rates3.size() == 1)
	  {
	    floatLeg3Idx = new KVPFloatLeg3Idx
	      (
	       name.c_str(),
	       DppDateArrayToTDateVector(resetDates),
	       DppDateArrayToTDateVector(resetEffDates),
	       DppDateArrayToTDateVector(accStDates),
	       DppDateArrayToTDateVector(accEndDates),
	       DppDateArrayToTDateVector(payDates),
	       DppVectorFromCMLIBArray(notionals),
	       DppVectorFromCMLIBArray(formula),
	       KDayCc(dayCc.c_str()),
	       KStubConv(stubConv.c_str()),
	       rates1[0],
	       rates2[0],
	       rates3[0],
	       discZcName.c_str()
	       );
	  }
	else
	  {
	    throw KFailure("%s: invalid inputs: following "
			   "floatLeg constructions are currently supported:\n"
			   "1. Single rate object and array of formula;\n"
			   );
	  }
	

	ASSERT_OR_THROW(floatLeg3Idx);
	
	// If name is empty, then use pointer address
	if (name.empty())
		name = "floatLeg3Idx_" + String(format("%p", floatLeg3Idx));

	floatLeg3Idx->SetName(name.c_str());

	return Raw2SharedPointer(floatLeg3Idx);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nFLOATLEG3_GENERAL: failed.\n");
	return 0;
  }
}




MAGNET_DESCRIPTION("Construct a protection leg")
MAGNET_PARAMETER("", name, "protLeg", "Name of the object.")
MAGNET_PARAMETER("", startDate,    "", "start date")
MAGNET_PARAMETER("", endDate,      "", "end date")
MAGNET_PARAMETER("", notional,     "", "notional")
MAGNET_PARAMETER("", recovery,     "nil", "recovery rate, nil for default.")
MAGNET_PARAMETER("", payType,       "", "pay type: upon default or at maturity")
MAGNET_PARAMETER("", discZcName,    "", "The discount zero curve name")
MAGNET_RESULT("A Protection Leg Handle")
//
MAGNET_X_FUNCTION7(
    SharedPointer<KVPProtLeg>, PROTLEG, // Wrap type and name

                                    // (I) 1. Name of the object
    String              name,       ="",

                                    // (I) 2. start date
    Date          	startDate,      MAGNET_MANDATORY,
                                    // (I) 3. end date
    Date          	endDate,    MAGNET_MANDATORY,,
                                    // (I) 4. notional
    double        	notional,   MAGNET_MANDATORY,
                                    // (I) 5. recovery
    String       	recovery,  ="nil",
                                    // (I) 6. pay type
    String       	pay_def_or_mat, ="PAY_DEF", 
                                    // (I) 7. discount curve
    String&             discZcName, MAGNET_MANDATORY)
{
static  char routine[] = "MAGNET::PROTLEG";

    KVPProtLeg  *protLeg  = NULL;

    KProtPayConv payType;

  try {

    if (strchr(pay_def_or_mat.c_str(), 'd') != NULL ||
        strchr(pay_def_or_mat.c_str(), 'D') != NULL)
        payType = PAY_DEF;
    else if (
        strchr(pay_def_or_mat.c_str(), 'm') != NULL ||
        strchr(pay_def_or_mat.c_str(), 'M') != NULL)
        payType = PAY_MAT;
    else
    {
        throw KFailure("%s: invalid input for pay type (%s).\n",
                       routine, pay_def_or_mat.c_str());
    }


    protLeg = new KVPProtLeg(
           name.c_str(),
           startDate,
           endDate,
           notional,
           recovery.c_str(),
           payType,
           discZcName.c_str());


    ASSERT_OR_THROW(protLeg);

    // If name is empty, then use pointer address
    if (name.empty())
        name = "protLeg_" + String(format("%p", protLeg));

    protLeg->SetName(name.c_str());

    return Raw2SharedPointer(protLeg);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nPROTLEG: failed.\n");
	return 0;
  }

}



MAGNET_DESCRIPTION("Construct a default exposure")
MAGNET_PARAMETER("", name, "defExppsure", "Name of the object.")
MAGNET_PARAMETER("", startDate,    "", "start date")
MAGNET_PARAMETER("", endDate,      "", "end date")
MAGNET_PARAMETER("", settleDate,   "", "settle date")
MAGNET_PARAMETER("", frequency,    "", "frequency")
MAGNET_PARAMETER("", recovery,     "", "Default recovery rate")
MAGNET_PARAMETER("", underlying,   "", "Underlying instrument")
MAGNET_PARAMETER("", discZcName,    "", "The discount zero curve name")
MAGNET_RESULT("A Default Protection Handle")
//
MAGNET_FUNCTION7(
    SharedPointer<KVPDefProtect>, DEFEXPOSURE_SIMPLE, // Wrap type and name

    String                  name,                // (I) 1. Name of the object
    Date          	    startDate,           // (I) 2. start date
    Date          	    endDate,             // (I) 3. end date 
    TDateInterval&          frequency,           // (I) 4. freqency 
    double	            recovery,		 // (I) 5. recovery 
    SharedPointer<KVPAtom>  underlying,          // (I) 6. underlying
    String&                 discZcName)          // (I) 7. discount curve
{
static  char routine[] = "MAGNET::DEFEXPOSURE_SIMPLE";

    KVPDefProtect  *defExpos = NULL;

  try {

    defExpos = new KVPDefProtect(
           name.c_str(),
           startDate,
           endDate,
           frequency,
           DEF_EXPOSURE,
           recovery,
           discZcName.c_str());


    ASSERT_OR_THROW(defExpos);

    // Add underlying
    defExpos->AddDep(underlying);

    return Raw2SharedPointer(defExpos);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nDEFEXPOSURE_SIMPLE: failed.\n");
	return 0;
  }


}



//----------------------------------------------------------------------
// Construct a default exposure with general schedule
//
MAGNET_DESCRIPTION("Construct a default exposure object with general schedule")
MAGNET_PARAMETER("", name, "DefProtect", "Name of the object.")
MAGNET_PARAMETER("", startDates,    "", "Array of start dates")
MAGNET_PARAMETER("", endDates,      "", "Array of end dates")
MAGNET_PARAMETER("", settleDates,   "", "Array of underlying settlement dates")
MAGNET_PARAMETER("", recovery,      "", "Default recovery rate")
MAGNET_PARAMETER("", underlying,    "", "Underlying instrument")
MAGNET_PARAMETER("", discZcName,    "", "The discount zero curve name")
MAGNET_RESULT("A DefProtect Handle")
MAGNET_SEE_ALSO("DEFEXPOSURE_SIMPLE", "See the defexposure with simple schedule")

MAGNET_X_FUNCTION7(
        SharedPointer<KVPDefProtect>, DEFEXPOSURE_GENERAL, // Wrap type and name

					// (I) 1. name of the object
	String		name,			="",		

					// (I) 2. start dates
  const Array<Date>&	startDates,		MAGNET_MANDATORY,
					// (I) 3. end dates
  const Array<Date>&	endDates,		MAGNET_MANDATORY,	
					// (I) 4. settlement dates
  const Array<Date>&	settleDates,		MAGNET_MANDATORY,	
					// (I) 5. default recovery
        double	        recovery,		MAGNET_MANDATORY,
					// (I) 6. underlying asset
	SharedPointer<KVPAtom>	underlying,	MAGNET_MANDATORY,
					// (I) 7. discount curve
	String&		discZcName,		MAGNET_MANDATORY)		
{
static	char routine[] = "MAGNET::DEFEXPOSURE_GENERAL";

	//int		idx;

        KVPDefProtect   *defExpos = NULL;

  try {


	// Create object
	defExpos = new KVPDefProtect(
		name.c_str(),
		DppDateArrayToTDateVector(startDates),
		DppDateArrayToTDateVector(endDates),
		DppDateArrayToTDateVector(settleDates),
                DEF_EXPOSURE,
		recovery,
		discZcName.c_str());

	ASSERT_OR_THROW(defExpos);

	// Add underlying
	defExpos->AddDep(underlying);

	// If name is empty, then use pointer address
	if (name.empty())
		name = "DefExposure" + String(format("%p", defExpos));

	defExpos->SetName(name.c_str());

	return Raw2SharedPointer(defExpos);


  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nDEFEXPOSURE_GENERAL: failed.\n");
	return 0;
  }

}






//----------------------------------------------------------------------
// Construct a default knock-in with general schedule and knock-in
//
MAGNET_DESCRIPTION("Construct a default knock-in object with general schedule")
MAGNET_PARAMETER("", name, "DefKnockIn", "Name of the object.")
MAGNET_PARAMETER("", startDates,    "", "Array of start dates")
MAGNET_PARAMETER("", endDates,      "", "Array of end dates")
MAGNET_PARAMETER("", settleDates,   "", "Array of underlying settlement dates")
MAGNET_PARAMETER("", rebates,       "empty", "Array of rebate values")
MAGNET_PARAMETER("", underlying,    "", "Underlying instrument")
MAGNET_PARAMETER("", discZcName,    "", "The discount zero curve name")
MAGNET_RESULT("A DefProtect Handle")

MAGNET_X_FUNCTION7(
        SharedPointer<KVPDefProtect>, DEFKNOCKIN_GENERAL, // Wrap type and name

					// (I) 1. name of the object
	String		name,			="",		

					// (I) 2. start dates
  const Array<Date>&	startDates,		MAGNET_MANDATORY,
					// (I) 3. end dates
  const Array<Date>&	endDates,		MAGNET_MANDATORY,	
					// (I) 4. settlement dates
  const Array<Date>&	settleDates,		MAGNET_MANDATORY,	
					// (I) 5. default rebate
  const Array<double>&	rebates,		=Array<double>(),
					// (I) 6. underlying asset
	SharedPointer<KVPAtom>	underlying,	MAGNET_MANDATORY,
					// (I) 7. discount curve
	String&		discZcName,		MAGNET_MANDATORY)		
{
static	char routine[] = "MAGNET::DEFKNOCKIN_GENERAL";

	int		idx;

        KVPDefProtect   *defKnockIn = NULL;
	Array<double>	tmpRebates;

  try {

	//
	// Initialize rebates
	//
	if (rebates.empty())
	{
		tmpRebates.resize(startDates.size());
		
		for (idx=0; idx < startDates.size(); idx++)
			tmpRebates[idx] = 0e0;	
	}
	else
		tmpRebates = rebates;


	// Create object
	defKnockIn = new KVPDefProtect(
		name.c_str(),
		DppDateArrayToTDateVector(startDates),
		DppDateArrayToTDateVector(endDates),
		DppDateArrayToTDateVector(settleDates),
		DppVectorFromCMLIBArray(tmpRebates),
                DEF_KNOCKIN,
                0e0,   // recovery, not used for knock-in.
		discZcName.c_str());

	ASSERT_OR_THROW(defKnockIn);

	// Add underlying
	defKnockIn->AddDep(underlying);

	// If name is empty, then use pointer address
	if (name.empty())
		name = "DefKnockIn" + String(format("%p", defKnockIn));

	defKnockIn->SetName(name.c_str());

	return Raw2SharedPointer(defKnockIn);


  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nDEFKNOCKIN_GENERAL: failed.\n");
	return 0;
  }

}








//**********************************************************************
//
// Option and KO constructors
//
//**********************************************************************



//----------------------------------------------------------------------
// Construct a KVPOption object (simple)
//
MAGNET_DESCRIPTION("Construct an option with simple schedule")
MAGNET_PARAMETER("", name, "", "Name of the object.")
MAGNET_PARAMETER("", call_or_put,  "", "(C)all or (P)ut")
MAGNET_PARAMETER("", amer_or_euro, "", "(A)merican or (E)uropean exercises")
MAGNET_PARAMETER("", startDate,    "", "Option exercise start date")
MAGNET_PARAMETER("", maturityDate, "", "Option exercise end date")
MAGNET_PARAMETER("", frequency,    "", "Exercise frequency")
MAGNET_PARAMETER("", strike,       "", "Strike")
MAGNET_PARAMETER("", underlying,   "", "Underlying instrument")
MAGNET_PARAMETER("", discZcName,   "", "The discount zero curve name")
MAGNET_RESULT("An Option Handle")
MAGNET_SEE_ALSO("OPTION_GENERAL", "See the option with arbitrary schedule")

MAGNET_FUNCTION9(
	SharedPointer<KVPOption>, OPTION_SIMPLE,	// Wrapper type and name

	const String&	name,			// (I) 1. Name of the object

	String&		call_or_put,		// (I) 2. call/put type
	String&		amer_or_euro,		// (I) 3. American/European
	Date		startDate,		// (I) 4. start date
	Date		maturityDate,		// (I) 5. end date
  const TDateInterval&	frequency,		// (I) 6. freqency
	double		strike,			// (I) 7. strike
	SharedPointer<KVPAtom>	underlying,	// (I) 8. underlying asset
	String&		discZcName)		// (I) 9. discount curve
{
static	char routine[] = "MAGNET::OPTION_SIMPLE";

	int		nDays = 0;

	KDateInterval	notifDays = KDateInterval(nDays, FALSE);

	TBoolean	A_or_E = FALSE;
	TBoolean	stubAtEnd = FALSE;

	KVPOption	*option = NULL;

  try {

	ASSERT_OR_THROW(amer_or_euro.c_str());
	if (toupper(*(amer_or_euro.begin())) == 'A')
		A_or_E = TRUE;
	else if (toupper(*(amer_or_euro.begin())) == 'E')
		A_or_E = FALSE;
	else
		throw KFailure("%s: invalid option exercise type (%s)."
			       "Only American or European allowed.\n", 
				routine, amer_or_euro.c_str());
		

	// Create option
	option = new KVPOption(
		name.c_str(),
		KVPOptionType(call_or_put.c_str()),
		A_or_E,
		startDate,
		maturityDate,
		frequency,
		stubAtEnd,
		notifDays,
		strike,
		discZcName.c_str());
	ASSERT_OR_THROW(option);

	// Add underlying
	option->AddDep(underlying);

	return Raw2SharedPointer(option);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nOPTION_SIMPLE: failed.\n");
	return 0;
  }

}


//----------------------------------------------------------------------
// Construct a KVPOption object (arbitrary)
// Strike pay dates are assumed to be the same as
// settlement dates.
//
MAGNET_DESCRIPTION("Construct an option with simple schedule")
MAGNET_PARAMETER("", name, "option", "Name of the object.")
MAGNET_PARAMETER("", call_or_put,  "", "(C)all or (P)ut")
MAGNET_PARAMETER("", amer_or_euro, "", "(A)merican or (E)uropean exercises")
MAGNET_PARAMETER("", notifDates,   "", "Option exercise notification dates")
MAGNET_PARAMETER("", settleDates,  "", "Option exercise settlement dates")
MAGNET_PARAMETER("", strikes,      "", "Strikes: absolute values if the notionals below are not specified, or percentage of the notionals if the notionals below are specified")
MAGNET_PARAMETER("", notionals,    "empty", "Notionals (optional).  If the argument is used, then the strikes given above are percentages of the notionals")
MAGNET_PARAMETER("", underlying,   "", "Underlying instrument")
MAGNET_PARAMETER("", discZcName,   "", "The discount zero curve name")
MAGNET_RESULT("An Option Handle")
MAGNET_SEE_ALSO("OPTION_SIMPLE", "See the option with simple schedule")

MAGNET_X_FUNCTION9(
	SharedPointer<KVPOption>, OPTION_GENERAL, // Wrapper type and name

					// (I) 1. Name of the object	
	String		name,		="",

					// (I) 2. call/put type
	String&		call_or_put, 	MAGNET_MANDATORY,
					// (I) 3. American/European
	String&		amer_or_euro,	MAGNET_MANDATORY,
					// (I) 4. array of notification dates
	const Array<Date>& notifDates,	MAGNET_MANDATORY,
					// (I) 5. array of settlement dates
	const Array<Date>& settleDates,	MAGNET_MANDATORY,
					// (I) 6. array of strikes
	const Array<double>& strikes,	MAGNET_MANDATORY,
					// (I) 7. array of strikes
	const Array<double>& notionals,	=Array<double>(),
					// (I) 8. underlying asset
	SharedPointer<KVPAtom>	underlying, MAGNET_MANDATORY,
					// (I) 9. discount curve
	String&		discZcName,	MAGNET_MANDATORY)	
{
static	char routine[] = "MAGNET::OPTION_GENERAL";

	int		nDays = 0;

	TBoolean	A_or_E = FALSE;

	KDateInterval	notifDays = KDateInterval(nDays, FALSE);

	KVPOption	*option = NULL;

  try {
	
	ASSERT_OR_THROW(amer_or_euro.c_str());
	if (toupper(*(amer_or_euro.begin())) == 'A')
		A_or_E = TRUE;
	else if (toupper(*(amer_or_euro.begin())) == 'E')
		A_or_E = FALSE;
	else
		throw KFailure("%s: invalid option exercise type (%s)."
			       "Only American or European allowed.\n", 
				routine, amer_or_euro.c_str());


	// Create option

	//
	// 1. The strikes are absolute values when the notionals 
	//    are NOT specified
	//
	if (notionals.empty())
		option = new KVPOption(
			name.c_str(),
			KVPOptionType(call_or_put.c_str()),
			A_or_E,
			DppDateArrayToTDateVector(notifDates),
			DppDateArrayToTDateVector(settleDates),
			DppVectorFromCMLIBArray(strikes),
			notifDays,
			discZcName.c_str());
	else 
		option = new KVPOption(
			name.c_str(),
			KVPOptionType(call_or_put.c_str()),
			A_or_E,
			DppDateArrayToTDateVector(notifDates),
			DppDateArrayToTDateVector(settleDates),
			DppDateArrayToTDateVector(settleDates),
			DppVectorFromCMLIBArray(strikes),
			DppVectorFromCMLIBArray(notionals),
			notifDays,
			discZcName.c_str());

	ASSERT_OR_THROW(option);

	// Add underlying
	option->AddDep(underlying);

	// If name is empty, then use pointer address
	if (name.empty())
		name = "option_" + String(format("%p", option));

	option->SetName(name.c_str());

	return Raw2SharedPointer(option);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nOPTION_GENERAL: failed.\n");
	return 0;
  }

}




//----------------------------------------------------------------------
// Construct a KVPKnockIO (simple)
//
MAGNET_DESCRIPTION("Construct a knock-in/out object with simple schedule")
MAGNET_PARAMETER("", name, "", "Name of the object.")
MAGNET_PARAMETER("", knock_in_or_out,  "", "(I): Knock-in or (O): Knock-out")
MAGNET_PARAMETER("", in_or_out_barrier,"", "(I): Inside or (O): outside the barriers")
MAGNET_PARAMETER("", rateObs,  	       "", "Triggered rate index")
MAGNET_PARAMETER("", smooth, 	       "", "Nod smoothing method: (N)one, (S)ingle, or (D)ouble")
MAGNET_PARAMETER("", startDate,        "", "Start date of observations")
MAGNET_PARAMETER("", maturityDate,     "", "End date of observations")
MAGNET_PARAMETER("", frequency,        "", "Observational frequency")
MAGNET_PARAMETER("", notifDays,        "", "Number of notification days")
MAGNET_PARAMETER("", barrierLow,       "", "Lower range of barrier")
MAGNET_PARAMETER("", barrierHigh,      "", "Higher range of barrier")
MAGNET_PARAMETER("", rebate,           "", "Rebate value (knock-out only)")
MAGNET_PARAMETER("", underlying,       "", "Underlying instrument")
MAGNET_PARAMETER("", discZcName,       "", "The discount zero curve name")
MAGNET_RESULT("An KnockIO Handle")
MAGNET_SEE_ALSO("KNOCKIO_GENERAL", "See the knock-IO with arbitrary schedule")

MAGNET_FUNCTION13(
	SharedPointer<KVPKnockIO>, KNOCKIO_SIMPLE,	// Wrapper type and name

	const String&	name,			// (I) 1. Name of the object	

	String&		knock_in_or_out,	// (I) 2. Knock-in or Knock-out
	String&		in_or_out_barrier,	// (I) 3. Inside or outside 
						//	  the barriers
	SharedPointer<KRate> rateObs,		// (I) 4. Triggered raet index 
	String&		smooth,			// (I) 5. Node smoothing method 
	Date		startDate,		// (I) 6. start date
	Date		maturityDate,		// (I) 7. end date
  const TDateInterval&	frequency,		// (I) 8. freqency
	double		barrierLow,		// (I) 9. low range of barrier
	double		barrierHigh,		// (I) 10. high range of barrier
	double		rebate,			// (I) 11. knock-out rebate
	SharedPointer<KVPAtom>	underlying,	// (I) 12. underlying asset
	String&		discZcName)		// (I) 13. discount curve
{
static	char routine[] = "MAGNET::KNOCKIO_SIMPLE";

	KKnockIO        kioType, kioWindow;
	KSmooth         kSmooth;
 
	int		nDays = 0;

	KDateInterval	notifDays = KDateInterval(nDays, FALSE);

	TBoolean	stubAtEnd = FALSE;

  try {

	//
	// Generate constant barrier and rebate schedules
	//
	TDate		barrierD[2] = {startDate, maturityDate};
	KVector(TDate)	barrierDates(barrierD, barrierD+2);
 
	KVector(double) barrierLos(2, barrierLow);
	KVector(double) barrierHis(2, barrierHigh);
	KVector(double) rebates(2, rebate);
 
        KVPKnockIO      *knockIO = NULL;


	// Knock in/out type
        //
	switch (toupper(knock_in_or_out[0])) {
	case 'I':	
		kioType = CRX_KNOCK_IN;
		break;
	case 'O':	
		kioType = CRX_KNOCK_OUT;
		break;
	case 'N':	
		kioType = CRX_NONE;
		break;
	default:
		throw KFailure("%s: invalid knock in/out type (%s).\n",
				routine,
				knock_in_or_out.c_str());
	}
 
	// Knock in/out window type
	//
	switch (toupper(in_or_out_barrier[0])) {
	case 'I':	
                kioWindow = CRX_KNOCK_IN;
		break;
	case 'O':	
                kioWindow = CRX_KNOCK_OUT;
		break;
	default:
                throw KFailure("%s: invalid knock in/out type (%s).\n",
                                routine,
				in_or_out_barrier.c_str());
	}
 
        // Smoothing type
        //
	switch (toupper(smooth[0])) {
	case 'D':	
                kSmooth = DOUBLE_SMOOTH;
		break;
	case 'S':	
                kSmooth = SINGLE_SMOOTH;
		break;
	case 'N':	
                kSmooth = NO_SMOOTH;
		break;
	default:
                throw KFailure("%s: invalid smoothing type (%s).\n",
                                routine,
                                smooth.c_str());
	}


	// Create option
	knockIO = new KVPKnockIO(
		name.c_str(),
		kioType,
		kioWindow,
		rateObs,
		kSmooth,
		startDate,
		maturityDate,
		frequency,
		stubAtEnd,
		notifDays,
		barrierDates,
		barrierLos,
		barrierHis,
		rebates,
		discZcName.c_str());

	ASSERT_OR_THROW(knockIO);

	// Add underlying
	knockIO->AddDep(underlying);

	return Raw2SharedPointer(knockIO);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nKNOCKIO_SIMPLE: failed.\n");
	return 0;
  }

}



//----------------------------------------------------------------------
// Construct a KVPKnockIO (general)
//
MAGNET_DESCRIPTION("Construct a knock-in/out object with general schedule")
MAGNET_PARAMETER("", name, "knockIO", "Name of the object.")
MAGNET_PARAMETER("", knock_in_or_out,  "", "(I): Knock-in or (O): Knock-out")
MAGNET_PARAMETER("", in_or_out_barrier,"", "(I): Inside or (O): Outside the barriers")
MAGNET_PARAMETER("", smooth,           "", "Nod smoothing method: (N)one, (S)ingle, or (D)ouble")
MAGNET_PARAMETER("", rateObs,          "", "Triggered rate index. Can be a single rate object or an array of rate objects")
MAGNET_PARAMETER("", observDates,      "", "Array of trigger rate observation dates")
MAGNET_PARAMETER("", observEffDates,   "", "Array of trigger rate effective dates")
MAGNET_PARAMETER("", settleDates,      "", "Array of underlying settlement dates")
MAGNET_PARAMETER("", barrierLows,      "", "Array of low barriers")
MAGNET_PARAMETER("", barrierHighs,     "", "Array of higher barriers")
MAGNET_PARAMETER("", rebates,          "empty", "Array of rebate values (knock-out only)")
MAGNET_PARAMETER("", underlying,       "", "Underlying instrument")
MAGNET_PARAMETER("", discZcName,       "", "The discount zero curve name")
MAGNET_RESULT("An KnockIO Handle")
MAGNET_SEE_ALSO("KNOCKIO_SIMPLE", "See the knock-IO with simple schedule")

MAGNET_X_FUNCTION13(
	SharedPointer<KVPKnockIO>, KNOCKIO_GENERAL,// Wrapper type and name

					// (I) 1. name of the object
	String		name,			="",		

					// (I) 2. Knock-in or Knock-out
	String&		knock_in_or_out,	MAGNET_MANDATORY,
					// (I) 3. Inside or outside 
					//	  the barriers
	String&		in_or_out_barrier,	MAGNET_MANDATORY,	
					// (I) 4. Node smoothing method 
	String&		smooth,			MAGNET_MANDATORY,
					// (I) 5. Triggered raet index 
  const	Array<SharedPointer<KRate> >& rateObs,	MAGNET_MANDATORY,
					// (I) 6. observation dates
  const Array<Date>&	observDates,		MAGNET_MANDATORY,
					// (I) 7. observation dates
  const Array<Date>&	observEffDates,		MAGNET_MANDATORY,	
					// (I) 8. settlement dates
  const Array<Date>&	settleDates,		MAGNET_MANDATORY,	
					// (I) 9. low barriers
  const Array<double>&	barrierLows,		MAGNET_MANDATORY,	
					// (I) 10. high barriers
  const Array<double>&	barrierHighs,		MAGNET_MANDATORY,	
					// (I) 11. knock-out rebate
  const Array<double>&	rebates,		=Array<double>(),
					// (I) 12. underlying asset
	SharedPointer<KVPAtom>	underlying,	MAGNET_MANDATORY,
					// (I) 13. discount curve
	String&		discZcName,		MAGNET_MANDATORY)		
{
static	char routine[] = "MAGNET::KNOCKIO_GENERAL";

	KKnockIO        kioType, kioWindow;
	KSmooth         kSmooth;
 
	int		idx;

	Array<double>	tmpRebates;
	KVector(SharedPointer<KRate>)	vRates;

        KVPKnockIO      *knockIO = NULL;

  try {

	// Knock in/out type
        //
	switch (toupper(knock_in_or_out[0])) {
	case 'I':	
		kioType = CRX_KNOCK_IN;
		break;
	case 'O':	
		kioType = CRX_KNOCK_OUT;
		break;
	case 'N':	
		kioType = CRX_NONE;
		break;
	default:
		throw KFailure("%s: invalid knock in/out type (%s).\n",
				routine,
				knock_in_or_out.c_str());
	}
 
	// Knock in/out window type
	//
	switch (toupper(in_or_out_barrier[0])) {
	case 'I':	
                kioWindow = CRX_KNOCK_IN;
		break;
	case 'O':	
                kioWindow = CRX_KNOCK_OUT;
		break;
	default:
                throw KFailure("%s: invalid knock in/out type (%s).\n",
                                routine,
				in_or_out_barrier.c_str());
	}
 
        // Smoothing type
        //
	switch (toupper(smooth[0])) {
	case 'D':	
                kSmooth = DOUBLE_SMOOTH;
		break;
	case 'S':	
                kSmooth = SINGLE_SMOOTH;
		break;
	case 'N':	
                kSmooth = NO_SMOOTH;
		break;
	default:
                throw KFailure("%s: invalid smoothing type (%s).\n",
                                routine,
                                smooth.c_str());
	}

	
	//
	// Initialize rebates
	//
	if (rebates.empty())
	{
		tmpRebates.resize(observDates.size());
		
		for (idx=0; idx < observDates.size(); idx++)
			tmpRebates[idx] = 0e0;	
	}
	else
		tmpRebates = rebates;

	// Create option
	// 1. One trigger rate
	if (rateObs.size() == 1)
	    knockIO = new KVPKnockIO(
		name.c_str(),
		kioType,
		kioWindow,
		rateObs[0],
		kSmooth,
		DppDateArrayToTDateVector(observDates),
		DppDateArrayToTDateVector(observEffDates),
		DppDateArrayToTDateVector(settleDates),
		DppVectorFromCMLIBArray(barrierLows),
		DppVectorFromCMLIBArray(barrierHighs),
		DppVectorFromCMLIBArray(tmpRebates),
		discZcName.c_str());
	else 
	{
	    for (idx=0; idx < observDates.size(); idx++)
		vRates.insert(vRates.end(), rateObs[idx]);

	    knockIO = new KVPKnockIO(
		name.c_str(),
		kioType,
		kioWindow,
		kSmooth,
		DppDateArrayToTDateVector(observDates),
		DppDateArrayToTDateVector(observEffDates),
		DppDateArrayToTDateVector(settleDates),
		vRates,
		DppVectorFromCMLIBArray(barrierLows),
		DppVectorFromCMLIBArray(barrierHighs),
		DppVectorFromCMLIBArray(tmpRebates),
		discZcName.c_str());
	}

	ASSERT_OR_THROW(knockIO);

	// Add underlying
	knockIO->AddDep(underlying);

	// If name is empty, then use pointer address
	if (name.empty())
		name = "knockIO" + String(format("%p", knockIO));

	knockIO->SetName(name.c_str());

	return Raw2SharedPointer(knockIO);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nKNOCKIO_GENERAL: failed.\n");
	return 0;
  }

}



//----------------------------------------------------------------------
// Construct a KVPKnockIO2Idx (simple)
//
MAGNET_DESCRIPTION("Construct a dual knock-in/out object with simple schedule")
MAGNET_PARAMETER("", name, "", "Name of the object.")
MAGNET_PARAMETER("", knock_in_or_out,  "", "(I): Knock-in or (O): Knock-out")
MAGNET_PARAMETER("", in_or_out_bar1,   "", "(I): Inside or (O): outside the 1st barriers")
MAGNET_PARAMETER("", rateObs1,         "", "1st rriggered rate index")
MAGNET_PARAMETER("", in_or_out_bar2,   "", "(I): Inside or (O): outside the 2nd barriers")
MAGNET_PARAMETER("", rateObs2,         "", "1st rriggered rate index")
MAGNET_PARAMETER("", smooth,           "", "Node smoothing method: (N)one, (S)ingle, or (D)ouble")
MAGNET_PARAMETER("", startDate,        "", "Start date of observations")
MAGNET_PARAMETER("", maturityDate,     "", "End date of observations")
MAGNET_PARAMETER("", frequency,        "", "Observational frequency")
MAGNET_PARAMETER("", notifDays,        "", "Number of notification days")
MAGNET_PARAMETER("", barrierLow1,      "", "Lower range of 1st barrier")
MAGNET_PARAMETER("", barrierHigh1,     "", "Higher range of 1st barrier")
MAGNET_PARAMETER("", barrierLow2,      "", "Lower range of 2nd barrier")
MAGNET_PARAMETER("", barrierHigh2,     "", "Higher range of 2nd barrier")
MAGNET_PARAMETER("", rebate,           "", "Rebate value (knock-out only)")
MAGNET_PARAMETER("", underlying,       "", "Underlying instrument")
MAGNET_PARAMETER("", discZcName,       "", "The discount zero curve name")
MAGNET_RESULT("An Dual KnockIO Handle")
MAGNET_SEE_ALSO("KNOCKIO_DUAL GENERAL", "See the dual knock-IO with arbitrary schedule")

MAGNET_FUNCTION17(
    SharedPointer<KVPKnockIO2Idx>, KNOCKIO_DUAL_SIMPLE, // Wrapper type and name

    const String&   name,           // (I) 1. Name of the object    

    String&     knock_in_or_out,    // (I) 2. Knock-in or Knock-out
    String&     in_or_out_barrier1, // (I) 3. Inside or outside 
    SharedPointer<KRate> rateObs1,  // (I) 4. 1st triggered rate index 
    String&     in_or_out_barrier2, // (I) 5. Inside or outside 
    SharedPointer<KRate> rateObs2,  // (I) 6. 2nd triggered rate index 
    String&     smooth,         // (I) 7. Node smoothing method 
    Date        startDate,      // (I) 8. start date
    Date        maturityDate,       // (I) 9. end date
    const TDateInterval&    frequency,      // (I) 10. freqency
    double      barrierLow1,        // (I) 11. low range of 1st barrier
    double      barrierHigh1,       // (I) 12. high range of 1st barrier
    double      barrierLow2,        // (I) 13. low range of 2nd barrier
    double      barrierHigh2,       // (I) 14. high range of 2nd barrier
    double      rebate,         // (I) 15. knock-out rebate
    SharedPointer<KVPAtom>  underlying, // (I) 16. underlying asset
    String&     discZcName)     // (I) 17. discount curve
{
static  char routine[] = "MAGNET::KNOCKIO_SIMPLE";

    KKnockIO        kioType, kioWindow1, kioWindow2;
    KSmooth         kSmooth;
 
    int     nDays = 0;

    KDateInterval   notifDays = KDateInterval(nDays, FALSE);

    TBoolean    stubAtEnd = FALSE;

  try {

    //
    // Generate constant barrier and rebate schedules
    //
    TDate       barrierD[2] = {startDate, maturityDate};
    KVector(TDate)  barrierDates(barrierD, barrierD+2);
 
    KVector(double) barrierLos1(2, barrierLow1);
    KVector(double) barrierHis1(2, barrierHigh1);
    KVector(double) barrierLos2(2, barrierLow2);
    KVector(double) barrierHis2(2, barrierHigh2);
    KVector(double) rebates(2, rebate);
 
    KVPKnockIO2Idx   *knockIO2Idx = NULL;


    // Knock in/out type
        //
    switch (toupper(knock_in_or_out[0])) {
    case 'I':   
        kioType = CRX_KNOCK_IN;
        break;
    case 'O':   
        kioType = CRX_KNOCK_OUT;
        break;
    case 'N':   
        kioType = CRX_NONE;
        break;
    default:
        throw KFailure("%s: invalid knock in/out type (%s).\n",
                routine,
                knock_in_or_out.c_str());
    }
 
    // Knock in/out window type
    //
    switch (toupper(in_or_out_barrier1[0])) {
    case 'I':   
                kioWindow1 = CRX_KNOCK_IN;
        break;
    case 'O':   
                kioWindow1 = CRX_KNOCK_OUT;
        break;
    default:
                throw KFailure("%s: invalid knock in/out type (%s).\n",
                                routine,
                in_or_out_barrier1.c_str());
    }
    switch (toupper(in_or_out_barrier2[0])) {
    case 'I':   
                kioWindow2 = CRX_KNOCK_IN;
        break;
    case 'O':   
                kioWindow2 = CRX_KNOCK_OUT;
        break;
    default:
                throw KFailure("%s: invalid knock in/out type (%s).\n",
                                routine,
                in_or_out_barrier2.c_str());
    }
 
        // Smoothing type
        //
    switch (toupper(smooth[0])) {
    case 'D':   
                kSmooth = DOUBLE_SMOOTH;
        break;
    case 'S':   
                kSmooth = SINGLE_SMOOTH;
        break;
    case 'N':   
                kSmooth = NO_SMOOTH;
        break;
    default:
                throw KFailure("%s: invalid smoothing type (%s).\n",
                                routine,
                                smooth.c_str());
    }


    // Create option
    knockIO2Idx = new KVPKnockIO2Idx(
        name.c_str(),
        kioType,
        kioWindow1,
        rateObs1,
        kioWindow2,
        rateObs2,
        kSmooth,
        startDate,
        maturityDate,
        frequency,
        stubAtEnd,
        notifDays,
        barrierDates,
        barrierLos1,
        barrierHis1,
        barrierLos2,
        barrierHis2,
        rebates,
        discZcName.c_str());

    ASSERT_OR_THROW(knockIO2Idx);

    // Add underlying
    knockIO2Idx->AddDep(underlying);

    return Raw2SharedPointer(knockIO2Idx);

  }
  catch (KFailure) {
    CM_THROW RuntimeError ("\nKNOCKIO_DUAL_SIMPLE: failed.\n");
    return 0;
  }

}



//----------------------------------------------------------------------
// Construct a KVPKnockIO2Idx (general)
//
MAGNET_DESCRIPTION("Construct a dual knock-in/out object with general schedule")
MAGNET_PARAMETER("", name, "knockIO2Idx", "Name of the object.")
MAGNET_PARAMETER("", knock_in_or_out,  "", "(I): Knock-in or (O): Knock-out")
MAGNET_PARAMETER("", in_or_out_bar1,   "", "(I): Inside or (O): Outside the 1st barriers")
MAGNET_PARAMETER("", rateObs1,         "", "1st triggered rate index. Can be a single rate object or an array of rate objects")
MAGNET_PARAMETER("", in_or_out_bar2,   "", "(I): Inside or (O): Outside the 2nd barriers")
MAGNET_PARAMETER("", rateObs2,         "", "2nd triggered rate index. Can be a single rate object or an array of rate objects")
MAGNET_PARAMETER("", smooth,           "", "Nod smoothing method: (N)one, (S)ingle, or (D)ouble")
MAGNET_PARAMETER("", observDates,      "", "Array of trigger rate observation dates")
MAGNET_PARAMETER("", observEffDates,   "", "Array of trigger rate effective dates")
MAGNET_PARAMETER("", settleDates,      "", "Array of underlying settlement dates")
MAGNET_PARAMETER("", barrierLows1,     "", "Array of 1st low barriers")
MAGNET_PARAMETER("", barrierHighs1,    "", "Array of 1st higher barriers")
MAGNET_PARAMETER("", barrierLows2,     "", "Array of 2nd low barriers")
MAGNET_PARAMETER("", barrierHighs2,    "", "Array of 2nd higher barriers")
MAGNET_PARAMETER("", rebates,          "empty", "Array of rebate values (knock-out only)")
MAGNET_PARAMETER("", underlying,       "", "Underlying instrument")
MAGNET_PARAMETER("", discZcName,       "", "The discount zero curve name")
MAGNET_RESULT("An KnockIO Handle")
MAGNET_SEE_ALSO("KNOCKIO_GENERAL", "See the knock-IO with general schedule")

MAGNET_X_FUNCTION17(
    SharedPointer<KVPKnockIO2Idx>, KNOCKIO_DUAL_GENERAL,// Wrapper type and name

                    // (I) 1. name of the object
    String      name,           ="",        

                    // (I) 2. Knock-in or Knock-out
    String&     knock_in_or_out,    MAGNET_MANDATORY,
                    // (I) 3.-6. Inside or outside 
                    //    the triggered rate indices and barriers
    String&     in_or_out_barrier1, MAGNET_MANDATORY,
    const Array<SharedPointer<KRate> >& rateObs1, MAGNET_MANDATORY,
    String&     in_or_out_barrier2, MAGNET_MANDATORY,
    const Array<SharedPointer<KRate> >& rateObs2, MAGNET_MANDATORY,
                    // (I) 7. Node smoothing method 
    String&     smooth,         MAGNET_MANDATORY,
                    // (I) 8. observation dates
  const Array<Date>&    observDates,        MAGNET_MANDATORY,
                    // (I) 9. observation dates
  const Array<Date>&    observEffDates,     MAGNET_MANDATORY,   
                    // (I) 10. settlement dates
  const Array<Date>&    settleDates,        MAGNET_MANDATORY,   
                    // (I) 11. low barriers 1
  const Array<double>&  barrierLows1,       MAGNET_MANDATORY,   
                    // (I) 12. high barriers 1
  const Array<double>&  barrierHighs1,      MAGNET_MANDATORY,   
                    // (I) 13. low barriers 2
  const Array<double>&  barrierLows2,       MAGNET_MANDATORY,   
                    // (I) 14. high barriers 2
  const Array<double>&  barrierHighs2,      MAGNET_MANDATORY,   
                    // (I) 15. knock-out rebate
  const Array<double>&  rebates,        =Array<double>(),
                    // (I) 16. underlying asset
    SharedPointer<KVPAtom>  underlying, MAGNET_MANDATORY,
                    // (I) 17. discount curve
    String&     discZcName,     MAGNET_MANDATORY)       
{
static  char routine[] = "MAGNET::KNOCKIO_DUAL_GENERAL";

    KKnockIO        kioType, kioWindow1, kioWindow2;
    KSmooth         kSmooth;
 
    int     idx;

    Array<double>   tmpRebates;
    KVector(SharedPointer<KRate>)   vRates1;
    KVector(SharedPointer<KRate>)   vRates2;

    KVPKnockIO2Idx  *knockIO2Idx = NULL;

  try {

    // Knock in/out type
        //
    switch (toupper(knock_in_or_out[0])) {
    case 'I':   
        kioType = CRX_KNOCK_IN;
        break;
    case 'O':   
        kioType = CRX_KNOCK_OUT;
        break;
    case 'N':   
        kioType = CRX_NONE;
        break;
    default:
        throw KFailure("%s: invalid knock in/out type (%s).\n",
                routine,
                knock_in_or_out.c_str());
    }
 
    // Knock in/out window type
    //
    switch (toupper(in_or_out_barrier1[0])) {
    case 'I':   
                kioWindow1 = CRX_KNOCK_IN;
        break;
    case 'O':   
                kioWindow1 = CRX_KNOCK_OUT;
        break;
    default:
                throw KFailure("%s: invalid knock in/out type (%s).\n",
                                routine,
                in_or_out_barrier1.c_str());
    }
    switch (toupper(in_or_out_barrier2[0])) {
    case 'I':   
                kioWindow2 = CRX_KNOCK_IN;
        break;
    case 'O':   
                kioWindow2 = CRX_KNOCK_OUT;
        break;
    default:
                throw KFailure("%s: invalid knock in/out type (%s).\n",
                                routine,
                in_or_out_barrier2.c_str());
    }

    // Smoothing type
    //
    switch (toupper(smooth[0])) {
    case 'D':   
                kSmooth = DOUBLE_SMOOTH;
        break;
    case 'S':   
                kSmooth = SINGLE_SMOOTH;
        break;
    case 'N':   
                kSmooth = NO_SMOOTH;
        break;
    default:
                throw KFailure("%s: invalid smoothing type (%s).\n",
                                routine,
                                smooth.c_str());
    }

    
    //
    // Initialize rebates
    //
    if (rebates.empty())
    {
        tmpRebates.resize(observDates.size());
        
        for (idx=0; idx < observDates.size(); idx++)
            tmpRebates[idx] = 0e0;  
    }
    else
        tmpRebates = rebates;

    // Create option
    ASSERT_OR_THROW (rateObs1.size() == rateObs2.size());
    // 1. One trigger rate
    if (rateObs1.size() == 1)
        knockIO2Idx = new KVPKnockIO2Idx(
        name.c_str(),
        kioType,
        kioWindow1,
        rateObs1[0],
        kioWindow2,
        rateObs2[0],
        kSmooth,
        DppDateArrayToTDateVector(observDates),
        DppDateArrayToTDateVector(observEffDates),
        DppDateArrayToTDateVector(settleDates),
        DppVectorFromCMLIBArray(barrierLows1),
        DppVectorFromCMLIBArray(barrierHighs1),
        DppVectorFromCMLIBArray(barrierLows2),
        DppVectorFromCMLIBArray(barrierHighs2),
        DppVectorFromCMLIBArray(tmpRebates),
        discZcName.c_str());
    else 
    {
        for (idx=0; idx < observDates.size(); idx++)
        vRates1.insert(vRates1.end(), rateObs1[idx]);
        vRates2.insert(vRates2.end(), rateObs2[idx]);

        knockIO2Idx = new KVPKnockIO2Idx(
        name.c_str(),
        kioType,
        kioWindow1,
        vRates1,
        kioWindow2,
        vRates2,
        kSmooth,
        DppDateArrayToTDateVector(observDates),
        DppDateArrayToTDateVector(observEffDates),
        DppDateArrayToTDateVector(settleDates),
        DppVectorFromCMLIBArray(barrierLows1),
        DppVectorFromCMLIBArray(barrierHighs1),
        DppVectorFromCMLIBArray(barrierLows2),
        DppVectorFromCMLIBArray(barrierHighs2),
        DppVectorFromCMLIBArray(tmpRebates),
        discZcName.c_str());
    }

    ASSERT_OR_THROW(knockIO2Idx);

    // Add underlying
    knockIO2Idx->AddDep(underlying);

    // If name is empty, then use pointer address
    if (name.empty())
        name = "knockIO2Idx" + String(format("%p", knockIO2Idx));

    knockIO2Idx->SetName(name.c_str());

    return Raw2SharedPointer(knockIO2Idx);

  }
  catch (KFailure) {
    CM_THROW RuntimeError ("\nKNOCKIO_DUAL_GENERAL: failed.\n");
    return 0;
  }

}



//**********************************************************************
//
// Construct pure IR market environment
//
//**********************************************************************
MAGNET_DESCRIPTION("Construct pure IR market environment")
MAGNET_PARAMETER("", today,        "", "Today's date")
MAGNET_PARAMETER("", zcCurves,     "", "Array of zero curve objects")
MAGNET_PARAMETER("", zcTypes,      "", "Array of zero curve types: 0=diffuse, 1,2 = deterministic spread, 3=basis")
MAGNET_PARAMETER("", zcNames,      "", "Array of user defined zero curve Names")
MAGNET_PARAMETER("", volExpDates,  "", "Array of vol expiration dates")
MAGNET_PARAMETER("", volMatDates,  "", "Array of vol underlying maturity dates")
MAGNET_PARAMETER("", volFreqs,     "", "Array of underlying rate frequencies")
MAGNET_PARAMETER("", volRates,     "", "Array of input volatilities")
MAGNET_PARAMETER("", volTypes,     "", "Input vol type: LOGNORMAL or NORMAL")
MAGNET_PARAMETER("", modParams,    "", "Model parameters: numFact, betas,, alphas, rhos, backbone")
MAGNET_PARAMETER("", tmxFlag,      "", "Y: TMX, N: FIX")
MAGNET_PARAMETER("", smileParams,  "", "q1, q2, fwd shift, Black iterations or MultiQ Smile Parameters")
MAGNET_PARAMETER("", liqBmkDates,  "", "Liquid benchmark dates")
MAGNET_PARAMETER("", treeParams,   "", "Tree parameters: ppy, smooth, numStd")
MAGNET_PARAMETER("", refixRates,   "empty", "Array of refix rate objects")
MAGNET_PARAMETER("", refixDates,   "empty", "Array of refix dates")
MAGNET_PARAMETER("", refixValues,  "empty", "Array of refix values")
MAGNET_RESULT("A Market Environment Handle")

MAGNET_X_FUNCTION17(
    SharedPointer<KMarketEnv>,	MARKET_ENV_IR, // Wrapper type and name

					// KMarketCurves
					      // (I) 1. today's date
	Date			today,	      MAGNET_MANDATORY,
					      // (I) 2. ALIB zero curve objs
 	const Array<const TZeroCurve*>& 	zcCurves,    MAGNET_MANDATORY,
					      // (I) 3. zc types (KV_DIFF, ..)
  	const Array<int>&	zcTypes,     MAGNET_MANDATORY,
					      // (I) 4. array of zc names
  	const Array<String>&	zcNames,     MAGNET_MANDATORY,

					// KVolDiag
					      // (I) 5. option expiration dates
	const Array<Date>&	volExpDates, MAGNET_MANDATORY,
					      // (I) 6. underlying mat dates
	const Array<Date>&	volMatDates, MAGNET_MANDATORY,
					      // (I) 7. underlying rate freqs
	const Array<int>&	volFreqs,    MAGNET_MANDATORY,
					      // (I) 8. volatilities
	const Array<double>&	volRates,    MAGNET_MANDATORY,
					      // (I) 9. percentage vol or bp vol
	String&			volType,     MAGNET_MANDATORY,

					// KMrParam
					      // (I) 10. modela parameters
	const Array<double>&	modParams,   MAGNET_MANDATORY,
                          // (I) 11. TMX or FIX tree
    String&                 tmxFlag,     MAGNET_MANDATORY,
					// KSmileParam
					      // (I) 12. q1, q2, qfsh, iteration or MultiQ Smile
	const Array<double>&	smileParams, MAGNET_MANDATORY,
					      // (I) 13. Liquid Benchmark dates
    const Array<Date>&      liqBmkDates, MAGNET_MANDATORY,
					// Tree Params
					      // (I) 14. ppy, smooth, numStd
	const Array<int>&	treeParams,  MAGNET_MANDATORY,
					// KResetBank
					      // (I) 15. array of float rates
const Array<SharedPointer<KRate> >& refixRates, =Array<SharedPointer<KRate> >(),
					      // (I) 16. array of reset dates
	const Array<Date>&	refixDates,      =Array<Date> (),
					      // (I) 17. array of reset values
	const Array<double>&	refixValues,     =Array<double> () )
{
static	char routine[] = "MAGNET::MARKET_ENV_IR";

	KMarketEnv *mkt= NULL;

  try {

	mkt = new KMarketEnv(
		today,
	  	zcCurves, 
  	  	zcTypes,
  	  	zcNames,   
	  	volExpDates, 
		volMatDates, 
		volFreqs,   
		volRates,  
		volType,  
		modParams,    
        tmxFlag,
		smileParams, 
        liqBmkDates,
		treeParams,
  		refixRates,
		refixDates, 
		refixValues);


	ASSERT_OR_THROW(mkt != NULL);

	return Raw2SharedPointer (mkt);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nMARKET_ENV_IR: failed.\n");
	return 0;
  }
}




//**********************************************************************
//
// Construct Basis market environment
//
//**********************************************************************
MAGNET_DESCRIPTION("Construct basis market environment")
MAGNET_PARAMETER("", today,        "", "Today's date")
MAGNET_PARAMETER("", zcCurves,     "", "Array of zero curve objects")
MAGNET_PARAMETER("", zcTypes,      "", "Array of zero curve types: 0=diffuse, 1,2 = deterministic spread, 3=basis")
MAGNET_PARAMETER("", zcNames,      "", "Array of user defined zero curve Names")
MAGNET_PARAMETER("", irVolExpDates,"", "Array of IR vol expiration dates")
MAGNET_PARAMETER("", irVolMatDates,"", "Array of IR vol underlying maturity dates")
MAGNET_PARAMETER("", irVolFreqs,   "", "Array of IR underlying rate frequencies")
MAGNET_PARAMETER("", irVolRates,   "", "Array of IR input volatilities")
MAGNET_PARAMETER("", irVolTypes,   "", "Input IR vol type: LOGNORMAL or NORMAL")
MAGNET_PARAMETER("", irModParams,  "", "IR model parameters: numFact, betas, alphas, rhos, backbone")
MAGNET_PARAMETER("", tmxFlag,      "", "Y: TMX; N: FIX")
MAGNET_PARAMETER("", irSmileParams,"", "IR smile parameters: q1, q2, fwd shift, Black iterations or MultiQ smile parameters")
MAGNET_PARAMETER("", irLiqBmkDates,"", "Liquid benchmark dates")
MAGNET_PARAMETER("", bsInfo,"", "Extra basis curve info: reference libor curve name, reference discount curve name, DCC for basis rate, DCC for Libor rate, spread type: (S)pread, (P)ercentage.")
MAGNET_PARAMETER("", bsVolExpDates,"", "Array of basis vol expiration dates")
MAGNET_PARAMETER("", bsVolRates,   "", "Array of basis input volatilities")
MAGNET_PARAMETER("", bsVolFreqs,   "", "Array of basis underlying rate frequencies")
MAGNET_PARAMETER("", bsVolType,    "", "Input basis vol type: LOGNORMAL or NORMAL")
MAGNET_PARAMETER("", bsModParams,  "", "Basis model parameters: numFact, betas, alphas, rhos, backbone")
MAGNET_PARAMETER("", bsSmileParams,"", "Basis smile parameters: q1, q2, fwd shift, Black iterations")
MAGNET_PARAMETER("", corrIRBS,     "", "Correlation between IR and basis")
MAGNET_PARAMETER("", treeParams,   "", "Tree parameters: ppy, smooth, numStd")
MAGNET_PARAMETER("", refixRates,   "empty", "Array of refix rate objects")
MAGNET_PARAMETER("", refixDates,   "empty", "Array of refix dates")
MAGNET_PARAMETER("", refixValues,  "empty", "Array of refix values")
MAGNET_RESULT("A Market Environment Handle")

MAGNET_X_FUNCTION25(
	SharedPointer<KMarketEnv>, MARKET_ENV_BASIS,  // Wrapper type and name

					// KMarketCurves
					      // (I) 1. today's date
	Date			today,	      MAGNET_MANDATORY,
					      // (I) 2. ALIB zero curve objs
 	const Array<const TZeroCurve*>& 	zcCurves,    MAGNET_MANDATORY,
					      // (I) 3. zc types (KV_DIFF, ..)
  	const Array<int>&	zcTypes,     MAGNET_MANDATORY,
					      // (I) 4. array of zc names
  	const Array<String>&	zcNames,     MAGNET_MANDATORY,

					// IR KVolDiag
					      // (I) 5. option expiration dates
	const Array<Date>&	irVolExpDates, MAGNET_MANDATORY,
					      // (I) 6. underlying mat dates
	const Array<Date>&	irVolMatDates, MAGNET_MANDATORY,
					      // (I) 7. underlying rate freqs
	const Array<int>&	irVolFreqs,    MAGNET_MANDATORY,
					      // (I) 8. volatilities
	const Array<double>&	irVolRates,    MAGNET_MANDATORY,
					      // (I) 9. percentage vol or bp vol
	String&			irVolType,     MAGNET_MANDATORY,

					// IR KMrParam
					      // (I) 10. model parameters
	const Array<double>&	irModParams,   MAGNET_MANDATORY,
                          // (I) 11. TMX or FIX tree
    String&                 tmxFlag,     MAGNET_MANDATORY,
					// IR KSmileParam
					      // (I) 12. q1, q2, qfsh, iteration or MultiQ Smile
	const Array<double>&	irSmileParams, MAGNET_MANDATORY,
					      // (I) 13. Liquid Benchmark dates
    const Array<Date>&      irLiqBmkDates, MAGNET_MANDATORY,

					// Basis curve
						/* (I) 14. Extra param for basis
                                                 * [1] Ref basis Libor curve
                                                 * [2] Ref basis Disc curve
                                                 * [3] DCC for basis rate
                                                 * [4] DCC for Libor rate
                                                 * [5] Basis type (Spread,
                                                 *      Percentage)  */
	const Array<String>&	bsInfo,		MAGNET_MANDATORY,
					// Basis KVolDiag
					      // (I) 15. option expiration dates
	const Array<Date>&	bsVolExpDates, MAGNET_MANDATORY,
					      // (I) 16. volatilities
	const Array<double>&	bsVolRates,    MAGNET_MANDATORY,
					      // (I) 17. underlying rate freqs
	const Array<int>&	bsVolFreqs,    MAGNET_MANDATORY,
					      // (I) 18. percentage or bp vol
	String&			bsVolType,     MAGNET_MANDATORY,

					// Basis KMrParam
					      // (I) 19. modela parameters
	const Array<double>&	bsModParams,   MAGNET_MANDATORY,
					// KSmileParam
					      // (I) 20. q1, q2, qfsh, iteration
	const Array<double>&	bsSmileParams, MAGNET_MANDATORY,

					// IR-BS correlation
					      // (I) 21. IR-BS correlation
	double			corrIRBS,       MAGNET_MANDATORY,
					// Tree Params
					      // (I) 22. ppy, smooth, numStdev
	const Array<int>&	treeParams,  MAGNET_MANDATORY,
					// KResetBank
					      // (I) 23. array of float rates
const Array<SharedPointer<KRate> >& refixRates, =Array<SharedPointer<KRate> >(),
					      // (I) 24. array of refix dates
	const Array<Date>&	refixDates,      =Array<Date> (),
					      // (I) 25. array of refix values
	const Array<double>&	refixValues,     =Array<double> () )
{
static	char routine[] = "MAGNET::MARKET_ENV_BASIS";

	KMarketEnv *mkt= NULL;

  try {
	mkt = new KMarketEnv(
		today,
	  	zcCurves, 
  	  	zcTypes,
  	  	zcNames,   
	  	irVolExpDates, 
		irVolMatDates, 
		irVolFreqs,   
		irVolRates,  
		irVolType,  
		irModParams,	
        tmxFlag,
		irSmileParams, 
        irLiqBmkDates,
		bsInfo,
	  	bsVolExpDates, 
		bsVolRates,  
		bsVolFreqs,   
		bsVolType,  
		bsModParams,	
		bsSmileParams, 
		corrIRBS,
		treeParams,
  		refixRates,
		refixDates, 
		refixValues);

	ASSERT_OR_THROW(mkt != NULL);

	return Raw2SharedPointer (mkt);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nMARKET_ENV_BASIS: failed.\n");
	return 0;
  }
}


//**********************************************************************
//
// Construct Credit market environment
//
//**********************************************************************
MAGNET_DESCRIPTION("Construct credit market environment")
MAGNET_PARAMETER("", today,        "", "Today's date")
MAGNET_PARAMETER("", zcCurves,     "", "Array of zero curve objects")
MAGNET_PARAMETER("", zcTypes,      "", "Array of zero curve types: 0=diffuse, 1,2 = deterministic spread, 6=credit")
MAGNET_PARAMETER("", zcNames,      "", "Array of user defined zero curve Names")
MAGNET_PARAMETER("", zcInterps,    "", "Array of user defined zero curve interp")
MAGNET_PARAMETER("", irVolExpDates,"", "Array of IR vol expiration dates")
MAGNET_PARAMETER("", irVolMatDates,"", "Array of IR vol underlying maturity dates")
MAGNET_PARAMETER("", irVolFreqs,   "", "Array of IR underlying rate frequencies")
MAGNET_PARAMETER("", irVolRates,   "", "Array of IR input volatilities")
MAGNET_PARAMETER("", irVolTypes,   "", "Input IR vol type: LOGNORMAL or NORMAL")
MAGNET_PARAMETER("", irModParams,  "", "IR model parameters: numFact, betas, alphas, rhos, backbone")
MAGNET_PARAMETER("", tmxFlag,      "", "Y: TMX; N: FIX")
MAGNET_PARAMETER("", irSmileParams,"", "IR smile parameters: q1, q2, fwd shift, Black iterations or MultiQ smile parameters")
MAGNET_PARAMETER("", irLiqBmkDates,"", "Liquid benchmark dates")
MAGNET_PARAMETER("", crIRDiscName, "", "IR discount curve name")
MAGNET_PARAMETER("", recovery, 	   "", "recovery rate")
MAGNET_PARAMETER("", crVolExpDates,"", "Array of credit vol expiration dates")
MAGNET_PARAMETER("", crVolMatDates,"", "Array of credit vol maturity dates")
MAGNET_PARAMETER("", crVolFreqs,   "", "Array of crdit underlying rate frequencies")
MAGNET_PARAMETER("", crVolRates,   "", "Array of credit input volatilities")
MAGNET_PARAMETER("", crVolType,    "", "Input credit vol type: LOGNORMAL or NORMAL")
MAGNET_PARAMETER("", crModParams,  "", "credit model parameters: numFact, betas, alphas, rhos, backbone")
MAGNET_PARAMETER("", crSmileParams,"", "credit smile parameters: q1, q2, fwd shift, Black iterations")
MAGNET_PARAMETER("", corrIRCR,     "", "Correlation between IR and credit")
MAGNET_PARAMETER("", treeParams,   "", "Tree parameters: ppy, smooth, numStd")
MAGNET_PARAMETER("", refixRates,   "empty", "Array of refix rate objects")
MAGNET_PARAMETER("", refixDates,   "empty", "Array of refix dates")
MAGNET_PARAMETER("", refixValues,  "empty", "Array of refix values")
MAGNET_RESULT("A Market Environment Handle")

MAGNET_X_FUNCTION28(
	SharedPointer<KMarketEnv>, MARKET_ENV_CREDIT, // Wrapper type and name

					// KMarketCurves
					      // (I) 1. today's date
	Date			today,	      MAGNET_MANDATORY,
					      // (I) 2. ALIB zero curve objs
 	const Array<const TZeroCurve*>& 	zcCurves,    MAGNET_MANDATORY,
					      // (I) 3. zc types (KV_DIFF, ..)
  	const Array<int>&	zcTypes,     MAGNET_MANDATORY,
					      // (I) 4. array of zc names
  	const Array<String>&	zcNames,     MAGNET_MANDATORY,
					      // (I) 5. array of zc interps
  	const Array<int>&	zcInterps,   MAGNET_MANDATORY,

					// IR KVolDiag
					      // (I) 6. option expiration dates
	const Array<Date>&	irVolExpDates, MAGNET_MANDATORY,
					      // (I) 7. underlying mat dates
	const Array<Date>&	irVolMatDates, MAGNET_MANDATORY,
					      // (I) 8. underlying rate freqs
	const Array<int>&	irVolFreqs,    MAGNET_MANDATORY,
					      // (I) 9. volatilities
	const Array<double>&	irVolRates,    MAGNET_MANDATORY,
					      // (I) 10.percentage vol or bp vol
	String&			irVolType,     MAGNET_MANDATORY,

					// IR KMrParam
					      // (I) 11. model parameters
	const Array<double>&	irModParams,   MAGNET_MANDATORY,
                          // (I) 12. TMX or FIX tree
    String&                 tmxFlag,     MAGNET_MANDATORY,
					// IR KSmileParam
					      // (I) 13. q1, q2, qfsh, iteration or MultiQ Smile
	const Array<double>&	irSmileParams, MAGNET_MANDATORY,
					      // (I) 14. Liquid Benchmark dates
    const Array<Date>&      irLiqBmkDates, MAGNET_MANDATORY,

                    // Credit curve
					       // (I) 15. IR discount curve
    	String                 &crIRDiscName,  MAGNET_MANDATORY,
						// (I) 16. recovery rate
    	double                 recovery,       MAGNET_MANDATORY,
                    // Credit vol
						// (I) 17. option exp dates
    	const Array<Date>      &crVolExpDates, MAGNET_MANDATORY,
						// (I) 18. underlying mat dates
    	const Array<Date>      &crVolMatDates, MAGNET_MANDATORY,
						// (I) 19. underlying rate freqs
    	const Array<int>       &crVolFreqs,    MAGNET_MANDATORY,
						// (I) 20. volatilities
    	const Array<double>    &crVolRates,    MAGNET_MANDATORY,
						// (I) 21. percentage or bp vols
    	String                 &crVolType,     MAGNET_MANDATORY,
                    // KMrParam
						// (I) 22. model parameters
    	const Array<double>    &crModParams,   MAGNET_MANDATORY,
                    // KSmileParam
						// (I) 23. q1, q2, qfsh, iter
    	const Array<double>    &crSmileParams, MAGNET_MANDATORY,

						// (I) 24. IR-CR correlation
    	double                 corrIRCR,       MAGNET_MANDATORY,

					// Tree Params
					      // (I) 25. ppy, smooth, numStdev
	const Array<int>&	treeParams,  MAGNET_MANDATORY,
					// KResetBank
					      // (I) 26. array of float rates
const Array<SharedPointer<KRate> >& refixRates, =Array<SharedPointer<KRate> >(),
					      // (I) 27. array of refix dates
	const Array<Date>&	refixDates,      =Array<Date> (),
					      // (I) 28. array of refix values
	const Array<double>&	refixValues,     =Array<double> () )
{
static	char routine[] = "MAGNET::MARKET_ENV_CREDIT";

	KMarketEnv *mkt= NULL;

  try {
    mkt = new KMarketEnv(
        today,
        zcCurves, 
        zcTypes,
        zcNames,   
        zcInterps,   
        irVolExpDates, 
        irVolMatDates, 
        irVolFreqs,   
        irVolRates,  
        irVolType,  
        irModParams,	
        tmxFlag,
        irSmileParams, 
        irLiqBmkDates,
        crIRDiscName,
        recovery, 
        crVolExpDates, 
        crVolMatDates, 
        crVolFreqs,   
        crVolRates,  
        crVolType,  
        crModParams,	
        crSmileParams, 
        corrIRCR,
        treeParams,
        refixRates,
        refixDates, 
        refixValues);

    ASSERT_OR_THROW(mkt != NULL);

    return Raw2SharedPointer (mkt);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nMARKET_ENV_CREDIT: failed.\n");
	return 0;
  }

}




//**********************************************************************
//
// Pricing KVPAtom
//
//**********************************************************************
MAGNET_DESCRIPTION("Price an instrument")
MAGNET_PARAMETER("", product,      "", "Product to be priced")
MAGNET_PARAMETER("", mktEnv,       "", "Market environment")
MAGNET_PARAMETER("", logLevel,     "0","Logging level")
MAGNET_RESULT("Price of the instrument")

MAGNET_X_FUNCTION3(
	double,		PRICE,		     // Wrapper type and name
	
					// (I) 1. Product to be evaluated
  	SharedPointer<KVPInstr>	  product,	MAGNET_MANDATORY, 
					// (I) 2. Market environments
  	SharedPointer<KMarketEnv> mktEnv,	MAGNET_MANDATORY,
					// (I) 3. Logging level
	int			  logLevel,	=0)    
{
static	char routine[] = "MAGNET::PRICE";


	double			discZeroShift;
	string			discCurveName;	// product discount curve name 

	String	fileName = String("C:/run.log");

	KMap(String, double) 	results;

    DR_ErrorSetCallBack(DR_ErrorCallback);
  try {

	// Enable logging, and define the error messege output file.
	//
	debugLevel = logLevel;

	if (debugLevel > 0)
	{
		DppLogMsgSetFile(fileName.c_str());
		dppLog << *mktEnv << endl;
	}


	

	//
	// Test for NULL pointer
	//
	if( !product )
		throw KFailure("%s: no instrument to be priced.\n", routine);

	discCurveName = product->GetDiscName();

    //
    // Call the main pricing routine
    //

    if (!mktEnv->IsCredit())    // Basis/IR tree
        if (!mktEnv->IsTmx())
        {
            KVPRootPrice(
                results,

                product,
                mktEnv->mMktCurves,
                mktEnv->mIRVolDiag,
                mktEnv->mIRMrParam,
                mktEnv->mIRSmileParam,
                mktEnv->mBSVolDiag,
                mktEnv->mBSMrParam,
                mktEnv->mBSSmileParam,
                mktEnv->mIRBSCorr,
                mktEnv->mVPResetBank,
                logLevel);
        }
        else
        {
            KVPRootPrice_TMX(
                results,

                product,
                mktEnv->mMktCurves,
                mktEnv->mIRVolDiag,
                mktEnv->mIRMrParam,
                mktEnv->mIRSmileParam,
                mktEnv->mBSVolDiag,
                mktEnv->mBSMrParam,
                mktEnv->mBSSmileParam,
                mktEnv->mIRBSCorr,
                mktEnv->mVPResetBank,
                logLevel);
        }

    else        // Credit tree
        KVPRootCreditPrice(
            results,

            product,
            mktEnv->mMktCurves,
            mktEnv->mIRVolDiag,
            mktEnv->mIRMrParam,
            mktEnv->mIRSmileParam,
            mktEnv->mCRVolDiag,
            mktEnv->mCRMrParam,
            mktEnv->mCRSmileParam,
            mktEnv->mIRCRCorr,
            mktEnv->mVPResetBank,
            logLevel);

	//
	// Output price.
	// Forward value at value date.
	//
	discZeroShift = mktEnv->mMktCurves.ZeroShift(discCurveName);
	double  npv = results["PV"] / discZeroShift;
	

	return npv;


  }
  catch (KFailure) {
	//CM_THROW RuntimeError ("\nSY failed. Please see Aladdin error log.\n");
	CM_THROW RuntimeError (DR_ErrorRetrieve());
	return 0;
  }

}

//**********************************************************************
//
// Gte factor vols and dates for KVPAtom
//
//**********************************************************************
MAGNET_DESCRIPTION("Return dates and factor vols for an instrument's pricing tree.")
MAGNET_PARAMETER("", product,      "", "Product to be priced")
MAGNET_PARAMETER("", mktEnv,       "", "Market environment")
MAGNET_PARAMETER("", logLevel,     "0","Logging level")
MAGNET_RESULT("Factor vols and dates")

MAGNET_X_FUNCTION3(
	Matrix<double>,		GET_FACTOR_VOLS,		     // Wrapper type and name
	
					// (I) 1. Product to be evaluated
  	SharedPointer<KVPInstr>	  product,	MAGNET_MANDATORY, 
					// (I) 2. Market environments
  	SharedPointer<KMarketEnv> mktEnv,	MAGNET_MANDATORY,
					// (I) 3. Logging level
	int			  logLevel,	=0)    
{
    static	char routine[] = "MAGNET::GET_FACTOR_VOLS";
    //double			discZeroShift;
	string			discCurveName;	// product discount curve name 
	String	fileName = String("C:/run.log");
    Matrix<double>          resultsXL;
    KVector(TDate)          volDates;
    KVector(KVector(double)) volRates;
    int i=0;

  try {

	// Enable logging, and define the error messege output file.
	//
	debugLevel = logLevel;

	if (debugLevel > 0)
	{
		DppLogMsgSetFile(fileName.c_str());
		dppLog << *mktEnv << endl;
	}


	

	//
	// Test for NULL pointer
	//
	if( !product )
		throw KFailure("%s: no instrument to be priced.\n", routine);

	discCurveName = product->GetDiscName();

    //
    // Call the main pricing routine
    //

    if (!mktEnv->IsCredit())    // Basis/IR tree
        KVPRootFactorVols(
            volDates,
            volRates,
            product,
            mktEnv->mMktCurves,
            mktEnv->mIRVolDiag,
            mktEnv->mIRMrParam,
            mktEnv->mIRSmileParam,
            mktEnv->mBSVolDiag,
            mktEnv->mBSMrParam,
            mktEnv->mBSSmileParam,
            mktEnv->mIRBSCorr,
            mktEnv->mVPResetBank,
            logLevel);
    else        // Credit tree
        KVPRootCreditFactorVols(
            volDates,
            volRates,
            product,
            mktEnv->mMktCurves,
            mktEnv->mIRVolDiag,
            mktEnv->mIRMrParam,
            mktEnv->mIRSmileParam,
            mktEnv->mCRVolDiag,
            mktEnv->mCRMrParam,
            mktEnv->mCRSmileParam,
            mktEnv->mIRCRCorr,
            mktEnv->mVPResetBank,
            logLevel);

	//
	// Output price.
	// Forward value at value date.
	//
	//discZeroShift = mktEnv->mMktCurves.ZeroShift(discCurveName);
	//double  npv = results["PV"] / discZeroShift;
	
    //FILE* el = fopen("C:\\vol.log","a");
    //fprintf(el, "Resizing results to: [%d, 1+%d]\n", volDates.size(), volRates.size());
    //fflush(el);
    resultsXL.resize(volDates.size(),1+volRates.size());
    TDate today = mktEnv->mMktCurves.mToday;
    for (i=0; i<volDates.size(); i++) {
        resultsXL[i][0] = (double)(volDates[i]-today);
        //fprintf(el, "\tvol rates size %d is %d\n", i, volRates[i].size());
        //fflush(el);
        int j;
        for (j=0; j<volRates.size(); j++) {
            //fprintf(el, "Asking for volRates[%d][%d]\n", i, j);
            resultsXL[i][j+1] = volRates[j][i];
        }
    }
    //fflush(el);
    //fclose(el);

	return resultsXL;


  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nSY failed. Please see Aladdin error log.\n");
	return resultsXL;
  }

}



//**********************************************************************
//
// Print KVPAtom to specified file
//
//**********************************************************************
MAGNET_DESCRIPTION("Print the details of an instrument to a file")
MAGNET_PARAMETER("", product,      "", "Product to be printed")
MAGNET_PARAMETER("", fileName,     "C:\run.log", "File name")
MAGNET_RESULT("File name to be printed")

MAGNET_X_FUNCTION2(
	String,		PRINT,			// Wrapper type and name
	
						//(I) 1. Product to be printed
  	SharedPointer<KVPAtom>	product,  	MAGNET_MANDATORY,
						// (I) 2. Output file name
	String	fileName, ="C:/run.log")
{
static	char routine[] = "MAGNET::PRINT";
	ofstream 	oFile;

	
  try {
	 // Open output file
	 //
	 oFile.open(fileName.c_str(), ios_base::out | ios_base::app);
	 if( !oFile )
		throw KFailure("%s: failed to open output file %s.\n",
						routine, fileName);

	 oFile << *product << endl;
	 
	 // Close the file
	 //
	 oFile.close();

	 return fileName;
  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nPRINT: failed.\n");
	return 0;
  }

}


//**********************************************************************
//
// Print KMarketEvn to specified file
//
//**********************************************************************
MAGNET_DESCRIPTION("Print the details of an instrument to a file")
MAGNET_PARAMETER("", mktEnv,      "", "Market environment to be printed")
MAGNET_PARAMETER("", fileName,     "C:\run.log", "File name")
MAGNET_RESULT("File name to be printed")

MAGNET_X_FUNCTION2(
	String,		PRINT_MARKET,		// Wrapper type and name
	
						//(I) 1. Product to be printed
  	SharedPointer<KMarketEnv> mktEnv,  	MAGNET_MANDATORY,
						// (I) 2. Output file name
	String	fileName, ="C:/run.log")
{
static	char routine[] = "MAGNET::PRINT_MARKET";
	ofstream 	oFile;

  try {
	
	 // Open output file
	 //
	 oFile.open(fileName.c_str(), ios_base::out | ios_base::app);
	 if( !oFile )
		throw KFailure("%s: failed to open output file %s.\n",
						routine, fileName);

	 oFile << *mktEnv << endl;
	 
	 // Close the file
	 //
	 oFile.close();

	 return fileName;
  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nPRINT_MARKET: failed.\n");
	return 0;
  }

}



//**********************************************************************
//
// Print Object in Yacction format to specified file
//
//**********************************************************************
MAGNET_DESCRIPTION("Print the tree structured, super yacc formatted deal description to a file")
MAGNET_PARAMETER("", rootProd,   "", "Root product to be printed")
MAGNET_PARAMETER("", mkt,      "", "Market environment")
MAGNET_PARAMETER("", fileName,     "C:/supyac_cr_t.dat", "File name")
MAGNET_RESULT("File name to be printed")

MAGNET_X_FUNCTION3(
	String,		YACC_PRINT,		// Wrapper type and name
	
						//(I) 1. Root product
  	SharedPointer<KVPInstr>	rootProd,  	MAGNET_MANDATORY,
						//(I) 2. Market environment
  	SharedPointer<KMarketEnv> mkt, 		MAGNET_MANDATORY,
						// (I) 3. Output file name
	String	fileName, ="nil")
{
static	char routine[] = "MAGNET::YACC_PRINT";
	ofstream 	oFile;

  try {
	
        // Default wrapper file name depends on trade type: IR or CR
        //
        if (fileName == "nil")
        {
            if (!mkt->IsCredit())    // Basis/IR tree
                fileName = String("C:/supyac_x.dat");
            else  // Credit
                fileName = String("C:/supyac_cr_x.dat");
        }

	// Open output file
	//
	oFile.open(fileName.c_str(), ios_base::out);
	if( !oFile )
	throw KFailure("%s: failed to open output file %s.\n",
			routine, fileName);

	//
	// 1. Start of deal file
	//
	oFile << "# DealFile" << endl;

	//
	// 2. Product descriptions
	//	
	rootProd->SetWriteFlag(true); 
	rootProd->YacctionWriteRecursive(oFile); 
	
	//
	// 3. Reset Bank
	//	
	(mkt->mVPResetBank).YacctionWrite(oFile);

	
	//
	// 4. End of deal file
	//	
	oFile << "DISC = \"" << rootProd->GetDiscName() << "\";" << endl;
	oFile << "EVAL(" << rootProd->GetName() << ");" << endl;
	oFile << "END" << endl;
	
	//
	// 5. Market and model section
	//	
	mkt->YacctionWrite(oFile);
	
	 
	// Close the file
	//
	oFile.close();

	return fileName;

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nYACC_PRINT: failed.\n");
	return 0;
  }

}




//**********************************************************************
//
// Resets
//
//**********************************************************************

//----------------------------------------------------------------------
// Adds all rate resets
//
/*
MAGNET_FUNCTION3(
	SharedPointer<KResetBank>,	RESETBANK,	// Wrapper type and name

  const Array<SharedPointer<KRate> >& rsRates, // (I) 1. array of floating rates
	const Array<Date>&	rsDates,      // (I) 2. array of reset dates
	const Array<double>&	rsValues)     // (I) 3. array of reset values
{
static	char routine[] = "MAGNET::RESETBANK";


	int		idx;
	KResetBank	vpResetBank;	// Reset bank
	KRate		rsRateNS;

  try {
	// Check size consistency
	ASSERT_OR_THROW(rsRates.size() == rsDates.size());
	ASSERT_OR_THROW(rsRates.size() == rsValues.size());

	// Add all reset rates
	for (idx = 0; idx < rsRates.size(); idx++)
	{
		rsRateNS = (*rsRates[idx]);
		rsRateNS.SetSpread(0e0);

		vpResetBank.Insert(rsRateNS, rsDates[idx], rsValues[idx]);
	}


	return Raw2SharedPointer(&vpResetBank);

  }
  catch (KFailure) {
	CM_THROW RuntimeError ("\nRESETBANK: failed.\n");
	return 0;
  }

}






//**********************************************************************
//
// Construct market curves
//
//**********************************************************************
MAGNET_FUNCTION4(
	SharedPointer<KMarketCurves>,	MARKET_ZC,	// Wrapper type and name

	Date			today,	   // (I) 1. today's date
	const Array<const TZeroCurve*> 	&zcCurves, // (I) 2. ALIB zero curve objects
	const Array<int>	&zcTypes,  // (I) 3. zc types (KV_DIFF, etc.)
	const Array<String>	&zcNames)  // (I) 4. array of zc names
{
static	char routine[] = "MAGNET::MARKET_ZC";


	KMarketCurves	mktCurves;	// curves and curve types
	int		idx;

	// Check size consistency
	ASSERT_OR_THROW(zcCurves.size() == zcTypes.size());
	ASSERT_OR_THROW(zcCurves.size() == zcNames.size());

	//
	// Today's date
	//
	mktCurves.mToday = today;

	//
	// Add all the curves
	//
	for (idx = 0; idx < zcCurves.size(); idx++)
	{
		ASSERT_OR_THROW(zcCurves[idx]->curve);	

		mktCurves.Insert(zcCurves[idx]->curve, 
				 zcTypes[idx], 
				 zcCurves[idx]->curve->fBaseDate, 
				 zcNames[idx]);
	}

	return Raw2SharedPointer(&mktCurves);

}
*/


//**********************************************************************
//
// Construct vol data
//
//**********************************************************************
/*
MAGNET_FUNCTION5(
	SharedPointer<KVolDiag>,	MARKET_VOL,	// Wrapper type and name

	const Array<Date> &volExpDates,	// (I) 1. option expiration dates
	const Array<Date> &volMatDates,	// (I) 2. underlying maturity dates
	const Array<int>  &volFreqs,	// (I) 3. underlying rate frequencies
	const Array<double> &volRates,	// (I) 4. array of volatilities
	String		  &volType)	// (I) 5. percentage vol of bp vol
{
static	char routine[] = "MAGNET::MARKET_VOL";

 try {

	KVolDiag	irVolDiag;	// IR volatility data.

	// Check size consistency
	ASSERT_OR_THROW(volExpDates.size() == volMatDates.size());
	ASSERT_OR_THROW(volExpDates.size() == volFreqs.size());
	ASSERT_OR_THROW(volExpDates.size() == volRates.size());

	irVolDiag = KVolDiag(DppDateArrayToTDateVector(volExpDates),
			     DppDateArrayToTDateVector(volMatDates),
			     DppVectorFromCMLIBArray(volFreqs),
			     DppVectorFromCMLIBArray(volRates));

	if (toupper(volType[0]) == 'N')
		irVolDiag.mVolType = NORMVOL;
	else if (toupper(volType[0]) == 'L')
		irVolDiag.mVolType = LOGVOL;
	else
		throw KFailure("%s: invalid vol type (%s).\n", 
				routine, volType);


	return Raw2SharedPointer(&irVolDiag);

 }
 catch (KFailure) {
	DppErrMsg("%s: failed.\n", routine);
	//return;
 }
}
*/



//**********************************************************************
//
// Construct mr data
//
//**********************************************************************
/*
MAGNET_FUNCTION8(
	SharedPointer<KMrParam>,  MODEL_PARAM,	// Wrapper type and name

	int			numFact,	// (I) 1. number of factors
	const Array<double>	&betas,		// (I) 2. mean reversion 
	const Array<double>	&alphas,	// (I) 3. weighting parameters
	const Array<double>	&rhos,		// (I) 4. correlation parameters
	double			backBone, 	// (I) 5. backbone
	int			ppy,		// (I) 6. ppy
	double			smoothFactor,	// (I) 7. tree smoothing factor
	int			numStdevCut)	// (I) 8. num stdev to cut tree
{
static	char routine[] = "MAGNET::MODEL_PARAM";

 try {
	KMrParam	irMrParam;	// IR mr data.

	int		idx;

	//
	// Error checking
	//
	if (numFact > 3 || numFact < 1)
		throw KFailure("%s: invalid number of factors (%d).\n", 
				routine, numFact);

	if (betas.size()  <  numFact ||
	    alphas.size() <  numFact ||		
	    rhos.size()   <  numFact*(numFact-1)/2 ||
	    betas.size () != alphas.size())
		throw KFailure("%s: size of input model parameters "
			       "inconsistent with factor number.\n.", 
				routine);

	irMrParam.mNumFact = numFact;

	for (idx = 0; idx < numFact; idx++)
	{
		irMrParam.mBeta[idx]  = betas[idx];
		irMrParam.mAlpha[idx] = alphas[idx];
	}

	for (idx = 0; idx < numFact*(numFact-1)/2; idx++)
		irMrParam.mRho[idx] = rhos[idx];

	irMrParam.mPpy          = ppy;
	irMrParam.mSmoothFact   = smoothFactor;
	irMrParam.mNumStdevCut  = numStdevCut;

	return Raw2SharedPointer(&irMrParam);

 }
 catch (KFailure) {
	DppErrMsg("%s: failed.\n", routine);
	//return;
 }
}

*/


//**********************************************************************
//
// Construct smile data
//
//**********************************************************************
/*
MAGNET_OVERLOAD(1, "normal")
MAGNET_FUNCTION1(
	SharedPointer<KSmileParam>,	SMILE_PARAM,	// Wrapper type and name

	String		&type)			// (I) 1. Constructor type
{
static	char routine[] = "MAGNET::SMILE_PARAM(normal)";

 try {

	KSmileParam		irSmileParam;	// IR skew data.

	//
	// User smile conventions:  1=Normal, 0=Lognormal
	// Model smile conventions: 0=Normal, 1=Lognormal
	//
	irSmileParam->mQ1      = 0.0;
	irSmileParam->mQ2      = 0.0;
	irSmileParam->mQF      = 0.0;
	irSmileParam->mNumIter = 0;
	
	return Raw2SharedPointer(&irSmileParam);

 }
 catch (KFailure) {
	DppErrMsg("%s: failed.\n", routine);
	//return;
 }
}




MAGNET_OVERLOAD(1, "lognormal")
MAGNET_FUNCTION1(
	SharedPointer<KSmileParam>,	SMILE_PARAM,	// Wrapper type and name

	String		&type)		// (I) 1. Constructor type
{
static	char routine[] = "MAGNET::SMILE_PARAM(lognormal)";

 try {

	KSmileParam		irSmileParam;	// IR skew data.

	//
	// User smile conventions:  1=Normal, 0=Lognormal
	// Model smile conventions: 0=Normal, 1=Lognormal
	//
	irSmileParam->mQ1      = 1.0;
	irSmileParam->mQ2      = 1.0;
	irSmileParam->mQF      = 0.0;
	irSmileParam->mNumIter = 0;
	
	return Raw2SharedPointer(&irSmileParam);

 }
 catch (KFailure) {
	DppErrMsg("%s: failed.\n", routine);
	//return;
 }
}




MAGNET_OVERLOAD(1, "general")
MAGNET_FUNCTION2(
	SharedPointer<KSmileParam>,	SMILE_PARAM,	// Wrapper type and name

	String			&type,	     // (I) 1. Constructor type
	const Array<double>	&smileParams)// (I) 2. q1, q2, qfsh, iteration
{
static	char routine[] = "MAGNET::SMILE_PARAM(general)";

 try {

	KSmileParam		irSmileParam;	// IR skew data.

	//
	// Error checking
	//
	if (smileParams.size() != 4)
		throw KFailure("%s: need 4 smile inputs (%d).\n", 
				routine, smileParams.size());

	//
	// User smile conventions:  1=Normal, 0=Lognormal
	// Model smile conventions: 0=Normal, 1=Lognormal
	//
	irSmileParam->mQ1      = 1.0 - smileParams[0];
	irSmileParam->mQ2      = 1.0 - smileParams[1];
	irSmileParam->mQF      = smileParams[2];
	irSmileParam->mNumIter = (int)smileParams[3];
	
	return Raw2SharedPointer(&irSmileParam);

 }
 catch (KFailure) {
	DppErrMsg("%s: failed.\n", routine);
	//return;
 }
}
*/
