/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	
 * Function:	
 * Author:	David Liu
 ***************************************************************/
#include "kstdinc.h"    /* Standard definitions & error hadling */
#include "kutilios.h"
#include "krate.h"

#define	_kmrntree_SRC
#include "kmrntree.h"

extern "C" {
#include "cgeneral.h"                   /* Has stdlib.h */
#include "bastypes.h"                   /* TDateList */
#include "cerror.h"                     /* GtoErrMsg */
#include "cmemory.h"                    /* MALLOC, FREE */
#include "extract.h"                    /* GtoExtractArray */
#include "macros.h"                     /* MAX */

#include "datelist.h"                   /* TDateList routines */
#include "ldate.h"                      /* GtoDtFwdAny */
#include "convert.h"                    /* GtoFormatDate */
#include "zerodate.h"                   /* TZeroDates */
#include "termtype.h"                   /* TFloatRateArray */

#include "drlio.h"			/* IO */
#include "drlmem.h"			/* memory allocation */

};


void 	CheckValidDates(TDate	EvDate,
			TDate	ErDate);


//--------------------------------------------------------------
// Check date validity
//
void
CheckValidDates(TDate   EvDate, 
		TDate   ErDate)
{
static	char	routine[] = "CheckValidDates";

	// Check the validity of curve name
	// to be done

        if ((EvDate <= 1L) || (ErDate >= 500000L))
                throw KFailure("%s: invalid evaluation date %ld.\n", 
			      routine, EvDate);

        if ((EvDate <= 1L) || (ErDate >= 500000L))
                throw KFailure("%s: invalid early usage date %ld.\n", 
			      routine, ErDate);

	if (EvDate < ErDate)
                throw KFailure("%s: EvDate (%s) < ErDate (%s).\n",
			      routine, GtoFormatDate(EvDate), 
			      GtoFormatDate(ErDate));


}



//--------------------------------------------------------------
// Copy constructor
// NOT A TRUE COPY !!!.  SLICES ARE SHALLOW COPIED
// 
KCBank::KCBank(const KCBank &cb)
{
	// copy every member
	mName      = cb.mName;
	mVTree     = cb.mVTree;
	mTpIdx     = cb.mTpIdx;
        mMaxErDate = cb.mMaxErDate;
	mLocked    = cb.mLocked;
	
	mEvDates   = cb.mEvDates;
	mErDates   = cb.mErDates;
	mSlices    = cb.mSlices;
	
}


//--------------------------------------------------------------
//
KCBank::~KCBank()
{

	for(KMap(TDate,KTSlice*)::iterator iterSlice=mSlices.begin();
		iterSlice != mSlices.end(); ++iterSlice) {
                delete (*iterSlice).second;
        }

	mSlices.clear();

	return;
}



//--------------------------------------------------------------
// Assignment operator
KCBank&
KCBank::operator=(const KCBank &cb)
{
	// copy every member
	mName      = cb.mName;
	mVTree     = cb.mVTree;
	mTpIdx     = cb.mTpIdx;
        mMaxErDate = cb.mMaxErDate;
	mLocked    = cb.mLocked;
	
	mEvDates   = cb.mEvDates;
	mErDates   = cb.mErDates;
	mSlices    = cb.mSlices;
	
	return (*this);
}



//--------------------------------------------------------------
// InsertDates does the following:
// 1. Insert evaluation and earliest usage dates in ascending order.
// 2. If EvDate already exists, update ErDate with earliest one.
// 3. Update the maximum earliest usage date mMaxErDate.
//
void
KCBank::InsertDates(TDate 	EvDate,		// (I) Evaluation date
		    TDate	ErDate)		// (I) Usage date
{
static	char	routine[] = "KCBank::InsertDate";

	// Check EvDate >= ErDate
	CheckValidDates(EvDate, ErDate);

	// Find the lower bound and insert in ascending order
	if(mEvDates.empty() && mErDates.empty())  // first element
	{
		mEvDates.insert(mEvDates.begin(), EvDate);
		mErDates.insert(TDate_TDate(EvDate, ErDate));
	} else {	// non-empty. find a proper position to insert
		KVector(TDate)::iterator iterEvDate =
			lower_bound(mEvDates.begin(), mEvDates.end(), EvDate);
 
		if (iterEvDate == mEvDates.end()){	//longer than existings
			mEvDates.insert(mEvDates.end(), EvDate);
			mErDates.insert(TDate_TDate(EvDate, ErDate));
		}	
		else if (*iterEvDate == EvDate) {	// exists
		    KMap(TDate, TDate)::iterator itDate=mErDates.find(EvDate);
		    if (itDate != mErDates.end())
			(*itDate).second = MIN(ErDate, (*itDate).second);
		    else
	    		throw KFailure("%s: no ErDate corresponding "
				       "to EvDate (%s) in the bank.\n", 
					routine,
					GtoFormatDate(EvDate));
		}
		else {		// non-exist, insert in ascending order
			mEvDates.insert(iterEvDate, EvDate);
			mErDates.insert(TDate_TDate(EvDate, ErDate));
		}
	}

	// Update mMaxErDate	
	mMaxErDate = MAX(mMaxErDate, ErDate);	

}




//--------------------------------------------------------------
// Get a slice on EvDate.
//
void
KCBank::GetSlice(KTSlice &ts,  TDate EvDate)
{
static	char	routine[] = "KCBank::GetSlice";

	KMap(TDate, KTSlice*)::iterator itSlice = mSlices.find(EvDate);
 
        if (itSlice == mSlices.end())
		throw KFailure("%s: no slice available on EvDate (%s).\n",
				routine,
				GtoFormatDate(EvDate));

	// Copy the slice
	ts = *((*itSlice).second); 

}



//--------------------------------------------------------------
//
void
KCBank::Update(int tpIdx)
{
static	char	routine[] = "KCBank::Update";

	TDate	tmpEvDate;
	TDate	tmpErDate;
	TDate	currDate;

	String	discCurveName;
	
	KVector(TDate) 	obsoleteSlices;

 try{
	
	currDate = (mVTree->TPDate)(tpIdx);

	if (mVTree == NULL)
	    throw KFailure("%s: non-initialized mVTree associated "
			   "with the bank .\n", routine);
		
	if (IsEmpty())
	    return; 	// does nothing


	// compact the bank if necessary 
	if (currDate < mMaxErDate)
	{
	    mMaxErDate = -1;  // reset max EarliestUseDate

	    // compact the bank and reduce numActiveSlices 
	    // Loop over the existing slices
	    for (
	        KMap(TDate, KTSlice*)::iterator iterSlice = mSlices.begin();
		iterSlice != mSlices.end(); ++iterSlice) 
	    {
		tmpEvDate = (*iterSlice).first;
		tmpErDate = mErDates[tmpEvDate];

		if (currDate < tmpErDate)
		{
		    // Store the keys of obsolate slices in a 
		    // temporary array, and remove them later.
		    //
		    // Avoid delete the key here to mess up
		    // with the map.
		    obsoleteSlices.push_back(tmpEvDate);

		    mErDates.erase(tmpEvDate);

		    KVector(TDate)::iterator itEvDate = 
			    find(mEvDates.begin(), mEvDates.end(), tmpEvDate);
		    mEvDates.erase(itEvDate);
		} 
		else
		    // iterSlice is an active slice
		    // update the mMaxErDate 
		    mMaxErDate = MAX(mMaxErDate, tmpErDate);

	    } 

	    // Free the memory of the obsolate slices 
	    // and remove them from the bank now
	    for (KVector(TDate)::iterator iterDate = obsoleteSlices.begin();
		 iterDate != obsoleteSlices.end(); ++iterDate)
	    {
		delete mSlices[*iterDate];;
		mSlices.erase(*iterDate);
	    }

	    obsoleteSlices.clear();

	} // (currDate < mMaxErDate) 


	// now dev the active slices 
	for (
	        KMap(TDate, KTSlice*)::iterator iterSlice  = mSlices.begin();
		iterSlice != mSlices.end(); 
		++iterSlice)
	{
	    if (debugLevel > DEBUG_LEVEL_CBANK)
	    	dppLog << "================ " 
		       << ((*iterSlice).second)->GetSliceName() 
		       << " ================" << endl;

	    discCurveName = ((*iterSlice).second)->GetCurveName();

	    ((*iterSlice).second)->Dev(discCurveName);
	
	}

	// Update the time index
	mTpIdx = tpIdx;

    }
    catch(KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }


}


//--------------------------------------------------------------
//

ostream& 
operator <<(ostream& os, KCBank& cb)
{
	os << "BANK: " << cb.GetCBankName() << endl;
	
	os << "TP = " 
	   << format("%s", GtoFormatDate(cb.TPDateCurrent())) << endl;
	for(KMap(TDate, KTSlice*)::iterator iter = cb.mSlices.begin(); 
			iter != cb.mSlices.end(); ++iter)
	{
		os << format("(Ev=%s; Er=%s)\t%s", 
			     GtoFormatDate((*iter).first), 
		   	     GtoFormatDate(cb.mErDates[(*iter).first]),
			     ((*iter).second)->GetSliceName().c_str()) 
		   << endl;

		/*
		os << *((*iter).second) 
		   << endl;
		*/
	} 
	return(os);
}









/**************************************************************

//--------------------------------------------------------------
//
KTSlice&
KCBank::PopSlice()
{
static	char	routine[] = "KCBank::PopSlice";

	KTSlice tmpSlice;

 try{
	if (IsLocked() ||
	    IsFull() )
	    throw KFailure("%s: bank is out of slices or locked.\n", routine);
	
	mLocked = true;
	tmpSlice = mSlices[mNTotSlice-1];

	// set the slice pointer to NULL to maintain safe FREE 
	mSlices[mNTotSlice-1] = NULL;

	return tmpSlice;

    }
    catch(KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }

}




//--------------------------------------------------------------
//
void
KCBank::PushSlice(KTSlice &slice,
		  TDate    EvDate,
		  TDate    ErDate)
{
static	char	routine[] = "KCBank::PushSlice";

	int     i;

 try{
	if (!mLocked ||
	    (mSlices[mNTotSlice-1]).mData != NULL) 
	    throw KFailure("%s: invalid bank condition.\n", routine);

	if (ErDate > EvDate) 
 	    throw KFailure("%s: ErDate > EvDate.\n", routine);	

	// make sure the EvDate is in strictly decreasing order 
	if (mNActiveSlices > 0)
	{
	    if (mEvDates[0]<=EvDate)
		throw KFailure("%s: New slice is evaluated later "
			       " than elements in the bank.\n", routine);
    	}

	// shift the lowest obsolete slice to the top 
	if (mNActiveSlices < mNTotSlice-1)
	{
	    mSlices[mNTotSlice-1] = mSlices[mNActiveSlices];
	}

	// nudge up the elements one by one 
	for (i=mNActiveSlices-1;i>=0;i--)
	{
	   mSlices[i+1]  = mSlices[i];
	   mEvDates[i+1] = mEvDates[i];
	   mErDates[i+1] = mErDates[i];
	}

	// add the element to 0 offset 
	mSlices[0]  = slice;
	mEvDates[0] = EvDate;
	mErDates[0] = ErDate;
	(mNActiveSlices)++;

	// update the internal parameters 
	mMaxErDate = MAX(mMaxErDate, ErDate);

	mLocked = FALSE;

    }
    catch(KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }

}

******************************************************************/




