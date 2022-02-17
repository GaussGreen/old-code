/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	kcbank.h
 * Function:	
 * Author:	David Liu
 ***************************************************************/
#ifndef	_kcbank_H
#define	_kcbank_H

#include "kstdinc.h"
#include <vector>

#include "kvtree.h"

class KVTree;
class KTSlice;


//--------------------------------------------------------------
/**
 * Base class definition for date indexed claim bank,
 * which consists of an array of slices indexed by evaluation date
 * (last date where the value of the slice is defined).
 * To each evaluation date corresponds and earlier usage date
 * in the tree.
 */

class	KCBank {
public:

	/** Default constructor. */
			KCBank() 
			{
			 	mVTree = NULL; 
			 	mName = K_DEFAULT_NAME;
				mMaxErDate = -1;
				mTpIdx = -1;
				mLocked = false;
			}

	/** Constructor. */
			KCBank(KVTree &vt,
			       const String &name = K_DEFAULT_NAME)
			{
				mVTree = &vt;
				mName = name;
				mMaxErDate = -1;
				mTpIdx = -1;
				mLocked = false;
			}

        /** Copy constructor. */
			KCBank(const KCBank &cb);


	/** Destructor  */
virtual			~KCBank();


	/**
	 * Updates the bank at timepoint "tpIdx", i.e. 
	 * performs dev, discard unused bank.
	 */
virtual void		Update(int tpIdx);


	/** Set the name of the claim bank for debugging. */
	void		SetCBankName(const String &name) {mName = name;}

	/** Get the curve name of claim bank. */
	const String&	GetCBankName()	const {return mName;}

	/**
	 * Insert a pair of evaluation and usage dates in the bank
	 *  (the allocation of slice will be done later after date 
	 * optimization). It performs the following:
	 * 1. Insert evaluation and earliest usage dates in ascending order.
	 * 2. If EvDate already exists, update ErDate with earliest one.
	 * 3. Update the maximum earliest usage date mMaxErDate.
	 */
virtual	void		InsertDates(TDate       EvDate,
				    TDate       ErDate);

	/**
	 * Insert a (evaluation date, slice) pair in the bank.
	 * WARNING: the pointer to the slice passed be NOT be freed
	 * manually (this is done by the destructor).
	 */
virtual void		InsertSlice(TDate EvDate, KTSlice *ts)
			{ mSlices.insert(KMap(TDate, KTSlice*)::value_type(
						EvDate, ts));}

	/** Get claim values on a time slice. */
virtual void		GetSlice(KTSlice &ts, TDate EvDate); 

	/** Current time point. */
virtual int		TPIdxCurrent() { return mVTree->TPIdxCurrent();}

	/** Current date point. */
virtual TDate		TPDateCurrent() { return mVTree->TPDateCurrent();}

	/** Is the bank full. */

virtual bool		IsFull() 
			{ return mEvDates.size() == mSlices.size() ||
				 mEvDates.size() == 0;}

	/** Is the bank empty. */
virtual bool		IsEmpty()
			{ return mSlices.size()  == 0 ||
				 mEvDates.size() == 0;}

	/** Is the bank locked. */
virtual bool		IsLocked() 
			{return mLocked ||
				mEvDates.size() == mSlices.size();}

	/** Assignment operator. */
KCBank& 		operator=(const KCBank &cb);

	/** Prints summary of a bank. */
friend	ostream&	operator<<(ostream& os, KCBank &cb);


friend class KVTree;
friend class KBirTree;
friend class KTMXirTree;

protected:
				/** Bank name (for debugging). */
	String	 mName;
				/** Virtual tree attached to. */
	KVTree   *mVTree;
				/** Current time node. */
	int	 mTpIdx;
				/** Maximum earliest use date. */
	TDate	 mMaxErDate;
				/** If TRUE, disable all activities but add */
	bool	 mLocked;

				/** Evaluation dates. */
	KVector(TDate)        mEvDates;
				/** Earliest usage dates for mEvDates. */
	KMap(TDate, TDate)    mErDates;
				/** Slices corresponding to mEvDates. */
	KMap(TDate, KTSlice*) mSlices;


};

#endif	/* kcbank_H */
