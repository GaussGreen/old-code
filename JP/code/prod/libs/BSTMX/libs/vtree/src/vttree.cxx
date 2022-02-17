/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	
 * Function:	Virtual tree and time slice functions
 * Author:	E. Ben-Artzi - C. Daher (from D. Gallager)
 * Revision:	$Header$
 ***************************************************************/
#include "kvtree.h"	/* Prototype Consistency */
#include "kutilios.h"


//--------------------------------------------------------------
//

void
KVTree::InsertDateList(int numDates, TDate *dates)
{
	for (int idx=0; idx<numDates; idx++) 
		Insert(dates[idx]);
}



//--------------------------------------------------------------
//

void
KVTree::MapCurveName(const String &cName, int idx)
{
	static char routine[] = "KVTree::MapCurveName";

	KMap(String,int)::iterator p = mCurveIdxTable.find(cName);

	if(p==mCurveIdxTable.end()) {
		mCurveIdxTable[cName] = idx;
	}
	else if ((*p).second != idx) {
	     	throw KFailure("%s: Attempted to map curve %s "
			       "into index %d.\n" 
			       "Pair <%s, %d> already exists.\n",
				routine, cName.c_str(), idx,
				cName.c_str(), idx);
	}
}


//--------------------------------------------------------------
//

int
KVTree::GetCurveIdx(const String &cName)
{
	static char routine[] = "KVTree::GetCurveIdx";

	KMap(String,int)::iterator p = mCurveIdxTable.find(cName);

	if(p != mCurveIdxTable.end()) 
		return ((*p).second);
	else {
		throw KFailure("%s: Curve name %s is not found "
			       "in the tree.\n", 
			       routine, cName.c_str());
	}
	return (K_DEFAULT_IDX);
}


//--------------------------------------------------------------
//

const String&
KVTree::GetCurveName(int idx)
{
	static char routine[] = "KVTree::GetCurveName";

	KMap(String,int)::iterator p = mCurveIdxTable.begin();

	for( ; p != mCurveIdxTable.end(); ++p) {
		if((*p).second == idx)
			return ((*p).first);
	}
	throw KFailure("%s: curve index %d is not set in the tree.\n", 
			routine, idx);
	
	return (K_DEFAULT_NAME);
}



//--------------------------------------------------------------
//

KTSlice::KTSlice(
	KVTree &vt,
	const String &tsName,
	const String &curveName)
{
	static char routine[] = "KTSlice::KTSlice";

	try {
		mVTree = &vt;
		SetSliceName(tsName);
		SetCurveIdx(curveName);
		SetTpIdx(-1);
		SetSliceDim(-1);
		mData = NULL;
		mVTree->TSliceCreate(*this);
	}
	catch (KFailure) {
		throw KFailure("%s: failed to construct time slice %s "
			       "with curve %s.\n",
			       routine, tsName.c_str(), curveName.c_str());
	}

}



//--------------------------------------------------------------
// Copy constructor

KTSlice::KTSlice(const KTSlice &ts)
{
	// copy every member
	mName     = ts.mName;
	mCurveIdx = ts.mCurveIdx;
	mVTree    = ts.mVTree;
	mTpIdx    = ts.mTpIdx;
	mSliceDim = ts.mSliceDim;
	
	this->mData = NULL;

	*this = ts;	// Using "=" operator to allocate slice
			// and copy mData. 
}



//--------------------------------------------------------------
//

void
KTSlice::CheckNonEmpty(const String &debugName) const
{
	static char routine[] = "KTSlice::CheckNonEmpty";

	if (mData == NULL) {
		if (! debugName.empty()) 
			throw KFailure("%s: Time slice %s is empty.\n",
					routine, GetSliceName().c_str());
	}
}


//--------------------------------------------------------------
//
//

KTSlice*
KTSliceNewVector(
	int size,
	KVTree &vt,
	const String &curveName,
	const String &tsName)
{
	static char routine[] = "KTSlice::KTSliceNewVector";

	KTSlice	*p = new KTSlice [size];
	ASSERT_OR_THROW(p != NULL);

	try {
		for (int idx = 0; idx < size; idx++) {
			KTSlice *ts = &p[idx];
			ts->mVTree = &vt;
			ts->SetSliceName(tsName);
			ts->SetCurveIdx(curveName);
			ts->mVTree->TSliceCreate(*ts);
		}
	}
	catch (KFailure) {
		throw KFailure("%s: failed to construct vector of slices %s curve %s.\n",
		routine, tsName.c_str(), curveName.c_str());
	}

	return(p);
}



//--------------------------------------------------------------
//

KTSliceNode::KTSliceNode()
{
	mIndex = NULL;
	mSlice = NULL;
	
}


//--------------------------------------------------------------
//

KTSliceNode::~KTSliceNode()
{
	delete [] mIndex;
}

