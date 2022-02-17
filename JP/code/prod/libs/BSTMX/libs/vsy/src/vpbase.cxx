/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "vpbase.h"
#include "kutilios.h"

//===============================================================
// VPTree
//===============================================================


//---------------------------------------------------------------
/*
 * Recursive destructor: it goes down the dependency tree
 * and erase all nodes by recursive calls to delete.
 */

KVPInstr::~KVPInstr()
{

/*
	//
	// Recursively delete dependencies
	//

	// for all successors
	KVector(KVPInstr*)::iterator	p,q;

	while ((p = mSuc.begin()) != mSuc.end()) {
		// For each successor, we remove it from all its predecessors
		// to avoid freeing twice as we recurse down a possibly
		// recomnmbining tree.
		KVPInstr	*sucNode = *p;
		for (q  = sucNode->mPre.begin();
		     q != sucNode->mPre.end();
		     q++) {
			KVPInstr *pSucNode = *q;
			if (pSucNode == (KVPInstr*)this) {
			} else {
			    pSucNode->RemoveSuc(sucNode);
			}
		}

		// Now delete successor
		delete sucNode;

		// Erase from list of successors
		mSuc.erase(p);
	}
	//dppLog << GetName() << ": KVPInstr deleted." << endl;

*/

}



//---------------------------------------------------------------


void
KVPInstr::SetDiscName(const char *zcName)
{
	if (zcName == NULL) {
		mDiscZcName = String(format("%p", this));
	} else {
		mDiscZcName = String(zcName);
	}
}


