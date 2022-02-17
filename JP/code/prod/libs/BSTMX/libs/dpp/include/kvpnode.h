/****************************************************************
 * Module:	
 * Submodule:	
 * File:	kvpnode.h
 * Function:	
 * Author:	Christian Daher, David Liu
 *****************************************************************/
#ifndef	_kvpnode_H
#define	_kvpnode_H

#include <vector>

#include "kstdinc.h"
#include "ktypes.h"


//--------------------------------------------------------------
/**
 * Template class for a node of a oriented tree structure.
 * 
 * It is derived from CMLib's Object to utilize the smart
 * pointer implementation, which allows memory management using
 * the reference-counting scheme.
 * 
 * It used as base class facility to build tree structures
 * where each node is a class and contains a vector
 * of pointers to its predecessors and successors.
 * In order to define a class X to also contain
 * vectors to other (predecessors and successors) X classes,
 * one MUST declare X by <BR>
 * <TT>class X : public KVPNode<X> { ... }; </TT>
 */


template <class T> class KVPNode : public Object {
public:

	
	/**
	 * Returns root xnode.
	 */
	KVPNode&	Root()
	{
		return (mPre.size() != 0 ? mPre[0]->Root() : *this);
	}


	/**
	 * Link two nodes.
	 * Sets "preNode" (resp. "sucNode") to be a precessor
	 * of "sucNode" (resp. successor of "preNode").
	 */
friend	T&	Link(T& preNode, T& sucNode)
	{
		preNode.InsertSuc(&sucNode);
		sucNode.InsertPre(&preNode);
		return(preNode);
	}

	/**
	 * Unlink two nodes.
	 */
friend	T&	Unlink(T& preNode, T& sucNode)
	{
		preNode.RemoveSuc(&sucNode);
		sucNode.RemovePre(&preNode);
		return(preNode);
	}


#ifdef	_SKIP
	/*
	 * Find a given node "xnode" in the successors.
	 */
	KVPNode* Find(KVPNode *xnode)
	{
		int idx;
		if (this == xnode) {
			return &this;
		}

		for (idx=0; idx<mSuc.size(); idx++) {
		}
		Suc[idx]->LoopReset();
	}

	/*
	 * Find a node "xnode" starting from the root.
	 */
	int	FindFronRoot(KVPNode *xnode)
			{ return Root().Find(xnode);}
#endif


	/**
	 * Resets recursively down in tree.
	 */
	void	SetFlag(long value)
	{
		int	idx;
		mUpdateFlag = value;
		for (idx=0; idx<mSuc.size(); idx++)
			mSuc[idx]->SetFlag(value);
	}


	/**
	 * Gets the value of the flag.
	 */
	long	GetFlag()
	{
		return mUpdateFlag;
	}

#ifdef	_SKIP
	/*
	 * Resets a loop to recurse in the tree.
	 */
	void	LoopReset()
	{
		int	idx;
		mUpdateFlag = FALSE;
		for (idx=0; idx<mSuc.size(); idx++)
			Suc[idx]->LoopReset();
	}

	/*
	 * Advances to next xnode in a loop.
	 */
	KVPNode*	LoopNext()
	{
		int	idx;
		KVPNode	*next;

		if (mSuc.size() > 0) {
			for (idx=0; idx<=mSuc-1; idx++) {
				next = mSuc[idx]->LoopNext();
				if (next != NULL) {
					return (next);
				}
			}
		}
		if (resetFlag == FALSE) {
			mUpdateFlag = TRUE;
			return (this);
		} else {
			return(NULL);
		}
	}



#endif



protected:

	/**
	 * Inserts a successor.
	 */
	KVPNode&	InsertSuc(SharedPointer<T> xnode)
	{
		int	idx;
		for (idx=0; idx<mSuc.size(); idx++)
			if (mSuc[idx] == xnode) return(*this);
		mSuc.insert(mSuc.end(), xnode);
		return(*this);
	}

	/**
	 * Removes a successor.
	 */
	KVPNode&	RemoveSuc(SharedPointer<T> xnode)
	{
		for (KVector(SharedPointer<T>)::iterator p = mSuc.begin();
		     p != mSuc.end();
		     p++) {
			if (*p == xnode) {
				mSuc.erase(p);
				break;
			}
		}
		return(*this);
	}



	/**
	 * Inserts a predecessor.
	 */
	KVPNode&        InsertPre(SharedPointer<T> xnode)
	{
	int     idx;
		for (idx=0; idx<mPre.size(); idx++)
			if (mPre[idx] == xnode) return(*this);
		mPre.insert(mPre.end(), xnode);
		return(*this);
	}
 
	/**
	 * Removes a predecessor "x" (or nothing if "x" is not a predecessor).
	 */
	KVPNode&        RemovePre(SharedPointer<T> xnode)
	{
		for (KVector(T*)::iterator p = mPre.begin();
			p != mPre.end();
			p++) {
			if (*p == xnode) {
				mPre.erase(p);
				break;
			}
		}
		return(*this);
	}



				/** Array of ptr to predecessors. */
	KVector(SharedPointer<T>) mPre;
				/** Array of ptr to successors. */
	KVector(SharedPointer<T>) mSuc;

				/** flag used to recurse in the tree. */
	long			  mUpdateFlag;

};




#endif	/*_kvpnode_H*/




