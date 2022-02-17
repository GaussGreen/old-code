/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtlbase.h
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_vtlbase_H
#define	_vtlbase_H

#include "vpbase.h"
#include "kvtree.h"
#include "krstbank.h"	// Past rate resets information (KResetBank)



//--------------------------------------------------------------
/**
 * Base class for an instrument tree tool instrument.
 */

class KVPToolAtom : public KVPNode<KVPToolAtom> {
public:

	/**
	 * Creates recursively a new KVPToolAtom and all corresponding
	 * dependencies given a instrument and
	 * a virtual tree.
	 */
friend	SharedPointer<KVPToolAtom>	
NewToolRecursive(
	  const	SharedPointer<KVPAtom> ins,	// (I) instrument to value
		KVTree &vt);			// (I) virtual tree


	/**
	 * Destructor.
	 */
virtual	~KVPToolAtom() {};
	

	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPToolAtom");}

	/**
         * Type name.
         */
virtual int     IsType(const char *name) const
	{ return (!strcmp(name, TypeName()));}

	/**
	 * Object name.
	 */
	const char*	GetName(); 

	/**
	 * Returns the curve name assocaited to the tool.
	 * This name is used by the KVtree to allocate the proper
	 * geometry for the slice (i.e. dimemsion in a model
	 * that supports different dimensions for different 
	 * assets living in the tree).
	 */
virtual	const String	&GetCurveName() = 0;

	/**
	 * This routine internally initializes the tool
	 * ONCE THE TREE IS BUILD.
	 * It is to be called AFTER the tree calibration, but
	 * BEFORE the rollback loop.
	 */
	void	Initialize();


	/**
	 * Updates the tool at each step in the tree.
	 */
virtual	void	Update() = 0;


	/**
	 * Returns the current value of the stream.
	 */
virtual	KTSlice&	GetValue() = 0;


	/**
	 * Calculates and returns the results at time 0.
	 * Default method is to calculate just the field "PV"
	 * by method GetValue().GetCenter().
	 */
virtual	KMap(String,double)  GetResults();


	/**
	 * Returns the KVPAtom to which the KVPToolAtom is associated.
	 */
virtual	SharedPointer<KVPAtom>	Atom() = 0;


	/**
	 * Find if an VPAtom is associated to a tool or its successors.
	 * The routine uses the Atom() method to retrieve the KVPAtom
	 * associated to a KVPToolAtom.
	 * Returns the corresponding tool or NULL.
	 */
	SharedPointer<KVPToolAtom> FindAtom(SharedPointer<KVPAtom> ins);


	/**
	 * Adds another instrument "instr" to the list
	 * of dependencies of the instrument.
	 */
virtual	KVPToolAtom&	AddDep(SharedPointer<KVPToolAtom> instr)
	{
		InsertSuc(instr);
		return(*this);
	}

	/**
	 * Returns the number of dependencies.
	 */
	int	NumDep()
	{
		return (mSuc.size());
	}

	/**
	 * Returns the number of dependencies.
	 */
	int	NumDep() const
	{
		return (mSuc.size());
	}

	/**
	 * Returns instrument dependency of index idx.
	 */
	SharedPointer<KVPToolAtom> 	Dep(int idx)
	{
		return (mSuc[idx]);
	}

	/**
	 * Returns instrument dependency of index idx.
	 */
	SharedPointer<KVPToolAtom>	Dep(int idx) const
	{
		return (mSuc[idx]);
	}


	/**
	 * Writes to a stream.
	 */
virtual	ostream& Put(ostream& os, int ident = 0) const;

	/**
	 * Writes in yacction format.
	 */
virtual ostream& YacctionWrite( ostream& os, int indent=FALSE)
	{Put(os, indent); return os;}


	/**
	 * Writes to a stream.
	 */
friend	ostream& operator<<(ostream& os, KVPToolAtom& object)
		{ object.Put(os); return (os);}



protected:

	/** Private initializator */
	KVPToolAtom(KVTree &vt)
	{
		mVTree = &vt;
	}


	/**
	 * This routine is to be called at each time step in the rollback
	 * loop. It is to be used INTERNALLY in each derived tool
	 * to check whether the tool needs update or has already
	 * been updated.
	 * It returns FALSE is the tool (and all its dependencies)
	 * have been updated. Otherwise, the routine updates
	 * all dependencies (with Update) and returns TRUE.
	 */
	int	NeedsUpdate();


	/**
	 * This routine is to be called at each time step in the rollback
	 * loop. It is to be used INTERNALLY in each derived tool
	 * to set the internal flag that keeps track of the current
	 * update status of the tool.
	 * 
	 */
	void	UpdateDone();


				/** Virtual tree linked to */
	KVTree		*mVTree;

private:
				/** Name (for debugging) */
	string          mName;

};


ostream& 
PutPV(SharedPointer<KVPToolAtom> atom, ostream& os, int ident);


/**
 * Performs the valuation of a product on the tree and
 * returnd the mtm.
 */

void	KVPAtomKVTreeEvalNpv(
	KVPAtom &vpRoot,	// (I) instrument to value
	KVTree &vt,		// (I) virtual tree for evaluation
	double *npv);		// (O) value





#endif




