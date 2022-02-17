/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	kvpatom.h
 * Function:	
 * Author:	David Liu
 ***************************************************************/
#ifndef	_kvpatom_H
#define	_kvpatom_H
#include "kstdinc.h"
#include "kvpnode.h"


//--------------------------------------------------------------
/**
 * Class for product component object
 * @version 2.0
 */

class KVPAtom : public KVPNode<KVPAtom> {
public:

	/**
	 * Default constructor.
	 */
	KVPAtom() : KVPNode<KVPAtom> () {mWriteFlag = false;}

	/**
	 * Default constructor.
	 */
	KVPAtom(const char *name) : KVPNode<KVPAtom> ()
	{
		SetName(name);
		mWriteFlag = false;
	};
 
	/**
	 * Copy constructor
	 */
	KVPAtom(const KVPAtom& atom)
	{
		mName = atom.GetName();
		mWriteFlag = atom.GetWriteFlag();
	}

	/**
	 * Sets the name of the instrument.
	 */
	void    SetName(const char *name);
 
	/**
	 * Gets the name of the instrument.
	 */
virtual const char*     GetName() const
	{
		return mName.c_str();
	}

	/**
	 * Type name.
	 */
virtual const char*	TypeName() const {return("KVPAtom");};
 

	/**
	 * Type name.
	 */
virtual int	IsType(const char *) const;


	/**
	 * Adds another instrument "instr" to the list
	 * of dependencies of the instrument.
	 * Reference the dependent only, not the precedent
	 * to avoid cyclic references in shared pointer.
	 */
virtual	KVPAtom&	AddDep(SharedPointer<KVPAtom> instr)
	{
		InsertSuc(instr);
		return(*this);
	}
 
	/**
	 * Returns the number of dependencies.
	 */
	int     NumDep()
	{
		return (mSuc.size());
	}
 
	/**
	 * Returns the number of dependencies.
	 */
	int     NumDep() const
	{
		return (mSuc.size());
	}
 
	/**
	 * Returns instrument dependency of index idx.
	 */
	SharedPointer<KVPAtom>  Dep(int idx)
	{
		return (mSuc[idx]);
	}
 
	/**
	 * Returns instrument dependency of index idx.
	 */
	SharedPointer<KVPAtom>  Dep(int idx) const
	{
		return (mSuc[idx]);
	}


	
	/**
	 * Returns the number of precedents.
	 */
	int     NumPre()
	{
		return (mPre.size());
	}
 
	/**
	 * Returns the number of precedents.
	 */
	int     NumPre() const
	{
		return (mPre.size());
	}
 
	/**
	 * Returns instrument precedent of index idx.
	 */
	SharedPointer<KVPAtom>  Pre(int idx)
	{
		return (mPre[idx]);
	}
 
	/**
	 * Returns instrument precedent of index idx.
	 */
	SharedPointer<KVPAtom>  Pre(int idx) const
	{
		return (mPre[idx]);
	}



	/**
	 * Writes to a stream.
	 */
virtual ostream& Put(ostream& os, int indent=FALSE) const 
		{return os;} 


	/**
	 * Writes details of node and all its dependencies to a stream.
	 */
	ostream& PutRecursive(ostream& os) const;
	

	/**
	 * Writes the dependency tree (in a "readable" format).
	 */
	ostream&	PutTree(ostream& os, int ident = 0);

	/**
	 * Write to stream in Yacction format
	 */
virtual ostream& YacctionWrite(ostream& os, int indent=FALSE)
		{Put(os, indent); return os;}

	/**
	 * Write recursively from bottom-up in yacction format
	 */
virtual	ostream& YacctionWriteRecursive(ostream& os);

	
	/**
	 * Set the yacction write status
	 * true: writable; false: no write
	 */
virtual void	SetWriteFlag(bool flag);

	/**
	 * Set the flag to false to prevent duplicating writing
	 */
void	WriteDone()
	{
		mWriteFlag = false;
	}

	/**
	 * Get the write status
	 */
bool	GetWriteFlag() const
	{
		return mWriteFlag;
	}



friend	ostream& operator<<(ostream& os, const KVPAtom& kObj)
	{ kObj.Put(os); return (os);}


	
protected:
			/** For debugging and yacction */
	String	mName;
			/** Yacction write status. true: writable;
			 *  false: no write */
	bool	mWriteFlag;
};


#endif  // _kvpatom_H

