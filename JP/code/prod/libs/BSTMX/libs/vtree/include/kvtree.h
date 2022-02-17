/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	kvtree.h
 * Function:	
 * Author:	
 * Revision:	$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/crxvtree/include/kvtree.h,v 1.2 2005/07/01 17:03:58 dliu Exp $
 ***************************************************************/
#ifndef	_kvtree_H
#define	_kvtree_H

#include "kstdinc.h"
#include "krate.h"
#include "kcplxrate.h"

class KTSlice;
class KTSliceNode;

extern	int      debugLevel;

//--------------------------------------------------------------
/**
 * Enumeration of the slice operations.
 */

enum KOper { COPY, ADD, SUB, MULT, POW, DIV, LOG, MAX, MIN, GEQ, LEQ, STVAR };
enum KTSComp { SAME, NEXT };


//--------------------------------------------------------------
/**
 * Pure virtual base class definition for tree.
 * <BR>
 * It provides a facility to internally map curves names
 * names into indices
 * through the methods MapCurveName, GetCurveName and GetCurveIdx.
 */

class	KVTree {
public:

	/** Destructor. */
virtual			~KVTree() {} ;

	/** Inserts a critical date. */
virtual	void		Insert(TDate critDate) = 0;

	/** Inserts a KZeroReset */
virtual void		Insert(const KZeroReset&, bool isCrit=TRUE) = 0;

	/** Inserts a KRateReset, and return the modelReset date */
virtual TDate		Insert(const KRateReset&, bool isCrit=TRUE) = 0;

	/** Inserts a KCplxRateReset, and return the modelReset date */
virtual TDate		Insert(const KCplxRateReset&, bool isCrit=TRUE) = 0;

	/** Inserts a KRateReset between reset date and end date, 
	 *  and return the modelReset date 
	 */
virtual TDate		Insert(const KRateReset&, TDate, bool isCrit=TRUE) = 0;

virtual TDate		Insert(const KCplxRateReset&, TDate, bool isCrit=TRUE) = 0;

	/** Gets a KZeroReset */
virtual void		Get(KTSlice& ts, const KZeroReset&) = 0;

	/** Gets a KRateReset */
virtual void		Get(KTSlice& ts, const KRateReset&) = 0;

	/** Gets a KCplxRateReset */
virtual void		Get(KTSlice& ts, const KCplxRateReset&) = 0;

	/** Inserts an array of critical dates. */
	void		InsertDateList(int numDates, TDate *critDates);

	/**
	 * Setup tree timeline.
	 */
virtual void            SetUpTimeline() = 0;

	/** Performs tree calibration. */
virtual	void		Calibrate() = 0;

	/** Updates tree at time node "tpIdx". */
virtual void		Update(int tpIdx) = 0;



	/* Dev a slice. */
virtual KTSlice&	TSliceDev(KTSlice&, 
				  const String &discCurveName = K_DEFAULT_NAME) = 0;

	/* Ev a slice. */
virtual KTSlice&	TSliceEv(KTSlice&) = 0;

	/** 
	 * Called by KSlice constructor.
	 * Computes and sets the slice dimension corresponding  
	 * to slice curveIdx, creates slice data etc.
	 * Alternatively  - is empty, and all the above functions are
	 * performed by a static (tree specific) routine that is
	 * only called by TS operations on "as needed" basis. */
virtual	KTSlice&	TSliceCreate(KTSlice&) = 0;

	/** Destroys a slice. */
virtual	void		TSliceDestroy(KTSlice&) = 0;

	/** Get the center value of a slice at time 0. */
virtual	double		TSliceGetCenter(KTSlice&) = 0;

	/**
	 * Sets the centar value ot a (co)-slice at time 0
	 * (adjoint of TSliceGetCenter).
	 */
virtual	void		TSliceSetCenter(KTSlice &ts, double value) = 0;



	/** Compare two slices based on comparision type. 
	 *  if SAME, two slices have same dimension and time point
	 *  if NEXT, two slices have same dimension and one time point
	 *  apart, ts1(t+1) and ts2(t). 
	 */
virtual	bool		TSliceCompare(KTSlice& ts1, KTSlice& ts2,
				      KTSComp compType = SAME) = 0;

	/**
	 * Performs a scalar operation on a slice
	 * with a double argument (slice is modified
	 * on exit. 
	 * Argument "what" determines the type of the
	 * operation: <BR>
	 * COPY sets the slice to the argument,<BR>
	 * ADD (SUB, MULT, DIV)  adds (subtracts, mutliplies,
	 * divides) the argument to the slice, <BR>
	 * MIN (MAX) takes the minimum (maximum) 
	 * of the argument and the slice and puts it in
	 * the slice.
	 * Special operations may be available
	 * in derived classes.
	 */
virtual	KTSlice&	TSliceScalarOper(KTSlice&, double, KOper) = 0;

	/**
	 * Performs a arithmetic operation on a slice
	 * (modified on exit) with another slice argument
	 * (unchanged on exit).
	 * Argument "what" determines the type of the
	 * operation: <BR>
	 * COPY sets the slice to the argument,<BR>
	 * ADD (SUB, MULT, DIV)  adds (subtracts, mutliplies,
	 * divides) the argument to the slice, <BR>
	 * MIN (MAX) takes the minimum (maximum) 
	 * of the argument and the slice and puts it in
	 * the slice.
	 * Special operations may be available
	 * in derived classes.
	 */
virtual	KTSlice&	TSliceUnaryOper(KTSlice&, const KTSlice&, KOper) = 0;


	/** Prints a slice on a stream. */
virtual	void		TSlicePut(KTSlice&, ostream& os, int minMaxOnly) = 0;

	/**
	 * Performs a special operation defined
	 * by the derived class.
	 */
virtual	void		TSliceSpecialOper(KTSlice&, char* what, ...) = 0;




	/**
	 * Create and set the begining of slice node 
	 */
virtual	KTSliceNode&	TSliceNodeBegin(KTSlice&, KTSliceNode&) = 0;

	/**
	 * Test the end of of slice node 
	 */
virtual	bool		TSliceNodeEnd(KTSliceNode&) = 0;

	/**
	 * Move to next slice node 
	 */
virtual	KTSliceNode&	TSliceNodeNext(KTSliceNode&) = 0;

	/**
	 * Get the value of slice node 
	 */
virtual	double&		TSliceAccessNodeValue(KTSlice&, KTSliceNode&) = 0;

	/**
	 *  Get the maximum difference of node values of slice 
	 *  in nearest neighbors of specified node 
	 */
virtual	double		TSliceGetNodeStepMax(KTSlice&, KTSliceNode&) = 0;

	/** Prints a slice node on a stream. */
virtual	void		TSliceNodePut(KTSliceNode&, ostream& os) = 0;



	/**
	 * Get the last time point index: the main tree
	 * loop is 
	 * for(tpIdx=vt.TPNum(); tpIdx>=; tpIdx--)...
	 */
virtual	int		TPNum() = 0;

	/** Returns today's date. */
virtual	TDate		TPToday() = 0;

	/** Returns the current time point index. */
virtual	int		TPIdxCurrent() = 0;

	/** Returns the date corresponding to the current timepoint. */
virtual	TDate		TPDateCurrent() = 0;

	/** Returns the year fraction to the current timepoint. */
virtual	double		TPTimeCurrent() = 0;

	/** Returns the time point index corresponding to the input date */
virtual int             TPIdx(TDate date) = 0;
			// use GetDLOffset function. need to be done

	/** Returns the date corresponding to the current timepoint. */
virtual	TDate		TPDate(int tpIdx) = 0;

	/** Return TRUE/FALSE whether the time point  is critical. */
virtual	int		TPIsCriticalDate(int tpIdx) = 0;

	/**
	 * Get IR discount curve name.
	 * a) for IR curve, just return itself,
	 * b) for Basis curve, return the IR reference discount curve,
	 * c) for CDS curve, return the IR reference discount curve,
	 */
virtual	const String&   GetIRDiscCurveName(const String& curveName)
                           {return curveName;}

	/** 
	 * Insert a curve into the index table.
	 * If the curve name already exists with the different index 
	 * it fails, otherwise with the same index it is ok.
	 */
	void             MapCurveName(const String& curveName, int idx);

	/**
	 * Get curve name from its index.
	 * Returns K_DEFAULT_NAME if index is not in the table.
	 */
	const String&    GetCurveName(int idx);             

	/**
	 * Get curve index from its name.
	 * Returns K_DEFAULT_IDX if name is not in the table,
	 * or if name=="Default".
	 */
	int              GetCurveIdx(const String& curveName);


	/** Returns value date of specified curve. */
virtual	TDate		 GetValueDate(int idx) = 0;

private:
	/**
	 * Internal map between names and indices.
	 */
	KMap(String,int)	mCurveIdxTable;

};


//--------------------------------------------------------------
/**
 * Class definition for an tree time slice.
 * This class does NOT contain any virtual functions.
 * At the level of the KVTree base class, the slices contain the
 * follwing data:<BR>
 * 1. a debugging name (tsName) used for error messages, <BR>
 * 2. a curve index to specify model dependent characteristics
 *    (such as discounting to be applied, etc.). The curve index
 *    is mapped to a curve name throught the curve name/index mapping
 *    facility of the KVTree. <BR>
 * 3. a slice dimension number (also model dependent).<BR>
 * 4. a model dependent data representation.
 */

class	KTSlice {
public:
	/**
	 * Creates a new empty slice. The curve name
	 * is used to specify the dimension, when applicable.
	 * Calls tree constructor TSliceCreate. 
	 */
			KTSlice(KVTree& vt,
				const String& tsName = K_DEFAULT_NAME,
				const String& curveName = K_DEFAULT_NAME);

	/** Copy constructor. */
			KTSlice(const KTSlice&);

	/** Slice destructor. */
			~KTSlice()
			{mVTree->TSliceDestroy(*this);}

	/** Sets the name of the slice. */
        void		SetSliceName(const String& name) {mName = name;}

	/** Gets the name of the slice. */
	const String&	GetSliceName() const { return mName;}

	/** Sets the curve index of the slice, given a curve NAME */
	void		SetCurveIdx(const String& curveName) 
	                {mCurveIdx = mVTree->GetCurveIdx(curveName);}

	/** Sets the curve index of the slice, given a curve index */
	void		SetCurveIdx(int curveIdx) 
	                {mCurveIdx = curveIdx;}


	/** Gets the curve index of the slice. */
        int     	GetCurveIdx() const { return mCurveIdx;}

	/** Gets the curve name of the slice - for error messages. */
        const String&   GetCurveName() const { return mVTree->GetCurveName(mCurveIdx);}

	/** Sets the dimension of the slice. */
	void		SetSliceDim(int sliceDim) { mSliceDim = sliceDim;}

	/** Gets the slice dimension */
	int     	GetSliceDim() const { return mSliceDim; }

	/** Gets the current timepoint of the slice. */
	int	        GetTpIdx() const { return mTpIdx;}

	/** Sets the current timepoint of the slice. */
	void	        SetTpIdx(int tpIdx) { mTpIdx = tpIdx;}

	/** Slice DEV operation. */
	KTSlice&	Dev(const String &discCurveName)
		    {return mVTree->TSliceDev(*this, discCurveName);}

	/** Slice EV operation. */
	KTSlice&	Ev()
		    {return mVTree->TSliceEv(*this);}

	/** Returns the center value of a slice at time zero. */
	double		GetCenter()
		    {return mVTree->TSliceGetCenter(*this);}

	/** Clears (empties) the time slice. */
	void		Clear()
		    {mVTree->TSliceDestroy(*this); mData = NULL;}

	/** Checks that two slices are consistent. */
	int		CheckConsistent(KTSlice &ts, KTSComp compType)
		    {return mVTree->TSliceCompare(*this, ts, compType);}

	/** Arithmetic operation with double argument.*/
	KTSlice&	operator= (const double value)
		    {return mVTree->TSliceScalarOper(*this, value, COPY);}

	/** Arithmetic operation with double argument.*/
	KTSlice&	operator+= (const double value)
		    {return mVTree->TSliceScalarOper(*this, value, ADD);}

	/** Arithmetic operation with double argument. */
	KTSlice&	operator-= (const double value)
		    {return mVTree->TSliceScalarOper(*this, value, SUB);}

	/** Arithmetic operation with double argument. */
	KTSlice&	operator*= (const double value)
		    {return mVTree->TSliceScalarOper(*this, value, MULT);}

	/** Arithmetic operation with double argument. */
	KTSlice&	operator/= (const double value)
		    {return mVTree->TSliceScalarOper(*this, value, DIV);}

	/**
	 * Takes the maximum of the argument and the slice
	 */
	KTSlice&	max(const double value)
		    {return mVTree->TSliceScalarOper(*this, value, MAX);}

	/**
	 * Takes the power of the slice
	 */
	KTSlice&	operator^= (const long value)
		    {return mVTree->TSliceScalarOper(*this, value, POW);}

	/**
	 * Takes the natural log of the slice
	 */
	KTSlice&	Log()
		    {return mVTree->TSliceScalarOper(*this, 1e0, LOG);}


	/**
	 * Takes the minimum of the argument and the slice
	 * and puts it in the slice.
	 */
	KTSlice&	min(const double value)
		    {return mVTree->TSliceScalarOper(*this, value, MIN);}

	/** Arithmetic operation with slice argument. */
	KTSlice&	operator= (const KTSlice& ts)
		    {return mVTree->TSliceUnaryOper(*this, ts, COPY);}

	/** Arithmetic operation with slice argument. */
	KTSlice&	operator+= (const KTSlice& ts)
		    {return mVTree->TSliceUnaryOper(*this, ts, ADD);}

	/** Arithmetic operation with slice argument. */
	KTSlice&	operator-= (const KTSlice& ts)
		    {return mVTree->TSliceUnaryOper(*this, ts, SUB);}

	/** Arithmetic operation with slice argument. */
	KTSlice&	operator*= (const KTSlice& ts)
		    {return mVTree->TSliceUnaryOper(*this, ts, MULT);}

	/** Arithmetic operation with slice argument. */
	KTSlice&	operator/= (const KTSlice& ts)
		    {return mVTree->TSliceUnaryOper(*this, ts, DIV);}


	/** Takes the maximum of the argument and the slice 
	 * and puts it in the slice.
	 */
	KTSlice&	max(const KTSlice& ts)
			    {return mVTree->TSliceUnaryOper(*this, ts, MAX);}


	/** Takes the minimum of the argument and the slice
	 * and puts it in the slice.
	 */
	KTSlice&	min(const KTSlice& ts)
			    {return mVTree->TSliceUnaryOper(*this, ts, MIN);}




	/** 
	 *  Slice node operator 
	 *
	 * [] operator 
	 */
	double&		operator[](KTSliceNode &tsNode)
		    	{return mVTree->TSliceAccessNodeValue(*this, tsNode);}

	/**
	 *  Get the maximum difference of node values of slice 
	 *  in nearest neighbors of specified node 
	 */
virtual	double		GetNodeStepMax(KTSliceNode &tsNode)
			{return mVTree->TSliceGetNodeStepMax(*this, tsNode);}



	/** 
	 *  A set of binary operations of slices 
	 *  Should be used with care because of overhead of
	 *  slice copy operation.
	 */
friend  KTSlice 	TSliceBinaryOper(const KTSlice &ts1, 
					 const KTSlice &ts2, 
					 KOper oper);

friend	KTSlice	operator+(const KTSlice &ts1, const KTSlice &ts2)
			{return TSliceBinaryOper(ts1, ts2, ADD);}

friend	KTSlice	operator-(const KTSlice &ts1, const KTSlice &ts2)
			{return TSliceBinaryOper(ts1, ts2, SUB);}

friend	KTSlice	operator*(const KTSlice &ts1, const KTSlice &ts2)
			{return TSliceBinaryOper(ts1, ts2, MULT);}

friend	KTSlice	operator/(const KTSlice &ts1, const KTSlice &ts2)
			{return TSliceBinaryOper(ts1, ts2, DIV);}


	/** Checks that a slice is nonempty and returns TRUE/FALSE.  */
	bool		IsEmpty() const { return (mData == NULL); }

	/** Checks that a slice is nonempty and creates 
	 * an exception if it is.
	 */
	void 		CheckNonEmpty(const String& debugName = "") const;


	/** Prints a time slice */
friend	ostream&	operator<<(ostream& os, KTSlice& ts)
			    { ts.mVTree->TSlicePut(ts, os, FALSE); return os;}


	/** Creates an array of slices. */
friend	KTSlice*	KTSliceNewVector(
				int size,
				KVTree& vt,
				const String& curveName = K_DEFAULT_NAME,
				const String& name = K_DEFAULT_NAME);

friend class KVTree;
friend class KTSliceNode;

public:
				/** Pointer to the virtual tree attached to. */
	KVTree	*mVTree;
				/** Slice dimension (model dependent). */
	int	mSliceDim;
				/** Slice data (model dependent.) */
	void	*mData;	


private:
				/** Slice name for debugging. */
	String	mName;
				/** Curve index (model dependent). */
	int	mCurveIdx;
				/** Current time node (model dependent). */
	int	mTpIdx;
				/** Private void slice constructor
				 * (used only internally).  */
		KTSlice()
		{ mVTree = NULL; SetTpIdx(-1); SetSliceDim(-1); mData = NULL; }


};


//--------------------------------------------------------------
/**
 * Class definition for node indices of tree time slice.
 */

class	KTSliceNode {
public:

	/** Default contructor */
			KTSliceNode();

	/** Slice destructor. */
			~KTSliceNode();

	/** Pointer to the first node in the slice */
	KTSliceNode&	begin(KTSlice &ts)
		{return ts.mVTree->TSliceNodeBegin(ts, *this);}

	/** True/False if it is the last node in the slice */
	bool		end()
		{return mSlice->mVTree->TSliceNodeEnd(*this);}

	/** Pointer to the next node in the slice */
	KTSliceNode&	next()
		{return mSlice->mVTree->TSliceNodeNext(*this);}

	/** Increment operation */
	KTSliceNode&	operator++()
		    {return next();}

	/** Prints a slice node */
friend  ostream&        operator<<(ostream& os, KTSliceNode& tsNode)
		   { tsNode.mSlice->mVTree->TSliceNodePut(tsNode, os); 
		     return os;}

friend class KVTree;


public:
				/** Slice node indices (model dependent). */
	void	*mIndex;
				/** Virtual slice on which to operate. */
	KTSlice	*mSlice;
};


#endif /* _kvtree_H */
