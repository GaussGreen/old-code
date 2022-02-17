/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pricingstates.h
 *
 *  \brief pricing states summarizes the state of the  world
 *
 *	\author  E. Benhamou JM Prie
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_PRICINGSTATES_H
#define _INGPINFRA_PRICINGSTATES_H

#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/port.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/rootobject.h"
#include "gpbase/mempoolmatrix.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// General Forward Declarations
class ARM_PricingStatesContext;

/// Forward Declarations for downcasts. 
class ARM_SFRMPricingStatesContext;
class ARM_Q1FPricingStatesContext;

//////////////////////////////////////////////
/// \class ARM_PayoffStates
/// \brief
/// This class contains the payoff
/// value for each future states
//////////////////////////////////////////////
class ARM_PayoffStates : public ARM_RootObject
{
private:
	// Payoff value for each states
	std::vector<double>							itsPayoffs;
	// Payoff value for each states and each flows
	ARM_GP_Matrix							itsIntermediatePayoffs;

    /// To save snapshots of payoff states
	ARM_VectorPtrVector itsPayoffSnapshots;

    /// To activate or not IntermediatePayoffs & Snapshots computations
    bool itsOtherPayoffsFlag;

	// To activate or not the Intermediate Values computations
	bool itsIVFlag;
public:
	ARM_PayoffStates(size_t nbStates =0,size_t nbIntermediatePayoffs=0,size_t nbSnapshots=0,bool otherPayoffsFlag=true,bool ivFlag=true) ;
	ARM_PayoffStates( const ARM_PayoffStates& rhs );
	ARM_PayoffStates& operator=( const ARM_PayoffStates& rhs );
	virtual ~ARM_PayoffStates();

	inline void SetPayoffs(const std::vector<double>& payoffs) {itsPayoffs=payoffs;}
	inline void SetPayoffs(const ARM_GP_Vector& payoffs) {itsPayoffs = payoffs.GetValues();}
	inline const std::vector<double>& GetPayoffs() { return itsPayoffs; }
	inline const ARM_GP_Vector& GetPayoffs(int i) { return ARM_GP_Vector(itsPayoffs); }
	/// stateIdx is the index for the state, while payoffIdx corresponds to the payoff nb
    inline double GetPayoff(size_t stateIdx) const {return itsPayoffs[stateIdx];}
	inline void SetPayoff(size_t stateIdx, double payoff) {itsPayoffs[stateIdx]=payoff;}

	/// Accessors (get, set matrix and double type) for the payoffs snapshot
	/// Accessors for payoff snapshots
    inline const ARM_GP_VectorPtr& GetPayoffSnapshot(size_t idx) const {return itsPayoffSnapshots[idx];}
    inline void SetPayoffSnapshot(size_t idx,const ARM_GP_VectorPtr& payoff) {itsPayoffSnapshots[idx]=payoff;}
    inline void SetPayoffSnapshot(size_t idx,const ARM_VectorPtr& payoff) {itsPayoffSnapshots[idx]= ARM_GP_VectorPtr( new ARM_GP_Vector(*payoff));}
    inline void push_backPayoffSnapshot(const ARM_GP_VectorPtr& payoff) {itsPayoffSnapshots.push_back(payoff);}
	inline size_t GetPayoffSnapshotsSize() const {return itsPayoffSnapshots.size();}
    inline const ARM_GP_Matrix& GetIntermediatePayoffs() const {return itsIntermediatePayoffs;}
    inline void SetIntermediatePayoffs(const ARM_GP_Matrix& intermediatePayoffs) {itsIntermediatePayoffs=intermediatePayoffs;}
    inline double GetIntermediatePayoff(size_t stateIdx,size_t payoffIdx) const {return itsIntermediatePayoffs(payoffIdx,stateIdx);}
    inline void SetIntermediatePayoff(size_t stateIdx,size_t payoffIdx,double payoff) {itsIntermediatePayoffs(payoffIdx,stateIdx)=payoff;}

    inline void resizePayoff(size_t nbStates) {itsPayoffs.resize(nbStates);}
	inline void resizeIntermediatePayoffs(size_t nbPayoffs,size_t stateSize = 0) {itsIntermediatePayoffs.resize(nbPayoffs,stateSize? stateSize : itsIntermediatePayoffs.GetColsNb());}
	inline void resizePayoffSnapshots(size_t nbPayoffSnapshots) {itsPayoffSnapshots.resize(nbPayoffSnapshots);}

	inline bool GetOtherPayoffsFlag() const { return itsOtherPayoffsFlag; }
    inline void SetOtherPayoffsFlag(bool otherPayoffsFlag) { itsOtherPayoffsFlag = otherPayoffsFlag; }
	inline bool GetIVFlag() const { return itsIVFlag; }
    inline void SetIVFlag(bool IVFlag) { itsIVFlag = IVFlag; }

	inline size_t GetIntermediatePayoffsSize() const {return itsIntermediatePayoffs.GetRowsNb();}
	inline size_t size() const {return itsPayoffs.size();}
    inline size_t GetIntermediatePayoffStatesSize() const {return itsIntermediatePayoffs.GetColsNb();}

	/// Standard ARM Object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const { return "ARM_PayoffStates"; }
};


//////////////////////////////////////////////
/// \class ARM_PricingStates
/// \brief
/// ARM_PricingStates class defines a global
/// future state with model state variables
/// and associated payoffs
//////////////////////////////////////////////
class ARM_PricingStates : public ARM_RootObject
{
private:
    /// Rows = Number of ModelStates or Payoffs
    /// Cols = Number of States
	ARM_GP_MatrixPtr						itsModelStates;
	ARM_GP_MatrixPtr						itsNumMethodStates;
	ARM_PricingStatesContextPtrVectorPtr	itsPricingStatesContextVector;
	ARM_MemPool_Matrix						itsPayoffs;
	ARM_MemPool_Matrix						itsProbaChanges;
	vector<ARM_PayoffStates>				itsPayoffStatesVector;

public:
	ARM_PricingStates(
		size_t nbStates =0, 
		size_t nbModelStates =0, 
		size_t nbPayoffs=0,
        size_t nbNumMethodStates=0,
		size_t nbIntermediatePayoffs=0,
		size_t nbPayoffStates=0,
		size_t nbSnapshots=0,
        bool otherPayoffsFlag=true,
		size_t nbProbaChanges=0);
	ARM_PricingStates( const ARM_PricingStates& rhs );
	ARM_PricingStates& operator=( const ARM_PricingStates& rhs );
	virtual ~ARM_PricingStates();

    /// Accessors (get, set matrix and double type) for the modelstates,numethodstates & payoff
    inline const ARM_GP_MatrixPtr& GetModelStates() const {return itsModelStates;}
    inline ARM_GP_MatrixPtr& GetModelStates() {return itsModelStates;}
	inline const ARM_GP_MatrixPtr& GetNumMethodStates() const {return itsNumMethodStates;}
    inline ARM_GP_MatrixPtr& GetNumMethodStates() {return itsNumMethodStates;}

	inline void SetModelStates(const ARM_GP_MatrixPtr& modelStates) {itsModelStates=modelStates;}
	inline void SetNumMethodStates(const ARM_GP_MatrixPtr& modelStates) {itsNumMethodStates=modelStates;}
    inline double GetModelState(size_t stateIdx,size_t modelStateIdx) const {return (*itsModelStates)(modelStateIdx,stateIdx);}
    inline void SetModelState(size_t stateIdx,size_t modelStateIdx,double modelState) {(*itsModelStates)(modelStateIdx,stateIdx)=modelState;}
    inline double GetNumMethodState(size_t stateIdx,size_t modelStateIdx) const {return (*itsNumMethodStates)(modelStateIdx,stateIdx);}
    inline void SetNumMethodState(size_t stateIdx,size_t modelStateIdx,double modelState) { (*itsNumMethodStates)(modelStateIdx,stateIdx)=modelState;}
    
	inline void SetPayoffs(const ARM_MemPool_Matrix& payoffs) {itsPayoffs=payoffs;}
	inline void SetProbaChanges(const ARM_MemPool_Matrix& probaChanges) {itsProbaChanges = probaChanges;}
	/// stateIdx is the index for the state, while payoffIdx corresponds to the payoff nb
    inline double GetPayoff(size_t stateIdx,size_t payoffIdx) const {return itsPayoffs(payoffIdx,stateIdx);}
	inline void SetPayoff(size_t stateIdx,size_t payoffIdx,double payoff) {itsPayoffs(payoffIdx,stateIdx)=payoff;}
	inline double GetProbaChange(size_t stateIdx,size_t probaChangeIdx) const {return itsProbaChanges(probaChangeIdx,stateIdx);}
// FIXMEFRED: mig.vc8 (22/05/2007 18:06:21):cast
	inline ARM_GP_VectorPtr GetProbaChangeVec(size_t probaChangeIdx) { return static_cast<ARM_GP_VectorPtr>(new ARM_GP_T_Vector<double>(itsProbaChanges.GetRow(probaChangeIdx))); }
	inline void SetProbaChange(size_t stateIdx,size_t probaChangeIdx,double probaChange) {itsProbaChanges(probaChangeIdx,stateIdx)=probaChange;}

	inline ARM_PayoffStates& GetPayoffStates(size_t idx) { return itsPayoffStatesVector[idx]; }

// FIXMEFRED: mig.vc8 (22/05/2007 18:06:18):cast
	inline ARM_GP_VectorPtr GetPayoffVec(size_t payoffIdx) {return static_cast<ARM_GP_VectorPtr>(new ARM_GP_Vector(itsPayoffs.GetRow(payoffIdx)));}
	inline ARM_MemPool_Matrix& GetPayoffs() { return itsPayoffs; }
	/// the indexVec is const because it is assumed to be already sorted!
    inline void RemovePayoffVec( const vector<size_t>& indexVec) {itsPayoffs.Remove_NRowsSorted(indexVec);}

	/// Accessors (get, set matrix and double type) for the payoffs snapshot
    inline void resizePayoffs(size_t nbPayoffs, size_t stateSize = 0) {itsPayoffs.resize(nbPayoffs, stateSize? stateSize: itsPayoffs.GetColsNb());}
	inline void resizeProbaChanges(size_t nbProbaChanges, size_t stateSize = 0) {itsProbaChanges.resize(nbProbaChanges, stateSize? stateSize: itsProbaChanges.GetColsNb());}
    inline void resizeModelStates(size_t nbModelStates,size_t stateSize = 0) {itsModelStates->resize(nbModelStates,stateSize ? stateSize  : itsModelStates->GetColsNb());}
    inline void resizeNumMethodStates(size_t nbModelStates, size_t stateSize = 0) {itsNumMethodStates->resize(nbModelStates, stateSize ? stateSize : itsNumMethodStates->GetColsNb());}
	void resizePayoffStatesVector(size_t nbPayoffStates);

    inline size_t ModelStatesSize() const {return itsModelStates->GetRowsNb();}
	inline size_t NumMethodStatesSize() const {return itsNumMethodStates->GetRowsNb();}
    inline size_t GetPayoffsSize() const {return itsPayoffs.GetRowsNb();}
	inline size_t GetProbaChangesSize() const {return itsProbaChanges.GetRowsNb(); }
	inline size_t GetStatesOfPayoffsSize() const {return itsPayoffs.GetColsNb();}
	inline size_t GetStatesOfProbaChangesSize() const {return itsProbaChanges.GetColsNb();}
	inline size_t GetPayoffStatesSize() const {return itsPayoffStatesVector.size();}
    inline size_t size() const {return itsNumMethodStates->GetRowsNb()>0 ? itsNumMethodStates->GetColsNb() : (itsModelStates->GetRowsNb()>0 ? itsModelStates->GetColsNb() : 0);}
    
	inline ARM_MemPool_Matrix::iterator payoffsBeginIterator( size_t payoffIdx ) { return itsPayoffs.GetIthRowIterator(payoffIdx);}

	inline void FreeMemoryPool() { itsPayoffs.FreeMemoryPool(); }

	/// Accessors for the pricing states context
    inline const ARM_PricingStatesContextPtrVectorPtr& GetPricingStatesContextVector() const {return itsPricingStatesContextVector; }
    inline void SetPricingStatesContextVector(const ARM_PricingStatesContextPtrVectorPtr& pricingStatesContextVector) { itsPricingStatesContextVector = pricingStatesContextVector;}
	inline const ARM_PricingStatesContextPtr& GetPricingStatesContext( int pricingStatesContextIdx ) const { return (*itsPricingStatesContextVector)[pricingStatesContextIdx]; }
	inline void EnqueuePricingStatesContext( const ARM_PricingStatesContextPtr& pricingStatesContext ) const { itsPricingStatesContextVector->push_back( pricingStatesContext ); }
	void AddPricingStatesContextVector( const ARM_PricingStatesContextPtrVectorPtr& pricingStatesContextVector ) const;

	/// Copy of the PricingStatesContextPtrVector (each PricingStatesContext is duplicated)
	const ARM_PricingStatesContextPtrVectorPtr GetCopyOfPricingStatesContextVector() const;

    void SetOtherPayoffsFlag(bool otherPayoffsFlag);
	void SetIVFlag(bool ivFlag);

	/// Standard ARM Object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const { return "ARM_PricingStates"; }
};

//////////////////////////////////////////////
/// \class ARM_PricingStatesContext
/// \brief
/// This class contains some specific data
/// linked to a model
//////////////////////////////////////////////
class ARM_PricingStatesContext : public ARM_RootObject
{
public: 
	/// Downcasting
	virtual ARM_SFRMPricingStatesContext * ToSFRMPricingStatesContext() { ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "This is not an SFRMPRicingStatesContext Object" ); }
	virtual ARM_Q1FPricingStatesContext * ToQ1FPricingStatesContext() { ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "This is not an Q1FPRicingStatesContext Object" ); }

	/// destructor
	virtual ~ARM_PricingStatesContext() {} ;

	/// standard ARM_RootObject Support
	virtual ARM_Object * Clone() const { return new ARM_PricingStatesContext();}
	virtual string toString(const string& indent="", const string& nextIndent="") const { return "ARM_SFRMPricingStates"; };
};


//////////////////////////////////////////////
/// \class PayoffFunc, GetAuxiliaryPayoffFunc, 
/// GetPayoff, GetIntermediatePayoff
/// \brief
/// These class are functor to access the 
/// the payoff of the pricing states
/// They will be used by the ComputeExpectation
/// method of the slice.
//////////////////////////////////////////////

class PayoffFunc
{
public:
	virtual ~PayoffFunc() {}; 

	virtual double operator()(size_t stateIdx) const = 0;
};

class GetAuxiliaryPayoffFunc : public PayoffFunc
{
public:
	GetAuxiliaryPayoffFunc(
		const ARM_PricingStatesPtr& states)
	: 
	itsStates(states)
	{
	}

	virtual ~GetAuxiliaryPayoffFunc() {};

private:
	GetAuxiliaryPayoffFunc(
		const GetAuxiliaryPayoffFunc& rhs);
	GetAuxiliaryPayoffFunc& operator=(const GetAuxiliaryPayoffFunc& rhs);

public:
	
	void SetPayoffIdx(int payoffIdx) { itsPayoffIdx = payoffIdx; }

	virtual double operator()(size_t stateIdx) const
	{
		return itsStates->GetPayoff(stateIdx, itsPayoffIdx);
	}

private:
	ARM_PricingStatesPtr itsStates;
	int itsPayoffIdx;
};

class GetProbaChangeFunc : public PayoffFunc
{
public:
	GetProbaChangeFunc(
		const ARM_PricingStatesPtr& states)
	: 
	itsStates(states)
	{
	}

	virtual ~GetProbaChangeFunc() {};

private:
	GetProbaChangeFunc(
		const GetProbaChangeFunc& rhs);
	GetProbaChangeFunc& operator=(const GetProbaChangeFunc& rhs);

public:
	
	void SetProbaChangeIdx(int probaChangeIdx) { itsProbaChangeIdx = probaChangeIdx; }

	virtual double operator()(size_t stateIdx) const
	{
		return itsStates->GetProbaChange(stateIdx, itsProbaChangeIdx);
	}

private:
	ARM_PricingStatesPtr itsStates;
	int itsProbaChangeIdx;
};

class GetPayoffFunc : public PayoffFunc
{
public:
	GetPayoffFunc(
		const ARM_PricingStatesPtr& states)
	: 
	itsStates(states)
	{ 
	}

	virtual ~GetPayoffFunc() {};
	
private:
	GetPayoffFunc(const GetPayoffFunc&);
	GetPayoffFunc& operator=(const GetPayoffFunc&);

public:
	
	void SetPayoffIdx(int payoffIdx) { itsPayoffIdx = payoffIdx; }

	virtual double operator()(size_t stateIdx) const
	{
		return itsStates->GetPayoffStates(itsPayoffIdx).GetPayoff(stateIdx);
	}

private:
	ARM_PricingStatesPtr itsStates;
	int itsPayoffIdx;
};


class GetIntermediatePayoffFunc : public PayoffFunc
{
public:
	GetIntermediatePayoffFunc(
		const ARM_PricingStatesPtr& states)
	: 
	itsStates(states)
	{
	}

	virtual ~GetIntermediatePayoffFunc() {};

	void SetPayoffIdx(int payoffIdx) { itsPayoffIdx = payoffIdx; }
	void SetInterIdx(int interIdx) { itsInterIdx = interIdx; }


private:
	GetIntermediatePayoffFunc(const GetIntermediatePayoffFunc&);
	GetIntermediatePayoffFunc& operator=(const GetIntermediatePayoffFunc&);

public:
	virtual double operator()(size_t stateIdx) const
	{
		return itsStates->GetPayoffStates(itsPayoffIdx).GetIntermediatePayoff(stateIdx,itsInterIdx);
	}

private:
	ARM_PricingStatesPtr itsStates;
	int itsPayoffIdx;
	int itsInterIdx;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

