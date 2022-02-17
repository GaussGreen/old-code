/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file slice.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPNUMMETHODS_SLICE_H
#define _INGPNUMMETHODS_SLICE_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpbase/rootobject.h"
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"

/// gpinfra
#include "gpinfra/pricingstates.h"

/// gpnumlib
#include "gpnumlib/numfunction.h"

/// gpnummethods
#include "gpnummethods/typedef.h"
#include "gpnummethods/transitor.h"
#include "gpnummethods/treeindex.h"

#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////
/// SLICE PART
/// a slice of a tree.
////////////////////////////////////

/// forward declaration
class ARM_Node1D;
class ARM_NodeND;
struct ARM_SamplerBase;
struct ARM_TruncatorBase;

class ARM_Slice1DBase;
class ARM_Slice1D;
class ARM_Slice1DCstSymProba;
class ARM_SliceNDBase;
class ARM_SliceND;
class ARM_SliceNDCstSymProba;
class ARM_PricingStates;
class ARM_TransitionStates;

///---------------------------
/// ARM_SliceBase
///---------------------------
struct ARM_SliceBase : public ARM_RootObject
{
private:
    size_t itsIndex;
	ARM_GP_VectorPtr itsSpotProbas;
	ARM_GP_VectorPtr itsArrowDebreuPrices;

public:
    ARM_SliceBase( size_t index = 0) : itsIndex(index) {}
    ARM_SliceBase(const ARM_SliceBase& rhs);

    /// downcast operators
    virtual ARM_Slice1DBase* ToSlice1DBase();
    virtual const ARM_Slice1DBase* ToSlice1DBase() const;

    virtual ARM_Slice1D* ToSlice1D();
    virtual const ARM_Slice1D* ToSlice1D() const;

    virtual ARM_Slice1DCstSymProba* ToSlice1DCstSymProba();
    virtual const ARM_Slice1DCstSymProba* ToSlice1DCstSymProba() const;

    virtual ARM_SliceNDBase* ToSliceNDBase();
    virtual const ARM_SliceNDBase* ToSliceNDBase() const;

    virtual ARM_SliceND* ToSliceND();
    virtual const ARM_SliceND* ToSliceND() const;

    virtual ARM_SliceNDCstSymProba* ToSliceNDCstSymProba();
    virtual const ARM_SliceNDCstSymProba* ToSliceNDCstSymProba() const;

	/// accessors to spot and Arrow Debreu prices
	ARM_GP_VectorPtr GetSpotProbas() const { return itsSpotProbas; }
	double GetSpotProbas(size_t i) const { return (*itsSpotProbas)[i]; }
	ARM_GP_VectorPtr GetArrowDebreuPrices() const { return itsArrowDebreuPrices; }
	double GetArrowDebreuPrices(size_t i) const { return (*itsArrowDebreuPrices)[i]; }
	void SetSpotProbas( const ARM_GP_VectorPtr& spotProbas ) { itsSpotProbas=spotProbas; }
	void SetArrowDebreuPrices( const ARM_GP_VectorPtr& arrowDebreuPrices ) { itsArrowDebreuPrices=arrowDebreuPrices; }

    /// accessor
    size_t GetIndex() const { return itsIndex; }
    virtual ARM_GP_VectorPtr GetDriftCorrectionVect() const=0;
    virtual void SetDriftCorrectionVect(const ARM_GP_VectorPtr& driftCorrection)=0;

    virtual size_t size() const = 0;

    virtual void LinkedWithNextSlice( ARM_SliceBase* nextSlice, ARM_SamplerBase* sampler, ARM_TruncatorBase* truncator, bool needProbas, ARM_TransitionStates& transStates, const ARM_GP_VectorPtr& driftCorrection, size_t driftOffset, double targetPayoff) = 0;
    virtual double ComputeExpectation( size_t stateIdx, const ARM_SliceBase* nextSlice, const PayoffFunc& payoffFunc, ARM_TransitionStates& transStates) const = 0;

    virtual ARM_GP_VectorPtr ComputeDriftCorrection(const ARM_PricingModel* model, ARM_SamplerBase* sampler, double nextTime, size_t driftOffset, double targetPayoff, bool isLocalDf=true) = 0;
};



///---------------------------
/// Slice1DFunct
///---------------------------
class Slice1DFunct : public CC_NS( ARM_GP, UnaryFunc)<double,double>
{
private:
    size_t itsTimeIdx;
    double itsDt;
    const ARM_PricingModel* itsModel;
    ARM_GP_VectorPtr itsXStates;
    ARM_GP_VectorPtr itsADPrices;
    mutable ARM_PricingStatesPtr itsStates;

public:
    Slice1DFunct(size_t timeIdx, double dt, const ARM_PricingModel* model, const ARM_GP_VectorPtr& xStates, const ARM_GP_VectorPtr& adPrices);
    virtual ~Slice1DFunct() {}

    virtual double operator () (double x) const;

    ARM_GP_VectorPtr GetLocalDf(double x) const;

};


///---------------------------
/// ARM_Slice1DBase
///---------------------------
class ARM_Slice1DBase : public ARM_SliceBase
{
private:
    int itsMin;
	int itsMax;
	double itsSpaceStep;

    double itsDriftCorrection;

    /// Initial process states, X states, are save for efficency
    ARM_GP_VectorPtr itsXStates;

protected:
    void UpdateProbasAndArrowDebreu( bool needProbas, bool needADPrices, bool isLocalDf, ARM_GP_VectorPtr localDf, ARM_Slice1DBase* nextSlice1D);

public:
    ARM_Slice1DBase( size_t index, int min, int max, double spaceStep)
    :   ARM_SliceBase(index),  itsMin( min ), itsMax( max ), itsSpaceStep( spaceStep ), itsDriftCorrection(0.0)
    {}
    ARM_Slice1DBase(const ARM_Slice1DBase& rhs);

    virtual ~ARM_Slice1DBase() {};
    virtual ARM_Object* Clone() const = 0;

    virtual ARM_Slice1DBase* ToSlice1DBase() { return this; }
    virtual const ARM_Slice1DBase* ToSlice1DBase() const { return this; }

    /// Accessors
    int GetMin() const { return itsMin; }
    void SetMin(int min) { itsMin = min; }
    int GetMax() const { return itsMax; }
    void SetMax(int max) { itsMax = max; }

    virtual ARM_GP_VectorPtr GetDriftCorrectionVect() const { return ARM_GP_VectorPtr( new ARM_GP_Vector(1,itsDriftCorrection) ); }
    virtual void SetDriftCorrectionVect(const ARM_GP_VectorPtr& driftCorrection) { itsDriftCorrection = (*driftCorrection)[0]; }
    double GetDriftCorrection() const { return itsDriftCorrection; }
    void SetDriftCorrection(double driftCorrection) { itsDriftCorrection=driftCorrection; }

    double GetSpaceStep() const { return itsSpaceStep; }
    void SetSpaceStep(double spaceStep) { itsSpaceStep=spaceStep; }

    const ARM_GP_VectorPtr& GetXStates() const { return itsXStates; }
    void SetXStates(const ARM_GP_VectorPtr& xStates) { itsXStates=xStates; }

	virtual double GetProbaUp(size_t stateIdx) const = 0;
	virtual double GetProbaDown(size_t stateIdx) const = 0;
	virtual double GetNextNodeIndex(size_t stateIdx) const = 0;

    /// Compute the slice size (=number of states)
    size_t size() const { return itsMax-itsMin+1; }

    /// Generate states of the actual process sampled by the slice (Z states)
    inline ARM_GP_VectorPtr GetZStates() const;

    /// Compute the drift correction to fufill AOA
    virtual ARM_GP_VectorPtr ComputeDriftCorrection(const ARM_PricingModel* model, ARM_SamplerBase* sampler, double nextTime, size_t driftOffset, double targetPayoff, bool isLocalDf=true);

};

inline ARM_GP_VectorPtr ARM_Slice1DBase::GetZStates() const
{
    size_t nbStates = size();
    ARM_GP_VectorPtr zStates( new ARM_GP_Vector( nbStates, 0.0 ) );
    int i; /// take care : size_t i => strange behaviour !
    for( i=0; i<nbStates; ++i )
        (*zStates)[i] = (i+GetMin()) * GetSpaceStep();

    return zStates;
}


///---------------------------
/// ARM_Slice1D
///---------------------------
class ARM_Slice1D : public ARM_Slice1DBase
{
private:
    /// For efficiency needs in creation/destruction no more vector of 1D node objects
	std::vector<double> itsProbaUp;
	std::vector<double> itsProbaDown;
	ARM_IntVector itsNextNodeIndex;

public:
    ARM_Slice1D( size_t index, int min, int max, double spaceStep)
        :   ARM_Slice1DBase(index,min,max,spaceStep), itsProbaUp(max-min+1), itsProbaDown(max-min+1), itsNextNodeIndex(max-min+1)
    {}
    ARM_Slice1D(const ARM_Slice1D& rhs)
        :   ARM_Slice1DBase(rhs), itsProbaUp(rhs.itsProbaUp), itsProbaDown(rhs.itsProbaDown), itsNextNodeIndex(rhs.itsNextNodeIndex)
    {}
    ARM_Slice1D& operator = (const ARM_Slice1D& rhs);

    virtual ~ARM_Slice1D() {}
    virtual ARM_Object* Clone() const { return new ARM_Slice1D(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const;

    virtual ARM_Slice1D* ToSlice1D() { return this; }
    virtual const ARM_Slice1D* ToSlice1D() const { return this; }

    /// Accessors
	virtual double GetProbaUp(size_t stateIdx) const {  return itsProbaUp[stateIdx]; }
	void SetProbaUp(size_t stateIdx, double proba) {  itsProbaUp[stateIdx]=proba; }
	virtual double GetProbaDown(size_t stateIdx) const {  return itsProbaDown[stateIdx]; }
	void SetProbaDown(size_t stateIdx, double proba) { itsProbaDown[stateIdx]=proba; }
	virtual double GetNextNodeIndex(size_t stateIdx) const {  return itsNextNodeIndex[stateIdx]; }
	void SetNextNodeIndex(size_t stateIdx, int nextNodeIndex) { itsNextNodeIndex[stateIdx]=nextNodeIndex; }

    /// Slice construction
    virtual void LinkedWithNextSlice( ARM_SliceBase* nextSlice, ARM_SamplerBase* sampler, ARM_TruncatorBase* truncator, bool needProbas, ARM_TransitionStates& transStates, const ARM_GP_VectorPtr& driftCorrection, size_t driftOffset, double targetPayoff);
    virtual double ComputeExpectation( size_t stateIdx, const ARM_SliceBase* nextSlice, const PayoffFunc& payoffFunc, ARM_TransitionStates& transStates) const;
};


///---------------------------
/// ARM_Slice1DCstSymProba
///---------------------------
class ARM_Slice1DCstSymProba : public ARM_Slice1DBase
{
private:
	double itsProba;                // cst proba to go up or down for the whole slice
    std::vector<double> itsProbaUpDown;   // to optimise computation time : may be 0 if truncated hedge

public:
    ARM_Slice1DCstSymProba( size_t index, int min, int max, double spaceStep, double proba )
        :   ARM_Slice1DBase(index,min,max,spaceStep), itsProba( proba ), itsProbaUpDown( max-min+1, proba )
    {}
    ARM_Slice1DCstSymProba(const ARM_Slice1DCstSymProba& rhs)
        :   ARM_Slice1DBase(rhs), itsProba(rhs.itsProba), itsProbaUpDown(rhs.itsProbaUpDown)
    {}
    ARM_Slice1DCstSymProba& operator = (const ARM_Slice1DCstSymProba& rhs);

    virtual ~ARM_Slice1DCstSymProba() {}

    virtual ARM_Object* Clone() const { return new ARM_Slice1DCstSymProba(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const;

    virtual ARM_Slice1DCstSymProba* ToSlice1DCstSymProba() { return this; }
    virtual const ARM_Slice1DCstSymProba* ToSlice1DCstSymProba() const { return this; }

    virtual double GetProbaUp(size_t stateIdx) const {  return itsProbaUpDown[stateIdx]; }
    virtual double GetProbaDown(size_t stateIdx) const {  return itsProbaUpDown[stateIdx]; }
    virtual double GetNextNodeIndex(size_t stateIdx) const {  return stateIdx + GetMin(); }

    /// Link with next slice is implicit but truncation may update slice range
    virtual void LinkedWithNextSlice( ARM_SliceBase* nextSlice, ARM_SamplerBase* sampler, ARM_TruncatorBase* truncator, bool needProbas, ARM_TransitionStates& transStates, const ARM_GP_VectorPtr& driftCorrection, size_t driftOffset, double targetPayoff);
    virtual double ComputeExpectation( size_t stateIdx, const ARM_SliceBase* nextSlice, const PayoffFunc& payoffFunc, ARM_TransitionStates& transStates) const;
};


///---------------------------
/// SliceNDFunct
///---------------------------
class SliceNDFunct : public CC_NS( ARM_GP, UnaryFunc)<double,double>
{
private:
    size_t itsTimeIdx;
    double itsDt;
    const ARM_PricingModel* itsModel;
    ARM_GP_MatrixPtr itsXStates;
    ARM_GP_VectorPtr itsADPrices;
    ARM_GP_VectorPtr itsDriftCorrection;
    size_t itsDriftOffset;
    mutable ARM_PricingStatesPtr itsStates;

public:
    SliceNDFunct(size_t timeIdx, double dt, const ARM_PricingModel* model, const ARM_GP_MatrixPtr& xStates, const ARM_GP_VectorPtr& adPrices,const ARM_GP_VectorPtr& driftCorrection,size_t driftOffset);
    virtual ~SliceNDFunct() {}

    virtual double operator () (double x) const;

    ARM_GP_VectorPtr GetLocalDf(double x) const;

    const ARM_PricingStatesPtr& GetStates() const { return itsStates; }


};


///---------------------------
/// ARM_SliceNDBase
///---------------------------
class ARM_SliceNDBase : public ARM_SliceBase
{
private:
	ARM_IntVector itsMin;
	ARM_IntVector itsMax;
	std::vector<double> itsSpaceStep;
    ARM_TransitorBase* itsTransitor;

    ARM_GP_VectorPtr itsDriftCorrection;

    /// Multi-indexes of connected nodes (explicitly or implicitly)
    /// No more ND node object for efficiency needs
	ARM_IntMatrix itsNextNodeIndex;

    /// Initial process states, X states, are save for efficency
    ARM_GP_MatrixPtr itsXStates;

protected:
    void UpdateProbasAndArrowDebreu( bool needProbas, bool needADPrices, bool isLocalDf, ARM_GP_VectorPtr localDf, ARM_SliceBase* nextSlice, ARM_TransitionStates& transStates);

public:
    ARM_SliceNDBase( size_t nbDims, size_t index, const ARM_IntVector& min, const ARM_IntVector& max, const std::vector<double>& spaceStep)
        :   ARM_SliceBase(index), itsMin( min ), itsMax( max ), itsSpaceStep( spaceStep ), itsNextNodeIndex(0,0), itsTransitor(NULL), itsDriftCorrection(NULL)
    { itsTransitor = ARM_TransitorFactory::CreateTransitor( nbDims ); }
    ARM_SliceNDBase(const ARM_SliceNDBase& rhs);

    virtual ~ARM_SliceNDBase() { delete itsTransitor; }

    virtual ARM_Object* Clone() const = 0;

    virtual ARM_SliceNDBase* ToSliceNDBase() { return this; }
    virtual const ARM_SliceNDBase* ToSliceNDBase() const { return this; }

    /// Accessors
    int GetMin(size_t i) const { return itsMin[i]; }
    void SetMin(size_t i,int min) { itsMin[i] = min; }
    const ARM_IntVector& GetMin() const { return itsMin; }
    void SetMin(const ARM_IntVector& min) { itsMin = min; }

    virtual ARM_GP_VectorPtr GetDriftCorrectionVect() const { return itsDriftCorrection; }
    virtual void SetDriftCorrectionVect(const ARM_GP_VectorPtr& driftCorrection) { itsDriftCorrection=driftCorrection; }

    int GetMax(size_t i) const { return itsMax[i]; }
    void SetMax(size_t i,int max) { itsMax[i] = max; }
    const ARM_IntVector& GetMax() const { return itsMax; }
    void SetMax(const ARM_IntVector& max) { itsMax = max; }

    const std::vector<double>& GetSpaceStep() const { return itsSpaceStep; }
    double GetSpaceStep(size_t i) const { return itsSpaceStep[i]; }
    void SetSpaceStep(const std::vector<double>& spaceStep) { itsSpaceStep=spaceStep; }
    void SetSpaceStep(size_t i, double spaceStep) { itsSpaceStep[i]=spaceStep; }

	void ResizeNextNodeIndex(size_t nbStates, size_t nbDims) { itsNextNodeIndex.resize(nbStates,nbDims); }
	int GetNextNodeIndex(size_t stateIdx, size_t dimIdx) const {  return itsNextNodeIndex(stateIdx,dimIdx); }
	void SetNextNodeIndex(size_t stateIdx, size_t dimIdx, int nextNodeIndex) { itsNextNodeIndex(stateIdx,dimIdx)=nextNodeIndex; }

    const ARM_TransitorBase& GetTransitor() const { return *itsTransitor; }

    const ARM_GP_MatrixPtr& GetXStates() const { return itsXStates; }
    void SetXStates(const ARM_GP_MatrixPtr& xStates) { itsXStates=xStates; }

	size_t dim() const { return itsSpaceStep.size(); }

    /// Compute the slice size (=number of states)
    inline virtual size_t size() const;

    /// Get proba up & down of the diffusion in the ith direction
    inline virtual double GetProbaUp(size_t stateIdx, size_t dimIdx) const = 0;
    inline virtual double GetProbaDown(size_t stateIdx, size_t dimIdx) const = 0;

    /// Generate states of the actual process sampled by the slice (Z states)
    inline virtual ARM_GP_MatrixPtr GetZStates() = 0;

    /// Compute the conditional expectation for a current state and payoff
    virtual double ComputeExpectation( size_t stateIdx, const ARM_SliceBase* nextSlice, const PayoffFunc& payoffFunc, ARM_TransitionStates& transStates) const;

    /// Compute the drift correction to fufill AOA
    virtual ARM_GP_VectorPtr ComputeDriftCorrection(const ARM_PricingModel* model, ARM_SamplerBase* sampler, double nextTime, size_t driftOffset, double targetPayoff, bool isLocalDf=true);
};

inline size_t ARM_SliceNDBase::size() const
{
    size_t size=1;
    for(size_t i=0;i<itsMin.size();++i)
        size *= itsMax[i]-itsMin[i]+1;
    return size;
}

///---------------------------
/// ARM_SliceND
///---------------------------
class ARM_SliceND : public ARM_SliceNDBase
{
private:
    /// For efficiency in creation/destruction no more vector of ND node objects
	ARM_GP_Matrix itsProbaUp;
	ARM_GP_Matrix itsProbaDown;

public:
    ARM_SliceND( size_t nbDims, size_t index, const ARM_IntVector& min, const ARM_IntVector& max, const std::vector<double>& spaceStep)
        :   ARM_SliceNDBase(nbDims,index,min,max,spaceStep), itsProbaUp(0,0), itsProbaDown(0,0)
    {}
    ARM_SliceND(const ARM_SliceND& rhs)
        : ARM_SliceNDBase(rhs), itsProbaUp(rhs.itsProbaUp), itsProbaDown(rhs.itsProbaDown)
    {}
    ARM_SliceND& operator = (const ARM_SliceND& rhs);

    virtual ~ARM_SliceND() {}

    virtual ARM_Object* Clone() const { return new ARM_SliceND(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const;

    virtual ARM_SliceND* ToSliceND() { return this; }
    virtual const ARM_SliceND* ToSliceND() const { return this; }

    /// Accessors
	inline virtual double GetProbaUp(size_t stateIdx, size_t dimIdx) const {  return itsProbaUp(stateIdx,dimIdx); }
	void SetProbaUp(size_t stateIdx, size_t dimIdx, double proba) {  itsProbaUp(stateIdx,dimIdx)=proba; }
	inline virtual double GetProbaDown(size_t stateIdx, size_t dimIdx) const {  return itsProbaDown(stateIdx,dimIdx); }
	void SetProbaDown(size_t stateIdx, size_t dimIdx, double proba) { itsProbaDown(stateIdx,dimIdx)=proba; }

    /// Generate states of the actual process sampled by the slice (Z states)
    inline virtual ARM_GP_MatrixPtr GetZStates();

    virtual void LinkedWithNextSlice( ARM_SliceBase* nextSlice, ARM_SamplerBase* sampler, ARM_TruncatorBase* truncator, bool needProbas, ARM_TransitionStates& transStates, const ARM_GP_VectorPtr& driftCorrection, size_t driftOffset, double targetPayoff);
};

inline ARM_GP_MatrixPtr ARM_SliceND::GetZStates()
{
    size_t nbStates = size();
    size_t nbDims = dim();
    ARM_GP_MatrixPtr zStates( new ARM_GP_Matrix(nbDims,nbStates,0.0) );

    ARM_TreeIndex index(GetMin(),GetMax());
    size_t stateIdx,i;
    for( index.Reset(); (stateIdx=index.GetPosition()) < nbStates; ++index )
    {
        for( i=0; i<nbDims; ++i )
            (*zStates)(i,stateIdx) = index[i] * GetSpaceStep(i);
    }

    return zStates;
}


///---------------------------
/// ARM_SliceNDCstSymProba
///---------------------------
class ARM_SliceNDCstSymProba : public ARM_SliceNDBase
{
private:
	std::vector<double> itsProba;             // cst proba to go up or down for the whole slice
//    ARM_BoolMatrix itsTruncatedHedge;   // to say if (state,dim) is truncated

public:
    ARM_SliceNDCstSymProba( size_t nbDims, size_t index, const ARM_IntVector& min,  const ARM_IntVector& max, const std::vector<double>& spaceStep, const std::vector<double>& proba)
        :   ARM_SliceNDBase(nbDims,index,min,max,spaceStep), itsProba( proba )/*, itsTruncatedHedge(0)*/
    {}
    ARM_SliceNDCstSymProba(const ARM_SliceNDCstSymProba& rhs)
        : ARM_SliceNDBase(rhs), itsProba(rhs.itsProba)/*, itsTruncatedHedge(rhs.itsTruncatedHedge)*/
    {}
    ARM_SliceNDCstSymProba& operator = (const ARM_SliceNDCstSymProba& rhs);

    virtual ~ARM_SliceNDCstSymProba() {}

    virtual ARM_Object* Clone() const { return new ARM_SliceNDCstSymProba(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const;

    virtual ARM_SliceNDCstSymProba* ToSliceNDCstSymProba() { return this; }
    virtual const ARM_SliceNDCstSymProba* ToSliceNDCstSymProba() const { return this; }

    inline virtual double GetProbaUp(size_t stateIdx, size_t dimIdx) const {  return /*itsTruncatedHedge(stateIdx,dimIdx) ? 0.0 : */itsProba[dimIdx]; }
	inline virtual double GetProbaDown(size_t stateIdx, size_t dimIdx) const {  return /*itsTruncatedHedge(stateIdx,dimIdx) ? 0.0 :*/ itsProba[dimIdx]; }

    /// Generate states of the actual process sampled by the slice (Z states)
    inline virtual ARM_GP_MatrixPtr GetZStates();

    /// Link with next slice is implicit but truncation may update slice range
    virtual void LinkedWithNextSlice( ARM_SliceBase* nextSlice, ARM_SamplerBase* sampler, ARM_TruncatorBase* truncator, bool needProbas, ARM_TransitionStates& transStates, const ARM_GP_VectorPtr& driftCorrection, size_t driftOffset, double targetPayoff);
};

inline ARM_GP_MatrixPtr ARM_SliceNDCstSymProba::GetZStates()
{
    size_t nbStates = size();
    size_t nbDims = dim();
    ARM_GP_MatrixPtr zStates( new ARM_GP_Matrix(nbDims,nbStates,0.0) );

    ResizeNextNodeIndex(nbStates,nbDims);

    ARM_TreeIndex index(GetMin(),GetMax());
    size_t stateIdx,i;
    for( index.Reset(); (stateIdx=index.GetPosition()) < nbStates; ++index )
    {
        for( i=0; i<nbDims; ++i )
        {
            (*zStates)(i,stateIdx) = index[i] * GetSpaceStep(i);

            /// Implicit straight node connection
            SetNextNodeIndex(stateIdx,i,index[i]);
        }
    }

    return zStates;
}


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

