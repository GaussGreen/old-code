/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file meanrevertingsampler.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPNUMMETHODS_MEANREVERTINGSAMPLER_H
#define _INGPNUMMETHODS_MEANREVERTINGSAMPLER_H

#include "gpbase/port.h"
#include "sampler.h"


CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////
/// ARM_MeanRevertingSamplerBase
/// discretize for a mean reverting process7
/////////////////////////////////////////////

class ARM_MeanRevertingSamplerBase: public ARM_SamplerBase
{
private:
    double itsMinStdDev;

public:
    ARM_MeanRevertingSamplerBase( const ARM_SchedulerBase* scheduler,double minStdDev )
    :   ARM_SamplerBase( ARM_SchedulerBasePtr(scheduler ? static_cast<ARM_SchedulerBase*>(scheduler->Clone()) : NULL ) ),
        itsMinStdDev(minStdDev) {}

    ARM_MeanRevertingSamplerBase( const ARM_MeanRevertingSamplerBase& rhs );
    ARM_MeanRevertingSamplerBase& operator=(const ARM_MeanRevertingSamplerBase& rhs );
    virtual ~ARM_MeanRevertingSamplerBase();

	virtual ARM_Object* Clone() const = 0;

    double GetMinStdDev() const { return itsMinStdDev; }
	/// 1D equivalent for hybrid model
	virtual ARM_SamplerBase* CorrespondingSampler1D(size_t dim=0) const;
};


class ARM_MeanRevertingSampler1D : public ARM_MeanRevertingSamplerBase, public ARM_Sampler1DBase
{
private:
    /// Relative & absolute local drifts of X
    /// X(t+dt) = [RelDrift(t).X(t) + AbsDrift(t)].dt + (...).dW
    /// vector<> size is the tree schedule one
    ARM_GP_MatrixPtr itsRelDriftX;
    ARM_GP_MatrixPtr itsAbsDriftX;

    /// Local and global variances of X
    ARM_GP_Vector itsLocalVarX;
    ARM_GP_Vector itsGlobalVarX;


    void CopyNoCleanUp(const ARM_MeanRevertingSampler1D& rhs);

public:
    ARM_MeanRevertingSampler1D( const ARM_SchedulerBase* scheduler,double minStdDev=1.0e-3 )
    :   ARM_MeanRevertingSamplerBase(scheduler,minStdDev) {};
    ARM_MeanRevertingSampler1D( const ARM_MeanRevertingSampler1D& rhs );
    ARM_MeanRevertingSampler1D& operator=(const ARM_MeanRevertingSampler1D& rhs );
    virtual ~ARM_MeanRevertingSampler1D();

	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_MeanRevertingSampler1D"; }

    virtual ARM_TimeStepsAndSlices* Init(const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices=true);

    /// Accessors
    inline void SetRelDriftX(const ARM_GP_MatrixPtr& relDriftX) {itsRelDriftX=relDriftX;}
    inline const ARM_GP_MatrixPtr& GetRelDriftX() const { return itsRelDriftX;}
    inline void SetAbsDriftX(const ARM_GP_MatrixPtr& absDriftX) {itsAbsDriftX=absDriftX;}
    inline const ARM_GP_MatrixPtr& GetAbsDriftX() const { return itsAbsDriftX;}

    void SetLocalVarX(const ARM_MatrixVector& localVarX)
        {itsLocalVarX.resize(localVarX.size()); for(size_t i=0;i<localVarX.size();++i)  itsLocalVarX[i]=(*(localVarX[i]))(0,0);}

    void SetGlobalVarX(const ARM_MatrixVector& globalVarX)
        {itsGlobalVarX.resize(globalVarX.size()); for(size_t i=0;i<globalVarX.size();++i)  itsGlobalVarX[i]=(*(globalVarX[i]))(0,0);}

    /// sampler base 1D specific code
    virtual ARM_Sampler1DBase* ToSampler1DBase() { return this; }
    virtual const ARM_Sampler1DBase* ToSampler1DBase() const { return this; }

    virtual ARM_GP_VectorPtr ComputeConditionalMean(size_t sliceIdx, const ARM_GP_VectorPtr& zStates, ARM_GP_VectorPtr& xStates) const;
    virtual double GetLocalVar(size_t sliceIdx) const { return itsLocalVarX[sliceIdx]; }
    virtual double GetGlobalVar(size_t sliceIdx) const { return itsGlobalVarX[sliceIdx]; }
	virtual void ApplyRelDriftToIntXStates(size_t sliceIdx, ARM_GP_VectorPtr& integXStates) const;
	virtual void ApplyAbsDriftToXStates(size_t sliceIdx, ARM_GP_VectorPtr& XStates) const;
    virtual ARM_GP_VectorPtr ComputeZtoXStates(size_t sliceIdx, const ARM_GP_VectorPtr& zStates) const { return zStates; }

    /// X to Z space global converter : nothing because Z=X
    virtual void ComputeXtoZCenters(ARM_GP_MatrixPtr& globalDrifts) const {};

};


class ARM_DriftedMeanRevertingSampler1D : public ARM_MeanRevertingSampler1D
{
public:
    ARM_DriftedMeanRevertingSampler1D( const ARM_SchedulerBase* scheduler,double minStdDev=1.0e-3 )
    :   ARM_MeanRevertingSampler1D(scheduler,minStdDev) {};
    ARM_DriftedMeanRevertingSampler1D( const ARM_DriftedMeanRevertingSampler1D& rhs )
    :   ARM_MeanRevertingSampler1D(rhs) {};
    ARM_DriftedMeanRevertingSampler1D& operator=(const ARM_DriftedMeanRevertingSampler1D& rhs )
    { if( this != &rhs ) ARM_MeanRevertingSampler1D::operator =( rhs ); return *this; }
    virtual ~ARM_DriftedMeanRevertingSampler1D() {}

    virtual ARM_Object* Clone() const { return new ARM_DriftedMeanRevertingSampler1D(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_DriftedMeanRevertingSampler1D"; }

    /// sampler base 1D specific code
    virtual ARM_Sampler1DBase* ToSampler1DBase() { return this; }
    virtual const ARM_Sampler1DBase* ToSampler1DBase() const { return this; }

    virtual ARM_GP_VectorPtr ComputeConditionalMean(size_t sliceIdx, const ARM_GP_VectorPtr& zStates, ARM_GP_VectorPtr& xStates) const;
};


class ARM_MeanRevertingSamplerND: public ARM_MeanRevertingSamplerBase, public ARM_SamplerNDBase
{
private:
    /// Relative & absolute local drifts of X
    /// Xi(t+dt) = [RelDrift(i,t).Xi(t) + AbsDrift(i,t)].dt + (...).dW
    /// vector<> size is the tree schedule one
    ARM_GP_MatrixPtr itsRelDriftX;
    ARM_GP_MatrixPtr itsAbsDriftX;

    /// Local rotation matrixes from Z (independent) to or from X (correlated) spaces
    ARM_MatrixVector itsXtoZ;  // [Z] = [XtoZ][X]
    ARM_MatrixVector itsZtoX;  // [X] = [ZtoX][Z]

    ARM_MatrixVector itsGlobalXtoZ; // [Z] = [XtoZ][X] globally

    /// Variances in Z space (eigen values)
    ARM_VectorVector itsLocalVarZ;
    ARM_VectorVector itsGlobalVarZ;

	// Maximum number of factors really simulated by the Monte Carlo method
	int itsNbRank;

    void CopyNoCleanUp(const ARM_MeanRevertingSamplerND& rhs);
    void CleanUp();

public:
    ARM_MeanRevertingSamplerND( const ARM_SchedulerBase* scheduler,double minStdDev=1.0e-3, double nbRank=-1 )
    :   ARM_MeanRevertingSamplerBase(scheduler,minStdDev), itsNbRank(nbRank) {};
    ARM_MeanRevertingSamplerND( const ARM_MeanRevertingSamplerND& rhs );
    ARM_MeanRevertingSamplerND& operator=(const ARM_MeanRevertingSamplerND& rhs );
    virtual ~ARM_MeanRevertingSamplerND();

	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_MeanRevertingSamplerND"; }

    virtual ARM_TimeStepsAndSlices* Init(const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices=true);

    /// Accessors (ptr not cloned for efficiency)
    inline void SetRelDriftX(const ARM_GP_MatrixPtr& relDriftX) {itsRelDriftX=relDriftX;}
    inline const ARM_GP_MatrixPtr& GetRelDriftX() const { return itsRelDriftX;}
    inline void SetAbsDriftX(const ARM_GP_MatrixPtr& absDriftX) {itsAbsDriftX=absDriftX;}
    inline const ARM_GP_MatrixPtr& GetAbsDriftX() const { return itsAbsDriftX;}
    inline const ARM_MatrixVector& GetXtoZ() const { return itsXtoZ;}
    inline const ARM_MatrixVector& GetZtoX() const { return itsZtoX;}

    /// sampler base ND specific code
    virtual ARM_SamplerNDBase* ToSamplerNDBase() { return this; }
    virtual const ARM_SamplerNDBase* ToSamplerNDBase() const { return this; }

    virtual ARM_GP_MatrixPtr ComputeConditionalMean(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, const ARM_GP_VectorPtr& driftCorrection, ARM_GP_MatrixPtr& xStates) const;
    virtual const ARM_GP_Vector& GetLocalVar(size_t sliceIdx) const { return *(itsLocalVarZ[sliceIdx]); }
    virtual const ARM_GP_Vector& GetGlobalVar(size_t sliceIdx) const { return *(itsGlobalVarZ[sliceIdx]); }
	virtual const ARM_GP_Matrix& GetRelativeDrifts() const { return *itsRelDriftX; }

	virtual void ApplyRelDriftToIntXStates(size_t sliceIdx, ARM_GP_MatrixPtr& integXStates) const;
	virtual void ApplyAbsDriftToXStates(size_t sliceIdx, ARM_GP_MatrixPtr& XStates) const;
    virtual ARM_GP_MatrixPtr ComputeZtoXStates(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, ARM_GP_MatrixPtr& yStates = ARM_GP_MatrixPtr(NULL)) const;

    /// X to Z space global converter
    virtual void ComputeXtoZCenters(ARM_GP_MatrixPtr& globalDrifts) const;

	/// dimension
	virtual size_t dim() const { return itsZtoX[0]->rows(); }

	virtual int TotalDimension() const;
};

inline int ARM_MeanRevertingSamplerND::TotalDimension() const
{
	int dim = 0;
	for(int k = 0; k < itsLocalVarZ.size(); k++) 
		dim += itsLocalVarZ[k]->size();

	return dim;
}

class ARM_DriftedMeanRevertingSamplerND: public ARM_MeanRevertingSamplerND
{
public:
    ARM_DriftedMeanRevertingSamplerND( const ARM_SchedulerBase* scheduler,double minStdDev=1.0e-3 )
    :   ARM_MeanRevertingSamplerND(scheduler,minStdDev) {};
    ARM_DriftedMeanRevertingSamplerND( const ARM_DriftedMeanRevertingSamplerND& rhs )
    :   ARM_MeanRevertingSamplerND(rhs) {};
    ARM_DriftedMeanRevertingSamplerND& operator=(const ARM_DriftedMeanRevertingSamplerND& rhs )
    { if( this != &rhs ) ARM_MeanRevertingSamplerND::operator =( rhs ); return *this; }
    virtual ~ARM_DriftedMeanRevertingSamplerND() {};

    virtual ARM_Object* Clone() const { return new ARM_DriftedMeanRevertingSamplerND(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_DriftedMeanRevertingSamplerND"; }

    /// sampler base ND specific code
    virtual ARM_SamplerNDBase* ToSamplerNDBase() { return this; }
    virtual const ARM_SamplerNDBase* ToSamplerNDBase() const { return this; }

    virtual ARM_GP_MatrixPtr ComputeConditionalMean(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, const ARM_GP_VectorPtr& driftCorrection, ARM_GP_MatrixPtr& xStates) const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

