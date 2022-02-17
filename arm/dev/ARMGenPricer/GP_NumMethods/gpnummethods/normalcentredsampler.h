/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file normalcentredsampler.h
 *
 *  \brief 
 *
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPNUMMETHODS_NORMALCENTREDSAMPLER_H
#define _INGPNUMMETHODS_NORMALCENTREDSAMPLER_H


#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "sampler.h"

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////
/// ARM_NormalCentredSamplerBase
/// discretize with a constant variance
/////////////////////////////////////////////

class ARM_NormalCentredSamplerBase: public ARM_SamplerBase
{
public:
    ARM_NormalCentredSamplerBase( const ARM_SchedulerBase* scheduler )
        :    ARM_SamplerBase( ARM_SchedulerBasePtr( scheduler ? static_cast<ARM_SchedulerBase*>(scheduler->Clone()) : NULL ) )
        {}
    ARM_NormalCentredSamplerBase( const ARM_NormalCentredSamplerBase& rhs )
        : ARM_SamplerBase( rhs ) {}

    virtual ~ARM_NormalCentredSamplerBase() {};

	/// 1D equivalent for hybrid model
	virtual ARM_SamplerBase* CorrespondingSampler1D(size_t dim=0) const;

    /// X to Z space global converter : nothing is done because always centred
    virtual void ComputeXtoZCenters(ARM_GP_MatrixPtr& globalDrifts) const {};

};



/// 1d case
class ARM_NormalCentredSampler1D: public ARM_NormalCentredSamplerBase, public ARM_Sampler1DBase
{
private :
    /// Local and global variances of X
    ARM_GP_Vector itsLocalVarX;
    ARM_GP_Vector itsGlobalVarX;

    void CopyNoCleanUp(const ARM_NormalCentredSampler1D& rhs);

public:
    ARM_NormalCentredSampler1D( const ARM_SchedulerBase* scheduler )
        :    ARM_NormalCentredSamplerBase( scheduler ) {}
    ARM_NormalCentredSampler1D( const ARM_NormalCentredSampler1D& rhs );
	ASSIGN_OPERATOR(ARM_NormalCentredSampler1D)
    virtual ~ARM_NormalCentredSampler1D();

    /// accessor
    void SetLocalVarX(const ARM_MatrixVector& localVarX)
        {itsLocalVarX.resize(localVarX.size()); for(size_t i=0;i<localVarX.size();++i)  itsLocalVarX[i]=(*(localVarX[i]))(0,0);}

    void SetGlobalVarX(const ARM_MatrixVector& globalVarX)
        {itsGlobalVarX.resize(globalVarX.size()); for(size_t i=0;i<globalVarX.size();++i)  itsGlobalVarX[i]=(*(globalVarX[i]))(0,0);}

    virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_NormalCentredSampler1D"; }

    virtual ARM_TimeStepsAndSlices* Init(const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices=true);

    /// sampler base 1D specific code
    virtual ARM_Sampler1DBase* ToSampler1DBase() { return this; }
    virtual const ARM_Sampler1DBase* ToSampler1DBase() const { return this; }

    virtual ARM_GP_VectorPtr ComputeConditionalMean(size_t sliceIdx, const ARM_GP_VectorPtr& zStates, ARM_GP_VectorPtr& xStates) const
    { xStates = ComputeZtoXStates(sliceIdx,zStates); return ARM_GP_VectorPtr(xStates); }

    virtual double GetLocalVar(size_t sliceIdx) const { return itsLocalVarX[sliceIdx]; }
    virtual double GetGlobalVar(size_t sliceIdx) const { return itsGlobalVarX[sliceIdx]; }

	virtual void ApplyRelDriftToIntXStates(size_t sliceIdx, ARM_GP_VectorPtr& integXStates) const {};
	virtual void ApplyAbsDriftToXStates(size_t sliceIdx, ARM_GP_VectorPtr& XStates) const {};
    virtual ARM_GP_VectorPtr ComputeZtoXStates(size_t sliceIdx, const ARM_GP_VectorPtr& zStates) const { return zStates; }
};



/// multi-dim case
class ARM_NormalCentredSamplerND: public ARM_NormalCentredSamplerBase, public ARM_SamplerNDBase
{
private:
    /// Rotation matrixes from Z (independent) to/from X (correlated) spaces
    ARM_MatrixVector itsZtoX;
    ARM_MatrixVector itsXtoZ;

    /// Variances in Z space (eigen values)
    ARM_VectorVector itsLocalVarZ;
    ARM_VectorVector itsGlobalVarZ;

	// Maximum number of factors really simulated by the Monte Carlo method
	int itsNbRank;

	// Just used for importance sampling
	ARM_GP_Matrix itsFakeRelativeDrift;

    void CopyNoCleanUp(const ARM_NormalCentredSamplerND& rhs);
    void CleanUp();

public:
    ARM_NormalCentredSamplerND( const ARM_SchedulerBase* scheduler, size_t nbRank = -1)
        :    ARM_NormalCentredSamplerBase( scheduler ), itsNbRank(nbRank) {}
	ARM_NormalCentredSamplerND(const ARM_NormalCentredSamplerND& rhs );
    ASSIGN_OPERATOR(ARM_NormalCentredSamplerND)
    virtual ~ARM_NormalCentredSamplerND();

	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_NormalCentredSamplerND"; }

    virtual ARM_TimeStepsAndSlices* Init(const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices=true);

    /// sampler base ND specific code
    virtual ARM_SamplerNDBase* ToSamplerNDBase() { return this; }
    virtual const ARM_SamplerNDBase* ToSamplerNDBase() const { return this; }

    inline virtual ARM_GP_MatrixPtr ComputeConditionalMean(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, const ARM_GP_VectorPtr& driftCorrection, ARM_GP_MatrixPtr& xStates) const
    { xStates = ComputeZtoXStates(sliceIdx,zStates); return ARM_GP_MatrixPtr(xStates); }

    virtual const ARM_GP_Vector& GetLocalVar(size_t sliceIdx) const { return *(itsLocalVarZ[sliceIdx]); }
    virtual const ARM_GP_Vector& GetGlobalVar(size_t sliceIdx) const { return *(itsGlobalVarZ[sliceIdx]); }
	virtual const ARM_GP_Matrix& GetRelativeDrifts() const { return itsFakeRelativeDrift; }
	virtual void ApplyRelDriftToIntXStates(size_t sliceIdx, ARM_GP_MatrixPtr& integXStates) const {};
	virtual void ApplyAbsDriftToXStates(size_t sliceIdx, ARM_GP_MatrixPtr& XStates) const {};
    virtual ARM_GP_MatrixPtr ComputeZtoXStates(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, ARM_GP_MatrixPtr& yStates = ARM_GP_MatrixPtr(NULL)) const;

	/// dimension
	virtual size_t dim() const { return itsZtoX[0]->rows(); }

	virtual int TotalDimension() const;
};

inline int ARM_NormalCentredSamplerND::TotalDimension() const
{
	int dim = 0;
	for(int k = 0; k < itsLocalVarZ.size(); k++) 
		dim += itsLocalVarZ[k]->size();

	return dim;
}

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

