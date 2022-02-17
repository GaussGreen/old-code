/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file markoviandriftsampler.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPNUMMETHODS_MARKOVIANDRIFTSAMPLER_H
#define _INGPNUMMETHODS_MARKOVIANDRIFTSAMPLER_H

#include "gpbase/port.h"
#include "sampler.h"

CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////
/// ARM_MarkovianDriftSamplerBase
/// discretize for a mean reverting process
/////////////////////////////////////////////

class ARM_MarkovianDriftSamplerBase: public ARM_SamplerBase
{
private:
    double itsMinStdDev;

public:
    ARM_MarkovianDriftSamplerBase( const ARM_SchedulerBase* scheduler,double minStdDev )
    :   ARM_SamplerBase( ARM_SchedulerBasePtr( scheduler ? static_cast<ARM_SchedulerBase*>(scheduler->Clone()) : NULL ) ),
        itsMinStdDev(minStdDev) {}

    ARM_MarkovianDriftSamplerBase( const ARM_MarkovianDriftSamplerBase& rhs );
    ARM_MarkovianDriftSamplerBase& operator=(const ARM_MarkovianDriftSamplerBase& rhs );
    virtual ~ARM_MarkovianDriftSamplerBase();

	virtual ARM_Object* Clone() const = 0;
    double GetMinStdDev() const { return itsMinStdDev; }

    virtual bool IsIntegratedSampling() const { return false; }

	/// 1D equivalent for hybrid model
	virtual ARM_SamplerBase* CorrespondingSampler1D(size_t dim=0) const;
};


class ARM_MarkovianDriftSampler1D: public ARM_MarkovianDriftSamplerBase, public ARM_Sampler1DBase
{
private:
    /// Instantaneous volatilities of initial process X
    /// sampled w.r.t. tree schedule (for convertion to or from constant
    /// volatility process Z)
    ARM_GP_VectorPtr itsVolX;
    double itsVolZ;

    double itsVolZ2; // for local/global variance purpose

    /// Relative drift of Z variable (it includes volatility change MRS)
    std::vector<double> itsRelDriftZ;

    /// Absolute drift of X
    std::vector<double> itsAbsDriftX;

    void CopyNoCleanUp(const ARM_MarkovianDriftSampler1D& rhs);

public:
    ARM_MarkovianDriftSampler1D( const ARM_SchedulerBase* scheduler,double minStdDev=1.0e-3 )
        : ARM_MarkovianDriftSamplerBase(scheduler,minStdDev) {}
    ARM_MarkovianDriftSampler1D( const ARM_MarkovianDriftSampler1D& rhs );
    ARM_MarkovianDriftSampler1D& operator=(const ARM_MarkovianDriftSampler1D& rhs );
    virtual ~ARM_MarkovianDriftSampler1D();

	virtual ARM_Object* Clone() const { return new ARM_MarkovianDriftSampler1D(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_MarkovianDriftSampler1D"; }

    virtual ARM_TimeStepsAndSlices* Init(const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices=true);

    /// sampler base 1D specific code
    virtual ARM_Sampler1DBase* ToSampler1DBase() { return this; }
    virtual const ARM_Sampler1DBase* ToSampler1DBase() const { return this; }

    virtual ARM_GP_VectorPtr ComputeConditionalMean(size_t sliceIdx, const ARM_GP_VectorPtr& zStates, ARM_GP_VectorPtr& xStates) const;
    virtual double GetLocalVar(size_t sliceIdx) const;
    virtual double GetGlobalVar(size_t sliceIdx) const;

	virtual void ApplyRelDriftToIntXStates(size_t sliceIdx, ARM_GP_VectorPtr& integXStates) const {};
	virtual void ApplyAbsDriftToXStates(size_t sliceIdx, ARM_GP_VectorPtr& XStates) const {};
    virtual ARM_GP_VectorPtr ComputeZtoXStates(size_t sliceIdx, const ARM_GP_VectorPtr& zStates) const;

    /// X to Z space global converter
    virtual void ComputeXtoZCenters(ARM_GP_MatrixPtr& globalDrifts) const;
};


class ARM_MarkovianDriftSamplerNDBase : public ARM_MarkovianDriftSamplerBase, public ARM_SamplerNDBase
{
private:
    void CopyNoCleanUp(const ARM_MarkovianDriftSamplerNDBase& rhs);
    void CleanUp();

protected:
    /// Instantaneous volatilities of initial processes Xi
    /// sampled w.r.t. tree schedule (for convertion to or from constant
    /// volatility processes Yi)
    ARM_GP_MatrixPtr itsVolX;  // VolX[factor,schedule]
    std::vector<double> itsVolY;     // VolY[factor]

    /// Rotation matrixes for conversion to or from
    /// independant brownians Zi.
    /// May change each time correlations change
    /// vector<> size is equal to the tree schedule one for fast access
    /// but vectors/matrixes are shared through ptrs
    vector< ARM_GP_MatrixPtr > itsYtoZ; // [Z] = [YtoZ][Y]
    vector< ARM_GP_MatrixPtr > itsZtoY; // [Y] = [ZtoY][Z]

    vector< ARM_GP_MatrixPtr > itsGlobalYtoZ; // [Y] = [ZtoY][Z] globally

    /// Variances in Z space (eigen values)
    ARM_VectorVector itsLocalVarZ;
    ARM_VectorVector itsGlobalVarZ;

    /// Relative drift of Y
    ARM_VectorVector itsRelDriftY;

    /// Absolute drift of X
    ARM_VectorVector itsAbsDriftX;

	// Just used for importance sampling
	ARM_GP_Matrix itsFakeRelativeDrift;

public :
    ARM_MarkovianDriftSamplerNDBase( const ARM_SchedulerBase* scheduler,double minStdDev );
    ARM_MarkovianDriftSamplerNDBase( const ARM_MarkovianDriftSamplerNDBase& rhs );
    ARM_MarkovianDriftSamplerNDBase& operator=(const ARM_MarkovianDriftSamplerNDBase& rhs );
    virtual ~ARM_MarkovianDriftSamplerNDBase();

    virtual ARM_Object* Clone() const = 0;

    /// sampler base ND specific code
    virtual ARM_SamplerNDBase* ToSamplerNDBase() { return this; }
    virtual const ARM_SamplerNDBase* ToSamplerNDBase() const { return this; }

	virtual void ApplyRelDriftToIntXStates(size_t sliceIdx, ARM_GP_MatrixPtr& integXStates) const {};
	virtual void ApplyAbsDriftToXStates(size_t sliceIdx, ARM_GP_MatrixPtr& XStates) const {};

	virtual const ARM_GP_Matrix& GetRelativeDrifts() const { return itsFakeRelativeDrift; }

	virtual int TotalDimension() const;
};

inline int ARM_MarkovianDriftSamplerNDBase::TotalDimension() const
{
	int dim = 0;
	for(int k = 0; k < itsLocalVarZ.size(); k++) 
		dim += itsLocalVarZ[k]->size();

	return dim;
}

class ARM_MarkovianDriftSampler2D : public ARM_MarkovianDriftSamplerNDBase
{
public:
    ARM_MarkovianDriftSampler2D( const ARM_SchedulerBase* scheduler,double minStdDev=1.0e-3 )
        : ARM_MarkovianDriftSamplerNDBase(scheduler,minStdDev) {}
    ARM_MarkovianDriftSampler2D( const ARM_MarkovianDriftSampler2D& rhs )
        : ARM_MarkovianDriftSamplerNDBase(rhs) {}
    ARM_MarkovianDriftSampler2D& operator=(const ARM_MarkovianDriftSampler2D& rhs );
    virtual ~ARM_MarkovianDriftSampler2D() {}

    virtual ARM_Object* Clone() const { return new ARM_MarkovianDriftSampler2D(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_MarkovianDriftSampler2D"; }

    virtual ARM_TimeStepsAndSlices* Init(const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices=true);

    virtual ARM_GP_MatrixPtr ComputeConditionalMean(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, const ARM_GP_VectorPtr& driftCorrection, ARM_GP_MatrixPtr& xStates) const;
	virtual const std::vector<double>& GetLocalVar(size_t sliceIdx) const { return itsLocalVarZ[sliceIdx]->GetValues(); }
    virtual const std::vector<double>& GetGlobalVar(size_t sliceIdx) const { return itsGlobalVarZ[sliceIdx]->GetValues(); }
    virtual ARM_GP_MatrixPtr ComputeZtoXStates(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, ARM_GP_MatrixPtr& yStates= ARM_GP_MatrixPtr(NULL)) const;

    /// X to Z space global converter
    virtual void ComputeXtoZCenters(ARM_GP_MatrixPtr& globalDrifts) const;

	/// dimension
	virtual size_t dim() const { return 2; }
};


class ARM_MarkovianDriftSampler3D : public ARM_MarkovianDriftSamplerNDBase
{
public:
    ARM_MarkovianDriftSampler3D( const ARM_SchedulerBase* scheduler,double minStdDev=1.0e-3 )
        : ARM_MarkovianDriftSamplerNDBase(scheduler,minStdDev) {}
    ARM_MarkovianDriftSampler3D( const ARM_MarkovianDriftSampler3D& rhs )
        : ARM_MarkovianDriftSamplerNDBase(rhs) {}
    ARM_MarkovianDriftSampler3D& operator=(const ARM_MarkovianDriftSampler3D& rhs );
    virtual ~ARM_MarkovianDriftSampler3D() {}

	virtual ARM_Object* Clone() const { return new ARM_MarkovianDriftSampler3D(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_MarkovianDriftSampler3D"; }

    virtual ARM_TimeStepsAndSlices* Init(const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices=true);

    virtual ARM_GP_MatrixPtr ComputeConditionalMean(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, const ARM_GP_VectorPtr& driftCorrection, ARM_GP_MatrixPtr& xStates) const;
	virtual const std::vector<double>& GetLocalVar(size_t sliceIdx) const { return itsLocalVarZ[sliceIdx]->GetValues(); }
    virtual const std::vector<double>& GetGlobalVar(size_t sliceIdx) const { return itsGlobalVarZ[sliceIdx]->GetValues(); }
    virtual ARM_GP_MatrixPtr ComputeZtoXStates(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, ARM_GP_MatrixPtr& yStates = ARM_GP_MatrixPtr(NULL)) const;

    /// X to Z space global converter
    virtual void ComputeXtoZCenters(ARM_GP_MatrixPtr& globalDrifts) const;

	/// dimension
	virtual size_t dim() const { return 3; }
};


class ARM_MarkovianDriftSamplerND : public ARM_MarkovianDriftSamplerNDBase
{
private:
	size_t itsDim;

public:
    ARM_MarkovianDriftSamplerND( const ARM_SchedulerBase* scheduler,double minStdDev=1.0e-3)
        : ARM_MarkovianDriftSamplerNDBase(scheduler,minStdDev), itsDim( 0 )
	{}
    ARM_MarkovianDriftSamplerND( const ARM_MarkovianDriftSamplerND& rhs )
        : ARM_MarkovianDriftSamplerNDBase(rhs), itsDim(rhs.itsDim) {}
    ARM_MarkovianDriftSamplerND& operator=(const ARM_MarkovianDriftSamplerND& rhs );
    virtual ~ARM_MarkovianDriftSamplerND() {}

	virtual ARM_Object* Clone() const { return new ARM_MarkovianDriftSamplerND(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_MarkovianDriftSamplerND"; }

    virtual ARM_TimeStepsAndSlices* Init(const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices=true);

    virtual ARM_GP_MatrixPtr ComputeConditionalMean(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, const ARM_GP_VectorPtr& driftCorrection, ARM_GP_MatrixPtr& xStates) const;
	virtual const std::vector<double>& GetLocalVar(size_t sliceIdx) const { return (*itsLocalVarZ[sliceIdx]).GetValues(); }
    virtual const std::vector<double>& GetGlobalVar(size_t sliceIdx) const { return (*itsGlobalVarZ[sliceIdx]).GetValues(); }
    virtual ARM_GP_MatrixPtr ComputeZtoXStates(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, ARM_GP_MatrixPtr& yStates = ARM_GP_MatrixPtr(NULL)) const;

    /// X to Z space global converter
    virtual void ComputeXtoZCenters(ARM_GP_MatrixPtr& globalDrifts) const;

	/// dimension
	virtual size_t dim() const { return itsDim; }
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

