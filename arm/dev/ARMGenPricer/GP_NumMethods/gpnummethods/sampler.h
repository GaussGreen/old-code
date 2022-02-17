/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file sampler.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPNUMMETHODS_SAMPLER_H
#define _INGPNUMMETHODS_SAMPLER_H

#include "gpbase/port.h"

#include "gpbase/rootobject.h"
#include "gpinfra/nummethod.h"

#include "typedef.h"
#include "scheduler.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_PricingModel;


//// structure returned by a process Sampler
struct ARM_TimeStepsAndSlices
{
    ARM_GP_VectorPtr        itsTimeSteps;
    ARM_SliceVectorPtr      itsSlices;
};

/////////////////////////////////////////////
/// code specific to 1D
/////////////////////////////////////////////
struct ARM_Sampler1DBase
{
    virtual double GetLocalVar(size_t sliceIdx) const = 0;
    virtual double GetGlobalVar(size_t sliceIdx) const = 0;
    virtual ARM_GP_VectorPtr ComputeConditionalMean(size_t sliceIdx, const ARM_GP_VectorPtr& zStates, ARM_GP_VectorPtr& xStates) const = 0;
	virtual void ApplyRelDriftToIntXStates(size_t sliceIdx, ARM_GP_VectorPtr& integXStates) const = 0;
	virtual void ApplyAbsDriftToXStates(size_t sliceIdx, ARM_GP_VectorPtr& XStates) const = 0;
    virtual ARM_GP_VectorPtr ComputeZtoXStates(size_t sliceIdx, const ARM_GP_VectorPtr& zStates) const = 0;
};


/////////////////////////////////////////////
/// code specific to ND
/////////////////////////////////////////////
struct ARM_SamplerNDBase
{
	ARM_SamplerNDBase(){};

	ARM_SamplerNDBase(const ARM_SamplerNDBase& rhs ) : itsNbFactors(rhs.itsNbFactors)
	{
	};

	virtual size_t dim() const = 0;
    virtual const ARM_GP_Vector& GetLocalVar(size_t sliceIdx) const = 0;
    virtual const ARM_GP_Vector& GetGlobalVar(size_t sliceIdx) const = 0;
    virtual const ARM_GP_Matrix& GetRelativeDrifts() const = 0;
	virtual ARM_GP_MatrixPtr ComputeConditionalMean(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, const ARM_GP_VectorPtr& driftCorrection, ARM_GP_MatrixPtr& xStates) const = 0;
	virtual void ApplyRelDriftToIntXStates(size_t sliceIdx, ARM_GP_MatrixPtr& integXStates) const = 0;
	virtual void ApplyAbsDriftToXStates(size_t sliceIdx, ARM_GP_MatrixPtr& XStates) const = 0;
    virtual ARM_GP_MatrixPtr ComputeZtoXStates(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, ARM_GP_MatrixPtr& yStates = ARM_GP_MatrixPtr(NULL)) const = 0;

	virtual int TotalDimension() const = 0;

	const ARM_GP_T_Vector<size_t>&	nbFactors() const {return itsNbFactors;};
protected:
	ARM_GP_T_Vector<size_t> itsNbFactors;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

