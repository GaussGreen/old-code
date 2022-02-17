/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file mcmethod.h
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date October 2005
 */


#ifndef _INGPINFRA_IMPSAMPLER_H
#define _INGPINFRA_IMPSAMPLER_H

#include "gpbase/rootobject.h"
#include "gpbase/typedef.h"
#include "gpbase/curve.h"
#include "gpbase/curvetypedef.h"

#include <map>
CC_USING_NS( std, map )

CC_BEGIN_NAMESPACE( ARM )

// This class is the general interface to define importance sampling

struct ARM_SamplerBase;

class ARM_ImpSampler : public ARM_RootObject
{
public:
	enum ARM_ImpSamplerType { DummyImpSampler, PropImpSampler };
	typedef map<double,ARM_GP_VectorPtr> ExpMap;

	ARM_ImpSampler() : ARM_RootObject() {};
	ARM_ImpSampler(const ARM_ImpSampler& rhs) : ARM_RootObject(rhs) {};
	~ARM_ImpSampler() {};


	// Compute the drift to add at the process
	virtual ARM_VectorPtrVector ComputeDrifts(
		int factorNb,
		const ARM_GP_Vector& timeSteps,
		ARM_SamplerBase* sampler) = 0;

	// Compute the exponential to multiply with the payoff
	virtual void ComputeExps(
		int nbStates,
		int factorNb,
		const ARM_GP_Vector& timeSteps,
		const ARM_MatrixPtrVector& processStates,
		ARM_SamplerBase* sampler) = 0;

	// Those functions are used to apply the exponential importance sampling
	virtual void ProcessPaidPayoffs(ARM_VectorPtr& payoffs, double evalTime) const = 0;
	virtual void ProcessUnPaidPayoffs(ARM_VectorPtr& payoffs, double evalTime) const = 0;


};

// This class defines the "Dummy" importance Sampler
// It does nothing

class ARM_DummyImpSampler : public ARM_ImpSampler
{
public:
	ARM_DummyImpSampler();

	ARM_DummyImpSampler(const ARM_DummyImpSampler& rhs);
	ARM_DummyImpSampler& operator=(const ARM_DummyImpSampler& rhs);
	~ARM_DummyImpSampler() {};

	virtual ARM_Object* Clone() const { return new ARM_DummyImpSampler(*this); };

	// Compute the drift to add at the process
	virtual ARM_VectorPtrVector ComputeDrifts(
		int factorNb,
		const ARM_GP_Vector& timeSteps,
		ARM_SamplerBase* sampler);

	// Compute the exponential to multiply with the payoff
	// It does nothing
	virtual void ComputeExps(
		int nbStates,
		int factorNb,
		const ARM_GP_Vector& timeSteps,
		const ARM_MatrixPtrVector& processStates,
		ARM_SamplerBase* sampler) {};

	// Those functions are used to apply the exponential importance sampling
	// It does nothing in the dummy case
	virtual void ProcessPaidPayoffs(ARM_VectorPtr& payoffs, double evalTime) const {};
	virtual void ProcessUnPaidPayoffs(ARM_VectorPtr& payoffs, double evalTime) const {};

	/// the two methods to reimplement
    virtual string toString(const string& indent="", const string& nextIndent="") const;

private:
	ARM_GP_MultiCurvePtr itsAlpha;
};

// This class defines the "Proportional" importance Sampler
// It uses the simulated variable to carry out the variance
// reduction

class ARM_PropImpSampler : public ARM_ImpSampler
{
private:
	ARM_ImpSampler::ExpMap itsExpProcessStates;
public:
	ARM_PropImpSampler(const ARM_MultiCurve* alpha);

	ARM_PropImpSampler(const ARM_PropImpSampler& rhs);
	ARM_PropImpSampler& operator=(const ARM_PropImpSampler& rhs);
	~ARM_PropImpSampler() {};

	virtual ARM_Object* Clone() const { return new ARM_PropImpSampler(*this); };

	//Accessor
	ARM_GP_MultiCurvePtr GetAlpha() const { return itsAlpha; };

	// Compute the drift to add at the process
	virtual ARM_VectorPtrVector ComputeDrifts(
		int factorNb,
		const ARM_GP_Vector& timeSteps,
		ARM_SamplerBase* sampler);

	// Compute the exponential to multiply with the payoff
	virtual void ComputeExps(
		int nbStates,
		int factorNb,
		const ARM_GP_Vector& timeSteps,
		const ARM_MatrixPtrVector& processStates,
		ARM_SamplerBase* sampler);

	virtual void BootstrapAlpha(
		int factorNb,
	const ARM_GP_Vector& timeSteps,
	const ARM_GP_Vector& rowSteps,
	const ARM_GP_T_Vector<ARM_GP_Vector>& globalAlpha,
	ARM_SamplerBase* sampler);

	// Those functions are used to apply the exponential importance sampling
	// It does nothing in the dummy case
	virtual void ProcessPaidPayoffs(ARM_VectorPtr& payoffs, double evalTime) const;
	virtual void ProcessUnPaidPayoffs(ARM_VectorPtr& payoffs, double evalTime) const;

	/// the two methods to reimplement
    virtual string toString(const string& indent="", const string& nextIndent="") const;

	//Accessor
	void SetAlpha(const ARM_GP_MultiCurvePtr& alpha) { itsAlpha = alpha; }

private:
	ARM_GP_MultiCurvePtr itsAlpha;
};

CC_END_NAMESPACE()

#endif