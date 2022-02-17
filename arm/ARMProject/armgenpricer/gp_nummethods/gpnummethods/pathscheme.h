/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pathscheme.h
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date November 2005
 */


#ifndef _INGPNUMMETHOD_PATHSCHEME_H
#define _INGPNUMMETHOD_PATHSCHEME_H

#include "gpbase/rootobject.h"
#include "gpbase/typedef.h"
#include "gpbase/assignop.h"
#include "gpbase/gpmatrix.h"

#include "gpnumlib/typedef.h"

#include <map>
CC_USING_NS( std, map )

CC_BEGIN_NAMESPACE( ARM )

// This class is the general interface to define path schemes

struct ARM_SamplerBase;

class ARM_PathScheme : public ARM_RootObject
{
public:
	enum ARM_PathSchemeType { Incremental, BrownianBridge, IncAdaptative };

	ARM_PathScheme() : ARM_RootObject() {};
	ARM_PathScheme(const ARM_PathScheme& rhs) : ARM_RootObject(rhs) {};
	~ARM_PathScheme() {};

	// Initialize the the Path Scheme
	virtual void Init(int nbTimeSteps, const ARM_SamplerBase* sampler) = 0;

	// Compute the process state based on the scheme
	virtual void ComputeProcessStates(
		size_t bucketSize,
		const ARM_RandomGeneratorPtr& randGen,
		ARM_MatrixPtrVector& processSates)  const = 0;
};

// This class defines the Incremental Path Scheme

class ARM_IncrementalPathScheme : public ARM_PathScheme
{
public:
	ARM_IncrementalPathScheme() {};
	ARM_IncrementalPathScheme(const ARM_IncrementalPathScheme& rhs);
	ASSIGN_OPERATOR(ARM_IncrementalPathScheme)
	~ARM_IncrementalPathScheme();

	virtual ARM_Object* Clone() const { return new ARM_IncrementalPathScheme(*this); };

	// Initialize the path scheme
	virtual void Init(int nbTimeSteps, const ARM_SamplerBase* sampler);

	// Compute the process state based on the scheme
	virtual void ComputeProcessStates(
		size_t bucketSize,
		const ARM_RandomGeneratorPtr& randGen,
		ARM_MatrixPtrVector& processSates)  const;

	// Display the content of the object
    virtual string toString(const string& indent="", const string& nextIndent="") const;
private:
	ARM_VectorVector itsWeights;
};

// This class defines the Brownian Bridge Path Scheme

class ARM_BrownianBridgePathScheme : public ARM_PathScheme
{
public:
	ARM_BrownianBridgePathScheme() {};
	ARM_BrownianBridgePathScheme(const ARM_BrownianBridgePathScheme& rhs);
	ASSIGN_OPERATOR(ARM_BrownianBridgePathScheme)
	~ARM_BrownianBridgePathScheme() {};

	virtual ARM_Object* Clone() const { return new ARM_BrownianBridgePathScheme(*this); };

	// Initialize the path scheme
	virtual void Init(int nbTimeSteps, const ARM_SamplerBase* sampler);

	// Compute the process state based on the scheme
	virtual void ComputeProcessStates(
		size_t bucketSize,
		const ARM_RandomGeneratorPtr& randGen,
		ARM_MatrixPtrVector& processStates)  const;

	// Display the content of the object
    virtual string toString(const string& indent="", const string& nextIndent="") const;

private:
	std::vector<double> itsBridgeIndex;
	std::vector<double> itsLeftIndex;
	std::vector<double> itsRightIndex;
	ARM_GP_Matrix itsLeftWeights;
	ARM_GP_Matrix itsRightWeights;
	ARM_GP_Matrix itsSigmas;
};

class ARM_IncrementAdaptativePathScheme : public ARM_PathScheme
{
public:
	ARM_IncrementAdaptativePathScheme() {};
	ARM_IncrementAdaptativePathScheme(const ARM_IncrementAdaptativePathScheme& rhs);
	ASSIGN_OPERATOR(ARM_IncrementAdaptativePathScheme)
	~ARM_IncrementAdaptativePathScheme();

	virtual ARM_Object* Clone() const { return new ARM_IncrementAdaptativePathScheme(*this); };

	// Initialize the path scheme
	virtual void Init(int nbTimeSteps, const ARM_SamplerBase* sampler);

	// Compute the process state based on the scheme
	virtual void ComputeProcessStates(
		size_t bucketSize,
		const ARM_RandomGeneratorPtr& randGen,
		ARM_MatrixPtrVector& processSates)  const;

	// Display the content of the object
    virtual string toString(const string& indent="", const string& nextIndent="") const;
private:
	ARM_VectorVector itsWeights;
};

CC_END_NAMESPACE()

#endif