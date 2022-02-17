/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pathscheme.cpp
 *
 *  \brief
 *
*	\author  R. Guillemot
 *	\version 1.0
 *	\date November 2005
 */


#include "gpnummethods/pathscheme.h"
#include "gpnummethods/sampler.h"

#include "gpnumlib/random.h"

#include <math.h>


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_IncrementalPathScheme
///	Routine: ARM_IncrementalPathScheme
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////

ARM_IncrementalPathScheme::ARM_IncrementalPathScheme(const ARM_IncrementalPathScheme& rhs) 
: itsWeights(rhs.itsWeights) 
{
}


////////////////////////////////////////////////////
///	Class  : ARM_IncrementalPathScheme
///	Routine: ~ARM_IncrementalPathScheme
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////

ARM_IncrementalPathScheme::~ARM_IncrementalPathScheme() 
{
	DeletePointorVector<ARM_GP_Vector>(itsWeights);
}



////////////////////////////////////////////////////
///	Class  : ARM_IncrementalPathScheme
///	Routine: Init
///	Returns: 
///	Action : Initialize the Path Scheme
////////////////////////////////////////////////////

void ARM_IncrementalPathScheme::Init(int nbTimeSteps, const ARM_SamplerBase* sampler)
{
	// There is just a N dimensional 

	const ARM_SamplerNDBase* samplerND = sampler->ToSamplerNDBase();
	
	size_t nbFactors = samplerND->dim();

	itsWeights.resize(nbTimeSteps);

	for (size_t timeIdx = 0; timeIdx < nbTimeSteps-1; ++timeIdx)
	{
		size_t idxSize = samplerND->GetLocalVar(timeIdx).size();
		
		itsWeights[timeIdx] = new ARM_GP_Vector(idxSize);

		for (size_t factorIdx =0; factorIdx < idxSize; ++factorIdx)
			(*itsWeights[timeIdx])[factorIdx] = sqrt(samplerND->GetLocalVar(timeIdx)[factorIdx]);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_IncrementalPathScheme
///	Routine: ComputeProcessStates
///	Returns: 
///	Action : Compute the process states with the scheme
////////////////////////////////////////////////////

void ARM_IncrementalPathScheme::ComputeProcessStates(
	size_t bucketSize,
	const ARM_RandomGeneratorPtr& randGen,
	ARM_MatrixPtrVector& processStates)  const
{
	int nbTimeSteps = itsWeights.size();

	processStates.resize(nbTimeSteps);

	// There are one more time steps than induct !
	for (size_t timeIdx = 0; timeIdx < nbTimeSteps-1; ++timeIdx)
	{
		int nbFactors = itsWeights[timeIdx]->size();

		ARM_GP_MatrixPtr gaussian( new ARM_GP_Matrix(nbFactors,bucketSize) );

		randGen->draw(*gaussian);
		processStates[timeIdx] = ARM_GP_MatrixPtr(new ARM_GP_Matrix(nbFactors, bucketSize));
			
		for (size_t stateIdx = 0; stateIdx < bucketSize; stateIdx++)
		{
			// We multiply the gaussian with the local volatility of each factor
			for (size_t factorIdx = 0; factorIdx < nbFactors; ++factorIdx)
				(*processStates[timeIdx])(factorIdx,stateIdx) = (*itsWeights[timeIdx])[factorIdx]*(*gaussian)(factorIdx,stateIdx);
		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_IncrementalPathScheme
///	Routine: toString
///	Returns:
///	Action : Display the contents of the object
////////////////////////////////////////////////////
string ARM_IncrementalPathScheme::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;

	os << "Incremental";

	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_BrownianBridgePathScheme
///	Routine: ARM_BrownianBridgePathScheme
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////

ARM_BrownianBridgePathScheme::ARM_BrownianBridgePathScheme(const ARM_BrownianBridgePathScheme& rhs) 
: itsBridgeIndex(rhs.itsBridgeIndex),
itsLeftIndex(rhs.itsLeftIndex),
itsRightIndex(rhs.itsRightIndex),
itsLeftWeights(rhs.itsLeftWeights),
itsRightWeights(rhs.itsRightWeights)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_BrownianBridgePathScheme
///	Routine: Init
///	Returns: 
///	Action : Compute the process states with the scheme
////////////////////////////////////////////////////

void ARM_BrownianBridgePathScheme::Init(int nbTimeSteps, const ARM_SamplerBase* sampler)
{
	unsigned long i,j,k,l,m;

	const ARM_SamplerNDBase* samplerND = sampler->ToSamplerNDBase();
	size_t nbFactors = samplerND->dim();

	ARM_GP_Matrix globalVars(nbTimeSteps, nbFactors);

	for (i = 0; i < nbTimeSteps; ++i)
		for (j = 0; j < nbFactors; ++j)
		{
			if (i > 0)
				globalVars(i,j) = globalVars(i-1,j) + samplerND->GetLocalVar(i-1)[j];
			else
				globalVars(i,j) = 0.0;
		}

	// This code has been adapted from the Monte Carlo Methods in Finance
	// By Peter Jackel - Wiley Finance

#if defined(__GP_STRICT_VALIDATION)
	if( nbTimeSteps <= 0 )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": the number of time steps has to be strictly positive.");
#endif
	vector <unsigned long> map(nbTimeSteps);

	itsBridgeIndex.resize(nbTimeSteps-1);
	itsLeftIndex.resize(nbTimeSteps-1);
	itsRightIndex.resize(nbTimeSteps-1);
	itsLeftWeights.resize(nbTimeSteps-1,nbFactors);
	itsRightWeights.resize(nbTimeSteps-1,nbFactors);
	itsSigmas.resize(nbTimeSteps-1, nbFactors);
	// map is used to indicate which points are already constructed. If map[i] is zero, path point i
	// is yet unconstructed.  map[i]-1 is the index of the variate that constructs the path point # i.
	map[nbTimeSteps-1] = 1;									//  The first point in the construction is the global step.
	itsBridgeIndex[0] = nbTimeSteps-1;						//  The global step is constructed from the first variate.
	for (m = 0; m < nbFactors; ++m)
	{
		itsSigmas(0,m) = sqrt(globalVars(nbTimeSteps-1,m));	//  The variance of the global step is numberOfSteps*1.0.
		itsLeftWeights(0,m) = itsRightWeights(0,m) = 0.;	//  The global step to the last point in time is special.
	}
	for (j=1,i=1;i<nbTimeSteps-1;++i){
		while (map[j]) ++j;									//  Find the next unpopulated entry in the map.
		k=j;
		while ((!map[k])) ++k;								//  Find the next populated entry in the map from there.
		l=j+((k-1-j)>>1);									//  l-1 is now the index of the point to be constructed next.
		map[l]=i;
		itsBridgeIndex[i] = l;								//  The i-th Gaussian variate will be used to set point l-1.
		itsLeftIndex[i]   = j;
		itsRightIndex[i]  = k;
		for (m = 0; m < nbFactors; ++m)
		{
			if ( fabs(globalVars(k,m)-globalVars(j-1,m)) > K_NEW_DOUBLE_TOL)
			{
				itsLeftWeights(i,m)  = (globalVars(k,m)-globalVars(l,m))/(globalVars(k,m)-globalVars(j-1,m));
				itsRightWeights(i,m) = (globalVars(l,m)-globalVars(j-1,m))/(globalVars(k,m)-globalVars(j-1,m));
				itsSigmas(i,m) = sqrt((globalVars(k,m)-globalVars(l,m))*(globalVars(l,m)-globalVars(j-1,m))/(globalVars(k,m)-globalVars(j-1,m)));
			}
			else
			{
				itsLeftWeights(i,m)  = 1.0;
				itsRightWeights(i,m) = 0.0;
				itsSigmas(i,m) = 0.0;
			}
		}
		j=k+1;
		if (j>=nbTimeSteps) j=1;							//	Wrap around.
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_BrownianBridgePathScheme
///	Routine: ComputeProcessStates
///	Returns: 
///	Action : Compute the process states with the scheme
////////////////////////////////////////////////////

void ARM_BrownianBridgePathScheme::ComputeProcessStates(
	size_t bucketSize,
	const ARM_RandomGeneratorPtr& randGen,
	ARM_MatrixPtrVector& processStates) const
{
	int nbTimeSteps = itsSigmas.rows();
	int nbFactors = itsSigmas.cols();

	processStates.resize(nbTimeSteps);

	unsigned long i,j,k,l,m,n;

	ARM_GP_MatrixPtr gaussian( new ARM_GP_Matrix(nbFactors,bucketSize) );

	randGen->draw(*gaussian);
	processStates[0] = ARM_GP_MatrixPtr(new ARM_GP_Matrix(nbFactors, bucketSize));
	processStates[nbTimeSteps-1] = ARM_GP_MatrixPtr(new ARM_GP_Matrix(nbFactors, bucketSize));
	for (m =0; m < nbFactors; ++m)
		for (n = 0; n < bucketSize; ++n)
		{
			(*processStates[nbTimeSteps-1])(m,n) = itsSigmas(0,m)*(*gaussian)(m,n);       //  The global step.
		}

	for ( i = 1; i < nbTimeSteps; ++i)
	{
		randGen->draw(*gaussian);
		j = itsLeftIndex[i];
		k = itsRightIndex[i];
		l = itsBridgeIndex[i];
		processStates[l-1] = ARM_GP_MatrixPtr(new ARM_GP_Matrix(nbFactors, bucketSize));
		// Apply the brownian bridge for the simulation
		for (m =0; m < nbFactors; ++m)
			for (n = 0; n < bucketSize; ++n)
				if (j-1) (*processStates[l-1])(m,n) = itsLeftWeights(i,m)*(*processStates[j-2])(m,n) + itsRightWeights(i,m)*(*processStates[k-1])(m,n) + itsSigmas(i,m)*(*gaussian)(m,n);
				else   (*processStates[l-1])(m,n) = itsRightWeights(i,m)*(*processStates[k-1])(m,n) + itsSigmas(i,m)*(*gaussian)(m,n);
	}

	// We just calculate the process increment
	for (i = nbTimeSteps-1; i > 0; --i)
		for (m=0; m < nbFactors; ++m)
			for (n = 0; n < bucketSize; ++n)
				(*processStates[i])(m,n) -= (*processStates[i-1])(m,n);
}

////////////////////////////////////////////////////
///	Class  : ARM_BrownianBridgePathScheme
///	Routine: toString
///	Returns:
///	Action : Display the contents of the object
////////////////////////////////////////////////////
string ARM_BrownianBridgePathScheme::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;

	os << "Brownian Bridge\n";

	int nbTimeSteps = itsBridgeIndex.size();
	int nbFactors = itsSigmas.cols();

	int i,j;

	int sliceCount = 0;

	if (nbTimeSteps)
	{
		os << "Slice " << sliceCount << ": ";

		for (i = 0; i < nbTimeSteps; ++i)
		{
			os << itsBridgeIndex[i];
			if (i < nbTimeSteps-1)
			{
				if (itsBridgeIndex[i+1] < itsBridgeIndex[i])
				{
					sliceCount++;
					os << "\nSlice " << sliceCount << ": ";
				}
				else
					os << " ";
			}
			else
			{
				os << "\n\n";
			}
		}
		
		os << std::left;

		os << std::setw(15) << "Bridge Index" << "\t";
		os << std::setw(15) << "Left Index" << "\t";
		os << std::setw(15) << "RightIndex" << "\t";

		for (i = 0; i < nbFactors; ++i)
		{
			CC_Ostringstream oslw, osrw, ossgm;
			oslw << "Left Weight " << i;
			osrw << "Right Weight " << i;
			ossgm << "Sigma " << i;;

			os << std::setw(15) << oslw.str() << "\t";
			os << std::setw(15) << osrw.str() << "\t";
			os << std::setw(15) << ossgm.str();

			if (i < nbFactors-1)
				os << "\t";
		}

		os << "\n";

		for (i = 0; i < nbTimeSteps; ++i)
		{
			os << std::fixed << std::setprecision(0) << std::setw(15) << itsBridgeIndex[i] << "\t";
			os << std::fixed << std::setprecision(0) << std::setw(15) << (itsLeftIndex[i]?itsLeftIndex[i]-1:0) << "\t";
			os << std::fixed << std::setprecision(0) <<  std::setw(15) << itsRightIndex[i] << "\t";

			for (j = 0; j < nbFactors; ++j)
			{
				CC_Ostringstream oslw, osrw, ossgm;
				oslw << std::fixed << std::setprecision(2) << itsLeftWeights(i,j)*100  << "%";
				osrw << std::fixed << std::setprecision(2) << itsRightWeights(i,j)*100  << "%";
				ossgm << std::fixed << std::setprecision(2) << itsSigmas(i,j)*100  << "%";

				os  << std::setw(15) <<  oslw.str() << "\t";
				os  << std::setw(15) <<  osrw.str() << "\t";
				os  << std::setw(15) <<  ossgm.str() << "\t";

				if (j < nbFactors-1)
					os << "\t";
			}

			os << "\n";
		}
	}

	return os.str();
}

//////////////////////////////////////////////////////

ARM_IncrementAdaptativePathScheme::ARM_IncrementAdaptativePathScheme(const ARM_IncrementAdaptativePathScheme& rhs) 
: itsWeights(rhs.itsWeights) 
{
}


////////////////////////////////////////////////////
///	Class  : ARM_IncrementAdaptativePathScheme
///	Routine: ~ARM_IncrementAdaptativePathScheme
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////

ARM_IncrementAdaptativePathScheme::~ARM_IncrementAdaptativePathScheme() 
{
	DeletePointorVector<ARM_GP_Vector>(itsWeights);
}


////////////////////////////////////////////////////
///	Class  : ARM_IncrementalPathScheme
///	Routine: Init
///	Returns: 
///	Action : Initialize the Path Scheme
////////////////////////////////////////////////////

void ARM_IncrementAdaptativePathScheme::Init(int nbTimeSteps, const ARM_SamplerBase* sampler)
{
	// There is just a N dimensional 

	const ARM_SamplerNDBase* samplerND = sampler->ToSamplerNDBase();
	
	size_t nbFactors = samplerND->dim();

	itsWeights.resize(nbTimeSteps);

	for (size_t timeIdx = 0; timeIdx < nbTimeSteps-1; ++timeIdx)
	{
		size_t idxSize = samplerND->GetLocalVar(timeIdx).size();
		
		itsWeights[timeIdx] = new ARM_GP_Vector(idxSize);

		for (size_t factorIdx =0; factorIdx < idxSize; ++factorIdx)
			(*itsWeights[timeIdx])[factorIdx] = sqrt(samplerND->GetLocalVar(timeIdx)[factorIdx]);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_IncrementalPathScheme
///	Routine: ComputeProcessStates
///	Returns: 
///	Action : Compute the process states with the scheme
////////////////////////////////////////////////////

void ARM_IncrementAdaptativePathScheme::ComputeProcessStates(
	size_t bucketSize,
	const ARM_RandomGeneratorPtr& randGen,
	ARM_MatrixPtrVector& processStates)  const
{
	int nbTimeSteps = itsWeights.size();

	processStates.resize(nbTimeSteps);

	// There are one more time steps than induct !
	for (size_t timeIdx = 0; timeIdx < nbTimeSteps-1; ++timeIdx)
	{
		int nbFactors = itsWeights[timeIdx]->size();

		ARM_GP_MatrixPtr gaussian( new ARM_GP_Matrix(nbFactors,bucketSize) );

		randGen->draw(*gaussian);
		processStates[timeIdx] = ARM_GP_MatrixPtr(new ARM_GP_Matrix(nbFactors, bucketSize));
			
		for (size_t stateIdx = 0; stateIdx < bucketSize; stateIdx++)
		{
			// We multiply the gaussian with the local volatility of each factor
			for (size_t factorIdx = 0; factorIdx < nbFactors; ++factorIdx)
				(*processStates[timeIdx])(factorIdx,stateIdx) = (*itsWeights[timeIdx])[factorIdx]*(*gaussian)(factorIdx,stateIdx);
		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_IncrementalPathScheme
///	Routine: toString
///	Returns:
///	Action : Display the contents of the object
////////////////////////////////////////////////////
string ARM_IncrementAdaptativePathScheme::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;

	os << "Incremental and Adaptative";

	return os.str();
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/