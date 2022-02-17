/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file quasirandom.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */


#include "gpbase/removeidentifiedwarning.h"
#include "gpnumlib/quasirandom.h"
#include "gpbase/ostringstream.h"
#include "gpbase/ostringstream.h"
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"

#include <glob/expt.h>

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_QuasiRandom
///	Routine: constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_QuasiRandom::ARM_QuasiRandom(int firstSimulations)
:
itsDim(1),
itsFactorNb(1),
itsNbOfPointsList(0),
itsCurrentBucketPos(0),
itsCurrentPos(0),
itsCurrentDimPos(0),
itsCurrentTimeIndex(0),
itsCurrentFactorPos(0),
itsCurrentValues(0),
itsFirstSimulations(firstSimulations)
{
	IsThisQuasiRandom = true;
	itsStepFactorNb.resize(1,1);
}



////////////////////////////////////////////////////
///	Class  : ARM_QuasiRandom
///	Routine: ValidateDim (validation routine for dimension)
///	Returns: 
///	Action : throw an exception if the dim is too low
////////////////////////////////////////////////////
void ARM_QuasiRandom::ValidateDim( size_t vecDim )
{
	if( vecDim < dim()) 
	{
		CC_Ostringstream os;
		os  << "quasi random sequence with dim " << dim()
			<< " but used to draw!";
	    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str());
	}
};

////////////////////////////////////////////////////
///	Class  : ARM_QuasiRandom
///	Routine: drawOne
///	Returns: 
///	Action : function to return the current double random number
////////////////////////////////////////////////////
double ARM_QuasiRandom::DrawOne()
{
	double value = 0.0;
	int nbfactor = itsStepFactorNb.size() == 0 ? itsFactorNb : itsStepFactorNb[itsCurrentTimeIndex];

	if( (itsCurrentBucketPos < itsNbOfPointsList.size()) && 
		(itsCurrentTimeIndex < itsStepFactorNb.size()) && 
		(itsCurrentPos < itsNbOfPointsList[itsCurrentBucketPos]) &&
		(itsCurrentFactorPos >= nbfactor))
	{
		itsCurrentFactorPos = 0,
		itsCurrentPos++;
	}

	if( (itsCurrentBucketPos < itsNbOfPointsList.size()) && 
		(itsCurrentTimeIndex < itsStepFactorNb.size()) && 
		(itsCurrentPos >= itsNbOfPointsList[itsCurrentBucketPos] ))
	{
		itsCurrentFactorPos = 0,
		itsCurrentPos = 0;
		itsCurrentTimeIndex+=1;
		itsCurrentDimPos+=nbfactor;
	}

	if( (itsCurrentBucketPos < itsNbOfPointsList.size()) && 
		(itsCurrentTimeIndex >= itsStepFactorNb.size()) )
	{
		itsCurrentFactorPos = 0,
		itsCurrentPos = 0;
		itsCurrentDimPos = 0;
		itsCurrentTimeIndex = 0;
		itsCurrentBucketPos++;
	}

	if ((itsCurrentBucketPos < itsNbOfPointsList.size()) && 
		(itsCurrentTimeIndex < itsStepFactorNb.size()) && 
		(itsCurrentPos < itsNbOfPointsList[itsCurrentBucketPos]) &&
		(itsCurrentFactorPos < nbfactor))
	{
		// We output the random number  dimension by dimension
		value = itsCurrentValues[itsCurrentBucketPos][itsCurrentPos][itsCurrentDimPos+itsCurrentFactorPos];
		itsCurrentFactorPos++;
	}

	if (itsCurrentBucketPos >= itsNbOfPointsList.size())
	{
		CC_Ostringstream os;
		os  << "We try to use too much bucket: bucket pos = " << itsCurrentBucketPos << " >= " << itsNbOfPointsList.size();
	    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str());
	}

	return value;
}


////////////////////////////////////////////////////
///	Class  : ARM_QuasiRandom
///	Routine: reset
///	Returns: 
///	Action : reset the generator
////////////////////////////////////////////////////
void ARM_QuasiRandom::reset( size_t dim, const ARM_GP_T_Vector<size_t>& NbOfPointsList, size_t factorsNb )
{
	ARM_GP_T_Vector<size_t> factors(1,factorsNb);
	reset(NbOfPointsList,factors);
	return;

	itsNbOfPointsList = NbOfPointsList;
	itsDim = dim;
	itsFactorNb = factorsNb;
	SetDim( dim );

	// We simulate the first blank random numbers with fake bucket
	if (itsFirstSimulations)
	{
		itsNbOfPointsList.resize(NbOfPointsList.size()+1);
		itsNbOfPointsList[0] = itsFirstSimulations;
		for (size_t i = 1; i < NbOfPointsList.size()+1; ++i)
			itsNbOfPointsList[i] = NbOfPointsList[i-1];
		itsCurrentValues.resize(NbOfPointsList.size()+1);
	}
	else
	{
		itsNbOfPointsList = NbOfPointsList;
		itsCurrentValues.resize(NbOfPointsList.size());
	}

	// We compute all the random vector first
	itsCurrentValues.resize(itsNbOfPointsList.size());
	for (itsCurrentBucketPos = 0; itsCurrentBucketPos < itsNbOfPointsList.size(); ++itsCurrentBucketPos)
	{
		itsCurrentValues[itsCurrentBucketPos].resize(itsNbOfPointsList[itsCurrentBucketPos]);
		for (itsCurrentPos = 0; itsCurrentPos < itsNbOfPointsList[itsCurrentBucketPos]; ++itsCurrentPos)
		{
			itsCurrentValues[itsCurrentBucketPos][itsCurrentPos].resize(dim);
			DrawAll();
		}
	}

	// We remove the first simulations
	if (itsFirstSimulations)
	{
		itsNbOfPointsList.erase(itsNbOfPointsList.begin());
		itsCurrentValues.erase(itsCurrentValues.begin());
	}

	// We initialize the iterators 
	itsCurrentBucketPos = 0;
	itsCurrentPos = 0;
	itsCurrentDimPos = 0;
	itsCurrentFactorPos = 0;
}

void ARM_QuasiRandom::reset(const ARM_GP_T_Vector<size_t>& NbOfPointsList,const ARM_GP_T_Vector<size_t>& factorNb)
{
	itsNbOfPointsList = NbOfPointsList;
	itsDim = 0;
	for(int k = 0; k < factorNb.size(); k++) itsDim += factorNb[k];
	SetDim(itsDim);
	if(factorNb.size() == 0) return;
	itsFactorNb = factorNb[0];
	itsStepFactorNb = factorNb;

	if (itsFirstSimulations)
	{
		itsNbOfPointsList.resize(NbOfPointsList.size()+1);
		itsNbOfPointsList[0] = itsFirstSimulations;
		for (size_t i = 1; i < NbOfPointsList.size()+1; ++i)
			itsNbOfPointsList[i] = NbOfPointsList[i-1];
		itsCurrentValues.resize(NbOfPointsList.size()+1);
	}
	else
	{
		itsNbOfPointsList = NbOfPointsList;
		itsCurrentValues.resize(NbOfPointsList.size());
	}

	// We compute all the random vector first
	itsCurrentValues.resize(itsNbOfPointsList.size());
	for (itsCurrentBucketPos = 0; itsCurrentBucketPos < itsNbOfPointsList.size(); ++itsCurrentBucketPos)
	{
		itsCurrentValues[itsCurrentBucketPos].resize(itsNbOfPointsList[itsCurrentBucketPos]);
		for (itsCurrentPos = 0; itsCurrentPos < itsNbOfPointsList[itsCurrentBucketPos]; ++itsCurrentPos)
		{
			itsCurrentValues[itsCurrentBucketPos][itsCurrentPos].resize(itsDim);
			DrawAll();
		}
	}

	// We remove the first simulations
	if (itsFirstSimulations)
	{
		itsNbOfPointsList.erase(itsNbOfPointsList.begin());
		itsCurrentValues.erase(itsCurrentValues.begin());
	}

	// We initialize the iterators 
	itsCurrentBucketPos = 0;
	itsCurrentPos = 0;
	itsCurrentDimPos = 0;
	itsCurrentTimeIndex = 0;
	itsCurrentFactorPos = 0;
}

CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----