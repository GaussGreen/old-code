/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file transposer.cpp
 *  \brief 
 *	\author  R. GUILLEMOT
 *	\version 1.0
 *	\date February 2005
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpnumlib/argconvdefault.h"
#include "gpnumlib/transposer.h"
#include "gpnumlib/random.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_Transposer
///	Routine: ARM_Transposer
///	Returns: void
///	Action : Constructor
////////////////////////////////////////////////////

ARM_Transposer::ARM_Transposer(
const ARM_RandomGeneratorPtr& randomGen,
TransposeOrder order,
int firstSimulations) :
ARM_ManipulatorMethod(randomGen),
itsOrder(order),
itsFirstSimulations(firstSimulations)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_Transposer
///	Routine: ARM_Transposer
///	Returns: void
///	Action : Copy Constructor
////////////////////////////////////////////////////

ARM_Transposer::ARM_Transposer(const ARM_Transposer& rhs)
: ARM_ManipulatorMethod(rhs),
itsOrder(rhs.itsOrder),
itsFirstSimulations(rhs.itsFirstSimulations)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_Transposer
///	Routine: operator=
///	Returns: void
///	Action : operator =
////////////////////////////////////////////////////

ARM_Transposer& ARM_Transposer::operator=(const ARM_Transposer& rhs )
{
	if( this != & rhs )
	{
		ARM_ManipulatorMethod::operator =( rhs );
		itsRandomGen	= ARM_RandomGeneratorPtr( (ARM_RandomGenerator*) rhs.itsRandomGen->Clone() ); /// clone it for independence!
		itsOrder = rhs.itsOrder;
		itsFirstSimulations = rhs.itsFirstSimulations;
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_Transposer
///	Routine: Clone()
///	Returns: void
///	Action : Create a copy of the object
////////////////////////////////////////////////////

ARM_ManipulatorMethod* ARM_Transposer::Clone() const
{
	return new ARM_Transposer(*this);
}

////////////////////////////////////////////////////
///	Class  : ARM_Transposer
///	Routine: ~ARM_Transposer
///	Returns: void
///	Action : Destructor
////////////////////////////////////////////////////

ARM_Transposer::~ARM_Transposer()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_Transposer
///	Routine: reset
///	Returns: void
///	Action : Reset the transoser
////////////////////////////////////////////////////

void ARM_Transposer::resetwork(const ARM_GP_T_Vector<size_t>& nbOfPointsList, const ARM_GP_T_Vector<size_t>& factorNb )
{
	itsRandomGen->reset(nbOfPointsList,factorNb);
	itsNbFactor = factorNb[0];
	itsDim = factorNb.size();
	
	bool checkCstNbFactor = true;
	for(int i = 1; i < factorNb.size(); i++)
	{
		if(factorNb[i] != factorNb[0])
		{
			checkCstNbFactor = false;
			break;
		}
	}

	// We simulate the first blank random numbers with fake bucket
	if (itsFirstSimulations)
	{
		itsNbOfPointsList.resize(nbOfPointsList.size()+1);
		itsNbOfPointsList[0] = itsFirstSimulations*factorNb[0];
		itsCurrentValues.resize(nbOfPointsList.size()+1);
	}
	else
	{
		itsNbOfPointsList = nbOfPointsList;
		itsCurrentValues.resize(nbOfPointsList.size());
	}

	for (itsCurrentBucketPos = (itsFirstSimulations?1:0); itsCurrentBucketPos < itsNbOfPointsList.size(); ++itsCurrentBucketPos)
		itsCurrentValues[itsCurrentBucketPos].resize(itsDim);

	// We compute all the random vector first
	for (itsCurrentBucketPos = 0; itsCurrentBucketPos < itsNbOfPointsList.size(); ++itsCurrentBucketPos)
	{
		for (itsCurrentDimPos = 0; itsCurrentDimPos < itsDim; ++itsCurrentDimPos)
		{
			itsCurrentValues[itsCurrentBucketPos][itsCurrentDimPos].resize(factorNb[itsCurrentDimPos]*itsNbOfPointsList[itsCurrentBucketPos]);	
		}
	}

	if (itsOrder == BucketOrder)
	{
		for (itsCurrentDimPos = 0; itsCurrentDimPos < itsDim; ++itsCurrentDimPos)
		{
			for (itsCurrentBucketPos = 0; itsCurrentBucketPos < itsNbOfPointsList.size(); ++itsCurrentBucketPos)
			{
				for (itsCurrentPos = 0; itsCurrentPos < factorNb[itsCurrentDimPos]*itsNbOfPointsList[itsCurrentBucketPos]; ++itsCurrentPos)
				{		
					itsCurrentValues[itsCurrentBucketPos][itsCurrentDimPos][itsCurrentPos] = (*itsRandomGen)();
				}
			}
		}
	}
	else if (itsOrder == PathOrder)
	{
		if(checkCstNbFactor)
		{
			for (itsCurrentBucketPos = 0; itsCurrentBucketPos < itsNbOfPointsList.size(); ++itsCurrentBucketPos)
			{
				for (itsCurrentPos = 0; itsCurrentPos < itsNbOfPointsList[itsCurrentBucketPos]*factorNb[0]; ++itsCurrentPos)
				{
					for (itsCurrentDimPos = 0; itsCurrentDimPos < itsDim; ++itsCurrentDimPos)
					{
						itsCurrentValues[itsCurrentBucketPos][itsCurrentDimPos][itsCurrentPos] = (*itsRandomGen)();
					}
				}
			}
		}
		else
		{
			ARM_IntVector idx(itsDim,0);
			
			for (itsCurrentBucketPos = 0; itsCurrentBucketPos < itsNbOfPointsList.size(); ++itsCurrentBucketPos)
			{
				for(int i = 0; i < itsDim; i++) idx[i] = 0;

				for (itsCurrentPos = 0; itsCurrentPos < itsNbOfPointsList[itsCurrentBucketPos]; ++itsCurrentPos)
				{
					for (itsCurrentDimPos = 0; itsCurrentDimPos < itsDim; ++itsCurrentDimPos)
					{
						for (int currentFactorDimPos = 0; currentFactorDimPos < factorNb[itsCurrentDimPos]; currentFactorDimPos++, idx[itsCurrentDimPos]++)
						{
							itsCurrentValues[itsCurrentBucketPos][itsCurrentDimPos][idx[itsCurrentDimPos]] = (*itsRandomGen)();
						}
					}
				}
			}
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
	itsCurrentDimPos = 0;
	itsCurrentPos = 0;
}

////////////////////////////////////////////////////
///	Class  : ARM_Transposer
///	Routine: operator()
///	Returns: void
///	Action : Reset the transoser
////////////////////////////////////////////////////0

double ARM_Transposer::operator()() const
{
	double value = 0.0;

	if( (itsCurrentBucketPos < itsNbOfPointsList.size()) && 
		(itsCurrentDimPos < itsDim) && 
		(itsCurrentPos >= itsCurrentValues[itsCurrentBucketPos][itsCurrentDimPos].size() ))
	{
		CC_MUTABLE( ARM_Transposer, itsCurrentPos ) = 0;
		CC_MUTABLE( ARM_Transposer, itsCurrentDimPos )++;
	}

	if( (itsCurrentBucketPos < itsNbOfPointsList.size()) && 
		(itsCurrentDimPos >= itsDim) )
	{
		CC_MUTABLE( ARM_Transposer, itsCurrentPos ) = 0;
		CC_MUTABLE( ARM_Transposer, itsCurrentDimPos ) = 0;
		CC_MUTABLE( ARM_Transposer, itsCurrentBucketPos )++;
	}

	if ((itsCurrentBucketPos < itsNbOfPointsList.size()) && 
		(itsCurrentDimPos < itsDim) && 
		(itsCurrentPos < itsCurrentValues[itsCurrentBucketPos][itsCurrentDimPos].size() ))
	{
		// We output the random number  dimension by dimension
		value = itsCurrentValues[itsCurrentBucketPos][itsCurrentDimPos][itsCurrentPos];
		CC_MUTABLE( ARM_Transposer, itsCurrentPos )++;
	}

	return value;
}

////////////////////////////////////////////////////
///	Class  : ARM_Transposer
///	Routine: toString()
///	Returns: void
///	Action : Display the contents
////////////////////////////////////////////////////

string ARM_Transposer::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << " Transposer\n";
	os << " Time Dimension: " << itsDim << " Factor Dimension: " << itsNbFactor << "\n";

	os << " PathOrder : " << ARM_ArgConvReverse_RandGenOrder.GetString(itsOrder) << "\n";
	os << itsRandomGen->toString();

	return os.str();
}

CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----