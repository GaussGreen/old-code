/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file antitheticgen.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */


#include "gpnumlib/antitheticgen.h"
#include "gpbase/ostringstream.h"
#include "gpbase/singleton.h"
#include "gpbase/eventviewer.h"
#include "gpbase/gpvector.h"

#include "gpnumlib/quasirandom.h"
#include "gpnumlib/stddevfunc.h"

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///			 ARM_AntitheticOneGen			 ///////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_AntitheticOneGen
///	Routine: constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_AntitheticOneGen::ARM_AntitheticOneGen( const ARM_RandomGeneratorPtr& randomGen )
:	ARM_RandomGenerator(), 
	itsRandGen( randomGen), 
	itsCurrentIndex(2), 
	itsPreviousValues( new std::vector<double>(1,0.0) ),
	itsDim(-1),
	itsFactorDim(1)
{
	itsIsQuasiRandom = (dynamic_cast<ARM_QuasiRandom*>(&*randomGen)!=NULL);
}


////////////////////////////////////////////////////
///	Class  : ARM_AntitheticOneGen  
///	Routine: reset
///	Returns: void
///	Action : reset the random number generator
////////////////////////////////////////////////////
void ARM_AntitheticOneGen::reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorsNb )
{
	ARM_GP_T_Vector<size_t> nbOfPointsListCopy(nbOfPointsList);

	for(size_t i =0; i < nbOfPointsList.size(); ++i)
		nbOfPointsListCopy[i] /= 2;

	/// to make sure we have positive nbOfPoints;
	itsRandGen->reset( dim, nbOfPointsListCopy, factorsNb );
	ARM_TheEventViewer.Instance()->AddToMessage("Antithetic end reset 1");
	itsCurrentIndex = -1; 
	
	itsDim				= dim;
	ARM_TheEventViewer.Instance()->AddToMessage("Before delete");
	delete itsPreviousValues;
	ARM_TheEventViewer.Instance()->AddToMessage("After delete");
	itsPreviousValues	= new std::vector<double>(factorsNb,0.0);
	ARM_TheEventViewer.Instance()->AddToMessage("After new");
	itsFactorDim		= factorsNb;
	itsCurrentIndex		= 2*itsFactorDim; 

	ARM_TheEventViewer.Instance()->AddToMessage("Antithetic end reset 2");
}


////////////////////////////////////////////////////
///	Class  : ARM_AntitheticOneGen
///	Routine: constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
double ARM_AntitheticOneGen::DrawOne()
{
	if( ++itsCurrentIndex < itsFactorDim )
		return (*itsPreviousValues)[itsCurrentIndex];
	else
	{
		if( itsCurrentIndex < 2*itsFactorDim )
		{
			if (itsIsQuasiRandom)
				return (*itsPreviousValues)[itsCurrentIndex-itsFactorDim];
			else
				return -(*itsPreviousValues)[itsCurrentIndex-itsFactorDim];
		}
		else
		{
			itsCurrentIndex = 0;
			for(size_t i=0; i<itsFactorDim; ++i )
				(*itsPreviousValues)[i]=itsRandGen->DrawOne();
			return (*itsPreviousValues)[0];
		}
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_AntitheticOneGen
///	Routine: toString
///	Returns: 
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_AntitheticOneGen::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << "Antithetic One Generator " << CC_NS(std,endl);
	os << itsRandGen->toString();
	return os.str();
}




/////////////////////////////////////////////////////////////////
///	Class  : ARM_AntitheticOneGen
///	Routine: copy constructor, assignment operator
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_AntitheticOneGen::ARM_AntitheticOneGen( const ARM_AntitheticOneGen& rhs )
:	ARM_RandomGenerator( rhs ), 
	itsRandGen( (ARM_RandomGenerator*) rhs.itsRandGen->Clone() ),
	itsCurrentIndex(rhs.itsCurrentIndex ), 
	itsPreviousValues( new std::vector<double>(*rhs.itsPreviousValues) ),
	itsDim( rhs.itsDim ),
	itsFactorDim( rhs.itsFactorDim ),
	itsIsQuasiRandom(rhs.itsIsQuasiRandom)
{}


ARM_AntitheticOneGen& ARM_AntitheticOneGen::operator =(const ARM_AntitheticOneGen& rhs )
{
	if( this != & rhs )
	{
		ARM_RandomGenerator::operator =( rhs );
		itsRandGen			= ARM_RandomGeneratorPtr( (ARM_RandomGenerator*) rhs.itsRandGen->Clone() );
		itsCurrentIndex		= rhs.itsCurrentIndex;
		itsPreviousValues	= new std::vector<double>( *rhs.itsPreviousValues);
		itsDim				= rhs.itsDim;
		itsFactorDim		= rhs.itsFactorDim;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_AntitheticOneGen
///	Routine: desctructor (necessary even for pure virtual case)
///	Action : draw a vector of uniform random numbers
////////////////////////////////////////////////////
ARM_AntitheticOneGen::~ARM_AntitheticOneGen()
{
	delete itsPreviousValues;
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_AntitheticOneGen
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_AntitheticOneGen::Clone() const
{
	return new ARM_AntitheticOneGen(*this);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_AntitheticOneGen
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
void ARM_AntitheticOneGen::ValidateResult( double result )
{
#if defined(__GP_STRICT_VALIDATION)
	itsRandGen->ValidateResult( result );
#endif
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_AntitheticOneGen
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
void ARM_AntitheticOneGen::ValidateDim( size_t dim )
{
#if defined(__GP_STRICT_VALIDATION)
	itsRandGen->ValidateDim( dim);
#endif
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_AntitheticOneGen
///	Routine: StdDevComputationFunc
///	Returns: ARM_MeanAndStdDevFuncPtr
///	Action : functor to compute the std dev in antithetic scheme
/////////////////////////////////////////////////////////////////
ARM_MomentFuncPtr ARM_AntitheticOneGen::StdDevComputationFunc() const
{ 
// FIXMEFRED: mig.vc8 (22/05/2007 18:14:47):cast
	return static_cast<ARM_MomentFuncPtr>(new ARM_AntitheticMomentFunc); 
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_AntitheticOneGen
///	Routine: GetDistributionType
///	Returns: ARM_DistributionType 
///	Action : the distribution type
/////////////////////////////////////////////////////////////////
ARM_RandomGenerator::ARM_DistributionType ARM_AntitheticOneGen::GetDistributionType() const
{
	return itsRandGen->GetDistributionType();
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

