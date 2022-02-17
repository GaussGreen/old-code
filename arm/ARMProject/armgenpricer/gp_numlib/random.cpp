/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file random.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */


#include "gpnumlib/random.h"

/// gpbase
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/gpmatrix.h"

/// gpnumlib
#include "gpnumlib/stddevfunc.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_RandomGenerator
///	Routine: operator()
///	Returns: vector of random numbers
////////////////////////////////////////////////////
std::vector<double>& ARM_RandomGenerator::operator()( size_t Number )
{
	std::vector<double>* result = new std::vector<double>(Number);
	draw(*result);
	return *result;
};



////////////////////////////////////////////////////
///	Class  : ARM_RandomGenerator
///	Routine: daw(ARM_RC_Vector& Vec)
///	Action : draw a vector of random numbers
////////////////////////////////////////////////////
void ARM_RandomGenerator::draw( std::vector<double>& Vec)
{
	std::vector<double>::iterator 
		iter 	= Vec.begin(),
		End		= Vec.end();

#if defined(__GP_STRICT_VALIDATION)
	ValidateDim( Vec.size() );
#endif

	for( ; iter!=End; ++iter)
	{
#if defined(__GP_STRICT_VALIDATION)
		*iter = DrawOne();
		ValidateResult(*iter);
#else
		*iter = DrawOne();
#endif
	}
}

///////////////////////////////////////////////////
///	Class  : ARM_RandomGenerator
///	Routine: daw(ARM_RC_Vector& Vec)
///	Action : draw a vector of random numbers
////////////////////////////////////////////////////
void ARM_RandomGenerator::draw( ARM_GP_Matrix& Matrix)
{
	size_t row = Matrix.GetRowsNb();
	size_t col = Matrix.GetColsNb();
	double elem;

#if defined(__GP_STRICT_VALIDATION)
	ValidateDim( Matrix.size() );
#endif
	
	for( size_t j=0; j<col; ++j)	
	{
		for( size_t i=0; i<row; ++i)
		{
#if defined(__GP_STRICT_VALIDATION)
			elem = DrawOne(row);
			ValidateResult(elem);			
#else
			elem = DrawOne(row);
#endif
			Matrix.Elt(i,j) = elem;
		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_RandomGenerator
///	Routine: desctructor (necessary even for pure virtual case)
///	Action : draw a vector of random numbers
////////////////////////////////////////////////////
ARM_RandomGenerator::~ARM_RandomGenerator()
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandomGenerator
///	Routine: constructor, copy constructor, assignment operator
///	Returns: 
///	Action :
/////////////////////////////////////////////////////////////////

ARM_RandomGenerator::ARM_RandomGenerator()
:	ARM_RootObject(), itsCurrentValue(0)
{
	CC_ARM_SETNAME(ARM_RANDOMGEN);

	IsThisQuasiRandom = false;
}


ARM_RandomGenerator::ARM_RandomGenerator( const ARM_RandomGenerator& rhs )
:	ARM_RootObject( rhs ), itsCurrentValue( rhs.itsCurrentValue ), IsThisQuasiRandom(rhs.IsThisQuasiRandom)
{}


ARM_RandomGenerator& ARM_RandomGenerator::operator=(const ARM_RandomGenerator& rhs )
{
	if( this != & rhs )
	{
		ARM_RootObject::operator =( rhs );
		itsCurrentValue = rhs.itsCurrentValue;
		IsThisQuasiRandom = rhs.IsThisQuasiRandom;
	}
	return *this;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_RandomGenerator
///	Routine: ARM_ComputeStdCovFuncPtr
///	Returns: ARM_MomentFuncPtr
///	Action : functor to compute the std dev
/////////////////////////////////////////////////////////////////
ARM_MomentFuncPtr ARM_RandomGenerator::StdDevComputationFunc() const
{ 
// FIXMEFRED: mig.vc8 (22/05/2007 18:14:15):cast
	return static_cast<ARM_MomentFuncPtr>(new ARM_StdMomentFunc);
}

void ARM_RandomGenerator::reset(const ARM_GP_T_Vector<size_t>& NbOfPointsList, const ARM_GP_T_Vector<size_t>& factorNb)
{
	int totaldim = 0;
	for(int k = 0; k < factorNb.size(); k++) totaldim += factorNb[k];
	if(totaldim == 0) return;
	reset(totaldim, NbOfPointsList, factorNb.size() == 0 ? 0 : factorNb[0]);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
