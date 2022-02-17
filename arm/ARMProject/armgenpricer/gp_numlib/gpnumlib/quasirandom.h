/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file quasirandom.h
 *
 *  \brief General file for the quasi random nb generator
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _INGPNUMLIB_QUASIRANDOM_H
#define _INGPNUMLIB_QUASIRANDOM_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpbase/gpvector.h"
#include "random.h"
#include <vector>
CC_USING_NS(std,vector)


CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////
/// \base class ARM_QuasiRandom
///		abstract class for quasi random 
///		number generators
//////////////////////////////
class ARM_QuasiRandom : public ARM_RandomGenerator
{
public:
	ARM_QuasiRandom(int firstSimulations);
	virtual void SetDim( size_t dim ) = 0;
	inline virtual size_t dim() const { return itsDim; }
	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb );
	virtual void reset( const ARM_GP_T_Vector<size_t>& nbOfPointsList, const ARM_GP_T_Vector<size_t>& factorNb );
	virtual ARM_DistributionType GetDistributionType() const { return ARM_Uniform; }

protected:
	/// protected for easier and faster access from derived classes!
	/// save the current values and current position
	// Array of random values
	vector<vector<vector<double> > > itsCurrentValues;
	// Current position in the bucket
	int itsCurrentBucketPos;
	// Current position in the random vector
	int itsCurrentDimPos;
	// Index of the time
	int itsCurrentTimeIndex;
	// Current position of the current random vector
	int itsCurrentPos;
	// Current factor pos
	int itsCurrentFactorPos;
	// Number of random vector return by the generator
	ARM_GP_T_Vector<size_t> itsNbOfPointsList;
	// Dimentsion of the random vector
	int itsDim;
	// Nb of factors
	int itsFactorNb;
	// Skipped first simulations
	int itsFirstSimulations;

	ARM_GP_T_Vector<size_t> itsStepFactorNb;

private:
	virtual void DrawAll() = 0;

	/// function to do a range check and a dimension check!
	virtual double DrawOne();
	virtual void ValidateResult( double result ) = 0;
	virtual void ValidateDim( size_t dim );
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
