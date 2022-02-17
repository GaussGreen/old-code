/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file transposer.h
 *
 *  \brief Generator transposer
 *	\author  R. GUILLEMOT
 *	\version 1.0
 *	\date February 2005
 */

#ifndef _INGPNUMLIB_TRANSPOSER_H
#define _INGPNUMLIB_TRANSPOSER_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpbase/gpvector.h"
#include "manipulator.h"
#include "typedef.h"

#include <string>
CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////
/// Transpose the generator 
////////////////////////////////////////////////

struct ARM_Transposer : public ARM_ManipulatorMethod
{
public:
	enum TransposeOrder {BucketOrder, PathOrder};

	ARM_Transposer(const ARM_RandomGeneratorPtr&, TransposeOrder order, int firstSimulations);
	ARM_Transposer(const ARM_Transposer& rhs);
	ARM_Transposer& operator=(const ARM_Transposer& rhs);
	virtual ARM_ManipulatorMethod* Clone() const;
	virtual ~ARM_Transposer();

//	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb );
	virtual double operator()() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const;

//	ARM_RandomGeneratorPtr GetBaseGen() const {return itsRandomGen;}

protected:
	virtual void resetwork(const ARM_GP_T_Vector<size_t>& nbOfPointsList, const ARM_GP_T_Vector<size_t>& factorNb);

private:
	// Base random generator
//	ARM_RandomGeneratorPtr itsRandomGen;

	// Array of random values
	vector<vector<vector<double> > > itsCurrentValues;
	// Current position of the current bucket
	CC_IS_MUTABLE size_t itsCurrentBucketPos;
	// Current position in the random vector
	CC_IS_MUTABLE size_t itsCurrentDimPos;
	// Current position of the current random vector
	CC_IS_MUTABLE size_t itsCurrentPos;
	

	// Number of random vector return by the generator
	ARM_GP_T_Vector<size_t> itsNbOfPointsList;
	// Time Dimentsion of the random vector
	size_t itsDim;
	// Factor dimension of the random vector
	size_t itsNbFactor;
	// Transposition order
	TransposeOrder itsOrder;
	// Skipped first simulations
	int itsFirstSimulations;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
