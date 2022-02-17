/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file mersenetwister.h
 *  \brief General file for the Mersene Twister
 *	random nb generator
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */

#ifndef _INGPNUMLIB_MERSENETWISTER_H
#define _INGPNUMLIB_MERSENETWISTER_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpnumlib/uniform.h"
#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////////////////////
/// \class ARM_FastMersenneTwister 
/// fast implementation of the Mersenne Twister
/// the property of the mersenne twister is to be a one dimensional 
/// random number generator with a period of 2^19937 - 1
///	and gives a sequence that is 623-dimensionally equidistributed
/// for more details see the cpp file or
/// http://www.math.keio.ac.jp/~matumoto/MT2002/emt19937ar.html
////////////////////////////////////////////////////////////////////

#define N 624	/// length of itsState vector

class ARM_FastMersenneTwister : public ARM_UniformGenerator 
{
public:
	ARM_FastMersenneTwister(const ARM_FastMersenneTwister& rhs );
	ARM_FastMersenneTwister& operator=(const ARM_FastMersenneTwister& rhs );
	virtual ~ARM_FastMersenneTwister();

	ARM_FastMersenneTwister( unsigned long seed=4357U );
	virtual string toString(const string& indent="", const string& nextIndent="") const{ return "Fast Mersenne Twister Uniform Random Nb Generator"; }
	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb );

	/// standard ARM support
	virtual ARM_Object* Clone() const;

private:
	unsigned long itsSeed;
	vector<unsigned long> itsState;				///the array for the itsState vector
	int left;
	int initf;
	unsigned long *next;

	void next_state();
	void init_genrand(unsigned long s);	/// function to reload the generator!
	virtual double DrawOne();
};


/// the original Mersenne Twister algorithm
class ARM_OriginalMersenneTwister : public ARM_UniformGenerator 
{
public:
	ARM_OriginalMersenneTwister(const ARM_OriginalMersenneTwister& rhs );
	ARM_OriginalMersenneTwister& operator=(const ARM_OriginalMersenneTwister& rhs );
	virtual ~ARM_OriginalMersenneTwister();

	ARM_OriginalMersenneTwister( unsigned long seed=4357U );
	virtual string toString(const string& indent="", const string& nextIndent="") const{ return "Original Mersenne Twister Uniform Random Nb Generator"; }
	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb );

	/// standard ARM support
	virtual ARM_Object* Clone() const;

private:
	unsigned long itsSeed;
	vector<unsigned long> mt;				/// the array for the itsState vector
	int mti;							/// next seed
	void init_genrand(unsigned long s);	/// function to reload the generator!
	virtual double DrawOne();
};

#undef N


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
