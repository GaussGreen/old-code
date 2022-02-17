/*!
 *
 * Copyright (c) CDC IXIS CM March 2005 Paris
 *
 *	\file newsobol.h
 *
 *  \brief General file for the sobol sequence
 *	\author  A. Triki
 *	\version 1.0
 *	\date March 2005
 */

#ifndef _INGPNUMLIB_NEWSOBOL_H
#define _INGPNUMLIB_NEWSOBOL_H

#include "gpbase/port.h"
#include "quasirandom.h"

CC_BEGIN_NAMESPACE( ARM )


class ARM_NewSobol : public ARM_QuasiRandom
{
public:
	ARM_NewSobol(long itsSeed, int firstSimulations);
	ARM_NewSobol( const ARM_NewSobol& rhs );
	ARM_NewSobol& operator=(const ARM_NewSobol& rhs );
	virtual ~ARM_NewSobol();

	virtual void SetDim( size_t dim );
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual void reset( size_t dim, size_t nbOfPoints );

	/// standard ARM support
	virtual ARM_Object* Clone() const;

	/// validation function
	virtual void ValidateResult( double result ) {};

private:

/*new functions*/
	double d_uniform_01 ( int *seed );
	int i4_bit_hi1 ( int n );
	int i4_bit_lo0 ( int n );
	void i4_sobol ( int dim_num, int *seed, float quasi[ ] );

	virtual void DrawAll();
	void Init();

	// key variables for the computation of Sobol sequences
	/// index_ represents the number of draws already done
	unsigned long itsIndex; 
	vector<unsigned long> itsValues;
	long itsSeed;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
