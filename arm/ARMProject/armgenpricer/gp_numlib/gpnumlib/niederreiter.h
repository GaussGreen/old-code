/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file niederreiter.h
 *
 *  \brief Niederreiter low discrepancy sequence
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPNUMLIB_NIEDERREITER_H
#define _INGPNUMLIB_NIEDERREITER_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "quasirandom.h"

CC_BEGIN_NAMESPACE( ARM )


//////////////////////////////////////////////////////////////////////////////////
///  Niederreiter multidimentional low-discrepancy sequence generator
//////////////////////////////////////////////////////////////////////////////////
class ARM_Niederreiter : public ARM_QuasiRandom
{
public:
	ARM_Niederreiter(int firstSimulations);
	ARM_Niederreiter(const ARM_Niederreiter& rhs);
	ARM_Niederreiter& operator=(const ARM_Niederreiter& rhs );

	virtual ~ARM_Niederreiter();
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual void SetDim( size_t dim );

	/// standard ARM support
	virtual ARM_Object* Clone() const;

	/// validation function
	virtual void ValidateResult( double result ) {};

private:
	virtual void DrawAll();
	bool itsFirstUse;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
