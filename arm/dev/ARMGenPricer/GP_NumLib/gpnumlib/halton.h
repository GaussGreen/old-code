/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file halton.h
 *
 *  \brief Halton low discrepancy sequence
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPNUMLIB_HALTON_H
#define _INGPNUMLIB_HALTON_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "quasirandom.h"

CC_BEGIN_NAMESPACE( ARM )


//////////////////////////////////////////////////////////////////////////////////
///  Halton multidimentional low-discrepancy sequence generator
//////////////////////////////////////////////////////////////////////////////////
class ARM_Halton : public ARM_QuasiRandom
{
public:
	ARM_Halton(int firstSimulations);
	ARM_Halton(const ARM_Halton& rhs);
	ARM_Halton& operator=(const ARM_Halton& rhs );

	virtual ~ARM_Halton();
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual void SetDim( size_t dim );

	/// standard ARM support
	virtual ARM_Object* Clone() const;

	/// validation function
	virtual void ValidateResult( double result ) {};

private:
	virtual void DrawAll();
	int itsCountor;
	vector<int> itsBase;
	void initBase( size_t dim );
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
