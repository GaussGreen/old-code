/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file hammersley.h
 *
 *  \brief Hammersley low discrepancy sequence
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date April 2004
 */

#ifndef _INGPNUMLIB_HAMMERSLEY_H
#define _INGPNUMLIB_HAMMERSLEY_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "quasirandom.h"

CC_BEGIN_NAMESPACE( ARM )


//////////////////////////////////////////////////////////////////////////////////
///  hammersley multidimentional low-discrepancy sequence generator
//////////////////////////////////////////////////////////////////////////////////
class ARM_Hammersley : public ARM_QuasiRandom
{
public:
	ARM_Hammersley(int firstSimulations);
	ARM_Hammersley(const ARM_Hammersley& rhs);
	ARM_Hammersley& operator=(const ARM_Hammersley& rhs );

	virtual ~ARM_Hammersley();
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual void SetDim( size_t dim );
	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorNb );

	/// standard ARM support
	virtual ARM_Object* Clone() const;

	/// validation function
	virtual void ValidateResult( double result ) {};

private:
	int i_modp(int i, int j);
	virtual void DrawAll();

	void initBase( size_t dim );

	int itsCountor;
	int itsTotalNbOfPoints;
	vector<int> itsBase;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
